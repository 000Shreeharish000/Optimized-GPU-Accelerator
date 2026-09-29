#include "io/mps.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <cstring>
#include <unordered_map>

#include "io/inflate.h"
#include "util/json.h"
#include "util/log.h"

namespace pramana {

namespace {

enum class Section { None, Name, ObjSense, Rows, Columns, Rhs, Ranges, Bounds, Quad, End };

struct Tokens {
  const char* tok[8];
  int len[8];
  int n = 0;
  std::string str(int i) const { return std::string(tok[i], len[i]); }
};

void tokenize(const char* s, const char* e, Tokens& t) {
  t.n = 0;
  while (s < e && t.n < 8) {
    while (s < e && (*s == ' ' || *s == '\t')) ++s;
    if (s >= e) break;
    const char* b = s;
    while (s < e && *s != ' ' && *s != '\t') ++s;
    t.tok[t.n] = b;
    t.len[t.n] = static_cast<int>(s - b);
    ++t.n;
  }
}

std::string trimCopy(const std::string& s) {
  size_t b = s.find_first_not_of(" \t");
  if (b == std::string::npos) return "";
  size_t e = s.find_last_not_of(" \t");
  return s.substr(b, e - b + 1);
}

// Fixed-format fields (1-based columns): 2-3, 5-12, 15-22, 25-36, 40-47, 50-61.
void fixedFields(const char* s, const char* e, std::string f[6]) {
  static const int beg[6] = {1, 4, 14, 24, 39, 49};
  static const int len[6] = {2, 8, 8, 12, 8, 12};
  int L = static_cast<int>(e - s);
  for (int k = 0; k < 6; ++k) {
    if (beg[k] >= L) { f[k].clear(); continue; }
    int l = std::min(len[k], L - beg[k]);
    f[k] = trimCopy(std::string(s + beg[k], l));
  }
}

double parseNumber(const std::string& s, int line) {
  char* end = nullptr;
  double v = std::strtod(s.c_str(), &end);
  if (end == s.c_str()) throw PramanaError("MPS line " + std::to_string(line) + ": bad number '" + s + "'");
  // Fortran-style exponents like 1.0D+03 appear in some legacy files.
  if (*end == 'D' || *end == 'd') {
    std::string t = s;
    for (char& c : t)
      if (c == 'D' || c == 'd') c = 'E';
    v = std::strtod(t.c_str(), nullptr);
  }
  return v;
}

struct Reader {
  const MpsReadOptions& opts;
  Model m;
  std::unordered_map<std::string, int> rowIndex;  // -1 = objective, -2 = other free N rows
  std::unordered_map<std::string, int> colIndex;
  std::vector<char> rowType;
  std::vector<double> rhs;
  std::vector<int> triRow, triCol;
  std::vector<double> triVal;
  std::vector<int> qRow, qCol;
  std::vector<double> qVal;
  std::vector<char> hasLowerBound, hasUpperBound, colIsIntMarker;
  bool integerMarker = false;
  std::string objRowName;
  bool objSenseMax = false;
  int lineNo = 0;
  bool quadIsFullMatrix = false;
  bool fixedFormat = false;  // names contain blanks => parse data lines by columns

  explicit Reader(const MpsReadOptions& o) : opts(o) {}

  int getCol(const std::string& name, bool create) {
    auto it = colIndex.find(name);
    if (it != colIndex.end()) return it->second;
    if (!create) throw PramanaError("MPS line " + std::to_string(lineNo) + ": unknown column " + name);
    int j = static_cast<int>(m.colCost.size());
    colIndex.emplace(name, j);
    m.colNames.push_back(name);
    m.colCost.push_back(0);
    m.colLower.push_back(0);
    m.colUpper.push_back(kInf);
    m.colType.push_back(integerMarker ? VarType::Integer : VarType::Continuous);
    hasLowerBound.push_back(0);
    hasUpperBound.push_back(0);
    colIsIntMarker.push_back(integerMarker ? 1 : 0);
    return j;
  }

  int findRow(const std::string& name) {
    auto it = rowIndex.find(name);
    if (it == rowIndex.end())
      throw PramanaError("MPS line " + std::to_string(lineNo) + ": unknown row " + name);
    return it->second;
  }

  void addEntry(int col, const std::string& rowName, double v) {
    int r = findRow(rowName);
    if (r == -1) {
      m.colCost[col] += v;
    } else if (r >= 0) {
      triRow.push_back(r);
      triCol.push_back(col);
      triVal.push_back(v);
    }
  }

  void setRhs(const std::string& rowName, double v) {
    int r = findRow(rowName);
    if (r == -1) m.objOffset = -v;  // RHS on objective row: constant term is -rhs
    else if (r >= 0) rhs[r] = v;
  }

  void setRange(const std::string& rowName, double v) {
    int r = findRow(rowName);
    if (r < 0) return;
    double b = rhs[r];
    switch (rowType[r]) {
      case 'L': m.rowLower[r] = b - std::fabs(v); break;
      case 'G': m.rowUpper[r] = b + std::fabs(v); break;
      case 'E':
        if (v >= 0) m.rowUpper[r] = b + v;
        else m.rowLower[r] = b + v;
        break;
    }
    rangeSet[r] = 1;
  }
  std::vector<char> rangeSet;

  void setBound(const std::string& type, const std::string& colName, double v, bool hasValue) {
    int j = getCol(colName, false);
    std::string t = type;
    for (char& c : t) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    if (t == "UP") {
      m.colUpper[j] = v;
      hasUpperBound[j] = 1;
      // Classic MPS rule: negative upper bound with no explicit lower bound => lower = -inf.
      if (v < 0 && !hasLowerBound[j] && m.colLower[j] == 0) {
        m.colLower[j] = -kInf;
      }
    } else if (t == "LO") {
      m.colLower[j] = v;
      hasLowerBound[j] = 1;
    } else if (t == "FX") {
      m.colLower[j] = m.colUpper[j] = v;
      hasLowerBound[j] = hasUpperBound[j] = 1;
    } else if (t == "FR") {
      m.colLower[j] = -kInf;
      m.colUpper[j] = kInf;
      hasLowerBound[j] = hasUpperBound[j] = 1;
    } else if (t == "MI") {
      m.colLower[j] = -kInf;
      hasLowerBound[j] = 1;
    } else if (t == "PL") {
      m.colUpper[j] = kInf;
      hasUpperBound[j] = 1;
    } else if (t == "BV") {
      m.colType[j] = VarType::Integer;
      m.colLower[j] = 0;
      m.colUpper[j] = 1;
      hasLowerBound[j] = hasUpperBound[j] = 1;
    } else if (t == "LI") {
      m.colType[j] = VarType::Integer;
      m.colLower[j] = v;
      hasLowerBound[j] = 1;
    } else if (t == "UI") {
      m.colType[j] = VarType::Integer;
      m.colUpper[j] = v;
      hasUpperBound[j] = 1;
      if (v < 0 && !hasLowerBound[j] && m.colLower[j] == 0) m.colLower[j] = -kInf;
    } else if (t == "SC") {
      throw PramanaError("MPS: semi-continuous (SC) bounds are not supported");
    } else {
      throw PramanaError("MPS line " + std::to_string(lineNo) + ": unknown bound type " + type);
    }
    (void)hasValue;
  }

  void parse(const std::string& text) {
    const char* p = text.c_str();
    const char* end = p + text.size();
    Section sec = Section::None;
    Tokens t;
    while (p < end) {
      const char* ls = p;
      while (p < end && *p != '\n') ++p;
      const char* le = p;
      if (p < end) ++p;
      if (le > ls && le[-1] == '\r') --le;
      ++lineNo;
      if (le == ls || *ls == '*') continue;
      bool header = !(*ls == ' ' || *ls == '\t');
      tokenize(ls, le, t);
      if (t.n == 0) continue;
      if (header) {
        std::string h = t.str(0);
        for (char& c : h) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        if (h == "NAME") {
          sec = Section::Name;
          if (t.n > 1) m.name = std::string(t.tok[1], le - t.tok[1]);
          m.name = trimCopy(m.name);
        } else if (h == "OBJSENSE") {
          sec = Section::ObjSense;
          if (t.n > 1) {
            std::string s = t.str(1);
            if (s == "MAX" || s == "MAXIMIZE") objSenseMax = true;
          }
        } else if (h == "OBJSENSE" || h == "OBJSENCE") {
          sec = Section::ObjSense;
        } else if (h == "ROWS") {
          sec = Section::Rows;
        } else if (h == "COLUMNS") {
          sec = Section::Columns;
          rhs.assign(m.rowLower.size(), 0.0);
          rangeSet.assign(m.rowLower.size(), 0);
        } else if (h == "RHS") {
          sec = Section::Rhs;
        } else if (h == "RANGES") {
          sec = Section::Ranges;
          applyRhs();
        } else if (h == "BOUNDS") {
          sec = Section::Bounds;
          applyRhs();
        } else if (h == "QUADOBJ") {
          sec = Section::Quad;
          quadIsFullMatrix = false;
          applyRhs();
        } else if (h == "QMATRIX" || h == "QSECTION") {
          sec = Section::Quad;
          quadIsFullMatrix = true;
          applyRhs();
        } else if (h == "ENDATA") {
          sec = Section::End;
          break;
        } else if (h == "OBJNAME") {
          // ignore
        } else if (sec == Section::ObjSense) {
          if (h == "MAX" || h == "MAXIMIZE") objSenseMax = true;
        } else {
          throw PramanaError("MPS line " + std::to_string(lineNo) + ": unknown section " + h);
        }
        continue;
      }
      switch (sec) {
        case Section::ObjSense: {
          std::string s = t.str(0);
          for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
          if (s == "MAX" || s == "MAXIMIZE") objSenseMax = true;
          break;
        }
        case Section::Rows: {
          std::string type = t.str(0);
          std::string name;
          if (t.n > 2) fixedFormat = true;
          if (fixedFormat) {
            std::string f[6];
            fixedFields(ls, le, f);
            type = f[0].empty() ? type : f[0];
            name = f[1];
          } else {
            name = t.n >= 2 ? t.str(1) : "";
          }
          char c = static_cast<char>(std::toupper(static_cast<unsigned char>(type[0])));
          if (c == 'N') {
            if (objRowName.empty()) {
              objRowName = name;
              m.objName = name;
              rowIndex[name] = -1;
            } else {
              rowIndex[name] = -2;  // extra free rows are dropped
            }
          } else {
            int r = static_cast<int>(m.rowLower.size());
            rowIndex[name] = r;
            m.rowNames.push_back(name);
            rowType.push_back(c);
            m.rowLower.push_back(0);
            m.rowUpper.push_back(0);
            if (c != 'E' && c != 'L' && c != 'G')
              throw PramanaError("MPS line " + std::to_string(lineNo) + ": bad row type " + type);
          }
          break;
        }
        case Section::Columns: {
          if (t.n >= 3 && (t.str(1) == "'MARKER'" || t.str(1) == "MARKER")) {
            std::string what = t.str(2);
            if (what.find("INTORG") != std::string::npos) integerMarker = true;
            else if (what.find("INTEND") != std::string::npos) integerMarker = false;
            break;
          }
          if (!fixedFormat && (t.n == 3 || t.n == 5)) {
            int j = getCol(t.str(0), true);
            addEntry(j, t.str(1), parseNumber(t.str(2), lineNo));
            if (t.n == 5) addEntry(j, t.str(3), parseNumber(t.str(4), lineNo));
          } else {
            std::string f[6];
            fixedFields(ls, le, f);
            int j = getCol(f[1], true);
            if (!f[2].empty()) addEntry(j, f[2], parseNumber(f[3], lineNo));
            if (!f[4].empty()) addEntry(j, f[4], parseNumber(f[5], lineNo));
          }
          break;
        }
        case Section::Rhs:
        case Section::Ranges: {
          auto apply = [&](const std::string& row, double v) {
            if (sec == Section::Rhs) setRhs(row, v);
            else setRange(row, v);
          };
          if (fixedFormat) {
            std::string f[6];
            fixedFields(ls, le, f);
            if (!f[2].empty()) apply(f[2], parseNumber(f[3], lineNo));
            if (!f[4].empty()) apply(f[4], parseNumber(f[5], lineNo));
          } else if (t.n == 3 || t.n == 5) {
            apply(t.str(1), parseNumber(t.str(2), lineNo));
            if (t.n == 5) apply(t.str(3), parseNumber(t.str(4), lineNo));
          } else if (t.n == 2 || t.n == 4) {  // set name omitted
            apply(t.str(0), parseNumber(t.str(1), lineNo));
            if (t.n == 4) apply(t.str(2), parseNumber(t.str(3), lineNo));
          } else {
            std::string f[6];
            fixedFields(ls, le, f);
            if (!f[2].empty()) apply(f[2], parseNumber(f[3], lineNo));
            if (!f[4].empty()) apply(f[4], parseNumber(f[5], lineNo));
          }
          break;
        }
        case Section::Bounds: {
          std::string type = t.str(0);
          std::string tu = type;
          for (char& c : tu) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
          bool noValue = tu == "FR" || tu == "MI" || tu == "PL" || tu == "BV";
          if (fixedFormat) {
            std::string f[6];
            fixedFields(ls, le, f);
            setBound(f[0], f[2], noValue || f[3].empty() ? 0.0 : parseNumber(f[3], lineNo), !noValue);
          } else if (noValue) {
            // "FR BND X" or "FR X" (set name omitted); BV may carry a value.
            std::string col = t.n >= 3 ? t.str(2) : t.str(1);
            if (t.n >= 3 && colIndex.find(col) == colIndex.end() && colIndex.count(t.str(1))) col = t.str(1);
            setBound(type, col, 0, false);
          } else if (t.n == 4) {
            setBound(type, t.str(2), parseNumber(t.str(3), lineNo), true);
          } else if (t.n == 3) {
            // set name omitted
            setBound(type, t.str(1), parseNumber(t.str(2), lineNo), true);
          } else {
            std::string f[6];
            fixedFields(ls, le, f);
            setBound(f[0], f[2], parseNumber(f[3], lineNo), true);
          }
          break;
        }
        case Section::Quad: {
          if (t.n < 3) throw PramanaError("MPS line " + std::to_string(lineNo) + ": bad quadratic entry");
          int a = getCol(t.str(0), false);
          int b = getCol(t.str(1), false);
          double v = parseNumber(t.str(2), lineNo);
          if (quadIsFullMatrix) {
            qRow.push_back(a);
            qCol.push_back(b);
            qVal.push_back(v);
          } else {
            qRow.push_back(a);
            qCol.push_back(b);
            qVal.push_back(v);
            if (a != b) {
              qRow.push_back(b);
              qCol.push_back(a);
              qVal.push_back(v);
            }
          }
          break;
        }
        default:
          throw PramanaError("MPS line " + std::to_string(lineNo) + ": data outside a section");
      }
    }
    applyRhs();
    finish();
  }

  bool rhsApplied = false;
  void applyRhs() {
    if (rhsApplied || rhs.empty()) return;
    rhsApplied = true;
    for (size_t r = 0; r < rhs.size(); ++r) {
      switch (rowType[r]) {
        case 'E': m.rowLower[r] = m.rowUpper[r] = rhs[r]; break;
        case 'L': m.rowLower[r] = -kInf; m.rowUpper[r] = rhs[r]; break;
        case 'G': m.rowLower[r] = rhs[r]; m.rowUpper[r] = kInf; break;
      }
    }
  }

  void finish() {
    if (rhs.empty()) {  // no COLUMNS section at all
      rhs.assign(m.rowLower.size(), 0.0);
      rangeSet.assign(m.rowLower.size(), 0);
    }
    applyRhs();
    int n = static_cast<int>(m.colCost.size());
    for (int j = 0; j < n; ++j) {
      if (colIsIntMarker[j] && !hasUpperBound[j] && m.colType[j] == VarType::Integer)
        m.colUpper[j] = opts.defaultIntegerUpper;
      if (m.colLower[j] <= -kInfBoundThreshold) m.colLower[j] = -kInf;
      if (m.colUpper[j] >= kInfBoundThreshold) m.colUpper[j] = kInf;
    }
    for (size_t r = 0; r < m.rowLower.size(); ++r) {
      if (m.rowLower[r] <= -kInfBoundThreshold) m.rowLower[r] = -kInf;
      if (m.rowUpper[r] >= kInfBoundThreshold) m.rowUpper[r] = kInf;
    }
    m.A = SparseMatrix::fromTriplets(static_cast<int>(m.rowLower.size()), n, triRow, triCol, triVal);
    if (!qVal.empty()) m.Q = SparseMatrix::fromTriplets(n, n, qRow, qCol, qVal);
    else m.Q = SparseMatrix(0, 0);
    if (objSenseMax) m.sense = ObjSense::Maximize;
  }
};

}  // namespace

Model readMpsFromString(const std::string& text, const MpsReadOptions& opts) {
  Reader r(opts);
  if (looksGzipped(text)) r.parse(gunzip(text));
  else r.parse(text);
  r.m.validate();
  return std::move(r.m);
}

Model readMps(const std::string& path, const MpsReadOptions& opts) {
  std::string text = readTextFile(path);
  Model m = readMpsFromString(text, opts);
  if (m.name.empty() || m.name == "model") {
    size_t s = path.find_last_of("/\\");
    std::string base = path.substr(s == std::string::npos ? 0 : s + 1);
    m.name = base.substr(0, base.find('.'));
  }
  return m;
}

static std::string num(double v) {
  char buf[40];
  std::snprintf(buf, sizeof buf, "%.17g", v);
  return buf;
}

std::string writeMpsString(const Model& model) {
  Model m = model;
  m.ensureNames();
  std::string out;
  out += "NAME " + m.name + "\n";
  if (m.sense == ObjSense::Maximize) out += "OBJSENSE\n    MAX\n";
  out += "ROWS\n N  " + m.objName + "\n";
  std::vector<char> type(m.numRows());
  for (int i = 0; i < m.numRows(); ++i) {
    double lo = m.rowLower[i], up = m.rowUpper[i];
    char c;
    if (lo == up) c = 'E';
    else if (isFiniteBound(lo) && isFiniteBound(up)) c = 'L';  // L + range
    else if (isFiniteBound(lo)) c = 'G';
    else if (isFiniteBound(up)) c = 'L';
    else c = 'N';
    if (c == 'N') c = 'G';  // free row: G with -inf rhs is not representable; use range trick below
    type[i] = c;
    out += std::string(" ") + c + "  " + m.rowNames[i] + "\n";
  }
  out += "COLUMNS\n";
  bool inInt = false;
  int marker = 0;
  for (int j = 0; j < m.numCols(); ++j) {
    bool isInt = m.colType[j] == VarType::Integer;
    if (isInt != inInt) {
      out += "    MARKER" + std::to_string(marker++) + " 'MARKER' " + (isInt ? "'INTORG'" : "'INTEND'") + "\n";
      inInt = isInt;
    }
    bool any = false;
    if (m.colCost[j] != 0) {
      out += "    " + m.colNames[j] + " " + m.objName + " " + num(m.colCost[j]) + "\n";
      any = true;
    }
    for (int k = m.A.start[j]; k < m.A.start[j + 1]; ++k) {
      out += "    " + m.colNames[j] + " " + m.rowNames[m.A.index[k]] + " " + num(m.A.value[k]) + "\n";
      any = true;
    }
    if (!any) out += "    " + m.colNames[j] + " " + m.objName + " 0\n";
  }
  if (inInt) out += "    MARKER" + std::to_string(marker++) + " 'MARKER' 'INTEND'\n";
  out += "RHS\n";
  if (m.objOffset != 0) out += "    RHS " + m.objName + " " + num(-m.objOffset) + "\n";
  std::string ranges;
  for (int i = 0; i < m.numRows(); ++i) {
    double lo = m.rowLower[i], up = m.rowUpper[i];
    double r = 0;
    if (type[i] == 'E') r = lo;
    else if (type[i] == 'L') r = isFiniteBound(up) ? up : 0;
    else r = isFiniteBound(lo) ? lo : -1e30;
    if (!isFiniteBound(lo) && !isFiniteBound(up)) r = -1e30;  // free row written as G >= -1e30
    if (r != 0) out += "    RHS " + m.rowNames[i] + " " + num(r) + "\n";
    if (type[i] == 'L' && isFiniteBound(lo) && isFiniteBound(up) && lo != up)
      ranges += "    RNG " + m.rowNames[i] + " " + num(up - lo) + "\n";
  }
  if (!ranges.empty()) out += "RANGES\n" + ranges;
  out += "BOUNDS\n";
  for (int j = 0; j < m.numCols(); ++j) {
    double lo = m.colLower[j], up = m.colUpper[j];
    const std::string& c = m.colNames[j];
    bool isInt = m.colType[j] == VarType::Integer;
    if (lo == up) { out += " FX BND " + c + " " + num(lo) + "\n"; continue; }
    if (!isFiniteBound(lo) && !isFiniteBound(up)) { out += " FR BND " + c + "\n"; continue; }
    if (!isFiniteBound(lo)) out += " MI BND " + c + "\n";
    else if (lo != 0) out += " LO BND " + c + " " + num(lo) + "\n";
    if (isFiniteBound(up)) out += " UP BND " + c + " " + num(up) + "\n";
    else if (isInt) out += " PL BND " + c + "\n";
  }
  if (m.isQp()) {
    out += "QUADOBJ\n";
    for (int j = 0; j < m.Q.numCols; ++j)
      for (int k = m.Q.start[j]; k < m.Q.start[j + 1]; ++k) {
        int i = m.Q.index[k];
        if (i < j) continue;  // lower triangle (i >= j), each entry once
        out += "    " + m.colNames[i] + " " + m.colNames[j] + " " + num(m.Q.value[k]) + "\n";
      }
  }
  out += "ENDATA\n";
  return out;
}

void writeMps(const Model& model, const std::string& path) { writeTextFile(path, writeMpsString(model)); }

}  // namespace pramana
