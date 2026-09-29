#include "util/json.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

#include "util/common.h"

namespace pramana {

Json::Json(const std::vector<double>& v) : type_(Type::Array) {
  arr_.reserve(v.size());
  for (double d : v) arr_.emplace_back(d);
}

Json::Json(const std::vector<int>& v) : type_(Type::Array) {
  arr_.reserve(v.size());
  for (int d : v) arr_.emplace_back(d);
}

Json& Json::operator[](const std::string& key) {
  if (type_ == Type::Null) type_ = Type::Object;
  PRAMANA_CHECK(type_ == Type::Object, "json: not an object");
  for (auto& kv : obj_)
    if (kv.first == key) return kv.second;
  obj_.emplace_back(key, Json());
  return obj_.back().second;
}

const Json& Json::at(const std::string& key) const {
  static const Json kNull;
  for (auto& kv : obj_)
    if (kv.first == key) return kv.second;
  return kNull;
}

bool Json::has(const std::string& key) const {
  for (auto& kv : obj_)
    if (kv.first == key) return true;
  return false;
}

void Json::push(Json v) {
  if (type_ == Type::Null) type_ = Type::Array;
  PRAMANA_CHECK(type_ == Type::Array, "json: not an array");
  arr_.push_back(std::move(v));
}

size_t Json::size() const {
  if (type_ == Type::Array) return arr_.size();
  if (type_ == Type::Object) return obj_.size();
  return 0;
}

std::vector<double> Json::numVector() const {
  std::vector<double> out;
  for (auto& v : arr_) out.push_back(v.num());
  return out;
}

static void escapeString(std::string& out, const std::string& s) {
  out.push_back('"');
  for (unsigned char c : s) {
    switch (c) {
      case '"': out += "\\\""; break;
      case '\\': out += "\\\\"; break;
      case '\n': out += "\\n"; break;
      case '\r': out += "\\r"; break;
      case '\t': out += "\\t"; break;
      default:
        if (c < 0x20) {
          char buf[8];
          std::snprintf(buf, sizeof buf, "\\u%04x", c);
          out += buf;
        } else {
          out.push_back(static_cast<char>(c));
        }
    }
  }
  out.push_back('"');
}

static void dumpNumber(std::string& out, double v) {
  if (std::isnan(v)) { out += "null"; return; }
  if (std::isinf(v)) { out += v > 0 ? "1e308" : "-1e308"; return; }  // JSON has no inf
  if (v == std::floor(v) && std::fabs(v) < 1e15) {
    char buf[32];
    std::snprintf(buf, sizeof buf, "%.0f", v);
    out += buf;
    return;
  }
  char buf[40];
  std::snprintf(buf, sizeof buf, "%.17g", v);  // round-trip exact
  out += buf;
}

void Json::dumpTo(std::string& out, int indent, int depth) const {
  auto newline = [&](int d) {
    if (indent < 0) return;
    out.push_back('\n');
    out.append(static_cast<size_t>(indent * d), ' ');
  };
  switch (type_) {
    case Type::Null: out += "null"; break;
    case Type::Bool: out += bool_ ? "true" : "false"; break;
    case Type::Number: dumpNumber(out, num_); break;
    case Type::String: escapeString(out, str_); break;
    case Type::Array: {
      out.push_back('[');
      bool scalarArray = true;
      for (auto& v : arr_)
        if (v.type_ == Type::Array || v.type_ == Type::Object) scalarArray = false;
      for (size_t i = 0; i < arr_.size(); ++i) {
        if (i) out.push_back(',');
        if (!scalarArray) newline(depth + 1);
        arr_[i].dumpTo(out, indent, depth + 1);
      }
      if (!scalarArray && !arr_.empty()) newline(depth);
      out.push_back(']');
      break;
    }
    case Type::Object: {
      out.push_back('{');
      for (size_t i = 0; i < obj_.size(); ++i) {
        if (i) out.push_back(',');
        newline(depth + 1);
        escapeString(out, obj_[i].first);
        out += indent >= 0 ? ": " : ":";
        obj_[i].second.dumpTo(out, indent, depth + 1);
      }
      if (!obj_.empty()) newline(depth);
      out.push_back('}');
      break;
    }
  }
}

std::string Json::dump(int indent) const {
  std::string out;
  dumpTo(out, indent, 0);
  return out;
}

namespace {
struct Parser {
  const std::string& s;
  size_t p = 0;
  void ws() { while (p < s.size() && std::isspace(static_cast<unsigned char>(s[p]))) ++p; }
  [[noreturn]] void fail(const char* what) {
    throw PramanaError(std::string("json parse error: ") + what + " at offset " + std::to_string(p));
  }
  Json value() {
    ws();
    if (p >= s.size()) fail("unexpected end");
    char c = s[p];
    if (c == '{') {
      ++p;
      Json o = Json::object();
      ws();
      if (p < s.size() && s[p] == '}') { ++p; return o; }
      for (;;) {
        ws();
        if (s[p] != '"') fail("expected key");
        std::string k = str();
        ws();
        if (s[p] != ':') fail("expected ':'");
        ++p;
        o[k] = value();
        ws();
        if (s[p] == ',') { ++p; continue; }
        if (s[p] == '}') { ++p; return o; }
        fail("expected ',' or '}'");
      }
    }
    if (c == '[') {
      ++p;
      Json a = Json::array();
      ws();
      if (p < s.size() && s[p] == ']') { ++p; return a; }
      for (;;) {
        a.push(value());
        ws();
        if (s[p] == ',') { ++p; continue; }
        if (s[p] == ']') { ++p; return a; }
        fail("expected ',' or ']'");
      }
    }
    if (c == '"') return Json(str());
    if (s.compare(p, 4, "true") == 0) { p += 4; return Json(true); }
    if (s.compare(p, 5, "false") == 0) { p += 5; return Json(false); }
    if (s.compare(p, 4, "null") == 0) { p += 4; return Json(); }
    char* end = nullptr;
    double v = std::strtod(s.c_str() + p, &end);
    if (end == s.c_str() + p) fail("bad token");
    p = static_cast<size_t>(end - s.c_str());
    return Json(v);
  }
  std::string str() {
    ++p;  // opening quote
    std::string out;
    while (p < s.size() && s[p] != '"') {
      if (s[p] == '\\') {
        ++p;
        char e = s[p];
        if (e == 'n') out.push_back('\n');
        else if (e == 't') out.push_back('\t');
        else if (e == 'r') out.push_back('\r');
        else if (e == 'u') { out.push_back('?'); p += 4; }
        else out.push_back(e);
        ++p;
      } else {
        out.push_back(s[p++]);
      }
    }
    ++p;
    return out;
  }
};
}  // namespace

Json Json::parse(const std::string& text) {
  Parser ps{text};
  return ps.value();
}

std::string readTextFile(const std::string& path) {
  std::ifstream f(path, std::ios::binary);
  if (!f) throw PramanaError("cannot open " + path);
  std::ostringstream ss;
  ss << f.rdbuf();
  return ss.str();
}

void writeTextFile(const std::string& path, const std::string& text) {
  std::ofstream f(path, std::ios::binary);
  if (!f) throw PramanaError("cannot write " + path);
  f << text;
}

}  // namespace pramana
