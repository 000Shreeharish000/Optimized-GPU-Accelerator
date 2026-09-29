// Tests: MPS I/O, sparse LU, certifier / interval arithmetic, LDL^T, convexity.
#include <algorithm>
#include <cmath>
#include <string>

#include "cert/certifier.h"
#include "cert/interval.h"
#include "io/mps.h"
#include "ipm/ipm.h"
#include "ipm/ldl.h"
#include "lp/lu.h"
#include "testing.h"
#include "util/rng.h"

using namespace pramana;

static const char* kRangesMps = R"(NAME RNGTEST
OBJSENSE
    MAX
ROWS
 N  OBJ
 E  EPOS
 E  ENEG
 L  LR
 G  GR
COLUMNS
    MARKER  'MARKER'  'INTORG'
    X  OBJ  1  EPOS  1
    X  ENEG  1
    MARKER  'MARKER'  'INTEND'
    Y  OBJ  2  LR  1
    Y  GR  1
    Z  OBJ  -1  GR  1
RHS
    RHS  OBJ  -5  EPOS  4
    RHS  ENEG  4  LR  10
    RHS  GR  2
RANGES
    R  EPOS  3  ENEG  -3
    R  LR  4  GR  -6
BOUNDS
 UP BND  Z  -2
 UP BND  Y  8
QUADOBJ
    Y  Z  0.5
ENDATA
)";

PTEST(mps_ranges_bounds_markers) {
  Model m = readMpsFromString(kRangesMps);
  EXPECT(m.sense == ObjSense::Maximize);
  EXPECT(m.numRows() == 4 && m.numCols() == 3);
  EXPECT_NEAR(m.objOffset, 5.0, 0);                 // RHS on objective = -offset
  EXPECT_NEAR(m.rowLower[0], 4, 0);                 // E with R>0: [rhs, rhs+R]
  EXPECT_NEAR(m.rowUpper[0], 7, 0);
  EXPECT_NEAR(m.rowLower[1], 1, 0);                 // E with R<0: [rhs+R, rhs]
  EXPECT_NEAR(m.rowUpper[1], 4, 0);
  EXPECT_NEAR(m.rowLower[2], 6, 0);                 // L: [rhs-|R|, rhs]
  EXPECT_NEAR(m.rowUpper[3], 8, 0);                 // G: [rhs, rhs+|R|]
  EXPECT(m.colType[0] == VarType::Integer);
  EXPECT(std::isinf(m.colUpper[0]));                // INTORG default upper = +inf
  EXPECT(std::isinf(m.colLower[2]) && m.colLower[2] < 0);  // negative UP => lower = -inf
  EXPECT(m.Q.nnz() == 2);                           // off-diagonal mirrored
  // Round trip through the writer.
  Model r = readMpsFromString(writeMpsString(m));
  EXPECT(r.numRows() == m.numRows() && r.numCols() == m.numCols() && r.nnz() == m.nnz());
  for (int i = 0; i < m.numRows(); ++i) {
    EXPECT(r.rowLower[i] == m.rowLower[i]);
    EXPECT(r.rowUpper[i] == m.rowUpper[i]);
  }
  EXPECT(r.colLower == m.colLower && r.colUpper == m.colUpper && r.colCost == m.colCost);
  EXPECT(r.objOffset == m.objOffset && r.sense == m.sense && r.Q.nnz() == m.Q.nnz());
}

PTEST(mps_reads_gzip_netlib) {
  Model m = readMps(dataPath("data/netlib/afiro.mps.gz"));
  EXPECT(m.numRows() == 27 && m.numCols() == 32 && m.nnz() == 83);
}

PTEST(mps_fixed_format_names_with_spaces) {
  Model m = readMps(dataPath("data/netlib/forplan.mps.gz"));
  EXPECT(m.numRows() == 161 && m.numCols() == 421);
}

static SparseMatrix randomSparse(Rng& rng, int m, int n, double density) {
  std::vector<int> ri, ci;
  std::vector<double> v;
  for (int j = 0; j < n; ++j) {
    ri.push_back(j % m);  // keep structurally nonsingular (diagonal)
    ci.push_back(j);
    v.push_back(2.0 + rng.uniform());
    for (int i = 0; i < m; ++i)
      if (rng.uniform() < density) {
        ri.push_back(i);
        ci.push_back(j);
        v.push_back(rng.uniform() * 2 - 1);
      }
  }
  return SparseMatrix::fromTriplets(m, n, ri, ci, v);
}

PTEST(lu_ftran_btran_residuals_and_updates) {
  Rng rng(3);
  const int m = 120;
  SparseMatrix A = randomSparse(rng, m, m, 0.03);
  LuFactor lu;
  lu.setup(m, m, &A);
  std::vector<int> basis(m);
  for (int p = 0; p < m; ++p) basis[p] = p;
  EXPECT(lu.factorize(basis) == 0);
  // B x = b
  HVector v;
  v.setup(m);
  std::vector<double> b(m);
  for (int i = 0; i < m; ++i) v.array[i] = b[i] = rng.uniform();
  lu.ftran(v);
  std::vector<double> Bx(m, 0.0);
  A.multiply(v.array.data(), Bx.data());
  for (int i = 0; i < m; ++i) EXPECT_NEAR(Bx[i], b[i], 1e-10);
  // B^T y = e
  HVector w;
  w.setup(m);
  for (int p = 0; p < m; ++p) w.array[p] = b[p];
  lu.btran(w);
  std::vector<double> BTy(m, 0.0);
  A.multiplyTranspose(w.array.data(), BTy.data());
  for (int p = 0; p < m; ++p) EXPECT_NEAR(BTy[p], b[p], 1e-10);
  // Replace a column by a logical (e_i) through a PFI update and check again.
  HVector col;
  col.setup(m);
  col.array[5] = 1.0;
  lu.ftran(col);  // B^{-1} e_5
  int pos = 0;
  for (int p = 0; p < m; ++p)
    if (std::fabs(col.array[p]) > std::fabs(col.array[pos])) pos = p;
  lu.update(col, pos);
  basis[pos] = m + 5;
  HVector v2;
  v2.setup(m);
  for (int i = 0; i < m; ++i) v2.array[i] = b[i];
  lu.ftran(v2);
  LuFactor fresh;
  fresh.setup(m, m, &A);
  EXPECT(fresh.factorize(basis) == 0);
  HVector v3;
  v3.setup(m);
  for (int i = 0; i < m; ++i) v3.array[i] = b[i];
  fresh.ftran(v3);
  for (int p = 0; p < m; ++p) EXPECT_NEAR(v2.array[p], v3.array[p], 1e-9);
}

PTEST(lu_detects_rank_deficiency) {
  // Two identical columns => rank m-1.
  std::vector<int> ri{0, 1, 0, 1, 2}, ci{0, 0, 1, 1, 2};
  std::vector<double> v{1, 2, 1, 2, 3};
  SparseMatrix A = SparseMatrix::fromTriplets(3, 3, ri, ci, v);
  LuFactor lu;
  lu.setup(3, 3, &A);
  std::vector<int> basis{0, 1, 2};
  EXPECT(lu.factorize(basis) == 1);
  EXPECT(lu.deficientPositions.size() == 1 && lu.deficientRows.size() == 1);
}

PTEST(interval_arithmetic_is_rigorous) {
  // 0.1 + 0.2 is inexact in binary: bounds must straddle the exact decimal 0.3 sum of the doubles.
  double lo = addDown(0.1, 0.2), hi = addUp(0.1, 0.2);
  EXPECT(lo < hi);
  EXPECT(addDown(1.0, 2.0) == 3.0 && addUp(1.0, 2.0) == 3.0);  // exact => no widening
  EXPECT(mulDown(3.0, 0.1) <= 3.0 * 0.1 && mulUp(3.0, 0.1) >= 3.0 * 0.1);
  Interval z{-1, 2};
  EXPECT(std::isinf(minProduct(z, 0, kInf)));   // z can be negative, x -> +inf
  EXPECT(minProduct(Interval{1, 2}, 3, 5) == 3);
}

PTEST(safe_bound_valid_for_any_dual) {
  // min x + y  s.t. x + 2y >= 4, 3x + y >= 6, x,y in [0,10]  -> optimum 2.8 at (1.6, 1.2)
  Model m;
  m.addColumn(1, 0, 10, {}, {});
  m.addColumn(1, 0, 10, {}, {});
  m.addRow(4, kInf, {0, 1}, {1, 2});
  m.addRow(6, kInf, {0, 1}, {3, 1});
  // Optimal duals: y1 = 0.4, y2 = 0.2 -> bound equals the optimum.
  SafeBoundResult sb = safeDualBound(m.A, m.colCost, m.colLower, m.colUpper, m.rowLower, m.rowUpper, {0.4, 0.2});
  EXPECT_NEAR(sb.bound, 2.8, 1e-12);
  EXPECT(sb.bound <= 2.8);
  // Any other dual vector still yields a VALID (weaker) lower bound.
  Rng rng(1);
  for (int t = 0; t < 200; ++t) {
    std::vector<double> y{rng.uniform() * 2 - 0.5, rng.uniform() * 2 - 0.5};
    SafeBoundResult s2 = safeDualBound(m.A, m.colCost, m.colLower, m.colUpper, m.rowLower, m.rowUpper, y);
    EXPECT(s2.bound <= 2.8 + 1e-12);
  }
}

PTEST(certifier_rejects_false_claims) {
  Model m;
  m.addColumn(1, 0, 10, {}, {});
  m.addColumn(1, 0, 10, {}, {});
  m.addRow(4, kInf, {0, 1}, {1, 2});
  m.addRow(6, kInf, {0, 1}, {3, 1});
  SolutionClaim ok;
  ok.status = Status::Optimal;
  ok.x = {1.6, 1.2};
  ok.rowDual = {0.4, 0.2};
  EXPECT(certify(m, ok).accepted);
  SolutionClaim wrongX = ok;
  wrongX.x = {1.0, 1.0};  // infeasible
  EXPECT(!certify(m, wrongX).accepted);
  SolutionClaim subopt = ok;
  subopt.x = {2.0, 2.0};  // feasible but not optimal (4 vs 2.8)
  EXPECT(!certify(m, subopt).accepted);
  SolutionClaim nanX = ok;
  nanX.x = {std::nan(""), 1.0};
  EXPECT(!certify(m, nanX).accepted);
  SolutionClaim fakeInfeas;
  fakeInfeas.status = Status::Infeasible;
  fakeInfeas.farkas = {1.0, 1.0};
  EXPECT(!certify(m, fakeInfeas).accepted);  // the model is feasible: no ray can pass
}

PTEST(ldl_quasidefinite_solve) {
  Rng rng(9);
  const int n = 60, mm = 25;
  SparseMatrix A = randomSparse(rng, mm, n, 0.1);
  std::vector<int> ri, ci;
  std::vector<double> v;
  for (int j = 0; j < n; ++j) {
    ri.push_back(j);
    ci.push_back(j);
    v.push_back(-(1.0 + rng.uniform()));
    for (int k = A.start[j]; k < A.start[j + 1]; ++k) {
      ri.push_back(n + A.index[k]);
      ci.push_back(j);
      v.push_back(A.value[k]);
      ri.push_back(j);
      ci.push_back(n + A.index[k]);
      v.push_back(A.value[k]);
    }
  }
  for (int r = 0; r < mm; ++r) {
    ri.push_back(n + r);
    ci.push_back(n + r);
    v.push_back(1e-8);
  }
  SparseMatrix K = SparseMatrix::fromTriplets(n + mm, n + mm, ri, ci, v);
  std::vector<int8_t> signs(n + mm, 1);
  for (int j = 0; j < n; ++j) signs[j] = -1;
  LdlFactor f;
  f.analyze(K, minimumDegreeOrdering(K));
  f.factorize(K, signs, 1e-14);
  EXPECT(f.negativePivots() == n);  // inertia of a quasidefinite matrix
  std::vector<double> b(n + mm), x;
  for (auto& t : b) t = rng.uniform();
  x = b;
  f.solve(x);
  std::vector<double> Kx(n + mm, 0.0);
  K.multiply(x.data(), Kx.data());
  for (int i = 0; i < n + mm; ++i) EXPECT_NEAR(Kx[i], b[i], 1e-7);
}

PTEST(convexity_check_with_certificate) {
  // Q = [[1, 2],[2, 1]] is indefinite: v = (1,-1) has v'Qv = -2.
  SparseMatrix Q = SparseMatrix::fromTriplets(2, 2, {0, 1, 0, 1}, {0, 0, 1, 1}, {1, 2, 2, 1});
  ConvexityResult r = checkConvexity(Q);
  EXPECT(!r.convex);
  std::vector<double> Qv(2, 0.0);
  Q.multiply(r.certificate.data(), Qv.data());
  EXPECT(r.certificate[0] * Qv[0] + r.certificate[1] * Qv[1] < 0);  // verifiable evidence
  SparseMatrix P = SparseMatrix::fromTriplets(2, 2, {0, 1, 0, 1}, {0, 0, 1, 1}, {2, 1, 1, 2});
  EXPECT(checkConvexity(P).convex);
}

PTEST(lu_hypersparse_matches_dense) {
  Rng rng(21);
  const int m = 400;
  SparseMatrix A = randomSparse(rng, m, m, 0.004);
  LuFactor lu;
  lu.setup(m, m, &A);
  std::vector<int> basis(m);
  for (int p = 0; p < m; ++p) basis[p] = p;
  EXPECT(lu.factorize(basis) == 0);
  for (int trial = 0; trial < 40; ++trial) {
    if (trial == 20) {  // also exercise the PFI eta path
      HVector c;
      c.setup(m);
      c.array[3] = 1.0;
      c.index[0] = 3;
      c.count = 1;
      lu.ftran(c, true);
      int pos = c.index[0];
      for (int k = 0; k < c.count; ++k)
        if (std::fabs(c.array[c.index[k]]) > std::fabs(c.array[pos])) pos = c.index[k];
      lu.update(c, pos);
    }
    HVector s, d;
    s.setup(m);
    d.setup(m);
    int k1 = rng.below(m), k2 = rng.below(m);
    for (int k : {k1, k2}) {
      double v = rng.uniform() + 0.5;
      if (s.array[k] == 0.0) s.index[s.count++] = k;
      s.array[k] += v;
      d.array[k] += v;
    }
    lu.ftran(s, true);
    lu.ftran(d, false);
    for (int i = 0; i < m; ++i) EXPECT_NEAR(s.array[i], d.array[i], 1e-12);
    HVector bs, bd;
    bs.setup(m);
    bd.setup(m);
    bs.array[k1] = bd.array[k1] = 1.0;
    bs.index[0] = k1;
    bs.count = 1;
    lu.btran(bs, true);
    lu.btran(bd, false);
    for (int i = 0; i < m; ++i) EXPECT_NEAR(bs.array[i], bd.array[i], 1e-12);
  }
}
