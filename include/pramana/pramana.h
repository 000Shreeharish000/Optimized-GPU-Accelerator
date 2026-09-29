/* PRAMANA C API - certified LP / MILP / convex QP solver.
 *
 * Every string returned by the library is owned by the object it came from
 * (model or result) unless documented otherwise. Options are passed as a JSON
 * object, e.g. {"algorithm":"auto","time_limit":60,"presolve":true}.
 * See docs/API.md for the full option list.
 */
#ifndef PRAMANA_H
#define PRAMANA_H

#ifdef _WIN32
#ifdef PRAMANA_BUILDING_DLL
#define PRAMANA_API __declspec(dllexport)
#else
#define PRAMANA_API __declspec(dllimport)
#endif
#else
#define PRAMANA_API __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct pramana_model pramana_model;
typedef struct pramana_result pramana_result;

/* Status codes (pramana_result_status) */
enum {
  PRAMANA_NOT_SOLVED = 0,
  PRAMANA_OPTIMAL = 1,
  PRAMANA_INFEASIBLE = 2,
  PRAMANA_UNBOUNDED = 3,
  PRAMANA_INFEASIBLE_OR_UNBOUNDED = 4,
  PRAMANA_TIME_LIMIT = 5,
  PRAMANA_ITERATION_LIMIT = 6,
  PRAMANA_NODE_LIMIT = 7,
  PRAMANA_NUMERICAL_FAILURE = 8,
  PRAMANA_NONCONVEX = 9,
  PRAMANA_ERROR = 10
};

PRAMANA_API const char* pramana_version(void);
PRAMANA_API const char* pramana_last_error(void);

/* ---- models ---- */
PRAMANA_API pramana_model* pramana_model_create(const char* name, int maximize);
PRAMANA_API pramana_model* pramana_model_read(const char* path); /* MPS/QPS, .gz ok; NULL on error */
PRAMANA_API void pramana_model_free(pramana_model* m);
PRAMANA_API int pramana_add_col(pramana_model* m, double cost, double lower, double upper, int integer,
                                const char* name); /* returns column index */
PRAMANA_API int pramana_add_row(pramana_model* m, double lower, double upper, int nnz, const int* cols,
                                const double* vals, const char* name); /* returns row index */
PRAMANA_API int pramana_add_q(pramana_model* m, int i, int j, double v); /* adds v to Q_ij and Q_ji (once if i==j) */
PRAMANA_API void pramana_set_offset(pramana_model* m, double offset);
PRAMANA_API int pramana_num_cols(pramana_model* m);
PRAMANA_API int pramana_num_rows(pramana_model* m);
PRAMANA_API int pramana_write_mps(pramana_model* m, const char* path);

/* ---- solving ---- */
PRAMANA_API pramana_result* pramana_solve(pramana_model* m, const char* options_json);
PRAMANA_API int pramana_result_status(const pramana_result* r);
PRAMANA_API const char* pramana_result_status_name(const pramana_result* r);
PRAMANA_API double pramana_result_objective(const pramana_result* r);
PRAMANA_API double pramana_result_bound(const pramana_result* r);
PRAMANA_API int pramana_result_certified(const pramana_result* r); /* 1 if the certificate was accepted */
/* Copy vectors into caller buffers of the given length; return the number of entries available. */
PRAMANA_API int pramana_result_x(const pramana_result* r, double* out, int len);
PRAMANA_API int pramana_result_row_duals(const pramana_result* r, double* out, int len);
PRAMANA_API int pramana_result_reduced_costs(const pramana_result* r, double* out, int len);
PRAMANA_API int pramana_result_row_activity(const pramana_result* r, double* out, int len);
PRAMANA_API const char* pramana_result_json(pramana_result* r, int include_vectors);
PRAMANA_API void pramana_result_free(pramana_result* r);

/* ---- certified parametric analysis / case families (JSON in, JSON out; free with pramana_free_string) ---- */
PRAMANA_API char* pramana_parametric(pramana_model* m, const char* spec_json);
PRAMANA_API char* pramana_family(pramana_model* m, const char* spec_json);
PRAMANA_API char* pramana_gpu_info(void);
PRAMANA_API void pramana_free_string(char* s);

#ifdef __cplusplus
}
#endif

#endif /* PRAMANA_H */
