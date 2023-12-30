#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_vector.h>

extern gsl_vector *mu;
extern gsl_vector *z;
extern gsl_vector *mu_th;
extern gsl_vector *diff;
extern gsl_vector *solve;
extern gsl_matrix *cov;
extern gsl_matrix *cov_copy;
extern int N;

double mu_theory(double, float, float);
double mu_fitting(double, float, float);
double D_L(double, float, float);

double Likelihood(float h, float Omega_m) {
  gsl_matrix_memcpy(cov_copy, cov);
  gsl_vector_memcpy(solve, mu);
  gsl_vector_memcpy(diff, mu);

  for (int i = 0; i < N; ++i) {
    gsl_vector_set(mu_th, i, mu_theory(gsl_vector_get(z, i), h, Omega_m));
  }

  gsl_vector_sub(diff, mu_th);
  gsl_vector_memcpy(solve, diff);

  gsl_permutation *p = gsl_permutation_alloc(N);
  int signum;

  gsl_linalg_LU_decomp(cov_copy, p, &signum);
  gsl_linalg_LU_svx(cov_copy, p, solve);

  gsl_permutation_free(p);

  double dot;
  gsl_blas_ddot(diff, solve, &dot);

  return -0.5 * dot;
}  // End of Lkelihood ()
