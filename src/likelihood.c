#include <cblas.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mcmc.h"

/* Overwrite the lower triangle of cov with L, where cov = L L^T.
   Done once, since the covariance matrix never changes. */
void cholesky_cov(void) {
  for (int j = 0; j < N; ++j) {
    double d = cov[j * N + j] - cblas_ddot(j, &cov[j * N], 1, &cov[j * N], 1);
    if (d <= 0.) {
      fprintf(stderr, "Error! Covariance matrix is not positive definite\n");
      exit(EXIT_FAILURE);
    }
    cov[j * N + j] = sqrt(d);

    for (int i = j + 1; i < N; ++i) {
      cov[i * N + j] =
          (cov[i * N + j] - cblas_ddot(j, &cov[i * N], 1, &cov[j * N], 1)) /
          cov[j * N + j];
    }
  }
}  // End of cholesky_cov ()

/* Returns the log-likelihood, -chi^2 / 2 */
double Likelihood(double h, double Omega_m) {
  for (int i = 0; i < N; ++i) {
    mu_th[i] = mu_theory(z[i], h, Omega_m);
    diff[i] = mu[i] - mu_th[i];
  }

  /* diff^T C^-1 diff = |L^-1 diff|^2 */
  cblas_dtrsv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit, N, cov,
              N, diff, 1);

  return -0.5 * cblas_ddot(N, diff, 1, diff, 1);
}  // End of Likelihood ()
