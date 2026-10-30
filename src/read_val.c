#include <stdio.h>
#include <stdlib.h>

#include "mcmc.h"

static void alloc(void) {
  mu = (double *)malloc(N * sizeof(double));
  z = (double *)malloc(N * sizeof(double));
  mu_th = (double *)malloc(N * sizeof(double));
  diff = (double *)malloc(N * sizeof(double));
  cov = (double *)malloc(N * N * sizeof(double));

  if (!mu || !z || !mu_th || !diff || !cov) {
    fprintf(stderr, "Error! Cannot allocate memory for N = %d\n", N);
    exit(EXIT_FAILURE);
  }
}

void read_val(const char *fname_z_mu) {
  FILE *fp;

  if ((fp = fopen(fname_z_mu, "r")) == NULL) {
    fprintf(stderr, "Error! Cannot open file %s!\n", fname_z_mu);
    exit(EXIT_FAILURE);
  }

  /* skip the header line */
  int c;
  do {
    c = getc(fp);
  } while (c != '\n' && c != EOF);

  for (int i = 0; i < N; ++i) {
    if (fscanf(fp, "%lf%lf", &z[i], &mu[i]) != 2) {
      fprintf(stderr, "Error! Cannot read row %d of %s (expected %d rows)\n",
              i + 1, fname_z_mu, N);
      exit(EXIT_FAILURE);
    }
  }

  fclose(fp);
}

void read_cov(const char *fname_cov) {
  FILE *fp;

  if ((fp = fopen(fname_cov, "r")) == NULL) {
    fprintf(stderr, "Error! Cannot open file %s!\n", fname_cov);
    exit(EXIT_FAILURE);
  }

  if (fscanf(fp, "%d", &N) != 1 || N <= 0) {
    fprintf(stderr, "Error! Cannot read a valid size N from %s\n", fname_cov);
    exit(EXIT_FAILURE);
  }

  alloc();

  /* row-major: cov[i * N + j] */
  for (int i = 0; i < N * N; ++i) {
    if (fscanf(fp, "%lf", &cov[i]) != 1) {
      fprintf(stderr, "Error! Cannot read element %d of %s (expected %d)\n",
              i + 1, fname_cov, N * N);
      exit(EXIT_FAILURE);
    }
  }

  fclose(fp);
}

void dealloc(void) {
  free(mu);
  free(z);
  free(mu_th);
  free(diff);
  free(cov);
}
