#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <stdio.h>

extern gsl_vector *mu;
extern gsl_vector *z;
extern gsl_vector *mu_th;
extern gsl_vector *diff;
extern gsl_vector *solve;
extern gsl_matrix *cov;
extern gsl_matrix *cov_copy;
extern int N;

void alloc();

void alloc() {
  mu = gsl_vector_alloc(N);
  z = gsl_vector_alloc(N);
  mu_th = gsl_vector_alloc(N);
  diff = gsl_vector_alloc(N);
  solve = gsl_vector_alloc(N);
  cov = gsl_matrix_alloc(N, N);
  cov_copy = gsl_matrix_alloc(N, N);
}

void read_val(char *fname_z_mu_th) {
  mu = gsl_vector_alloc(N);
  z = gsl_vector_alloc(N);
  mu_th = gsl_vector_alloc(N);

  FILE *fp;

  if ((fp = fopen(fname_z_mu_th, "r")) == NULL) {
    fprintf(stderr, "Error!Cannot open file %s!\n", fname_z_mu_th);
    exit(1);
  }

  float x;
  float y;

  char c;
  do {
    c = getc(fp);
  } while (c != '\n');

  size_t Num;
  for (int i = 0; i < N; ++i) {
    Num = fscanf(fp, "%f%f", &x, &y);
    gsl_vector_set(z, i, x);
    gsl_vector_set(mu, i, y);
  }

  fclose(fp);
}

void read_cov(char *fname_cov) {
  FILE *fp = fopen(fname_cov, "r");
  size_t Num = fscanf(fp, "%d", &N);
  alloc();
  gsl_matrix_fscanf(fp, cov);
  fclose(fp);
}

void dealloc() {
  gsl_vector_free(mu);
  gsl_vector_free(z);
  gsl_vector_free(mu_th);
  gsl_vector_free(diff);
  gsl_vector_free(solve);
  gsl_matrix_free(cov);
  gsl_matrix_free(cov_copy);
}
