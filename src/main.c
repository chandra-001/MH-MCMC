#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <stdio.h>

gsl_vector *mu;
gsl_vector *z;
gsl_vector *mu_th;
gsl_vector *diff;
gsl_vector *solve;
gsl_matrix *cov;
gsl_matrix *cov_copy;
int N;

void read_val(char *);
void read_cov(char *);
double Likelihood(float, float);
void mh_mcmc(unsigned long, unsigned long, unsigned long, long, long,
             const char *);
void alloc();
void dealloc();

float Omega_m_init;
float h_init;

int main(int argc, char *argv[]) {
  read_cov("jla_mub_covmatrix.dat");
  read_val("jla_mub.txt");
  Omega_m_init = atof(argv[1]);
  h_init = atof(argv[2]);
  unsigned long seed1 = 3784556191;
  unsigned long seed2 = 2359912397;
  unsigned long seed3 = 7291679348;
  long burnin = 10000;
  long numsteps = 100000;
  gsl_set_error_handler_off();
  mh_mcmc(seed1, seed2, seed3, burnin, numsteps, "walker_chain.txt");
  dealloc();
}
