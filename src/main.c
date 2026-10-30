#include <stdio.h>
#include <stdlib.h>

#include "mcmc.h"

double *mu;
double *z;
double *mu_th;
double *diff;
double *cov;
int N;

int main(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s Omega_m_init h_init\n", argv[0]);
    return EXIT_FAILURE;
  }

  double Omega_m_init = atof(argv[1]);
  double h_init = atof(argv[2]);

  if (!(Omega_m_init > Om0_low && Omega_m_init < Om0_high) ||
      !(h_init > h_low && h_init < h_high)) {
    fprintf(stderr,
            "Error! Initial guess must satisfy %.2f < Omega_m < %.2f and "
            "%.2f < h < %.2f\n",
            Om0_low, Om0_high, h_low, h_high);
    return EXIT_FAILURE;
  }

  read_cov("jla_mub_covmatrix.dat");
  read_val("jla_mub.txt");
  cholesky_cov();
  unsigned long seed1 = 3784556191;
  unsigned long seed2 = 2359912397;
  unsigned long seed3 = 7291679348;
  long burnin = 10000;
  long numsteps = 100000;
  mh_mcmc(Omega_m_init, h_init, seed1, seed2, seed3, burnin, numsteps,
          "walker_chain.txt");
  dealloc();
}
