#ifndef MCMC_H
#define MCMC_H

/* Data, allocated in read_val.c */
extern double *mu;    /* observed distance moduli */
extern double *z;     /* redshifts */
extern double *mu_th; /* model distance moduli */
extern double *diff;  /* mu - mu_th, overwritten by the triangular solve */
extern double *cov;   /* N x N row-major covariance, then its Cholesky factor */
extern int N;

/* Flat prior box, defined in mh_mcmc.c */
extern const double Om0_low, Om0_high, h_low, h_high;

/* read_val.c */
void read_val(const char *fname_z_mu);
void read_cov(const char *fname_cov);
void dealloc(void);

/* func.c */
double D_L(double z, double h, double Omega_m);
double mu_theory(double z, double h, double Omega_m);
double mu_fitting(double z, double h, double Omega_m);

/* likelihood.c */
void cholesky_cov(void);
double Likelihood(double h, double Omega_m);

/* mh_mcmc.c */
void mh_mcmc(double Omega_m_init, double h_init, unsigned long seed1,
             unsigned long seed2, unsigned long seed3, long burnin,
             long numsteps, const char *filename);

#endif /* MCMC_H */
