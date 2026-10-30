#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "mcmc.h"

const double Om0_low = 0.20;
const double Om0_high = 0.40;
const double h_low = 0.60;
const double h_high = 0.80;

typedef struct {
  double Omega_m;
  double h;
  double logL; /* log-likelihood at (Omega_m, h), cached */
} Params;

/* Maximally equidistributed combined Tausworthe generator (L'Ecuyer 1996).
   Same algorithm and seeding as GSL's gsl_rng_taus, so a given seed gives
   the same sequence. */
typedef struct {
  unsigned long s1, s2, s3;
} rng_taus;

#define MASK 0xffffffffUL
#define TAUSWORTHE(s, a, b, c, d) \
  ((((s) & (c)) << (d)) & MASK) ^ (((((s) << (a)) & MASK) ^ (s)) >> (b))
#define LCG(n) ((69069 * (n)) & MASK)

static unsigned long rng_get(rng_taus *r) {
  r->s1 = TAUSWORTHE(r->s1, 13, 19, 4294967294UL, 12);
  r->s2 = TAUSWORTHE(r->s2, 2, 25, 4294967288UL, 4);
  r->s3 = TAUSWORTHE(r->s3, 3, 11, 4294967280UL, 17);
  return r->s1 ^ r->s2 ^ r->s3;
}

static void rng_set(rng_taus *r, unsigned long seed) {
  if (seed == 0) seed = 1;
  r->s1 = LCG(seed);
  if (r->s1 < 2) r->s1 += 2UL;
  r->s2 = LCG(r->s1);
  if (r->s2 < 8) r->s2 += 8UL;
  r->s3 = LCG(r->s2);
  if (r->s3 < 16) r->s3 += 16UL;
  for (int i = 0; i < 6; ++i) rng_get(r); /* warm up */
}

/* uniform on [0, 1) */
static double rng_uniform(rng_taus *r) { return rng_get(r) / 4294967296.0; }

/* uniform on [a, b) */
static double ran_flat(rng_taus *r, double a, double b) {
  double x = rng_uniform(r);
  return a * (1 - x) + b * x;
}

static void walker(rng_taus *, rng_taus *, rng_taus *, Params *);
static void proposal_dist_flat(double *, double *, rng_taus *, rng_taus *,
                               const Params *);

void mh_mcmc(double Omega_m_init, double h_init, unsigned long seed1,
             unsigned long seed2, unsigned long seed3, long burnin,
             long numsteps, const char *filename) {
  Params params = {Omega_m_init, h_init, Likelihood(h_init, Omega_m_init)};

  rng_taus r_prop_dist_h;
  rng_taus r_prop_dist_Omega_m;
  rng_taus r_accpt_rjct;

  rng_set(&r_prop_dist_h, seed1);
  rng_set(&r_prop_dist_Omega_m, seed2);
  rng_set(&r_accpt_rjct, seed3);

  for (long i = 0; i < burnin; ++i) {
    walker(&r_prop_dist_h, &r_prop_dist_Omega_m, &r_accpt_rjct, &params);
  }

  double *Omega_m_list = (double *)malloc(numsteps * sizeof(double));
  double *h_list = (double *)malloc(numsteps * sizeof(double));

  if (!Omega_m_list || !h_list) {
    fprintf(stderr, "Error! Cannot allocate memory for %ld steps\n", numsteps);
    exit(EXIT_FAILURE);
  }

  for (long i = 0; i < numsteps; ++i) {
    walker(&r_prop_dist_h, &r_prop_dist_Omega_m, &r_accpt_rjct, &params);
    Omega_m_list[i] = params.Omega_m;
    h_list[i] = params.h;
  }

  FILE *fp = fopen(filename, "w");

  if (fp == NULL) {
    perror("Error! Cannot open file");
    exit(EXIT_FAILURE);
  }

  for (long i = 0; i < numsteps; ++i) {
    fprintf(fp, "%f %f\n", Omega_m_list[i], h_list[i]);
  }

  fclose(fp);

  free(h_list);
  free(Omega_m_list);
} /* End of mh_mcmc() */

static void walker(rng_taus *r_prop_dist_h, rng_taus *r_prop_dist_Omega_m,
                   rng_taus *r_accpt_rjct, Params *params) {
  double Omega_m_new;
  double h_new;

  proposal_dist_flat(&Omega_m_new, &h_new, r_prop_dist_h, r_prop_dist_Omega_m,
                     params);

  bool condition = (Omega_m_new > Om0_low && Omega_m_new < Om0_high) &&
                   (h_new > h_low && h_new < h_high);

  if (condition == true) {
    /* Likelihood() takes (h, Omega_m) and returns the log-likelihood */
    double logL_new = Likelihood(h_new, Omega_m_new);
    double ratio = exp(logL_new - params->logL);

    if (ratio >= 1. || rng_uniform(r_accpt_rjct) <= ratio) {
      params->Omega_m = Omega_m_new;
      params->h = h_new;
      params->logL = logL_new;
    }
  }
} /* End of Walker() */

static void proposal_dist_flat(double *Om0, double *h, rng_taus *r1,
                               rng_taus *r2, const Params *params) {
  /* half-widths tuned for ~35% acceptance */
  double width1 = 0.07;
  double width2 = 0.016;

  *Om0 = ran_flat(r1, params->Omega_m - width1, params->Omega_m + width1);
  *h = ran_flat(r2, params->h - width2, params->h + width2);
} /* End of proposal_dist() */
