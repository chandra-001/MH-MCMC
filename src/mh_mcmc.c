#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define sqr(x) x *x
#define sqrt(x) gsl_sf_exp(gsl_sf_log(x) / 2.)

extern float Omega_m_init;
extern float h_init;

const float Om0_low = 0.20;
const float Om0_high = 0.40;
const float h_low = 0.60;
const float h_high = 0.80;

float Omega_m_new;
float h_new;
double ratio;
double u;

float *Omega_m_list;
float *h_list;

typedef struct {
  float Omega_m;
  float h;
} Params;

double Likelihood(float, float);
void walker(gsl_rng *, gsl_rng *, gsl_rng *, Params *);
void proposal_dist_flat(float *, float *, gsl_rng *, gsl_rng *, Params *);

void mh_mcmc(unsigned long seed1, unsigned long seed2, unsigned long seed3,
             long burnin, long numsteps, char *filename) {
  Params params = {Omega_m_init, h_init};

  gsl_rng *r_prop_dist_h = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng *r_prop_dist_Omega_m = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng *r_accpt_rjct = gsl_rng_alloc(gsl_rng_taus);

  gsl_rng_set(r_prop_dist_h, seed1);
  gsl_rng_set(r_prop_dist_Omega_m, seed2);
  gsl_rng_set(r_accpt_rjct, seed3);

  for (long i = 0; i < burnin; ++i) {
    walker(r_prop_dist_h, r_prop_dist_Omega_m, r_accpt_rjct, &params);
  }

  Omega_m_list = (float *)malloc(numsteps * sizeof(float));
  h_list = (float *)malloc(numsteps * sizeof(float));

  for (long i = 0; i < numsteps; ++i) {
    walker(r_prop_dist_h, r_prop_dist_Omega_m, r_accpt_rjct, &params);
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
  gsl_rng_free(r_accpt_rjct);
  gsl_rng_free(r_prop_dist_Omega_m);
  gsl_rng_free(r_prop_dist_h);
} /* End of mh_mcmc() */

void walker(gsl_rng *r_prop_dist_h, gsl_rng *r_prop_dist_Omega_m,
            gsl_rng *r_accpt_rjct, Params *params) {
  proposal_dist_flat(&Omega_m_new, &h_new, r_prop_dist_h, r_prop_dist_Omega_m,
                     params);

  bool condition = (Omega_m_new > Om0_low && Omega_m_new < Om0_high) &&
                   (h_new > h_low && h_new < h_high);

  if (condition == true) {
    ratio =
        Likelihood(Omega_m_new, h_new) / Likelihood(params->Omega_m, params->h);

    if (ratio >= 1.) {
      params->Omega_m = Omega_m_new;
      params->h = h_new;
    } else {
      u = gsl_rng_uniform(r_accpt_rjct);
      if (u <= ratio) {
        params->Omega_m = Omega_m_new;
        params->h = h_new;
      }
    }
  }
} /* End of Walker() */

void proposal_dist_flat(float *Om0, float *h, gsl_rng *r1, gsl_rng *r2,
                        Params *params) {
  double width1 = 5e-4;
  double width2 = 5e-4;

  *Om0 = gsl_ran_flat(r1, params->Omega_m - width1, params->Omega_m + width1);
  *h = gsl_ran_flat(r2, params->h - width2, params->h + width2);
} /* End of proposal_dist() */
