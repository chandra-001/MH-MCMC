#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <stddef.h>

double eta(double, float);

const float c = GSL_CONST_CGSM_SPEED_OF_LIGHT;
#define D_H(h) (c / (100. * h * 1e5))
#define Log10(x) (gsl_sf_log(x) / gsl_sf_log(10.))
#define nth_root(x, n) gsl_sf_exp(gsl_sf_log(x) / n)

static inline double E(double z, double Omega_m) {
  return sqrt(Omega_m * gsl_pow_3(1. + z) + (1. - Omega_m));
}

typedef struct {
  float Omega_m;
} params;

double integrand(double x, void *p) {
  params *params_loc = (params *)p;
  float Omega_m = (params_loc->Omega_m);
  return 1. / E(x, Omega_m);
}

double D_L(double z, float h, float Omega_m) {
  params params1 = {Omega_m};
  gsl_function F;
  F.function = &integrand;
  F.params = &params1;
  double epsabs = 1e-5;
  double epsrel = 1e-2;
  double result;
  double abserr;
  size_t neval;
  int val =
      gsl_integration_qng(&F, 0., z, epsabs, epsrel, &result, &abserr, &neval);
  return D_H(h) * (1. + z) * result;
}

double mu_theory(double z, float h, float Omega_m) {
  return (25. - 5. * Log10(h) + 5. * Log10(D_L(z, 1., Omega_m)));
}

double mu_fitting(double z, float h, float Omega_m) {
  double D_L1 =
      D_H(1.) * (1. + z) * (eta(1., Omega_m) - eta(1. / (1. + z), Omega_m));
  return 25. - 5 * Log10(h) + 5 * Log10(D_L1);
}

double eta(double a, float Omega_m) {
  double s = (1. - Omega_m) / Omega_m;
  double term1 = 2 * nth_root(gsl_pow_3(s) + 1., 2);
  double term2 = (1. / gsl_pow_4(a) - 0.1540 * s / gsl_pow_3(a) +
                  0.4304 * gsl_pow_2(s) / gsl_pow_2(a) +
                  0.19097 * gsl_pow_3(s) / a + 0.066941 * gsl_pow_4(s));

  return term1 * nth_root(term2, -8.);
}
