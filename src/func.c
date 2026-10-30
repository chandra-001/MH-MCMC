#include <math.h>

#include "mcmc.h"

static double eta(double, double);

static const double c = 2.99792458e10; /* speed of light in cm/s */
#define D_H(h) (c / (100. * (h) * 1e5))  /* Hubble distance in Mpc */

static inline double E(double z, double Omega_m) {
  return sqrt(Omega_m * (1. + z) * (1. + z) * (1. + z) + (1. - Omega_m));
}

/* Composite Simpson's rule for the integral of 1/E(x) from 0 to z */
double D_L(double z, double h, double Omega_m) {
  const int n = 32; /* number of subintervals, must be even */
  double dx = z / n;
  double sum = 1. / E(0., Omega_m) + 1. / E(z, Omega_m);
  for (int i = 1; i < n; ++i) {
    sum += (i % 2 ? 4. : 2.) / E(i * dx, Omega_m);
  }
  double result = sum * dx / 3.;
  return D_H(h) * (1. + z) * result;
}

double mu_theory(double z, double h, double Omega_m) {
  return (25. - 5. * log10(h) + 5. * log10(D_L(z, 1., Omega_m)));
}

/* Same as mu_theory(), using the Pen (1999) fitting formula for D_L */
double mu_fitting(double z, double h, double Omega_m) {
  double D_L1 =
      D_H(1.) * (1. + z) * (eta(1., Omega_m) - eta(1. / (1. + z), Omega_m));
  return 25. - 5 * log10(h) + 5 * log10(D_L1);
}

static double eta(double a, double Omega_m) {
  double s = cbrt((1. - Omega_m) / Omega_m); /* s^3 = (1 - Omega_m) / Omega_m */
  double term1 = 2 * sqrt(pow(s, 3) + 1.);
  double term2 = (1. / pow(a, 4) - 0.1540 * s / pow(a, 3) +
                  0.4304 * pow(s, 2) / pow(a, 2) + 0.19097 * pow(s, 3) / a +
                  0.066941 * pow(s, 4));

  return term1 * pow(term2, -1. / 8.);
}
