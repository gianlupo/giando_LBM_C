#include <math.h>
#include "param.h"

void zero(double rho[NX][NY][NZ],
          double u[NX][NY][NZ],
          double v[NX][NY][NZ],
          double w[NX][NY][NZ])
{

  for (int k = 0; k < NZ; ++k) {
    for (int j = 0; j < NY; ++j) {
      for (int i = 0; i < NX; ++i) {

        rho[i][j][k] = 1.0;
        u[i][j][k]   = 0.0;
        v[i][j][k]   = 0.0;
        w[i][j][k]   = 0.0;
      }
    }
  }

}

void uniform(double rho[NX][NY][NZ],
             double u[NX][NY][NZ],
             double v[NX][NY][NZ],
             double w[NX][NY][NZ])
{

  for (int k = 0; k < NZ; ++k) {
    for (int j = 0; j < NY; ++j) {
      for (int i = 0; i < NX; ++i) {

        rho[i][j][k] = 1.0;
        u[i][j][k]   = 0.66665*U0;
        v[i][j][k]   = 0.0;
        w[i][j][k]   = 0.0;
      }
    }
  }

}

void poiseuille(double rho[NX][NY][NZ],
                double u[NX][NY][NZ],
                double v[NX][NY][NZ],
                double w[NX][NY][NZ])
{

  for (int k = 0; k < NZ; ++k) {
    for (int j = 0; j < NY; ++j) {
      for (int i = 0; i < NX; ++i) {

        double x = (double)i;
        double z = (double)k / (NZ-1);
        double nu = (TAU - 0.5) / 3.0;
        double mu = 1.0 * nu;
        double G = 8 * mu * U0 / (NZ - 1.0) / (NZ - 1.0);

        rho[i][j][k] = 1.0 - G * (x - (NX-1));
        u[i][j][k] = 4.0 * U0 * z * (1.0 - z);
        v[i][j][k] = 0.0;
        w[i][j][k] = 0.0;
      }
    }
  }

}

void sinewave(double rho[NX][NY][NZ],
              double u[NX][NY][NZ],
              double v[NX][NY][NZ],
              double w[NX][NY][NZ])
{

  double lambda = 1.0 * NZ;
  double kappa = 2.0 * PI / lambda;

  for (int k = 0; k < NZ; ++k) {
    for (int j = 0; j < NY; ++j) {
      for (int i = 0; i < NX; ++i) {

        double z = k + 0.5;

        rho[i][j][k] = 1.0;
        u[i][j][k] = U0 * sin(kappa * z);
        v[i][j][k] = 0.0;
        w[i][j][k] = 0.0;
      }
    }
  }

}

void initialize(int ic,
                double rho[NX][NY][NZ],
                double u[NX][NY][NZ],
                double v[NX][NY][NZ],
                double w[NX][NY][NZ])
{
  switch (ic) {
    case 0:
      zero(rho, u, v, w);
      break;
    case 1:
      uniform(rho, u, v, w);
      break;
    case 2:
      poiseuille(rho, u, v, w);
      break;
    case 3:
      sinewave(rho, u, v, w);
      break;
  }
}
