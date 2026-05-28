#include <math.h>
#include "param.h"

void compute_force(double rho[NX][NY][NZ],
                   double u[NX][NY][NZ],
                   double v[NX][NY][NZ],
                   double w[NX][NY][NZ],
                   double fbx[NX][NY][NZ],
                   double fby[NX][NY][NZ],
                   double fbz[NX][NY][NZ])
{

  int fbulk = 0;

  double nu = (TAU - 0.5) / 3.0;

  switch (fbulk) {
    case 0: {
      for (int k = 0; k < NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
          for (int i = 0; i < NX; ++i) {

            fbx[i][j][k] = 0.0;
            fby[i][j][k] = 0.0;
            fbz[i][j][k] = 0.0;
          }
        }
      }
      break;
    }
    case 1: {
      double lambda = 1.0 * (NZ - 1);
      double kappa = 2.0 * PI / lambda;
      for (int k = 0; k < NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
          for (int i = 0; i < NX; ++i) {

            fbx[i][j][k] = nu * kappa * kappa * u[i][j][k];
            fby[i][j][k] = 0.0;
            fbz[i][j][k] = 0.0;
          }
        }
      }
      break;
    }
    case 2: {
      for (int k = 0; k < NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
          for (int i = 0; i < NX; ++i) {

            fbx[i][j][k] = 8.0 * rho[i][j][k] * nu * U0 / (NZ - 1) / (NZ - 1);
            fby[i][j][k] = 0.0;
            fbz[i][j][k] = 0.0;
          }
        }
      }
      break;
    }
  }

}
