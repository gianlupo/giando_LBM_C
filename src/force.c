#include <math.h>
#include "param.h"

void compute_force(double u[NX][NY][NZ],
                   double v[NX][NY][NZ],
                   double w[NX][NY][NZ],
                   double fx[NX][NY][NZ],
                   double fy[NX][NY][NZ],
                   double fz[NX][NY][NZ])
{
  double lambda = 1.0 * NZ;
  double kappa = 2.0 * PI / lambda;
  double nu = (TAU - 0.5) / 3.0;

  for (int k = 0; k < NZ; ++k) {
    for (int j = 0; j < NY; ++j) {
      for (int i = 0; i < NX; ++i) {

        fx[i][j][k] = 0.0;//nu * kappa * kappa * u[i][j][k];
        fy[i][j][k] = 0.0;
        fz[i][j][k] = 0.0;
      }
    }
  }
}
