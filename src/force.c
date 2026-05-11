#include <math.h>
#include "param.h"

void compute_force(double u[NX][NY],
                   double v[NX][NY],
                   double fx[NX][NY],
                   double fy[NX][NY])
{
  double lambda = 1.0 * NY;
  double kappa = 2.0 * PI / lambda;
  double nu = (TAU - 0.5) / 3.0;

  for (int j = 0; j < NY; ++j) {
    for (int i = 0; i < NX; ++i) {

      fx[i][j] = 0.0;//nu * kappa * kappa * u[i][j];
      fy[i][j] = 0.0;
    }
  }
}
