#include <math.h>
#include "param.h"

void initialize(double rho[NX][NY][NZ],
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
