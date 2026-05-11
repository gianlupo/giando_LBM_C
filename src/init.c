#include <math.h>
#include "param.h"

void initialize(double rho[NX][NY],
                double u[NX][NY],
                double v[NX][NY])
{
  double lambda = 1.0 * NY;
  double kappa = 2.0 * PI / lambda;

  for (int j = 0; j < NY; ++j) {
    for (int i = 0; i < NX; ++i) {

      double y = j + 0.5;

      rho[i][j] = 1.0;
      u[i][j] = U0 * sin(kappa * y);
      v[i][j] = 0.0;
    }
  }
}
