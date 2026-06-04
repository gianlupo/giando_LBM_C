#include <math.h>
#include "param.h"
#include "utils.h"

void compute_force(double rho[LNX_H][LNY_H][LNZ_H],
                   double u[LNX_H][LNY_H][LNZ_H],
                   double v[LNX_H][LNY_H][LNZ_H],
                   double w[LNX_H][LNY_H][LNZ_H],
                   double fbx[LNX_H][LNY_H][LNZ_H],
                   double fby[LNX_H][LNY_H][LNZ_H],
                   double fbz[LNX_H][LNY_H][LNZ_H])
{

  int fbulk = 0;

  double nu = (TAU - 0.5) / 3.0;

  switch (fbulk) {
    case 0: {
      for (int k = 1; k <= LNZ; ++k) {
        for (int j = 1; j <= LNY; ++j) {
          for (int i = 1; i <= LNX; ++i) {

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
      for (int k = 1; k <= LNZ; ++k) {
        for (int j = 1; j <= LNY; ++j) {
          for (int i = 1; i <= LNX; ++i) {

            fbx[i][j][k] = nu * kappa * kappa * u[i][j][k];
            fby[i][j][k] = 0.0;
            fbz[i][j][k] = 0.0;
          }
        }
      }
      break;
    }
    case 2: {
      for (int k = 1; k <= LNZ; ++k) {
        for (int j = 1; j <= LNY; ++j) {
          for (int i = 1; i <= LNX; ++i) {

            fbx[i][j][k] = 8.0 * rho[i][j][k] * nu * U0 / LZ / LZ;
            fby[i][j][k] = 0.0;
            fbz[i][j][k] = 0.0;
          }
        }
      }
      break;
    }
  }

}
