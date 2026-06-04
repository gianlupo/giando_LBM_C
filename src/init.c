#include <math.h>
#include "param.h"
#include "utils.h"

static void zero(double rho[LNX_H][LNY_H][LNZ_H],
                 double u[LNX_H][LNY_H][LNZ_H],
                 double v[LNX_H][LNY_H][LNZ_H],
                 double w[LNX_H][LNY_H][LNZ_H])
{

  for (int k = 1; k <= LNZ; ++k) {
    for (int j = 1; j <= LNY; ++j) {
      for (int i = 1; i <= LNX; ++i) {

        rho[i][j][k] = 1.0;
        u[i][j][k]   = 0.0;
        v[i][j][k]   = 0.0;
        w[i][j][k]   = 0.0;
      }
    }
  }

}

static void uniform(double rho[LNX_H][LNY_H][LNZ_H],
                    double u[LNX_H][LNY_H][LNZ_H],
                    double v[LNX_H][LNY_H][LNZ_H],
                    double w[LNX_H][LNY_H][LNZ_H])
{

  for (int k = 1; k <= LNZ; ++k) {
    for (int j = 1; j <= LNY; ++j) {
      for (int i = 1; i <= LNX; ++i) {

        rho[i][j][k] = 1.0;
        u[i][j][k]   = 0.66665*U0;
        v[i][j][k]   = 0.0;
        w[i][j][k]   = 0.0;
      }
    }
  }

}

static void poiseuille(double rho[LNX_H][LNY_H][LNZ_H],
                       double u[LNX_H][LNY_H][LNZ_H],
                       double v[LNX_H][LNY_H][LNZ_H],
                       double w[LNX_H][LNY_H][LNZ_H])
{

  double nu = (TAU - 0.5) / 3.0;
  double mu = 1.0 * nu;
  double G = 8 * mu * U0 / LZ / LZ;

  for (int k = 1; k <= LNZ; ++k) {
    for (int j = 1; j <= LNY; ++j) {
      for (int i = 1; i <= LNX; ++i) {

        double x = (double)gi(i) / LX;

        double upoi, vpoi, wpoi;
        poi(i, j, k, &upoi, &vpoi, &wpoi);

        rho[i][j][k] = (bc[0] == 'I' && bc[1] == 'O') ? (1.0 - G * (x - 1.0) * (LX)) : 1.0;
        u[i][j][k]   = upoi;
        v[i][j][k]   = vpoi;
        w[i][j][k]   = wpoi;
      }
    }
  }

}

static void sinewave(double rho[LNX_H][LNY_H][LNZ_H],
                     double u[LNX_H][LNY_H][LNZ_H],
                     double v[LNX_H][LNY_H][LNZ_H],
                     double w[LNX_H][LNY_H][LNZ_H])
{

  for (int k = 1; k <= LNZ; ++k) {
    for (int j = 1; j <= LNY; ++j) {
      for (int i = 1; i <= LNX; ++i) {

        double usin, vsin, wsin;
        sinvel(i, j, k, &usin, &vsin, &wsin);

        rho[i][j][k] = 1.0;
        u[i][j][k] = usin;
        v[i][j][k] = vsin;
        w[i][j][k] = wsin;
      }
    }
  }

}

void initialize(int ic,
                double rho[LNX_H][LNY_H][LNZ_H],
                double u[LNX_H][LNY_H][LNZ_H],
                double v[LNX_H][LNY_H][LNZ_H],
                double w[LNX_H][LNY_H][LNZ_H])
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
