#include "param.h"
#include "utils.h"

void equilibrium(double rho[LNX_H][LNY_H][LNZ_H],
                 double u[LNX_H][LNY_H][LNZ_H],
                 double v[LNX_H][LNY_H][LNZ_H],
                 double w[LNX_H][LNY_H][LNZ_H],
                 double feq[Q][LNX_H][LNY_H][LNZ_H])
{
  for (int q = 0; q < Q; ++q) {
    for (int k = 1; k <= LNZ; ++k) {
      for (int j = 1; j <= LNY; ++j) {
        for (int i = 1; i <= LNX; ++i) {

          double cu = c[q][0]    * u[i][j][k] + c[q][1]    * v[i][j][k] + c[q][2]    * w[i][j][k];

          double uu = u[i][j][k] * u[i][j][k] + v[i][j][k] * v[i][j][k] + w[i][j][k] * w[i][j][k];

          feq[q][i][j][k] = weights[q] * rho[i][j][k] * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu);
        }
      }
    }
  }
}

void stream(double f[Q][LNX_H][LNY_H][LNZ_H],
            double fstar[Q][LNX_H][LNY_H][LNZ_H])
{
  for (int q = 0; q < Q; ++q) {

    for (int k = 1; k <= LNZ; ++k) {
      for (int j = 1; j <= LNY; ++j) {
        for (int i = 1; i <= LNX; ++i) {

          int istar, jstar, kstar;
          istar = i - c[q][0];
          jstar = j - c[q][1];
          kstar = k - c[q][2];

          fstar[q][i][j][k] = f[q][istar][jstar][kstar];
        }
      }
    }
  }
}

void moments(double fstar[Q][LNX_H][LNY_H][LNZ_H],
             int    obst[LNX_H][LNY_H][LNZ_H],
             double fbx[LNX_H][LNY_H][LNZ_H],
             double fby[LNX_H][LNY_H][LNZ_H],
             double fbz[LNX_H][LNY_H][LNZ_H],
             double rho[LNX_H][LNY_H][LNZ_H],
             double u[LNX_H][LNY_H][LNZ_H],
             double v[LNX_H][LNY_H][LNZ_H],
             double w[LNX_H][LNY_H][LNZ_H])
{
  for (int k = 1; k <= LNZ; ++k) {
    for (int j = 1; j <= LNY; ++j) {
      for (int i = 1; i <= LNX; ++i) {

        double den = 0.0;
        double momx = 0.0;
        double momy = 0.0;
        double momz = 0.0;

        for (int q = 0; q < Q; ++q) {

          den  += fstar[q][i][j][k];
          momx += c[q][0] * fstar[q][i][j][k];
          momy += c[q][1] * fstar[q][i][j][k];
          momz += c[q][2] * fstar[q][i][j][k];
        }

        rho[i][j][k] = obst[i][j][k] + (1.0 - obst[i][j][k]) * den;
        u[i][j][k] = (1.0 - obst[i][j][k]) * (momx + 0.5 * fbx[i][j][k]) / den;
        v[i][j][k] = (1.0 - obst[i][j][k]) * (momy + 0.5 * fby[i][j][k]) / den;
        w[i][j][k] = (1.0 - obst[i][j][k]) * (momz + 0.5 * fbz[i][j][k]) / den;

      }
    }
  }
}

void guo(double fbx[LNX_H][LNY_H][LNZ_H],
         double fby[LNX_H][LNY_H][LNZ_H],
         double fbz[LNX_H][LNY_H][LNZ_H],
         double u[LNX_H][LNY_H][LNZ_H],
         double v[LNX_H][LNY_H][LNZ_H],
         double w[LNX_H][LNY_H][LNZ_H],
         double fguo[Q][LNX_H][LNY_H][LNZ_H])
{
  double taui = 1.0 / TAU;

  // Guo et al. Physical Review E, Vol. 65, 046308 (2002)

  for (int q = 0; q < Q; ++q) {

    for (int k = 1; k <= LNZ; ++k) {
      for (int j = 1; j <= LNY; ++j) {
        for (int i = 1; i <= LNX; ++i) {

          double cu = c[q][0]    * u[i][j][k]  + c[q][1]    * v[i][j][k]  + c[q][2]    * w[i][j][k];

          double cf = c[q][0]    * fbx[i][j][k] + c[q][1]    * fby[i][j][k] + c[q][2]    * fbz[i][j][k];

          double uf = u[i][j][k] * fbx[i][j][k] + v[i][j][k] * fby[i][j][k] + w[i][j][k] * fbz[i][j][k];

          fguo[q][i][j][k] = weights[q] * (1.0 - 0.5 * taui) * (3.0 * cf - 3.0 * uf + 9.0 * cu * cf);
        }
      }
    }
  }
}

void collide(int    obst[LNX_H][LNY_H][LNZ_H],
             double fstar[Q][LNX_H][LNY_H][LNZ_H],
             double feq[Q][LNX_H][LNY_H][LNZ_H],
             double fguo[Q][LNX_H][LNY_H][LNZ_H])
{
  double taui = 1.0 / TAU;

  for (int q = 0; q < Q; ++q) {

    for (int k = 1; k <= LNZ; ++k) {
      for (int j = 1; j <= LNY; ++j) {
        for (int i = 1; i <= LNX; ++i) {

          if (obst[i][j][k] == 0) {
            fstar[q][i][j][k] = fstar[q][i][j][k] - taui * (fstar[q][i][j][k] - feq[q][i][j][k]) + fguo[q][i][j][k];
          }
        }
      }
    }
  }
}
