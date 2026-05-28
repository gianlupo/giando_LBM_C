#include "param.h"

static int mod(int a, int b)
{
  int r = a % b;
  return (r < 0) ? r + b : r;
}

void equilibrium(double rho[NX][NY][NZ],
                 double u[NX][NY][NZ],
                 double v[NX][NY][NZ],
                 double w[NX][NY][NZ],
                 double feq[Q][NX][NY][NZ])
{
  for (int q = 0; q < Q; ++q) {
    for (int k = 0; k < NZ; ++k) {
      for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {

          double cu = c[q][0]    * u[i][j][k] + c[q][1]    * v[i][j][k] + c[q][2]    * w[i][j][k];

          double uu = u[i][j][k] * u[i][j][k] + v[i][j][k] * v[i][j][k] + w[i][j][k] * w[i][j][k];

          feq[q][i][j][k] = weights[q] * rho[i][j][k] * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu);
        }
      }
    }
  }
}

void stream(double f[Q][NX][NY][NZ],
            double fstar[Q][NX][NY][NZ])
{
  for (int q = 0; q < Q; ++q) {

    for (int k = 0; k < NZ; ++k) {
      for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {

          // includes periodic wrap-around in all directions
          int istar, jstar, kstar;
          istar = mod(i - c[q][0], NX);
          jstar = mod(j - c[q][1], NY);
          kstar = mod(k - c[q][2], NZ);

          fstar[q][i][j][k] = f[q][istar][jstar][kstar];
        }
      }
    }
  }
}

void moments(double fstar[Q][NX][NY][NZ],
             int obst[NX][NY][NZ],
             double fbx[NX][NY][NZ],
             double fby[NX][NY][NZ],
             double fbz[NX][NY][NZ],
             double rho[NX][NY][NZ],
             double u[NX][NY][NZ],
             double v[NX][NY][NZ],
             double w[NX][NY][NZ])
{
  for (int k = 0; k < NZ; ++k) {
    for (int j = 0; j < NY; ++j) {
      for (int i = 0; i < NX; ++i) {

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

void guo(double fbx[NX][NY][NZ],
         double fby[NX][NY][NZ],
         double fbz[NX][NY][NZ],
         double u[NX][NY][NZ],
         double v[NX][NY][NZ],
         double w[NX][NY][NZ],
         double fguo[Q][NX][NY][NZ])
{
  double taui = 1.0 / TAU;

  // Guo et al. Physical Review E, Vol. 65, 046308 (2002)

  for (int q = 0; q < Q; ++q) {
    for (int k = 0; k < NZ; ++k) {
      for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {

          double cu = c[q][0]    * u[i][j][k]  + c[q][1]    * v[i][j][k]  + c[q][2]    * w[i][j][k];

          double cf = c[q][0]    * fbx[i][j][k] + c[q][1]    * fby[i][j][k] + c[q][2]    * fbz[i][j][k];

          double uf = u[i][j][k] * fbx[i][j][k] + v[i][j][k] * fby[i][j][k] + w[i][j][k] * fbz[i][j][k];

          fguo[q][i][j][k] = weights[q] * (1.0 - 0.5 * taui) * (3.0 * cf - 3.0 * uf + 9.0 * cu * cf);
        }
      }
    }
  }
}

void collide(int obst[NX][NY][NZ],
             double fstar[Q][NX][NY][NZ],
             double feq[Q][NX][NY][NZ],
             double fguo[Q][NX][NY][NZ])
{
  double taui = 1.0 / TAU;

  for (int q = 0; q < Q; ++q) {

    for (int k = 0; k < NZ; ++k) {
      for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {

          if (obst[i][j][k] == 0) {
            fstar[q][i][j][k] = fstar[q][i][j][k] - taui * (fstar[q][i][j][k] - feq[q][i][j][k]) + fguo[q][i][j][k];
          }
        }
      }
    }
  }
}
