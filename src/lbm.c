#include "param.h"

void equilibrium(double rho[NX][NY],
                 double u[NX][NY],
                 double v[NX][NY],
                 double feq[Q][NX][NY])
{
  for (int q = 0; q < Q; ++q) {
    for (int j = 0; j < NY; ++j) {
      for (int i = 0; i < NX; ++i) {

        double cu = c[q][0] * u[i][j] + c[q][1] * v[i][j];

        double uu = u[i][j] * u[i][j] + v[i][j] * v[i][j];

        feq[q][i][j] = w[q] * rho[i][j] * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * uu);
      }
    }
  }
}

void stream(double f[Q][NX][NY],
            double fstar[Q][NX][NY])
{
  for (int q = 0; q < Q; ++q) {

    for (int j = 1; j < NY - 1; ++j) {
      for (int i = 1; i < NX - 1; ++i) {

        int istar = i - c[q][0];
        int jstar = j - c[q][1];

        fstar[q][i][j] = f[q][istar][jstar];
      }
    }
  }
}

void moments(double fstar[Q][NX][NY],
             double fx[NX][NY],
             double fy[NX][NY],
             double rho[NX][NY],
             double u[NX][NY],
             double v[NX][NY])
{
  for (int j = 0; j < NY; ++j) {
    for (int i = 0; i < NX; ++i) {

      double den = 0.0;
      double momx = 0.0;
      double momy = 0.0;

      for (int q = 0; q < Q; ++q) {

        den  += fstar[q][i][j];

        momx += c[q][0] * fstar[q][i][j];

        momy += c[q][1] * fstar[q][i][j];
      }

      rho[i][j] = den;
      u[i][j] = (momx + 0.5 * fx[i][j]) / den;
      v[i][j] = (momy + 0.5 * fy[i][j]) / den;
    }
  }
}

void guo(double fx[NX][NY],
         double fy[NX][NY],
         double u[NX][NY],
         double v[NX][NY],
         double fguo[Q][NX][NY])
{
  double taui = 1.0 / TAU;

  for (int q = 0; q < Q; ++q) {
    for (int j = 0; j < NY; ++j) {
      for (int i = 0; i < NX; ++i) {

        double cu = c[q][0] * u[i][j]  + c[q][1] * v[i][j];

        double cf = c[q][0] * fx[i][j] + c[q][1] * fy[i][j];

        double uf = u[i][j] * fx[i][j] + v[i][j] * fy[i][j];

        fguo[q][i][j] = w[q] * (1.0 - 0.5 * taui) * (3.0 * cf - 3.0 * uf + 9.0 * cu * cf);
      }
    }
  }
}

void collide(double fstar[Q][NX][NY],
             double feq[Q][NX][NY],
             double fguo[Q][NX][NY],
             double f[Q][NX][NY])
{
  double taui = 1.0 / TAU;

  for (int q = 0; q < Q; ++q) {

    for (int j = 0; j < NY; ++j) {
      for (int i = 0; i < NX; ++i) {

        f[q][i][j] = fstar[q][i][j] - taui * (fstar[q][i][j] - feq[q][i][j]) + fguo[q][i][j];
      }
    }
  }
}
