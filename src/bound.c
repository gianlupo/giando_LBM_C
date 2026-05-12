#include "param.h"

static int mod(int a, int b)
{
  int r = a % b;
  return (r < 0) ? r + b : r;
}

void impose_bc(int plane,
               char bctype,
               double f[Q][NX][NY][NZ],
               double fstar[Q][NX][NY][NZ])
{

  int imin;
  int imax;
  int jmin;
  int jmax;
  int kmin;
  int kmax;

  imin = 0;
  imax = NX;
  jmin = 0;
  jmax = NY;
  kmin = 0;
  kmax = NZ;

  switch (plane) {

    case 1:
      imin = 0;
      imax = 1;
      break;

    case 2:
      imin = NX - 1;
      imax = NX;
      break;

    case 3:
      jmin = 0;
      jmax = 1;
      break;

    case 4:
      jmin = NY - 1;
      jmax = NY;
      break;

    case 5:
      kmin = 0;
      kmax = 1;
      break;

    case 6:
      kmin = NZ - 1;
      kmax = NZ;
      break;

  }

  switch (bctype) {
    case 'P':
      // periodic b.c.
      for (int q = 0; q < Q; ++q) {
        for (int k = kmin; k < kmax; ++k) {
          for (int j = jmin; j < jmax; ++j) {
            for (int i = imin; i < imax; ++i) {

              int istar, jstar, kstar;

              istar = mod(i - c[q][0], NX);
              jstar = mod(j - c[q][1], NY);
              kstar = mod(k - c[q][2], NZ);

              fstar[q][i][j][k] = f[q][istar][jstar][kstar];
            }
          }
        }
      }
      break;
    case 'W':
      // wall b.c. (bounceback scheme)
      for (int k = kmin; k < kmax; ++k) {
        for (int j = jmin; j < jmax; ++j) {
          for (int i = imin; i < imax; ++i) {

            for (int q = 0; q < 5; ++q) {

              fstar[qo[plane-1][q]][i][j][k] = f[oppq[qo[plane-1][q]]][i][j][k];

            }

            for (int q = 0; q < 14; ++q) {

              int istar, jstar, kstar;

              istar = mod(i - c[qi[plane-1][q]][0], NX);
              jstar = mod(j - c[qi[plane-1][q]][1], NY);
              kstar = mod(k - c[qi[plane-1][q]][2], NZ);

              fstar[qi[plane-1][q]][i][j][k] = f[qi[plane-1][q]][istar][jstar][kstar];

            }
          }
        }
      }
      break;
  }

}

void boundary(double f[Q][NX][NY][NZ],
              double fstar[Q][NX][NY][NZ])
{
  int plane;
  char bctype;

  for (int b = 0; b < 6; b++) {
    plane = b + 1;
    bctype = bc[b];
    impose_bc(plane, bctype, f, fstar);

  }
}

void triperiodic_bc(double f[Q][NX][NY][NZ],
                    double fstar[Q][NX][NY][NZ])
{
  for (int q = 0; q < Q; ++q) {

    for (int k = 0; k < NZ; ++k) {
      for (int j = 0; j < NY; ++j) {

        int istar, jstar, kstar;

        // left
        istar = mod(0 - c[q][0], NX);
        jstar = mod(j - c[q][1], NY);
        kstar = mod(k - c[q][2], NZ);

        fstar[q][0][j][k] = f[q][istar][jstar][kstar];

        // right
        istar = mod((NX - 1) - c[q][0], NX);

        fstar[q][NX - 1][j][k] = f[q][istar][jstar][kstar];
      }
    }

    for (int k = 0; k < NZ; ++k) {
      for (int i = 0; i < NX; ++i) {

        int istar, jstar, kstar;

        // back
        istar = mod(i - c[q][0], NX);
        jstar = mod(0 - c[q][1], NY);
        kstar = mod(k - c[q][2], NZ);

        fstar[q][i][0][k] = f[q][istar][jstar][kstar];

        // front
        jstar = mod((NY - 1) - c[q][1], NY);

        fstar[q][i][NY - 1][k] = f[q][istar][jstar][kstar];
      }
    }

    for (int j = 0; j < NY; ++j) {
      for (int i = 0; i < NX; ++i) {

        int istar, jstar, kstar;

        // bottom
        istar = mod(i - c[q][0], NX);
        jstar = mod(j - c[q][1], NY);
        kstar = mod(0 - c[q][2], NZ);

        fstar[q][i][j][0] = f[q][istar][jstar][kstar];

        // top
        kstar = mod((NZ - 1) - c[q][2], NZ);

        fstar[q][i][j][NZ - 1] = f[q][istar][jstar][kstar];
      }
    }
  }
}
