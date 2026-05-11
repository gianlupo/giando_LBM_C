#include "param.h"

static int mod(int a, int b)
{
  int r = a % b;
  return (r < 0) ? r + b : r;
}

void periodic_bc(double f[Q][NX][NY][NZ],
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

    /*
    for (int i = 1; i < NX - 1; ++i) {

      int istar, kstar;

      // bottom
      istar = mod(i - c[q][0], NX);
      kstar = mod(0 - c[q][2], NZ);

      fstar[q][i][0] = f[q][istar][kstar];

      // top
      kstar = mod((NZ - 1) - c[q][2], NZ);

      fstar[q][i][NY - 1] = f[q][istar][kstar];
    }
    */
  }
}
