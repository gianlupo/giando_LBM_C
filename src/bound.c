#include "param.h"

static int mod(int a, int b)
{
  int r = a % b;
  return (r < 0) ? r + b : r;
}

void periodic_bc(double f[Q][NX][NY],
                 double fstar[Q][NX][NY])
{
  for (int q = 0; q < Q; ++q) {

    for (int j = 0; j < NY; ++j) {

      int istar, jstar;

      // left
      istar = mod(0 - c[q][0], NX);
      jstar = mod(j - c[q][1], NY);

      fstar[q][0][j] = f[q][istar][jstar];

      // right
      istar = mod((NX - 1) - c[q][0], NX);

      fstar[q][NX - 1][j] = f[q][istar][jstar];
    }

    for (int i = 1; i < NX - 1; ++i) {

      int istar, jstar;

      // bottom
      istar = mod(i - c[q][0], NX);
      jstar = mod(0 - c[q][1], NY);

      fstar[q][i][0] = f[q][istar][jstar];

      // top
      jstar = mod((NY - 1) - c[q][1], NY);

      fstar[q][i][NY - 1] = f[q][istar][jstar];
    }
  }
}
