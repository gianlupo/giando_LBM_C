#include "param.h"
#include <stdlib.h>
#include <stdio.h> //debug

static int mod(int a, int b)
{
  int r = a % b;
  return (r < 0) ? r + b : r;
}

void inlet(int i,
           int j,
           int k,
           double *uin,
           double *vin,
           double *win)
{

  double z = (k + 0.5) / NZ;
  *uin = 4.0 * U0 * z * (1.0 - z);
  *vin = 0.0;
  *win = 0.0;

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

    case 1: {
      imin = 0;
      imax = 1;
      break;
      }

    case 2: {
      imin = NX - 1;
      imax = NX;
      break;
      }

    case 3: {
      jmin = 0;
      jmax = 1;
      break;
      }

    case 4: {
      jmin = NY - 1;
      jmax = NY;
      break;
      }

    case 5: {
      kmin = 0;
      kmax = 1;
      break;
      }

    case 6: {
      kmin = NZ - 1;
      kmax = NZ;
      break;
      }

  }

  switch (bctype) {
    case 'P': {
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
      }
    case 'W': {
      // wall b.c. (bounceback scheme)
      for (int k = kmin; k < kmax; ++k) {
        for (int j = jmin; j < jmax; ++j) {
          for (int i = imin; i < imax; ++i) {

            for (int q = 0; q < 5; ++q) {

              fstar[qo[plane-1][q]][i][j][k] = f[opp[qo[plane-1][q]]][i][j][k];

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
    case 'I': {
      // inlet b.c.
      for (int k = kmin; k < kmax; ++k) {
        for (int j = jmin; j < jmax; ++j) {
          for (int i = imin; i < imax; ++i) {

            double uin;
            double vin;
            double win;

            inlet(i, j, k, &uin, &vin, &win);

            double rhoin = 0.0;

            for (int p = 0; p < 5; ++p) {
              int q = qo[plane-1][p];
              rhoin = rhoin + f[opp[q]][i][j][k];
            }

            for (int p = 0; p < 14; ++p) {
              int q = qi[plane-1][p];
              rhoin = rhoin + f[q][i][j][k];
            }

            rhoin = rhoin / (1.0 - uin);

            for (int p = 0; p < 5; ++p) {

              int q = qo[plane-1][p];

              int cn = c[q][0] * n[plane-1][0] + c[q][1] * n[plane-1][1] + c[q][2] * n[plane-1][2];
              int tx = c[q][0] - cn * n[plane-1][0];
              int ty = c[q][1] - cn * n[plane-1][1];
              int tz = c[q][2] - cn * n[plane-1][2];

              double cu = c[q][0] * uin + c[q][1] * vin + c[q][2] * win;
              double tu = tx      * uin + ty      * vin + tz      * win;

              double sumf  = 0.0;

              for (int r = 0; r < 14; ++r) {
                int s = qi[plane-1][r];
                int cst = c[s][0] * tx            + c[s][1] * ty            + c[s][2] * tz           ;
                int csn = c[s][0] * n[plane-1][0] + c[s][1] * n[plane-1][1] + c[s][2] * n[plane-1][2];
                sumf = sumf + cst * f[s][i][j][k] * (1 - abs(csn));
              }

              fstar[q][i][j][k] = f[opp[q]][i][j][k] + 1.0*(1.0 / 3.0) * tu * rhoin
                                                     + 1.0*(1.0 / 6.0) * cu * rhoin
                                                     - 0.0*(1.0 / 2.0) * sumf;
            }
            fstar[plane][i][j][k] = fstar[plane][i][j][k] + (1.0 / 6.0) * uin * rhoin;

            for (int p = 0; p < 14; ++p) {

              int q = qi[plane-1][p];

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
      }
    case 'O': {
      // outlet b.c.
      for (int k = kmin; k < kmax; ++k) {
        for (int j = jmin; j < jmax; ++j) {
          for (int i = imin; i < imax; ++i) {

            double rhoout = 1.0;
            double outlet = 0.0;

            for (int p = 0; p < 5; ++p) {
              int q = qo[plane-1][p];
              outlet = outlet + f[opp[q]][i][j][k];
            }

            for (int p = 0; p < 14; ++p) {
              int q = qi[plane-1][p];
              outlet = outlet + f[q][i][j][k];
            }

            outlet = -1.0 + (1.0 / rhoout) * outlet;


            //double uout = u[i-1][j][k];//outlet * n[plane-1][0];
            double uout ;//= outlet * n[plane-1][0];
            double vout ;//= 0.0;//outlet * n[plane-1][1];
            double wout ;//= 0.0;//outlet * n[plane-1][2];
            inlet(i, j, k, &uout, &vout, &wout);

            for (int p = 0; p < 5; ++p) {

              int q = qo[plane-1][p];

              int cn = c[q][0] * n[plane-1][0] + c[q][1] * n[plane-1][1] + c[q][2] * n[plane-1][2];
              int tx = c[q][0] - cn * n[plane-1][0];
              int ty = c[q][1] - cn * n[plane-1][1];
              int tz = c[q][2] - cn * n[plane-1][2];
              //if (j == 2 && k == 50) {printf("q = %d, tangent = %d %d %d\n",q, tx, ty, tz);}

              double cu = c[q][0] * uout + c[q][1] * vout + c[q][2] * wout;
              double tu = tx      * uout + ty      * vout + tz      * wout;

              double sumf  = 0.0;

              for (int r = 0; r < 14; ++r) {
                int s = qi[plane-1][r];
                int cst = c[s][0] * tx            + c[s][1] * ty            + c[s][2] * tz           ;
                int csn = c[s][0] * n[plane-1][0] + c[s][1] * n[plane-1][1] + c[s][2] * n[plane-1][2];
                ////if (j == 2 && k == 50) {printf("%d %d %d\n",s, csn, cst);}
                //if (j == 2 && k == 50) {printf("%d %d\n",s, (1-abs(csn)) * cst);}
                sumf = sumf + cst * f[s][i][j][k] * (1 - abs(csn));
              }

              fstar[q][i][j][k] = f[opp[q]][i][j][k] + 1.0*(1.0 / 3.0) * tu * rhoout
                                                     + 1.0*(1.0 / 6.0) * cu * rhoout
                                                     - 0.0*(1.0 / 2.0) * sumf;
            }
            fstar[plane][i][j][k] = fstar[plane][i][j][k] - (1.0 / 6.0) * uout * rhoout;

            for (int p = 0; p < 14; ++p) {

              int q = qi[plane-1][p];

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
      }
  }

}

void boundary(double f[Q][NX][NY][NZ],
              double fstar[Q][NX][NY][NZ])
{
  int plane;
  char bctype;

  //for (int b = 0; b < 6; b++) {
  //  plane = b + 1;
  //  bctype = bc[b];
  //  impose_bc(plane, bctype, f, fstar);
  //}

  for (int b = 0; b < 6; b++) {
    plane = b + 1;
    bctype = bc[b];
    if (bctype == 'I'){
      impose_bc(plane, bctype, f, fstar);
    }
  }
  for (int b = 0; b < 6; b++) {
    plane = b + 1;
    bctype = bc[b];
    if (bctype == 'O'){
      impose_bc(plane, bctype, f, fstar);
    }
  }
  for (int b = 0; b < 6; b++) {
    plane = b + 1;
    bctype = bc[b];
    if (bctype == 'W'){
      impose_bc(plane, 'W'   , f, fstar);
      impose_bc(plane, bctype, f, fstar);
    }
  }
  for (int b = 0; b < 6; b++) {
    plane = b + 1;
    bctype = bc[b];
    if (bctype == 'P'){
      impose_bc(plane, 'W'   , f, fstar);
      impose_bc(plane, bctype, f, fstar);
    }
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
