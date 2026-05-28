#include "param.h"
#include <stdlib.h>
#include <stdio.h>

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

  double z = (double)k / (NZ-1);
  *uin = 4.0 * U0 * z * (1.0 - z);
  *vin = 0.0;
  *win = 0.0;

}

void face_bc(int face,
             int face_n[3],
             char face_bctype,
             double fstar[Q][NX][NY][NZ])
{

  int imin;
  int imax;
  int jmin;
  int jmax;
  int kmin;
  int kmax;

  imin = 1;
  imax = NX - 1;
  jmin = 1;
  jmax = NY - 1;
  kmin = 1;
  kmax = NZ - 1;

  switch (face) {

    case 0: {
      imin = 0;
      imax = 1;
      break;
      }
    case 1: {
      imin = NX - 1;
      imax = NX;
      break;
      }
    case 2: {
      jmin = 0;
      jmax = 1;
      break;
      }
    case 3: {
      jmin = NY - 1;
      jmax = NY;
      break;
      }
    case 4: {
      kmin = 0;
      kmax = 1;
      break;
      }
    case 5: {
      kmin = NZ - 1;
      kmax = NZ;
      break;
      }
  }

  switch (face_bctype) {
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

              fstar[q][i][j][k] = fstar[q][istar][jstar][kstar];
            }
          }
        }
      }
      break;
      }
    case 'W': {
      // wall b.c. (bounceback scheme)
      //for (int k = kmin; k < kmax; ++k) {
      //  for (int j = jmin; j < jmax; ++j) {
      //    for (int i = imin; i < imax; ++i) {

      //      for (int q = 0; q < 5; ++q) {

      //        fstar[qo[face][q]][i][j][k] = fstar[opp[qo[face][q]]][i][j][k];

      //      }
      //    }
      //  }
      //}
      //
      // wall b.c. (Hecht-Harting scheme)
      for (int k = kmin; k < kmax; ++k) {
        for (int j = jmin; j < jmax; ++j) {
          for (int i = imin; i < imax; ++i) {

            double rhowall;
            double uwall;
            double vwall;
            double wwall;

            uwall = 0.0;
            vwall = 0.0;
            wwall = 0.0;

            double sumfqi = 0.0, sumfoppqo = 0.0;
            for (int p = 0; p < 14; ++p){
              sumfqi += fstar[qi[face][p]][i][j][k];
            }
            for (int p = 0; p < 5; ++p){
              sumfoppqo += fstar[opp[qo[face][p]]][i][j][k];
            }

            rhowall = (sumfqi + sumfoppqo) / (1.0 - uwall);

            for (int p = 0; p < 5; ++p) {

              int q = qo[face][p];

              int cn = c[q][0] * face_n[0] + c[q][1] * face_n[1] + c[q][2] * face_n[2];
              int tx = c[q][0] - cn * face_n[0];
              int ty = c[q][1] - cn * face_n[1];
              int tz = c[q][2] - cn * face_n[2];

              double cu = c[q][0] * uwall + c[q][1] * vwall + c[q][2] * wwall;
              double tu = tx      * uwall + ty      * vwall + tz      * wwall;

              double sumf  = 0.0;

              for (int r = 0; r < 14; ++r) {
                int s = qi[face][r];
                int cst = c[s][0] * tx        + c[s][1] * ty        + c[s][2] * tz       ;
                int csn = c[s][0] * face_n[0] + c[s][1] * face_n[1] + c[s][2] * face_n[2];
                sumf = sumf + cst * fstar[s][i][j][k] * (1 - abs(csn));
              }

              fstar[q][i][j][k] = fstar[opp[q]][i][j][k] + 1.0*(1.0 / 3.0) * tu * rhowall
                                                         + 1.0*(1.0 / 6.0) * cu * rhowall
                                                         - 1.0*(1.0 / 2.0) * sumf;
            }
            fstar[face+1][i][j][k] = fstar[face+1][i][j][k] + (1.0 / 6.0) * uwall * rhowall;
          }
        }
      }
      break;
      }
    case 'I': {
      // inlet b.c. (Hecht-Harting scheme)
      for (int k = kmin; k < kmax; ++k) {
        for (int j = jmin; j < jmax; ++j) {
          for (int i = imin; i < imax; ++i) {

            double rhoin;
            double uin;
            double vin;
            double win;

            inlet(i, j, k, &uin, &vin, &win);

            double sumfqi = 0.0, sumfoppqo = 0.0;
            for (int p = 0; p < 14; ++p){
              sumfqi += fstar[qi[face][p]][i][j][k];
            }
            for (int p = 0; p < 5; ++p){
              sumfoppqo += fstar[opp[qo[face][p]]][i][j][k];
            }

            rhoin = (sumfqi + sumfoppqo) / (1.0 - uin);

            for (int p = 0; p < 5; ++p) {

              int q = qo[face][p];

              int cn = c[q][0] * face_n[0] + c[q][1] * face_n[1] + c[q][2] * face_n[2];
              int tx = c[q][0] - cn * face_n[0];
              int ty = c[q][1] - cn * face_n[1];
              int tz = c[q][2] - cn * face_n[2];

              double cu = c[q][0] * uin + c[q][1] * vin + c[q][2] * win;
              double tu = tx      * uin + ty      * vin + tz      * win;

              double sumf  = 0.0;

              for (int r = 0; r < 14; ++r) {
                int s = qi[face][r];
                int cst = c[s][0] * tx        + c[s][1] * ty        + c[s][2] * tz       ;
                int csn = c[s][0] * face_n[0] + c[s][1] * face_n[1] + c[s][2] * face_n[2];
                sumf = sumf + cst * fstar[s][i][j][k] * (1 - abs(csn));
              }

              fstar[q][i][j][k] = fstar[opp[q]][i][j][k] + 1.0*(1.0 / 3.0) * tu * rhoin
                                                         + 1.0*(1.0 / 6.0) * cu * rhoin
                                                         - 1.0*(1.0 / 2.0) * sumf;
            }
            fstar[face+1][i][j][k] = fstar[face+1][i][j][k] + (1.0 / 6.0) * uin * rhoin;

          }
        }
      }
      break;
      }
    case 'O': {
      // outlet b.c. (Hecht-Harting scheme)
      for (int k = kmin; k < kmax; ++k) {
        for (int j = jmin; j < jmax; ++j) {
          for (int i = imin; i < imax; ++i) {

            double rhoout;
            double uout;
            double vout;
            double wout;

            rhoout = 1.0;

            double sumfqi = 0.0, sumfoppqo = 0.0;
            for (int p = 0; p < 14; ++p){
              sumfqi += fstar[qi[face][p]][i][j][k];
            }
            for (int p = 0; p < 5; ++p){
              sumfoppqo += fstar[opp[qo[face][p]]][i][j][k];
            }

            uout = -1.0 + (sumfqi + sumfoppqo) / rhoout;
            vout = 0.0;
            wout = 0.0;

            //inlet(i, j, k, &uout, &vout, &wout);

            for (int p = 0; p < 5; ++p) {

              int q = qo[face][p];

              int cn = c[q][0] * face_n[0] + c[q][1] * face_n[1] + c[q][2] * face_n[2];
              int tx = c[q][0] - cn * face_n[0];
              int ty = c[q][1] - cn * face_n[1];
              int tz = c[q][2] - cn * face_n[2];

              double cu = c[q][0] * uout + c[q][1] * vout + c[q][2] * wout;
              double tu = tx      * uout + ty      * vout + tz      * wout;

              double sumf  = 0.0;

              for (int r = 0; r < 14; ++r) {
                int s = qi[face][r];
                int cst = c[s][0] * tx        + c[s][1] * ty        + c[s][2] * tz;
                int csn = c[s][0] * face_n[0] + c[s][1] * face_n[1] + c[s][2] * face_n[2];
                sumf = sumf + cst * fstar[s][i][j][k] * (1 - abs(csn));
              }

              fstar[q][i][j][k] = fstar[opp[q]][i][j][k] + 1.0*(1.0 / 3.0) * tu * rhoout
                                                         + 1.0*(1.0 / 6.0) * cu * rhoout
                                                         - 1.0*(1.0 / 2.0) * sumf;
            }
            fstar[face+1][i][j][k] = fstar[face+1][i][j][k] - (1.0 / 6.0) * uout * rhoout;

          }
        }
      }
      break;
      }
  }

}

void edge_bc(int edge,
             int edge_face[2],
             int edge_n[2][3],
             char edge_bctype[2],
             double fstar[Q][NX][NY][NZ])
{

  int imin;
  int imax;
  int jmin;
  int jmax;
  int kmin;
  int kmax;

  imin = 1;
  imax = NX - 1;
  jmin = 0;
  jmax = NY;
  kmin = 1;
  kmax = NZ - 1;

  switch (edge) {

    case 0: {
      imin = 0;
      imax = 1;
      jmin = 0;
      jmax = 1;
      break;
      }
    case 1: {
      imin = 0;
      imax = 1;
      jmin = NY - 1;
      jmax = NY;
      break;
      }
    case 2: {
      imin = 0;
      imax = 1;
      kmin = 0;
      kmax = 1;
      break;
      }
    case 3: {
      imin = 0;
      imax = 1;
      kmin = NZ - 1;
      kmax = NZ;
      break;
      }
    case 4: {
      imin = NX - 1;
      imax = NX;
      jmin = 0;
      jmax = 1;
      break;
      }
    case 5: {
      imin = NX - 1;
      imax = NX;
      jmin = NY - 1;
      jmax = NY;
      break;
      }
    case 6: {
      imin = NX - 1;
      imax = NX;
      kmin = 0;
      kmax = 1;
      break;
      }
    case 7: {
      imin = NX - 1;
      imax = NX;
      kmin = NZ - 1;
      kmax = NZ;
      break;
      }
    case 8: {
      jmin = 0;
      jmax = 1;
      kmin = 0;
      kmax = 1;
      break;
      }
    case 9: {
      jmin = 0;
      jmax = 1;
      kmin = NZ - 1;
      kmax = NZ;
      break;
      }
    case 10: {
      jmin = NY - 1;
      jmax = NY;
      kmin = 0;
      kmax = 1;
      break;
      }
    case 11: {
      jmin = NY - 1;
      jmax = NY;
      kmin = NZ - 1;
      kmax = NZ;
      break;
      }

  }

  switch ((edge_bctype[0] << 8) | edge_bctype[1]) {
    case ('I' << 8) | 'P':{
      // inlet + periodic
      for (int k = kmin; k < kmax; ++k) {
        for (int j = jmin; j < jmax; ++j) {
          for (int i = imin; i < imax; ++i) {

            double rhoin;
            double uin;
            double vin;
            double win;

            inlet(i, j, k, &uin, &vin, &win);

            double sumfqi = 0.0, sumfoppqo = 0.0;
            for (int p = 0; p < 14; ++p){
              sumfqi += fstar[qi[edge_face[0]][p]][i][j][k];
            }
            for (int p = 0; p < 5; ++p){
              sumfoppqo += fstar[opp[qo[edge_face[0]][p]]][i][j][k];
            }

            rhoin = (sumfqi + sumfoppqo) / (1.0 - uin);

            for (int p = 0; p < 5; ++p) {

              int q = qo[edge_face[0]][p];

              int cn = c[q][0] * edge_n[0][0] + c[q][1] * edge_n[0][1] + c[q][2] * edge_n[0][2];
              int tx = c[q][0] - cn * edge_n[0][0];
              int ty = c[q][1] - cn * edge_n[0][1];
              int tz = c[q][2] - cn * edge_n[0][2];

              double cu = c[q][0] * uin + c[q][1] * vin + c[q][2] * win;
              double tu = tx      * uin + ty      * vin + tz      * win;

              double sumf  = 0.0;

              for (int r = 0; r < 14; ++r) {
                int s = qi[edge_face[0]][r];
                int cst = c[s][0] * tx           + c[s][1] * ty           + c[s][2] * tz          ;
                int csn = c[s][0] * edge_n[0][0] + c[s][1] * edge_n[0][1] + c[s][2] * edge_n[0][2];
                sumf = sumf + cst * fstar[s][i][j][k] * (1 - abs(csn));
              }

              fstar[q][i][j][k] = fstar[opp[q]][i][j][k] + 1.0*(1.0 / 3.0) * tu * rhoin
                                                         + 1.0*(1.0 / 6.0) * cu * rhoin
                                                         - 1.0*(1.0 / 2.0) * sumf;
            }
            fstar[edge_face[0]+1][i][j][k] = fstar[edge_face[0]+1][i][j][k] + (1.0 / 6.0) * uin * rhoin;

          }
        }
      }
      break;
      }
    case ('I' << 8) | 'W':{
      // inlet + wall
      int edge_par[3];
      int edge_nor[3];

      int is_buried_inlet[Q];
      int is_buried_wall[Q];
      int buried_inlet;
      int buried_wall;

      edge_par[0] = edge_n[0][1] * edge_n[1][2] - edge_n[0][2] * edge_n[1][1];
      edge_par[1] = edge_n[0][2] * edge_n[1][0] - edge_n[0][0] * edge_n[1][2];
      edge_par[2] = edge_n[0][0] * edge_n[1][1] - edge_n[0][1] * edge_n[1][0];

      edge_nor[0] = edge_n[0][0] + edge_n[1][0];
      edge_nor[1] = edge_n[0][1] + edge_n[1][1];
      edge_nor[2] = edge_n[0][2] + edge_n[1][2];

      for (int k = kmin; k < kmax; ++k) {
        for (int j = jmin; j < jmax; ++j) {
          for (int i = imin; i < imax; ++i) {

            for (int q = 0; q < Q; ++q) {
              is_buried_inlet[q] = 0;
              is_buried_wall[q]  = 0;
            }
            buried_inlet = 0;
            buried_wall  = 0;

            for (int p = 0; p < 5; ++p) {

              int q_inlet = qo[edge_face[0]][p];
              int q_wall  = qo[edge_face[1]][p];

              int cn_inlet = c[q_inlet][0] * edge_n[0][0] + c[q_inlet][1] * edge_n[0][1] + c[q_inlet][2] * edge_n[0][2];
              int cn_wall  = c[q_wall][0]  * edge_n[1][0] + c[q_wall][1]  * edge_n[1][1] + c[q_wall][2]  * edge_n[1][2];

              int tx_inlet = c[q_inlet][0] - cn_inlet * edge_n[0][0];
              int ty_inlet = c[q_inlet][1] - cn_inlet * edge_n[0][1];
              int tz_inlet = c[q_inlet][2] - cn_inlet * edge_n[0][2];
              int tx_wall  = c[q_wall][0]  - cn_wall  * edge_n[1][0];
              int ty_wall  = c[q_wall][1]  - cn_wall  * edge_n[1][1];
              int tz_wall  = c[q_wall][2]  - cn_wall  * edge_n[1][2];

              int cpar_inlet = c[q_inlet][0] * edge_par[0] + c[q_inlet][1] * edge_par[1] + c[q_inlet][2] * edge_par[2];
              int cpar_wall  = c[q_wall][0]  * edge_par[0] + c[q_wall][1]  * edge_par[1] + c[q_wall][2]  * edge_par[2];

              int cnor_inlet = (c[q_inlet][0] * edge_nor[0]   + c[q_inlet][1] * edge_nor[1]   + c[q_inlet][2] * edge_nor[2]) /
                               (c[q_inlet][0] * c[q_inlet][0] + c[q_inlet][1] * c[q_inlet][1] + c[q_inlet][2] * c[q_inlet][2]);
              int cnor_wall  = (c[q_wall][0]  * edge_nor[0]   + c[q_wall][1]  * edge_nor[1]   + c[q_wall][2]  * edge_nor[2]) /
                               (c[q_wall][0]  * c[q_wall][0]  + c[q_wall][1]  * c[q_wall][1]  + c[q_wall][2]  * c[q_wall][2]);

              is_buried_inlet[q_inlet] = (1 - abs(cpar_inlet)) * (1 - abs(cnor_inlet));
              is_buried_wall[q_wall]   = (1 - abs(cpar_wall))  * (1 - abs(cnor_wall));
              buried_inlet += q_inlet * is_buried_inlet[q_inlet];
              buried_wall  += q_wall  * is_buried_wall[q_wall];

              double sumf_inlet  = 0.0;
              double sumf_wall   = 0.0;

              for (int s = 0; s < Q; ++s) {

                int cst_inlet = c[s][0] * tx_inlet     + c[s][1] * ty_inlet     + c[s][2] * tz_inlet    ;
                int cst_wall  = c[s][0] * tx_wall      + c[s][1] * ty_wall      + c[s][2] * tz_wall     ;
                int csn_inlet = c[s][0] * edge_n[0][0] + c[s][1] * edge_n[0][1] + c[s][2] * edge_n[0][2];
                int csn_wall  = c[s][0] * edge_n[1][0] + c[s][1] * edge_n[1][1] + c[s][2] * edge_n[1][2];

                sumf_inlet += fstar[s][i][j][k] * cst_inlet * (1 - abs(csn_inlet)) * (1 - abs(csn_wall));
                sumf_wall  += fstar[s][i][j][k] * cst_wall  * (1 - abs(csn_inlet)) * (1 - abs(csn_wall));
              }

              fstar[q_inlet][i][j][k] = fstar[opp[q_inlet]][i][j][k] - 0.25 * sumf_inlet;
              fstar[q_wall][i][j][k]  = fstar[opp[q_wall]][i][j][k]  - 0.25 * sumf_wall;
            }

            double sum_buried = 0.0;
            for (int r = 1; r < Q; ++r) {
              sum_buried += fstar[r][i][j][k] * (1 - (is_buried_inlet[r] + is_buried_wall[r]));
            }
            fstar[buried_inlet][i][j][k] = (1.0 / 22.0) * sum_buried;
            fstar[buried_wall][i][j][k]  = (1.0 / 22.0) * sum_buried;
            fstar[0][i][j][k] = (12.0 / 22.0) * sum_buried;
          }
        }
      }
      break;
      }
    case ('O' << 8) | 'P':{
      // outlet + periodic
      for (int k = kmin; k < kmax; ++k) {
        for (int j = jmin; j < jmax; ++j) {
          for (int i = imin; i < imax; ++i) {

            double rhoout;
            double uout;
            double vout;
            double wout;

            rhoout = 1.0;

            double sumfqi = 0.0, sumfoppqo = 0.0;
            for (int p = 0; p < 14; ++p){
              sumfqi += fstar[qi[edge_face[0]][p]][i][j][k];
            }
            for (int p = 0; p < 5; ++p){
              sumfoppqo += fstar[opp[qo[edge_face[0]][p]]][i][j][k];
            }

            uout = -1.0 + (sumfqi + sumfoppqo) / rhoout;
            vout = 0.0;
            wout = 0.0;

            //inlet(i, j, k, &uout, &vout, &wout);

            for (int p = 0; p < 5; ++p) {

              int q = qo[edge_face[0]][p];

              int cn = c[q][0] * edge_n[0][0] + c[q][1] * edge_n[0][1] + c[q][2] * edge_n[0][2];
              int tx = c[q][0] - cn * edge_n[0][0];
              int ty = c[q][1] - cn * edge_n[0][1];
              int tz = c[q][2] - cn * edge_n[0][2];

              double cu = c[q][0] * uout + c[q][1] * vout + c[q][2] * wout;
              double tu = tx      * uout + ty      * vout + tz      * wout;

              double sumf  = 0.0;

              for (int r = 0; r < 14; ++r) {
                int s = qi[edge_face[0]][r];
                int cst = c[s][0] * tx           + c[s][1] * ty           + c[s][2] * tz          ;
                int csn = c[s][0] * edge_n[0][0] + c[s][1] * edge_n[0][1] + c[s][2] * edge_n[0][2];
                sumf = sumf + cst * fstar[s][i][j][k] * (1 - abs(csn));
              }

              fstar[q][i][j][k] = fstar[opp[q]][i][j][k] + 1.0*(1.0 / 3.0) * tu * rhoout
                                                         + 1.0*(1.0 / 6.0) * cu * rhoout
                                                         - 1.0*(1.0 / 2.0) * sumf;
            }
            fstar[edge_face[0]+1][i][j][k] = fstar[edge_face[0]+1][i][j][k] - (1.0 / 6.0) * uout * rhoout;

          }
        }
      }
      break;
      }
    case ('O' << 8) | 'W':{
      // outlet + wall
      int edge_par[3];
      int edge_nor[3];

      int is_buried_outlet[Q];
      int is_buried_wall[Q];
      int buried_outlet;
      int buried_wall;

      edge_par[0] = edge_n[0][1] * edge_n[1][2] - edge_n[0][2] * edge_n[1][1];
      edge_par[1] = edge_n[0][2] * edge_n[1][0] - edge_n[0][0] * edge_n[1][2];
      edge_par[2] = edge_n[0][0] * edge_n[1][1] - edge_n[0][1] * edge_n[1][0];

      edge_nor[0] = edge_n[0][0] + edge_n[1][0];
      edge_nor[1] = edge_n[0][1] + edge_n[1][1];
      edge_nor[2] = edge_n[0][2] + edge_n[1][2];

      for (int k = kmin; k < kmax; ++k) {
        for (int j = jmin; j < jmax; ++j) {
          for (int i = imin; i < imax; ++i) {

            for (int q = 0; q < Q; ++q) {
              is_buried_outlet[q] = 0;
              is_buried_wall[q]  = 0;
            }
            buried_outlet = 0;
            buried_wall   = 0;

            for (int p = 0; p < 5; ++p) {

              int q_outlet = qo[edge_face[0]][p];
              int q_wall   = qo[edge_face[1]][p];

              int cn_outlet = c[q_outlet][0] * edge_n[0][0] + c[q_outlet][1] * edge_n[0][1] + c[q_outlet][2] * edge_n[0][2];
              int cn_wall   = c[q_wall][0]   * edge_n[1][0] + c[q_wall][1]   * edge_n[1][1] + c[q_wall][2]   * edge_n[1][2];

              int tx_outlet = c[q_outlet][0] - cn_outlet * edge_n[0][0];
              int ty_outlet = c[q_outlet][1] - cn_outlet * edge_n[0][1];
              int tz_outlet = c[q_outlet][2] - cn_outlet * edge_n[0][2];
              int tx_wall   = c[q_wall][0]   - cn_wall   * edge_n[1][0];
              int ty_wall   = c[q_wall][1]   - cn_wall   * edge_n[1][1];
              int tz_wall   = c[q_wall][2]   - cn_wall   * edge_n[1][2];

              int cpar_outlet = c[q_outlet][0] * edge_par[0] + c[q_outlet][1] * edge_par[1] + c[q_outlet][2] * edge_par[2];
              int cpar_wall   = c[q_wall][0]   * edge_par[0] + c[q_wall][1]   * edge_par[1] + c[q_wall][2]   * edge_par[2];

              int cnor_outlet = (c[q_outlet][0] * edge_nor[0]    + c[q_outlet][1] * edge_nor[1]    + c[q_outlet][2] * edge_nor[2]) /
                                (c[q_outlet][0] * c[q_outlet][0] + c[q_outlet][1] * c[q_outlet][1] + c[q_outlet][2] * c[q_outlet][2]);
              int cnor_wall   = (c[q_wall][0]  * edge_nor[0]   + c[q_wall][1]  * edge_nor[1]   + c[q_wall][2]  * edge_nor[2]) /
                                (c[q_wall][0]  * c[q_wall][0]  + c[q_wall][1]  * c[q_wall][1]  + c[q_wall][2]  * c[q_wall][2]);

              is_buried_outlet[q_outlet] = (1 - abs(cpar_outlet)) * (1 - abs(cnor_outlet));
              is_buried_wall[q_wall]     = (1 - abs(cpar_wall))   * (1 - abs(cnor_wall));
              buried_outlet += q_outlet * is_buried_outlet[q_outlet];
              buried_wall   += q_wall   * is_buried_wall[q_wall];

              double sumf_outlet  = 0.0;
              double sumf_wall    = 0.0;

              for (int s = 0; s < Q; ++s) {

                int cst_outlet = c[s][0] * tx_outlet    + c[s][1] * ty_outlet    + c[s][2] * tz_outlet   ;
                int cst_wall   = c[s][0] * tx_wall      + c[s][1] * ty_wall      + c[s][2] * tz_wall     ;
                int csn_outlet = c[s][0] * edge_n[0][0] + c[s][1] * edge_n[0][1] + c[s][2] * edge_n[0][2];
                int csn_wall   = c[s][0] * edge_n[1][0] + c[s][1] * edge_n[1][1] + c[s][2] * edge_n[1][2];

                sumf_outlet += fstar[s][i][j][k] * cst_outlet * (1 - abs(csn_outlet)) * (1 - abs(csn_wall));
                sumf_wall   += fstar[s][i][j][k] * cst_wall   * (1 - abs(csn_outlet)) * (1 - abs(csn_wall));
              }

              fstar[q_outlet][i][j][k] = fstar[opp[q_outlet]][i][j][k] - 0.25 * sumf_outlet;
              fstar[q_wall][i][j][k]   = fstar[opp[q_wall]][i][j][k]   - 0.25 * sumf_wall;
            }

            double sum_buried = 0.0;
            for (int r = 1; r < Q; ++r) {
              sum_buried += fstar[r][i][j][k] * (1 - (is_buried_outlet[r] + is_buried_wall[r]));
            }
            double rhoout = 1.0;
            fstar[buried_outlet][i][j][k] = (1.0 / 22.0) * sum_buried;
            fstar[buried_wall][i][j][k]   = (1.0 / 22.0) * sum_buried;
            fstar[0][i][j][k] = (12.0 / 14.0) * (rhoout - sum_buried);
          }
        }
      }
      break;
      }
    case ('P' << 8) | 'W':{
      // periodic + wall
      for (int k = kmin; k < kmax; ++k) {
        for (int j = jmin; j < jmax; ++j) {
          for (int i = imin; i < imax; ++i) {

            double rhowall;
            double uwall;
            double vwall;
            double wwall;

            uwall = 0.0;
            vwall = 0.0;
            wwall = 0.0;

            double sumfqi = 0.0, sumfoppqo = 0.0;
            for (int p = 0; p < 14; ++p){
              sumfqi += fstar[qi[edge_face[1]][p]][i][j][k];
            }
            for (int p = 0; p < 5; ++p){
              sumfoppqo += fstar[opp[qo[edge_face[1]][p]]][i][j][k];
            }

            rhowall = (sumfqi + sumfoppqo) / (1.0 - uwall);

            for (int p = 0; p < 5; ++p) {

              int q = qo[edge_face[1]][p];

              int cn = c[q][0] * edge_n[1][0] + c[q][1] * edge_n[1][1] + c[q][2] * edge_n[1][2];
              int tx = c[q][0] - cn * edge_n[1][0];
              int ty = c[q][1] - cn * edge_n[1][1];
              int tz = c[q][2] - cn * edge_n[1][2];

              double cu = c[q][0] * uwall + c[q][1] * vwall + c[q][2] * wwall;
              double tu = tx      * uwall + ty      * vwall + tz      * wwall;

              double sumf  = 0.0;

              for (int r = 0; r < 14; ++r) {
                int s = qi[edge_face[1]][r];
                int cst = c[s][0] * tx           + c[s][1] * ty           + c[s][2] * tz          ;
                int csn = c[s][0] * edge_n[1][0] + c[s][1] * edge_n[1][1] + c[s][2] * edge_n[1][2];
                sumf = sumf + cst * fstar[s][i][j][k] * (1 - abs(csn));
              }

              fstar[q][i][j][k] = fstar[opp[q]][i][j][k] + 1.0*(1.0 / 3.0) * tu * rhowall
                                                         + 1.0*(1.0 / 6.0) * cu * rhowall
                                                         - 1.0*(1.0 / 2.0) * sumf;
            }
            fstar[edge_face[1]+1][i][j][k] = fstar[edge_face[1]+1][i][j][k] + (1.0 / 6.0) * uwall * rhowall;
          }
        }
      }
      break;
      }
  }
}

void obst_bounce(int obst[NX][NY][NZ],
                 double f[Q][NX][NY][NZ])
{

  for (int k = 0; k < NZ; ++k) {
    for (int j = 0; j < NY; ++j) {
      for (int i = 0; i < NX; ++i) {

        // simple implementation
        //
        //double ftmp[Q];
        //
        //for (int q = 0; q < Q; ++q){
        //  if (obst[i][j][k] == 1) {
        //    ftmp[q] = f[opp[q]][i][j][k];
        //  }
        //}
        //for (int q = 0; q < Q; ++q){
        //  if (obst[i][j][k] == 1) {
        //    f[q][i][j][k] = ftmp[q];
        //  }
        //}

        // Ladd and Verberg link-bounce-back
        //
        double ftmp[Q];
        double ftmp2[Q];

        for (int q = 0; q < Q; ++q) {

          if (obst[i][j][k] == 1) {

            int ip;
            int jp;
            int kp;

            ip = mod(i + c[q][0], NX);
            jp = mod(j + c[q][1], NY);
            kp = mod(k + c[q][2], NZ);

            if (obst[ip][jp][kp] == 0) {
              ftmp[q]  = f[opp[q]][i][j][k];
              ftmp2[q] = f[opp[q]][ip][jp][kp];
            }
          }
        }
        for (int q = 0; q < Q; ++q) {

          if (obst[i][j][k] == 1) {

            int ip;
            int jp;
            int kp;

            ip = mod(i + c[q][0], NX);
            jp = mod(j + c[q][1], NY);
            kp = mod(k + c[q][2], NZ);

            if (obst[ip][jp][kp] == 0) {
              f[q][i][j][k]    = ftmp[q];
              f[q][ip][jp][kp] = ftmp2[q];
            }
          }
        }

      }
    }
  }

}

void boundary(double fstar[Q][NX][NY][NZ])
{

  int edge_face[2];

  int face_n[3];
  int edge_n[2][3];

  char face_bctype;
  char edge_bctype[2];

  //for (int face = 0; face < 6; face++) {

  //  face_bctype = bc[face];
  //  if (face_bctype == 'P') {

  //    for (int m = 0; m < 3; m++){
  //        face_n[m] = n[face][m];
  //    }

  //    face_bc(face, face_n, face_bctype, fstar);
  //  }
  //}
  for (int face = 0; face < 6; face++) {

    face_bctype = bc[face];
    if (face_bctype == 'I') {

      for (int m = 0; m < 3; m++){
          face_n[m] = n[face][m];
      }

      face_bc(face, face_n, face_bctype, fstar);
    }
  }
  for (int face = 0; face < 6; face++) {

    face_bctype = bc[face];
    if (face_bctype == 'O') {

      for (int m = 0; m < 3; m++){
          face_n[m] = n[face][m];
      }

      face_bc(face, face_n, face_bctype, fstar);
    }
  }
  for (int face = 0; face < 6; face++) {

    face_bctype = bc[face];
    if (face_bctype == 'W') {

      for (int m = 0; m < 3; m++){
          face_n[m] = n[face][m];
      }

      face_bc(face, face_n, face_bctype, fstar);
    }
  }
  for (int edge = 0; edge < 12; edge++) {

    edge_face[0] = (edge / 4) + 0 + (edge / 10);
    edge_face[1] =  2 + (edge % 4)+ 2 * (edge / 8) - 2 * (edge / 10);

    edge_bctype[0] = bc[edge_face[0]];
    edge_bctype[1] = bc[edge_face[1]];

    for (int m = 0; m < 3; m++){
        edge_n[0][m] = n[edge_face[0]][m];
        edge_n[1][m] = n[edge_face[1]][m];
    }

    edge_bc(edge, edge_face, edge_n, edge_bctype, fstar);

  }

}
