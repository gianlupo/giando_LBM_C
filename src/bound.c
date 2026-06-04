#include <stdlib.h>
#include <stdio.h>
#include "param.h"
#include "utils.h"

void face_bc(int    face,
             int    face_n[3],
             char   face_bctype,
             double fstar[Q][LNX_H][LNY_H][LNZ_H])
{

  int imin = 1, imax = LNX;
  int jmin = 1, jmax = LNY;
  int kmin = 1, kmax = LNZ;

  switch (face) {

    case 0: {
      imin = 1;
      imax = 1;
      break;
      }
    case 1: {
      imin = LNX;
      imax = LNX;
      break;
      }
    case 2: {
      jmin = 1;
      jmax = 1;
      break;
      }
    case 3: {
      jmin = LNY;
      jmax = LNY;
      break;
      }
    case 4: {
      kmin = 1;
      kmax = 1;
      break;
      }
    case 5: {
      kmin = LNZ;
      kmax = LNZ;
      break;
      }
  }

  switch (face_bctype) {
    case 'P': {
      // periodic b.c. (handled by halo exchange)
      break;
    }
    case 'W': {
      // wall b.c. (bounceback scheme)
      //for (int k = kmin; k <= kmax; ++k) {
      //  for (int j = jmin; j <= jmax; ++j) {
      //    for (int i = imin; i <= imax; ++i) {

      //      for (int q = 0; q < 5; ++q) {

      //        fstar[qo[face][q]][i][j][k] = fstar[opp[qo[face][q]]][i][j][k];

      //      }
      //    }
      //  }
      //}
      //
      // wall b.c. (Hecht-Harting scheme)
      for (int k = kmin; k <= kmax; ++k) {
        for (int j = jmin; j <= jmax; ++j) {
          for (int i = imin; i <= imax; ++i) {

            double rhowall;
            double uwall, vwall, wwall;

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
                sumf += cst * fstar[s][i][j][k] * (1 - abs(csn));
              }

              fstar[q][i][j][k] = fstar[opp[q]][i][j][k] + (1.0 / 3.0) * tu * rhowall
                                                         + (1.0 / 6.0) * cu * rhowall
                                                         - (1.0 / 2.0) * sumf;
            }
            fstar[face+1][i][j][k] += (1.0 / 6.0) * uwall * rhowall;
          }
        }
      }
      break;
    }
    case 'I': {
      // inlet b.c. (Hecht-Harting scheme)
      for (int k = kmin; k <= kmax; ++k) {
        for (int j = jmin; j <= jmax; ++j) {
          for (int i = imin; i <= imax; ++i) {

            double rhoin;
            double uin, vin, win;

            poi(i, j, k, &uin, &vin, &win);

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
                sumf += cst * fstar[s][i][j][k] * (1 - abs(csn));
              }

              fstar[q][i][j][k] = fstar[opp[q]][i][j][k] + (1.0 / 3.0) * tu * rhoin
                                                         + (1.0 / 6.0) * cu * rhoin
                                                         - (1.0 / 2.0) * sumf;
            }
            fstar[face+1][i][j][k] += (1.0 / 6.0) * uin * rhoin;

          }
        }
      }
      break;
    }
    case 'O': {
      // outlet b.c. (Hecht-Harting scheme)
      for (int k = kmin; k <= kmax; ++k) {
        for (int j = jmin; j <= jmax; ++j) {
          for (int i = imin; i <= imax; ++i) {

            double rhoout;
            double uout, vout, wout;

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
                sumf += cst * fstar[s][i][j][k] * (1 - abs(csn));
              }

              fstar[q][i][j][k] = fstar[opp[q]][i][j][k] + (1.0 / 3.0) * tu * rhoout
                                                         + (1.0 / 6.0) * cu * rhoout
                                                         - (1.0 / 2.0) * sumf;
            }
            fstar[face+1][i][j][k] -= (1.0 / 6.0) * uout * rhoout;

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
             double fstar[Q][LNX_H][LNY_H][LNZ_H])
{
  // determine whether the current rank owns the global faces called
  int owns[2];
  owns[0] = (edge_face[0] == 0) ? (cart_coords[0] == 0)
          : (edge_face[0] == 1) ? (cart_coords[0] == NPROCX-1)
          : (edge_face[0] == 2) ? (cart_coords[1] == 0)
          : (edge_face[0] == 3) ? (cart_coords[1] == NPROCY-1)
          : (edge_face[0] == 4) ? (cart_coords[2] == 0)
          :                       (cart_coords[2] == NPROCZ-1);
  owns[1] = (edge_face[1] == 0) ? (cart_coords[0] == 0)
          : (edge_face[1] == 1) ? (cart_coords[0] == NPROCX-1)
          : (edge_face[1] == 2) ? (cart_coords[1] == 0)
          : (edge_face[1] == 3) ? (cart_coords[1] == NPROCY-1)
          : (edge_face[1] == 4) ? (cart_coords[2] == 0)
          :                       (cart_coords[2] == NPROCZ-1);
 
  //skip interior ranks
  if (!owns[0] || !owns[1]) return;

  int imin = 1, imax = LNX;
  int jmin = 1, jmax = LNY;
  int kmin = 1, kmax = LNZ;

  switch (edge) {

    case 0: {
      imin = 1;
      imax = 1;
      jmin = 1;
      jmax = 1;
      break;
      }
    case 1: {
      imin = 1;
      imax = 1;
      jmin = LNY;
      jmax = LNY;
      break;
      }
    case 2: {
      imin = 1;
      imax = 1;
      kmin = 1;
      kmax = 1;
      break;
      }
    case 3: {
      imin = 1;
      imax = 1;
      kmin = LNZ;
      kmax = LNZ;
      break;
      }
    case 4: {
      imin = LNX;
      imax = LNX;
      jmin = 1;
      jmax = 1;
      break;
      }
    case 5: {
      imin = LNX;
      imax = LNX;
      jmin = LNY;
      jmax = LNY;
      break;
      }
    case 6: {
      imin = LNX;
      imax = LNX;
      kmin = 1;
      kmax = 1;
      break;
      }
    case 7: {
      imin = LNX;
      imax = LNX;
      kmin = LNZ;
      kmax = LNZ;
      break;
      }
    case 8: {
      jmin = 1;
      jmax = 1;
      kmin = 1;
      kmax = 1;
      break;
      }
    case 9: {
      jmin = 1;
      jmax = 1;
      kmin = LNZ;
      kmax = LNZ;
      break;
      }
    case 10: {
      jmin = LNY;
      jmax = LNY;
      kmin = 1;
      kmax = 1;
      break;
      }
    case 11: {
      jmin = LNY;
      jmax = LNY;
      kmin = LNZ;
      kmax = LNZ;
      break;
      }

  }

  int  key = (edge_bctype[0] << 8) | edge_bctype[1];

  switch (key) {
    case ('I' << 8) | 'P':{
      // inlet + periodic
      for (int k = kmin; k <= kmax; ++k) {
        for (int j = jmin; j <= jmax; ++j) {
          for (int i = imin; i <= imax; ++i) {

            double rhoin;
            double uin, vin, win;

            poi(i, j, k, &uin, &vin, &win);

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
                sumf += cst * fstar[s][i][j][k] * (1 - abs(csn));
              }

              fstar[q][i][j][k] = fstar[opp[q]][i][j][k] + (1.0 / 3.0) * tu * rhoin
                                                         + (1.0 / 6.0) * cu * rhoin
                                                         - (1.0 / 2.0) * sumf;
            }
            fstar[edge_face[0]+1][i][j][k] += (1.0 / 6.0) * uin * rhoin;

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

      for (int k = kmin; k <= kmax; ++k) {
        for (int j = jmin; j <= jmax; ++j) {
          for (int i = imin; i <= imax; ++i) {

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
      for (int k = kmin; k <= kmax; ++k) {
        for (int j = jmin; j <= jmax; ++j) {
          for (int i = imin; i <= imax; ++i) {

            double rhoout;
            double uout, vout, wout;

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
                sumf += cst * fstar[s][i][j][k] * (1 - abs(csn));
              }

              fstar[q][i][j][k] = fstar[opp[q]][i][j][k] + 1.0*(1.0 / 3.0) * tu * rhoout
                                                         + 1.0*(1.0 / 6.0) * cu * rhoout
                                                         - 1.0*(1.0 / 2.0) * sumf;
            }
            fstar[edge_face[0]+1][i][j][k] -= (1.0 / 6.0) * uout * rhoout;

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

      for (int k = kmin; k <= kmax; ++k) {
        for (int j = jmin; j <= jmax; ++j) {
          for (int i = imin; i <= imax; ++i) {

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
      for (int k = kmin; k <= kmax; ++k) {
        for (int j = jmin; j <= jmax; ++j) {
          for (int i = imin; i <= imax; ++i) {

            double rhowall;
            double uwall, vwall, wwall;

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
                sumf += cst * fstar[s][i][j][k] * (1 - abs(csn));
              }

              fstar[q][i][j][k] = fstar[opp[q]][i][j][k] + (1.0 / 3.0) * tu * rhowall
                                                         + (1.0 / 6.0) * cu * rhowall
                                                         - (1.0 / 2.0) * sumf;
            }
            fstar[edge_face[1]+1][i][j][k] += (1.0 / 6.0) * uwall * rhowall;
          }
        }
      }
      break;
    }
  }
}

void obst_bounce(int    obst[LNX_H][LNY_H][LNZ_H],
                 double obst_omega[3],
                 double f[Q][LNX_H][LNY_H][LNZ_H])
{

  for (int k = 1; k <= LNZ; ++k) {
    for (int j = 1; j <= LNY; ++j) {
      for (int i = 1; i <= LNX; ++i) {
        if (obst[i][j][k] != 1) continue;
        // simple implementation
        //
        //double ftmp[Q];
        //
        //for (int q = 0; q < Q; ++q){
        //  ftmp[q] = f[opp[q]][i][j][k];
        //}
        //for (int q = 0; q < Q; ++q){
        //  f[q][i][j][k] = ftmp[q];
        //}

        // Ladd and Verberg link-bounce-back
        //
        double ftmp[Q];
        double ftmp2[Q];

        for (int q = 0; q < Q; ++q) {

          int ip = i + c[q][0];
          int jp = j + c[q][1];
          int kp = k + c[q][2];

          double rx = (double)gi(i) - OBST_XC;
          double ry = (double)gj(j) - OBST_YC;
          double rz = (double)gk(k) - OBST_ZC;

          double ubx = obst_omega[1]*rz - obst_omega[2]*ry;
          double uby = obst_omega[2]*rx - obst_omega[0]*rz;
          double ubz = obst_omega[0]*ry - obst_omega[1]*rx;
          double ubc = ubx*c[q][0] + uby*c[q][1] + ubz*c[q][2];

          if (obst[ip][jp][kp] == 0) {
            ftmp[q]  = f[opp[q]][i][j][k]   -6.0*weights[q]*ubc;
            ftmp2[q] = f[opp[q]][ip][jp][kp]+6.0*weights[q]*ubc;
          } else {
            ftmp [q]=f[q][i][j][k];
            ftmp2[q]=f[q][ip][jp][kp];
          }
        }
        for (int q = 0; q < Q; ++q) {

          int ip = i + c[q][0];
          int jp = j + c[q][1];
          int kp = k + c[q][2];

          if (obst[ip][jp][kp] == 0) {
            f[q][i][j][k]    = ftmp[q];
            f[q][ip][jp][kp] = ftmp2[q];
          }
        }

      }
    }
  }

}

void boundary(double fstar[Q][LNX_H][LNY_H][LNZ_H])
{
  // determine which faces are owned by the current rank
  int owns_face[6] = { 
    (cart_coords[0] == 0), 
    (cart_coords[0] == NPROCX-1), 
    (cart_coords[1] == 0), 
    (cart_coords[1] == NPROCY-1), 
    (cart_coords[2] == 0), 
    (cart_coords[2] == NPROCZ-1) 
  }; 

  // BCs on faces: I first, then O, then W
  for (int pass=0; pass<3; ++pass) {
    char target = (pass==0)?'I':(pass==1)?'O':'W';
    for (int face=0; face<6; ++face) {
      if (!owns_face[face]) continue;
      if (bc[face] != target) continue;
      int face_n[3] = {n[face][0], n[face][1], n[face][2]};
      face_bc(face, face_n, bc[face], fstar);
    }
  }

  // BCs on edges
  for (int edge=0; edge<12; ++edge) {
    int edge_face[2], edge_n[2][3];
    char edge_bctype[2];
    edge_face[0] = (edge/4) + 0 + (edge/10);
    edge_face[1] =  2 + (edge%4) + 2*(edge/8) - 2*(edge/10);
    edge_bctype[0] = bc[edge_face[0]];
    edge_bctype[1] = bc[edge_face[1]];
    for (int m = 0; m < 3; m++){
        edge_n[0][m] = n[edge_face[0]][m];
        edge_n[1][m] = n[edge_face[1]][m];
    }
    edge_bc(edge, edge_face, edge_n, edge_bctype, fstar);
  }

}
