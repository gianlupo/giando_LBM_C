#include <math.h>
#include <stdlib.h>
#include "param.h"

static int mod(int a, int b)
{
  int r = a % b;
  return (r < 0) ? r + b : r;
}

void build_obstacle(int obst[NX][NY][NZ])
{

  int shape = 0;

  switch (shape) {

    case 0: {
      for (int k = 0; k < NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
          for (int i = 0; i < NX; ++i) {

            double x = (double)i - XC;
            double y = (double)j - YC;
            double z = (double)k - ZC;

            if ((x*x + z*z) <= (OBST_R*OBST_R) ) {
              obst[i][j][k] = 1;
            } else {
              obst[i][j][k] = 0;
            }
          }
        }
      }
      break;
    }
    case 1: {
      double a = 1.0;
      double b = 5.0;
      for (int k = 0; k < NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
          for (int i = 0; i < NX; ++i) {

            double x = (double)i - XC;
            double y = (double)j - YC;
            double z = (double)k - ZC;

            if ( (abs(x) <= a * OBST_R) && (abs(z) <= b * OBST_R) ) {
              obst[i][j][k] = 1;
            } else {
              obst[i][j][k] = 0;
            }
          }
        }
      }
      break;
    }
  }

}

void compute_obst_torque(int obst[NX][NY][NZ],
                         double fstar[Q][NX][NY][NZ],
                         double obst_omega[3],
                         double obst_torq[3])
{

  obst_torq[0] = 0.0;
  obst_torq[1] = 0.0;
  obst_torq[2] = 0.0;

  for (int k = 0; k < NZ; ++k) {
    for (int j = 0; j < NY; ++j) {
      for (int i = 0; i < NX; ++i) {

        if (obst[i][j][k] == 1) {

          for (int q = 0; q < Q; ++q) {

            int ip;
            int jp;
            int kp;

            ip = mod(i + c[q][0], NX);
            jp = mod(j + c[q][1], NY);
            kp = mod(k + c[q][2], NZ);

            if (obst[ip][jp][kp] == 0) {

              double rx = (double)i - XC;
              double ry = (double)j - YC;
              double rz = (double)k - ZC;

              double ubx = obst_omega[1] * rz - obst_omega[2] * ry;
              double uby = obst_omega[2] * rx - obst_omega[0] * rz;
              double ubz = obst_omega[0] * ry - obst_omega[1] * rx;

              double ubc = ubx * c[q][0] + uby * c[q][1] + ubz * c[q][2];

              double fx = +1.0 * (fstar[q][i][j][k] - fstar[opp[q]][ip][jp][kp] - 6 * weights[q] * 1.0 * ubc) * c[q][0];
              double fy = +1.0 * (fstar[q][i][j][k] - fstar[opp[q]][ip][jp][kp] - 6 * weights[q] * 1.0 * ubc) * c[q][1];
              double fz = +1.0 * (fstar[q][i][j][k] - fstar[opp[q]][ip][jp][kp] - 6 * weights[q] * 1.0 * ubc) * c[q][2];

              obst_torq[0] += ry * fz - rz * fy;
              obst_torq[1] += rz * fx - rx * fz;
              obst_torq[2] += rx * fy - ry * fx;

            }

          }

        }

      }
    }
  }

}

void move_obstacle(int iter,
                   int obst[NX][NY][NZ],
                   double obst_theta[3],
                   double obst_torq[3],
                   double obst_omega[3])
{

  int movement = 0;

  switch (movement) {
    case 0: {
      break;
    }
    case 1: {
      double omega = 0.2 / (NZ - 1.0);
      double A = (NZ - 1.0) / 6.0;
      for (int k = 0; k < NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
          for (int i = 0; i < NX; ++i) {

            double x = (double)i - XC;
            double y = (double)j - YC;
            double z = (double)k - ZC + A * sin(omega*iter);

            if ((x*x + z*z) <= (OBST_R*OBST_R) ) {
              obst[i][j][k] = 1;
            } else {
              obst[i][j][k] = 0;
            }
          }
        }
      }
      break;
    }
    case 2: {
      double a = 1.0 * OBST_R;
      double b = 5.0 * OBST_R;
      double rho_s = 5.e6;
      double OBST_I = rho_s * a * b * (a*a + b*b) / 12.0;

      obst_omega[0] = obst_omega[0];
      obst_omega[1] = obst_omega[1] + obst_torq[1] / OBST_I; // dt = 1
      obst_omega[2] = obst_omega[2];

      for (int k = 0; k < NZ; ++k) {
        for (int j = 0; j < NY; ++j) {
          for (int i = 0; i < NX; ++i) {

            double x = (double)i - XC;
            double y = (double)j - YC;
            double z = (double)k - ZC;

            obst_theta[1] = obst_theta[1] - obst_omega[1]; // dt = 1
            double xrot = x*cos(obst_theta[1]) - z*sin(obst_theta[1]);
            double zrot = x*sin(obst_theta[1]) + z*cos(obst_theta[1]);
            if ( (abs(xrot) <= a) && (abs(zrot) <= b) ) {
              obst[i][j][k] = 1;
            } else {
              obst[i][j][k] = 0;
            }
          }
        }
      }
      break;
    }
  }

}
