#include <math.h>
#include <stdlib.h>
#include "param.h"
#include "utils.h"

void build_obstacle(int    obst[LNX_H][LNY_H][LNZ_H],
                    double obst_theta[3])
{

  int shape = 2;

  // initialise to 0 (with halos)
  for (int k = 0; k < LNZ_H; ++k) {
    for (int j = 0; j < LNY_H; ++j) {
      for (int i = 0; i < LNX_H; ++i) {
        obst[i][j][k] = 0;
      }
    }
  }


  switch (shape) {

    case 0: {
      // no obstacle
      break;
    }
    case 1: {
      // cylinder
      for (int k = 1; k <= LNZ; ++k) {
        for (int j = 1; j <= LNY; ++j) {
          for (int i = 1; i <= LNX; ++i) {

            double x = (double)gi(i) - OBST_XC;
            double y = (double)gj(j) - OBST_YC;
            double z = (double)gk(k) - OBST_ZC;

            obst[i][j][k] = ((x*x + z*z) <= OBST_R*OBST_R) ? 1 : 0;

          }
        }
      }
      break;
    }
    case 2: {
      // slab
      double a = 1.0*OBST_R;
      double b = 5.0*OBST_R;
      for (int k = 1; k <= LNZ; ++k) {
        for (int j = 1; j <= LNY; ++j) {
          for (int i = 1; i <= LNX; ++i) {

            double x = (double)gi(i) - OBST_XC;
            double y = (double)gj(j) - OBST_YC;
            double z = (double)gk(k) - OBST_ZC;
            double xrot =  x*cos(obst_theta[1]) - z*sin(obst_theta[1]);
            double zrot =  x*sin(obst_theta[1]) + z*cos(obst_theta[1]);

            obst[i][j][k] = ((abs(xrot) <= a) && (abs(zrot) <= b)) ? 1 : 0;

          }
        }
      }
      break;
    }
  }

}

void compute_obst_torque(int    obst[LNX_H][LNY_H][LNZ_H],
                         double fstar[Q][LNX_H][LNY_H][LNZ_H],
                         double obst_omega[3],
                         double obst_torq [3])
{
  double local_torq[3] = {0.0, 0.0, 0.0};

  for (int k = 1; k <= LNZ; ++k) {
    for (int j = 1; j <= LNY; ++j) {
      for (int i = 1; i <= LNX; ++i) {

        if (obst[i][j][k] != 1) continue;

        for (int q = 0; q < Q; ++q) {

          int ip = i + c[q][0];
          int jp = j + c[q][1];
          int kp = k + c[q][2];

          if (obst[ip][jp][kp] != 0) continue;

          double rx = (double)gi(i) - OBST_XC;
          double ry = (double)gj(j) - OBST_YC;
          double rz = (double)gk(k) - OBST_ZC;

          double ubx = obst_omega[1] * rz - obst_omega[2] * ry;
          double uby = obst_omega[2] * rx - obst_omega[0] * rz;
          double ubz = obst_omega[0] * ry - obst_omega[1] * rx;
          double ubc = ubx * c[q][0] + uby * c[q][1] + ubz * c[q][2];

          double fx = 1.0 * (fstar[q][i][j][k] - fstar[opp[q]][ip][jp][kp] - 6.0 * weights[q] * 1.0 * ubc) * c[q][0];
          double fy = 1.0 * (fstar[q][i][j][k] - fstar[opp[q]][ip][jp][kp] - 6.0 * weights[q] * 1.0 * ubc) * c[q][1];
          double fz = 1.0 * (fstar[q][i][j][k] - fstar[opp[q]][ip][jp][kp] - 6.0 * weights[q] * 1.0 * ubc) * c[q][2];

          local_torq[0] += ry * fz - rz * fy;
          local_torq[1] += rz * fx - rx * fz;
          local_torq[2] += rx * fy - ry * fx;

        }

      }
    }
  }

  MPI_Allreduce(local_torq, obst_torq, 3, MPI_DOUBLE, MPI_SUM, cart_comm);
}

void move_obstacle(int    iter,
                   int    obst[LNX_H][LNY_H][LNZ_H],
                   double obst_theta[3],
                   double obst_torq [3],
                   double obst_omega[3])

{

  int movement = 2;

  // clear halos
  for (int k = 0; k < LNZ_H; ++k) {
    for (int j = 0; j < LNY_H; ++j) {
      for (int i = 0; i < LNX_H; ++i) {
        obst[i][j][k] = 0;
      }
    }
  }

  switch (movement) {
    case 0: {
      // no movement
      break;
    }
    case 1: {
      // imposed movement
      double omega = 0.2 / (NZ - 1.0);
      double A = (NZ - 1.0) / 6.0;
      for (int k = 1; k <= LNZ; ++k) {
        for (int j = 1; j <= LNY; ++j) {
          for (int i = 1; i <= LNX; ++i) {

            double x = (double)gi(i) - OBST_XC;
            double y = (double)gj(j) - OBST_YC;
            double z = (double)gk(k) - OBST_ZC + A * sin(omega*iter);

            obst[i][j][k] = ((x*x + z*z) <= OBST_R*OBST_R) ? 1 : 0;

          }
        }
      }
      break;
    }
    case 2: {
      // coupled movement
      double a = 1.0 * OBST_R;
      double b = 5.0 * OBST_R;
      double rho_s = 50.0;
      double OBST_I = rho_s * a * b * (a*a + b*b) / 12.0;

      obst_omega[1] += obst_torq[1] / OBST_I; // dt = 1
      obst_theta[1] += obst_omega[1];         // dt = 1

      for (int k = 1; k <= LNZ; ++k) {
        for (int j = 1; j <= LNY; ++j) {
          for (int i = 1; i <= LNX; ++i) {

            double x = (double)gi(i) - OBST_XC;
            double y = (double)gj(j) - OBST_YC;
            double z = (double)gk(k) - OBST_ZC;
            double xrot = x*cos(obst_theta[1]) - z*sin(obst_theta[1]);
            double zrot = x*sin(obst_theta[1]) + z*cos(obst_theta[1]);

            obst[i][j][k] = ((abs(xrot) <= a) && (abs(zrot) <= b)) ? 1 : 0;

          }
        }
      }
      break;
    }
  }

}
