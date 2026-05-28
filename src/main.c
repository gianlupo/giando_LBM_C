/*
 * D3Q19 LBM
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "param.h"
#include "subroutines.h"

int main(void)
{
  static double rho[NX][NY][NZ];
  static double u[NX][NY][NZ];
  static double v[NX][NY][NZ];
  static double w[NX][NY][NZ];
  static double fbx[NX][NY][NZ];
  static double fby[NX][NY][NZ];
  static double fbz[NX][NY][NZ];

  static double Kin[NX][NY][NZ];

  static double f[Q][NX][NY][NZ];
  static double fstar[Q][NX][NY][NZ];
  static double feq[Q][NX][NY][NZ];
  static double fguo[Q][NX][NY][NZ];

  static int obst[NX][NY][NZ];

  static double obst_force[3];
  static double obst_torq[3];
  static double obst_vel[3];
  static double obst_omega[3];
  static double obst_theta[3];

  (void)system("mkdir -p data");

  // build obstacle geometry

  build_obstacle(obst);

  // initialize

  initialize(ic, rho, u, v, w);

  compute_force(rho, u, v, w, fbx, fby, fbz);

  equilibrium(rho, u, v, w, f);

  obst_force[0] = 0.0;
  obst_force[1] = 0.0;
  obst_force[2] = 0.0;

  obst_torq[0] = 0.0;
  obst_torq[1] = 0.0;
  obst_torq[2] = 0.0;

  obst_vel[0] = 0.0;
  obst_vel[1] = 0.0;
  obst_vel[2] = 0.0;

  obst_omega[0] = 0.0;
  obst_omega[1] = 0.0;
  obst_omega[2] = 0.0;

  obst_theta[0] = 0.0;
  obst_theta[1] = 0.0;
  obst_theta[2] = 0.0;

  // write initial condition

  compute_kinetic(rho, u, v, w, Kin);

  write_0d(0, 0.0,
           max_u(u),
           "u_peak");

  write_0d(0, 0.0,
           max_kinetic(rho, u, v, w),
           "K_peak");

  write_field_slice(0, rho, 2, NY/2, "rho_slice");
  write_field_slice(0, u,   2, NY/2, "u_slice");
  write_field_slice(0, v,   2, NY/2, "v_slice");
  write_field_slice(0, w,   2, NY/2, "w_slice");
  write_field_slice(0, Kin, 2, NY/2, "K_slice");

  write_field(0, rho, "rho");
  write_field(0, u, "u");
  write_field(0, v, "v");
  write_field(0, w, "w");
  write_field(0, Kin, "K");

  clock_t t0 = clock();

  // march time

  for (int iter = 1; iter <= MAXITER; ++iter) {

    obst_bounce(obst, f);

    stream(f, fstar);

    boundary(fstar);

    moments(fstar, obst, fbx, fby, fbz, rho, u, v, w);

    compute_force(rho, u, v, w, fbx, fby, fbz);

    guo(fbx, fby, fbz, u, v, w, fguo);

    equilibrium(rho, u, v, w, feq);

    collide(obst, fstar, feq, fguo);

    if (iter % 1 == 0) {

      compute_obst_drag(fstar, obst, obst_vel, obst_force);

      write_0d(iter,
               iter * DT_PHYS,
               (obst_force[0] / (1.0 * U0 * U0 * OBST_R * (NY))),
               "C_D");

      write_0d(iter,
               iter * DT_PHYS,
               (obst_force[2] / (1.0 * U0 * U0 * OBST_R * (NY))),
               "C_L");
    }


    compute_obst_torque(obst, fstar, obst_omega, obst_torq);

    move_obstacle(iter, obst, obst_theta, obst_torq, obst_omega);

    memcpy(f, fstar, sizeof(double) * Q * NX * NY * NZ);

    // write time series

    write_0d(iter,
             iter * DT_PHYS,
             max_u(u),
             "u_peak");

    write_0d(iter,
             iter * DT_PHYS,
             max_kinetic(rho, u, v, w),
             "K_peak");

    write_0d(iter,
             iter * DT_PHYS,
             obst_theta[1],
             "theta");

    // write fields

    if (iter % IOUT2D == 0) {

      compute_kinetic(rho, u, v, w, Kin);

      write_field_slice(iter, rho, 2, NY/2, "rho_slice");
      write_field_slice(iter, u,   2, NY/2, "u_slice");
      write_field_slice(iter, v,   2, NY/2, "v_slice");
      write_field_slice(iter, w,   2, NY/2, "w_slice");
      write_field_slice(iter, Kin, 2, NY/2, "K_slice");

      printf("Iteration = %d / %d\n",
             iter,
             MAXITER);
      printf("%f\n",max_u(u));
    }

    if (iter % IOUT3D == 0) {

      compute_kinetic(rho, u, v, w, Kin);

      write_field(iter, rho, "rho");
      write_field(iter, u, "u");
      write_field(iter, v, "v");
      write_field(iter, w, "w");
      write_field(iter, Kin, "K");

    }
  }

  clock_t t1 = clock();

  printf("\n");
  printf("FINE");
  printf("\n");

  printf("Elapsed time = %f s\n",
         (double)(t1 - t0) / CLOCKS_PER_SEC);

  return 0;
}
