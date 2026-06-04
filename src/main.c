/*
 * D3Q19 LBM
 * 3D domain decomposition with MPI and halo cells
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "param.h"
#include "subroutines.h"

int main(void)
{
  // initialize MPI and build geometry
  init_mpi();

  // get local MPI info
  local_mpi();

  // allocate arrays
  static double rho[LNX_H][LNY_H][LNZ_H];
  static double u[LNX_H][LNY_H][LNZ_H];
  static double v[LNX_H][LNY_H][LNZ_H];
  static double w[LNX_H][LNY_H][LNZ_H];
  static double fbx[LNX_H][LNY_H][LNZ_H];
  static double fby[LNX_H][LNY_H][LNZ_H];
  static double fbz[LNX_H][LNY_H][LNZ_H];

  static double f[Q][LNX_H][LNY_H][LNZ_H];
  static double fstar[Q][LNX_H][LNY_H][LNZ_H];
  static double feq[Q][LNX_H][LNY_H][LNZ_H];
  static double fguo[Q][LNX_H][LNY_H][LNZ_H];

  static int obst[LNX_H][LNY_H][LNZ_H];

  MPI_Barrier(cart_comm);

  // initialize
  initialize(ic, rho, u, v, w);

  compute_force(rho, u, v, w, fbx, fby, fbz);

  equilibrium(rho, u, v, w, f);

  // initialize obstacle
  double obst_force[3] = {0.0, 0.0, 0.0};
  double obst_torq[3]  = {0.0, 0.0, 0.0};
  double obst_vel[3]   = {0.0, 0.0, 0.0};
  double obst_omega[3] = {0.0, 0.0, 0.0};
  double obst_theta[3] = {0.0, 0.5*PI, 0.0};

  build_obstacle(obst, obst_theta);

  halo_exchange_int(obst);

  // write initial condition
  compute_obst_drag(f, obst, obst_vel, obst_force);

  write_0d(0, 0.0, obst_force[0] / (1.0 * U0 * U0 * OBST_R * NY), "C_D");
  write_0d(0, 0.0, obst_force[2] / (1.0 * U0 * U0 * OBST_R * NY), "C_L");

  write_0d(0, 0.0, max_u(u), "u_peak");
  write_0d(0, 0.0, obst_theta[1], "theta");
  write_0d(0, 0.0, obst_omega[1], "omega");
  write_0d(0, 0.0, obst_torq[1], "torq");


  write_field_slice(0, rho, 2, NY/2, "rho_slice");
  write_field_slice(0, u,   2, NY/2, "u_slice");
  write_field_slice(0, v,   2, NY/2, "v_slice");
  write_field_slice(0, w,   2, NY/2, "w_slice");

  write_field(0, rho, "rho");
  write_field(0, u, "u");
  write_field(0, v, "v");
  write_field(0, w, "w");

  clock_t t0 = clock();

  // march time
  for (int iter = 1; iter <= MAXITER; ++iter) {

    obst_bounce(obst, obst_omega, f);

    halo_exchange(f);

    stream(f, fstar);

    boundary(fstar);

    moments(fstar, obst, fbx, fby, fbz, rho, u, v, w);

    compute_force(rho, u, v, w, fbx, fby, fbz);

    guo(fbx, fby, fbz, u, v, w, fguo);

    equilibrium(rho, u, v, w, feq);

    collide(obst, fstar, feq, fguo);

    compute_obst_drag(fstar, obst, obst_vel, obst_force);

    compute_obst_torque(obst, fstar, obst_omega, obst_torq);

    move_obstacle(iter, obst, obst_theta, obst_torq, obst_omega);

    halo_exchange_int(obst);

    memcpy(f, fstar, sizeof(double) * Q * LNX_H * LNY_H * LNZ_H);

    // write 0D time series
    write_0d(iter, iter * DT, (obst_force[0] / (1.0 * U0 * U0 * OBST_R * (NY))), "C_D");
    write_0d(iter, iter * DT, (obst_force[2] / (1.0 * U0 * U0 * OBST_R * (NY))), "C_L");
    write_0d(iter, iter * DT, max_u(u), "u_peak");
    write_0d(iter, iter * DT, obst_torq[1], "torq");
    write_0d(iter, iter * DT, obst_omega[1], "omega");
    write_0d(iter, iter * DT, obst_theta[1], "theta");

    // write fields
    if (iter % IOUT2D == 0) {

      write_field_slice(iter, rho, 2, NY/2, "rho_slice");
      write_field_slice(iter, u,   2, NY/2, "u_slice");
      write_field_slice(iter, v,   2, NY/2, "v_slice");
      write_field_slice(iter, w,   2, NY/2, "w_slice");

      double maxu = max_u(u);
      if (cart_rank == 0) {
        printf("Iteration = %d / %d\n", iter, MAXITER);
        printf("max_u = %f\n", maxu);
      }

    }

    if (iter % IOUT3D == 0) {

      write_field(iter, rho, "rho");
      write_field(iter, u, "u");
      write_field(iter, v, "v");
      write_field(iter, w, "w");

    }
  }

  clock_t t1 = clock();

  if (cart_rank == 0) {
    printf("\nFINE\n");
    printf("Elapsed time = %f s\n",
           (double)(t1 - t0) / CLOCKS_PER_SEC);
  }

  // finalise MPI
  MPI_Comm_free(&cart_comm);
  MPI_Finalize();

  return 0;

}
