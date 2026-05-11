/*
 * D2Q9 LBM
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "param.h"
#include "subroutines.h"

int main(void)
{
  static double rho[NX][NY];
  static double u[NX][NY];
  static double v[NX][NY];
  static double fx[NX][NY];
  static double fy[NX][NY];

  static double K[NX][NY];

  static double f[Q][NX][NY];
  static double fstar[Q][NX][NY];
  static double feq[Q][NX][NY];
  static double fguo[Q][NX][NY];

  (void)system("mkdir -p data");

  // initialize

  initialize(rho, u, v);

  compute_force(u, v, fx, fy);

  equilibrium(rho, u, v, f);

  // write initial condition

  compute_kinetic(rho, u, v, K);

  write_0d(0, 0.0,
           max_u(u) * U0_PHYS,
           "u_peak");

  write_0d(0, 0.0,
           max_kinetic(rho, u, v),
           "K_peak");

  write_field(0, rho, "rho");
  write_field(0, u, "u");
  write_field(0, v, "v");
  write_field(0, K, "K");

  clock_t t0 = clock();

  // march time

  for (int iter = 1; iter <= MAXITER; ++iter) {

    periodic_bc(f, fstar);

    stream(f, fstar);

    moments(fstar, fx, fy, rho, u, v);

    compute_force(u, v, fx, fy);

    guo(fx, fy, u, v, fguo);

    equilibrium(rho, u, v, feq);

    collide(fstar, feq, fguo, f);

    // write time series

    write_0d(iter,
             iter * DT_PHYS,
             max_u(u) * U0_PHYS,
             "u_peak");

    write_0d(iter,
             iter * DT_PHYS,
             max_kinetic(rho, u, v),
             "K_peak");

    // write fields

    if (iter % IOUT == 0) {

      compute_kinetic(rho, u, v, K);

      write_field(iter, rho, "rho");
      write_field(iter, u, "u");
      write_field(iter, v, "v");
      write_field(iter, K, "K");

      printf("Iteration = %d / %d\n",
             iter,
             MAXITER);
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
