#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "param.h"

void initialize(double rho[NX][NY],
                double u[NX][NY],
                double v[NX][NY]);

void equilibrium(double rho[NX][NY],
                 double u[NX][NY],
                 double v[NX][NY],
                 double feq[Q][NX][NY]);

void stream(double f[Q][NX][NY],
            double fstar[Q][NX][NY]);

void moments(double fstar[Q][NX][NY],
             double rho[NX][NY],
             double u[NX][NY],
             double v[NX][NY]);

void collide(double fstar[Q][NX][NY],
             double feq[Q][NX][NY],
             double f[Q][NX][NY]);

void periodic_bc(double f[Q][NX][NY],
                 double fstar[Q][NX][NY]);

void write_field(int iter,
                 double field[NX][NY],
                 const char* name);

void write_0d(int iter,
              double time,
              double var,
              const char* name);

double max_u(double u[NX][NY]);

double max_kinetic(double rho[NX][NY],
                   double u[NX][NY],
                   double v[NX][NY]);

void compute_kinetic(double rho[NX][NY],
                     double u[NX][NY],
                     double v[NX][NY],
                     double K[NX][NY]);

int main(void)
{
    static double rho[NX][NY];
    static double u[NX][NY];
    static double v[NX][NY];

    static double K[NX][NY];

    static double f[Q][NX][NY];
    static double fstar[Q][NX][NY];
    static double feq[Q][NX][NY];

    (void)system("mkdir -p data");

    initialize(rho, u, v);

    equilibrium(rho, u, v, f);

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

    for (int iter = 1; iter <= MAXITER; ++iter) {

        periodic_bc(f, fstar);

        stream(f, fstar);

        moments(fstar, rho, u, v);

        equilibrium(rho, u, v, feq);

        collide(fstar, feq, f);

        write_0d(iter,
                 iter * DT_PHYS,
                 max_u(u) * U0_PHYS,
                 "u_peak");

        write_0d(iter,
                 iter * DT_PHYS,
                 max_kinetic(rho, u, v),
                 "K_peak");

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

    printf("Elapsed time = %f s\n",
           (double)(t1 - t0) / CLOCKS_PER_SEC);

    return 0;
}
