#ifndef SUBROUTINES_H
#define SUBROUTINES_H

#include "param.h"

// mpimod.c
void init_mpi  (void);

void local_mpi (void);


// obstacle.c
void build_obstacle      (int    obst[LNX_H][LNY_H][LNZ_H],
                          double obst_theta[3]);

void compute_obst_torque (int    obst[LNX_H][LNY_H][LNZ_H],
                          double fstar[Q][LNX_H][LNY_H][LNZ_H],
                          double obst_omega[3],
                          double obst_torq[3]);

void move_obstacle       (int    iter,
                          int    obst[LNX_H][LNY_H][LNZ_H],
                          double obst_theta[3], double obst_torq[3],
                          double obst_omega[3]);

// init.c
void initialize          (int    ic,
                          double rho[LNX_H][LNY_H][LNZ_H],
                          double u[LNX_H][LNY_H][LNZ_H],
                          double v[LNX_H][LNY_H][LNZ_H],
                          double w[LNX_H][LNY_H][LNZ_H]);

// force.c
void compute_force       (double rho[LNX_H][LNY_H][LNZ_H],
                          double u[LNX_H][LNY_H][LNZ_H],
                          double v[LNX_H][LNY_H][LNZ_H],
                          double w[LNX_H][LNY_H][LNZ_H],
                          double fbx[LNX_H][LNY_H][LNZ_H],
                          double fby[LNX_H][LNY_H][LNZ_H],
                          double fbz[LNX_H][LNY_H][LNZ_H]);

// lbm.c
void equilibrium         (double rho[LNX_H][LNY_H][LNZ_H],
                          double u[LNX_H][LNY_H][LNZ_H],
                          double v[LNX_H][LNY_H][LNZ_H],
                          double w[LNX_H][LNY_H][LNZ_H],
                          double feq[Q][LNX_H][LNY_H][LNZ_H]);

void stream              (double f[Q][LNX_H][LNY_H][LNZ_H],
                          double fstar[Q][LNX_H][LNY_H][LNZ_H]);

void moments             (double fstar[Q][LNX_H][LNY_H][LNZ_H],
                          int    obst[LNX_H][LNY_H][LNZ_H],
                          double fbx[LNX_H][LNY_H][LNZ_H],
                          double fby[LNX_H][LNY_H][LNZ_H],
                          double fbz[LNX_H][LNY_H][LNZ_H],
                          double rho[LNX_H][LNY_H][LNZ_H],
                          double u[LNX_H][LNY_H][LNZ_H],
                          double v[LNX_H][LNY_H][LNZ_H],
                          double w[LNX_H][LNY_H][LNZ_H]);

void guo                 (double fbx[LNX_H][LNY_H][LNZ_H],
                          double fby[LNX_H][LNY_H][LNZ_H],
                          double fbz[LNX_H][LNY_H][LNZ_H],
                          double u[LNX_H][LNY_H][LNZ_H],
                          double v[LNX_H][LNY_H][LNZ_H],
                          double w[LNX_H][LNY_H][LNZ_H],
                          double fguo[Q][LNX_H][LNY_H][LNZ_H]);

void collide             (int    obst [LNX_H][LNY_H][LNZ_H],
                          double fstar[Q][LNX_H][LNY_H][LNZ_H],
                          double feq  [Q][LNX_H][LNY_H][LNZ_H],
                          double fguo [Q][LNX_H][LNY_H][LNZ_H]);

// bound.c
void obst_bounce         (int    obst[LNX_H][LNY_H][LNZ_H],
                          double obst_omega[3],
                          double f[Q][LNX_H][LNY_H][LNZ_H]);

void boundary            (double fstar[Q][LNX_H][LNY_H][LNZ_H]);

// halo.c 
void halo_exchange       (double f[Q][LNX_H][LNY_H][LNZ_H]);

void halo_exchange_int   (int    a[LNX_H][LNY_H][LNZ_H]);

// output.c
void write_field         (int         iter,
                          double      field[LNX_H][LNY_H][LNZ_H],
                          const char *name);

void write_0d            (int         iter,
                          double      time,
                          double      var,
                          const char *name);

void write_field_slice   (int         iter,
                          double      slice[LNX_H][LNY_H][LNZ_H],
                          int         direction,
                          int         position,
                          const char *name);

double max_u             (double u[LNX_H][LNY_H][LNZ_H]);

void compute_obst_drag   (double fstar[Q][LNX_H][LNY_H][LNZ_H],
                          int    obst[LNX_H][LNY_H][LNZ_H],
                          double obst_vel[3],
                          double obst_force[3]);

#endif
