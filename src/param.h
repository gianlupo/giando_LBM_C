#ifndef PARAM_H
#define PARAM_H

#include <mpi.h>

// MPI processes
// NPROC = NPROCX * NPROCY * NPROCZ
// to be set at compilation e.g. make NPROCX=4 NPROCZ=2
#ifndef NPROCX
#define NPROCX 1
#endif
#ifndef NPROCY
#define NPROCY 1
#endif
#ifndef NPROCZ
#define NPROCZ 1
#endif

// global lattice parameters
#define NX 800
#define NY 1
#define NZ 400
#define Q 19

// local sub-domain size (no halo)
#define LNX  (NX / NPROCX)
#define LNY  (NY / NPROCY)
#define LNZ  (NZ / NPROCZ)

// local sub-domain size (with halo)
#define LNX_H (LNX + 2)
#define LNY_H (LNY + 2)
#define LNZ_H (LNZ + 2)

// maximum iterations
#define MAXITER 1000

// output every IOUT iterations
#define IOUT2D 100
#define IOUT3D 10000

// lattice constants (assigned in param.c)
extern const double weights[Q];
extern const int c[Q][3];
extern const int opp[Q];
extern const int qo[6][5];
extern const int qi[6][14];
extern const int n[6][3];

// global parameters (assigned in param.c)
extern const double PI;
extern const double TAU;
extern const double U0;

extern const double LX;
extern const double LY;
extern const double LZ;
extern const double DX;
extern const double DY;
extern const double DZ;
extern const double DT;

// obstacle parameters (assigned in param.c)
extern const double OBST_XC;
extern const double OBST_YC;
extern const double OBST_ZC;
extern const double OBST_R;

// boundary and initial conditions (assigned in param.c)
extern const char bc[6];
extern const int ic;

// MPI cartesian communicator and rank constants (assigned in mpimod.c)
extern MPI_Comm cart_comm;
extern int      world_size;
extern int      cart_rank;
extern int      cart_coords[3];

extern int neighbours[6];

extern int offset_x;
extern int offset_y;
extern int offset_z;

#endif
