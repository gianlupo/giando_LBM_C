#ifndef PARAM_H
#define PARAM_H

// lattice parameters
#define NX 200
#define NY 4
#define NZ 100
#define Q 19

// maximum iterations
#define MAXITER 2000

// output every IOUT iterations
#define IOUT2D 100
#define IOUT3D 10000

extern const double PI;
extern const double TAU;
extern const double U0;

extern const double LX_PHYS;
extern const double LY_PHYS;
extern const double LZ_PHYS;
extern const double DX_PHYS;
extern const double DY_PHYS;
extern const double DZ_PHYS;
extern const double U0_PHYS;
extern const double DT_PHYS;

extern const double weights[Q];
extern const int c[Q][3];
extern const int opp[Q];
extern const int qo[6][5];
extern const int qi[6][14];
extern const int n[6][3];

extern const char bc[6];

extern const int ic;

#endif
