#ifndef PARAM_H
#define PARAM_H

// lattice parameters
#define NX 200
#define NY 100
#define Q 9

// maximum iterations
#define MAXITER 5000

// output every IOUT iterations
#define IOUT 100

extern const double PI;
extern const double TAU;
extern const double U0;

extern const double LX_PHYS;
extern const double LY_PHYS;
extern const double DX_PHYS;
extern const double DY_PHYS;
extern const double U0_PHYS;
extern const double DT_PHYS;

extern const double w[Q];
extern const int c[Q][2];

#endif
