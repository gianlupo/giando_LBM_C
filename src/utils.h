#include <math.h>
#include "param.h"

// global indices (e.g. in x: from 0 to NX - 1, no halos)
static inline int gi(int i) { return offset_x + (i - 1); }
static inline int gj(int j) { return offset_y + (j - 1); }
static inline int gk(int k) { return offset_z + (k - 1); }

static int mod(int a, int b)
{
  int r = a % b;
  return (r < 0) ? r + b : r;
}

static void poi(int i,
                int j,
                int k,
                double *upoi,
                double *vpoi,
                double *wpoi)
{

  double z = (double)gk(k) / LZ;

  *upoi = 4.0 * U0 * z * (1.0 - z);
  *vpoi = 0.0;
  *wpoi = 0.0;
}


static void sinvel(int i,
                   int j,
                   int k,
                   double *usin,
                   double *vsin,
                   double *wsin)
{

  double lambda = 1.0 * LZ;
  double kappa = 2.0 * PI / lambda;
  double z = (double)gk(k);

  *usin = U0 * sin(kappa * z);
  *vsin = 0.0;
  *wsin = 0.0;
}
