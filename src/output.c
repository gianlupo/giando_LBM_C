#include <stdio.h>
#include "param.h"

void write_field(int iter,
                 double field[NX][NY],
                 const char* name)
{
  char filename[256];

  sprintf(filename, "data/%s_%06d.dat", name, iter);

  FILE* fp = fopen(filename, "w");

  if (!fp) {
    return;
  }

  for (int j = 0; j < NY; ++j) {
    for (int i = 0; i < NX; ++i) {

      fprintf(fp,
              "%d %d %.10f\n",
              i + 1,
              j + 1,
              field[i][j]);
    }
  }

  fclose(fp);
}

void write_0d(int iter,
              double time,
              double var,
              const char* name)
{
  char filename[256];

  sprintf(filename, "data/%s.time", name);

  FILE* fp = fopen(filename, "a");

  if (!fp) {
    return;
  }

  fprintf(fp,
          "%d %.10f %.10f\n",
          iter,
          time,
          var);

  fclose(fp);
}

double max_u(double u[NX][NY])
{
  double m = u[0][0];

  for (int j = 0; j < NY; ++j) {
    for (int i = 0; i < NX; ++i) {

      if (u[i][j] > m) {
        m = u[i][j];
      }
    }
  }

  return m;
}

double max_kinetic(double rho[NX][NY],
                   double u[NX][NY],
                   double v[NX][NY])
{
  double m = 0.0;

  for (int j = 0; j < NY; ++j) {
    for (int i = 0; i < NX; ++i) {

      double k = 0.5 * rho[i][j] * (u[i][j] * u[i][j] + v[i][j] * v[i][j]);

      if (k > m) {
        m = k;
      }
    }
  }

  return m * U0_PHYS * U0_PHYS;
}

void compute_kinetic(double rho[NX][NY],
                     double u[NX][NY],
                     double v[NX][NY],
                     double K[NX][NY])
{
  for (int j = 0; j < NY; ++j) {
    for (int i = 0; i < NX; ++i) {

      K[i][j] = 0.5 * rho[i][j] * (u[i][j] * u[i][j] + v[i][j] * v[i][j]) * U0_PHYS * U0_PHYS;
    }
  }
}
