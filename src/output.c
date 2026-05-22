#include <stdio.h>
#include "param.h"

void write_field(int iter,
                 double field[NX][NY][NZ],
                 const char* name)
{
  char filename[256];

  sprintf(filename, "data/%s_%06d.dat", name, iter);

  FILE* fp = fopen(filename, "w");

  if (!fp) {
    return;
  }

  for (int k = 0; k < NZ; ++k) {
    for (int j = 0; j < NY; ++j) {
      for (int i = 0; i < NX; ++i) {

        fprintf(fp,
                "%d %d %d %.10f\n",
                i + 1,
                j + 1,
                k + 1,
                field[i][j][k]);
      }
    }
  }

  fclose(fp);
}

void write_field_slice(int iter,
                       double slice[NX][NY][NZ],
                       int direction,
                       int position,
                       const char* name)
{
  char filename[256];

  sprintf(filename, "data/%s_%06d.dat", name, iter);

  FILE* fp = fopen(filename, "w");

  if (!fp) {
    return;
  }

  int imin;
  int imax;
  int jmin;
  int jmax;
  int kmin;
  int kmax;

  switch (direction) {
    // direction == 1 slice perpendicular to x
    case 1:
      imin = position;
      imax = position + 1;
      jmin = 0;
      jmax = NY;
      kmin = 0;
      kmax = NZ;
      break;
    // direction == 2 slice perpendicular to y
    case 2:
      imin = 0;
      imax = NX;
      jmin = position;
      jmax = position + 1;
      kmin = 0;
      kmax = NZ;
      break;
    // direction == 3 slice perpendicular to z
    case 3:
      imin = 0;
      imax = NX;
      jmin = 0;
      jmax = NY;
      kmin = position;
      kmax = position + 1;
      break;
    // default to direction == 2
    default:
      imin = 0;
      imax = NX;
      jmin = position;
      jmax = position + 1;
      kmin = 0;
      kmax = NZ;
      break;
  }

  for (int k = kmin; k < kmax; ++k) {
    for (int j = jmin; j < jmax; ++j) {
      for (int i = imin; i < imax; ++i) {

        fprintf(fp,
                "%d %d %d %.10f\n",
                i + 1,
                j + 1,
                k + 1,
                slice[i][j][k]);
      }
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

double max_u(double u[NX][NY][NZ])
{
  double m = u[0][0][0];

  for (int k = 0; k < NZ; ++k) {
    for (int j = 0; j < NY; ++j) {
      for (int i = 0; i < NX; ++i) {

        if (u[i][j][k] > m) {
          m = u[i][j][k];
        }
      }
    }
  }

  return m;
}

double max_kinetic(double rho[NX][NY][NZ],
                   double u[NX][NY][NZ],
                   double v[NX][NY][NZ],
                   double w[NX][NY][NZ])
{
  double m = 0.0;

  for (int k = 0; k < NZ; ++k) {
    for (int j = 0; j < NY; ++j) {
      for (int i = 0; i < NX; ++i) {

        double max_kin = 0.5 * rho[i][j][k] * (u[i][j][k] * u[i][j][k] +
                                               v[i][j][k] * v[i][j][k] +
                                               w[i][j][k] * w[i][j][k]);

        if (max_kin > m) {
          m = max_kin;
        }
      }
    }
  }

  return m;
}

void compute_kinetic(double rho[NX][NY][NZ],
                     double u[NX][NY][NZ],
                     double v[NX][NY][NZ],
                     double w[NX][NY][NZ],
                     double Kin[NX][NY][NZ])
{
  for (int k = 0; k < NZ; ++k) {
    for (int j = 0; j < NY; ++j) {
      for (int i = 0; i < NX; ++i) {

        Kin[i][j][k] = 0.5 * rho[i][j][k] * (u[i][j][k] * u[i][j][k] + 
                                             v[i][j][k] * v[i][j][k] +
                                             w[i][j][k] * w[i][j][k]);
      }
    }
  }
}
