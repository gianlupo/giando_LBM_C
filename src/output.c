#include <stdio.h>
#include <stdlib.h>
#include "param.h"
#include "utils.h"

void write_field(int iter,
                 double field[LNX_H][LNY_H][LNZ_H],
                 const char *name)
{
  // pack the local field into a flat buffer (no halos)
  int local_size = LNX * LNY * LNZ;
  double *send_buf = (double *)malloc(local_size * sizeof(double));
  if (!send_buf) return;

  int idx = 0;
  for (int i = 1; i <= LNX; ++i) {
    for (int j = 1; j <= LNY; ++j) {
      for (int k = 1; k <= LNZ; ++k) {
        send_buf[idx++] = field[i][j][k];
      }
    }
  }

  // gather on rank 0
  int nproc;
  MPI_Comm_size(cart_comm, &nproc);

  double *recv_buf    = NULL;
  int    *recv_counts = NULL;
  int    *displs      = NULL;
  int    *all_offsets = NULL;

  if (cart_rank == 0) {
    recv_buf    = (double *)malloc(NX * NY * NZ * sizeof(double));
    recv_counts = (int *)malloc(nproc * sizeof(int));
    displs      = (int *)malloc(nproc * sizeof(int));
    all_offsets = (int *)malloc(nproc * 3 * sizeof(int));
  }

  // gather recv counts
  MPI_Gather(&local_size, 1, MPI_INT, recv_counts, 1, MPI_INT, 0, cart_comm);

  if (cart_rank == 0) {
    displs[0] = 0;
    for (int r = 1; r < nproc; ++r) {
      displs[r] = displs[r-1] + recv_counts[r-1];
    }
  }

  // gather data
  MPI_Gatherv(send_buf, local_size, MPI_DOUBLE, recv_buf, recv_counts, displs, MPI_DOUBLE, 0, cart_comm);

  // gather offsets
  int my_off[3] = {offset_x, offset_y, offset_z};
  MPI_Gather(my_off, 3, MPI_INT, all_offsets, 3, MPI_INT, 0, cart_comm);
  
  // write global field
  if (cart_rank == 0) {
    // global array
    double *global = recv_buf;
    double *gfield = (double *)malloc(NX * NY * NZ * sizeof(double));
    int base = 0;
    for (int r = 0; r < nproc; ++r) {
      int ox = all_offsets[3*r + 0];
      int oy = all_offsets[3*r + 1];
      int oz = all_offsets[3*r + 2];
      int ri = 0;
      for (int i = 0; i < LNX; ++i) {
        for (int j = 0; j < LNY; ++j) {
          for (int k = 0; k < LNZ; ++k) {
            gfield[(ox+i)*NY*NZ + (oy+j)*NZ + (oz+k)] = recv_buf[displs[r] + ri++];
          }
        }
      }
    }

    char filename[256];
    sprintf(filename, "data/%s_%06d.dat", name, iter);
    FILE *fp = fopen(filename, "w");

    if (fp) {
      for (int i = 0; i < NX; ++i) {
        for (int j = 0; j < NY; ++j) {
          for (int k = 0; k < NZ; ++k) {

            fprintf(fp,
                    "%d %d %d %.10f\n",
                    i + 1,
                    j + 1,
                    k + 1,
                    gfield[i*NY*NZ + j*NZ + k]);
          }
        }
      }
      fclose(fp);
    }
    free(gfield);
    free(recv_buf);
    free(recv_counts);
    free(displs);
    free(all_offsets);
  }
  free(send_buf);
}

void write_field_slice(int iter,
                       double slice[LNX_H][LNY_H][LNZ_H],
                       int direction,
                       int position,
                       const char *name)
{
  // pack the local field into a flat buffer (no halos)
  int local_size = LNX * LNY * LNZ;
  double *send_buf = (double *)malloc(local_size * sizeof(double));
  if (!send_buf) return;

  int idx = 0;
  for (int i = 1; i <= LNX; ++i) {
    for (int j = 1; j <= LNY; ++j) {
      for (int k = 1; k <= LNZ; ++k) {
        send_buf[idx++] = slice[i][j][k];
      }
    }
  }

  // gather on rank 0
  int nproc;
  MPI_Comm_size(cart_comm, &nproc);

  double *recv_buf    = NULL;
  int    *recv_counts = NULL;
  int    *displs      = NULL;
  int    *all_offsets = NULL;

  if (cart_rank == 0) {
    recv_buf    = (double *)malloc(NX * NY * NZ * sizeof(double));
    recv_counts = (int *)malloc(nproc * sizeof(int));
    displs      = (int *)malloc(nproc * sizeof(int));
    all_offsets = (int *)malloc(nproc * 3 * sizeof(int));
  }

  // gather recv counts
  MPI_Gather(&local_size, 1, MPI_INT, recv_counts, 1, MPI_INT, 0, cart_comm);

  if (cart_rank == 0) {
    displs[0] = 0;
    for (int r = 1; r < nproc; ++r) {
      displs[r] = displs[r-1] + recv_counts[r-1];
    }
  }

  // gather data
  MPI_Gatherv(send_buf, local_size, MPI_DOUBLE, recv_buf, recv_counts, displs, MPI_DOUBLE, 0, cart_comm);

  // gather offsets
  int my_off[3] = {offset_x, offset_y, offset_z};
  MPI_Gather(my_off, 3, MPI_INT, all_offsets, 3, MPI_INT, 0, cart_comm);

  // write global slice
  if (cart_rank == 0) {
    // global array
    double *gfield = (double *)malloc(NX * NY * NZ * sizeof(double));
    for (int r = 0; r < nproc; ++r) {
      int ox = all_offsets[3*r + 0];
      int oy = all_offsets[3*r + 1];
      int oz = all_offsets[3*r + 2];
      int ri = 0;
      for (int i = 0; i < LNX; ++i) {
        for (int j = 0; j < LNY; ++j) {
          for (int k = 0; k < LNZ; ++k) {
            gfield[(ox+i)*NY*NZ + (oy+j)*NZ + (oz+k)] = recv_buf[displs[r] + ri++];
          }
        }
      }
    }

    char filename[256];
    sprintf(filename, "data/%s_%06d.dat", name, iter);
    FILE *fp = fopen(filename, "w");

    if (fp) {
      int imin, imax, jmin, jmax, kmin, kmax;
      imin = 0; imax = NX;
      jmin = 0; jmax = NY;
      kmin = 0; kmax = NZ;
      switch (direction) {
        // direction == 1 slice perpendicular to x
        case 1:
          imin = position;
          imax = position + 1;
          break;
        // direction == 2 slice perpendicular to y
        case 2:
          jmin = position;
          jmax = position + 1;
          break;
        // direction == 3 slice perpendicular to z
        case 3:
          kmin = position;
          kmax = position + 1;
          break;
        // default to direction == 2
        default:
          jmin = position;
          jmax = position + 1;
          break;
      }

      for (int i = imin; i < imax; ++i) {
        for (int j = jmin; j < jmax; ++j) {
          for (int k = kmin; k < kmax; ++k) {

            fprintf(fp,
                    "%d %d %d %.10f\n",
                    i + 1,
                    j + 1,
                    k + 1,
                    gfield[i*NY*NZ + j*NZ + k]);
          }
        }
      }
      fclose(fp);
    }
    free(gfield);
    free(recv_buf);
    free(recv_counts);
    free(displs);
    free(all_offsets);
  }
  free(send_buf);
}

void write_0d(int iter,
              double time,
              double var,
              const char* name)
{
  if (cart_rank != 0) return;

  char filename[256];
  sprintf(filename, "data/%s.time", name);
  FILE* fp = fopen(filename, "a");
  if (!fp) return;
  fprintf(fp, "%d %.10f %.10f\n", iter, time, var);
  fclose(fp);

}

double max_u(double u[LNX_H][LNY_H][LNZ_H])
{
  double local_m = u[1][1][1];

  for (int k = 1; k <= LNZ; ++k) {
    for (int j = 1; j <= LNY; ++j) {
      for (int i = 1; i <= LNX; ++i) {

        if (u[i][j][k] > local_m) {
          local_m = u[i][j][k];
        }
      }
    }
  }

  double global_m;
  MPI_Allreduce(&local_m, &global_m, 1, MPI_DOUBLE, MPI_MAX, cart_comm);
  return global_m;
}

void compute_obst_drag(double fstar[Q][LNX_H][LNY_H][LNZ_H],
                       int    obst[LNX_H][LNY_H][LNZ_H],
                       double obst_vel[3],
                       double obst_force[3])
{

  double local_obst_force[3] = {0.0, 0.0, 0.0};

  for (int k = 1; k <= LNZ; ++k) {
    for (int j = 1; j <= LNY; ++j) {
      for (int i = 1; i <= LNX; ++i) {

        if (obst[i][j][k] != 1) continue;

        for (int q = 0; q < Q; ++q) {

          int ip = i + c[q][0];
          int jp = j + c[q][1];
          int kp = k + c[q][2];

          if (obst[ip][jp][kp] == 0) {

            // Ladd Verberg link bounce back

            local_obst_force[0] += 1.0 * (fstar[q][i][j][k] - fstar[opp[q]][ip][jp][kp])* c[q][0];
            local_obst_force[1] += 1.0 * (fstar[q][i][j][k] - fstar[opp[q]][ip][jp][kp])* c[q][1];
            local_obst_force[2] += 1.0 * (fstar[q][i][j][k] - fstar[opp[q]][ip][jp][kp])* c[q][2];

          }
        }

      }
    }
  }

  MPI_Allreduce(local_obst_force, obst_force, 3, MPI_DOUBLE, MPI_SUM, cart_comm);
}
