#include <stdio.h>
#include "param.h"

void init_mpi(void) {

  MPI_Init(NULL, NULL);

  int world_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // sanity check
  if (world_size != NPROCX * NPROCY * NPROCZ) {
    if (world_rank == 0) {
      fprintf(stderr,
              "Error: launched with %d ranks but NPROCX*NPROCY*NPROCZ = %d\n",
              world_size, NPROCX * NPROCY * NPROCZ);
    }
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // build 3-D Cartesian communicator
  int dims[3]    = {NPROCX, NPROCY, NPROCZ};
  int periods[3];

  // periodic in all directions, overridden later by BCs
  periods[0] = 1;
  periods[1] = 1;
  periods[2] = 1;

  MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 1, &cart_comm);

}

void local_mpi(void) {

  MPI_Comm_rank(cart_comm, &cart_rank);
  MPI_Cart_coords(cart_comm, cart_rank, 3, cart_coords);

  // ranks of neighbour
  MPI_Cart_shift(cart_comm, 0, 1, &neighbours[0], &neighbours[1]);
  MPI_Cart_shift(cart_comm, 1, 1, &neighbours[2], &neighbours[3]);
  MPI_Cart_shift(cart_comm, 2, 1, &neighbours[4], &neighbours[5]);

  // global index offsets for this rank's sub-domain
  offset_x = cart_coords[0] * LNX;
  offset_y = cart_coords[1] * LNY;
  offset_z = cart_coords[2] * LNZ;

  if (cart_rank == 0) {
    printf("D3Q19 LBM – MPI parallel (%d×%d×%d = %d ranks)\n",
           NPROCX, NPROCY, NPROCZ, world_size);
    printf("Global domain: %d × %d × %d\n", NX, NY, NZ);
    printf("Local  domain: %d × %d × %d  (+ halo)\n", LNX, LNY, LNZ);
  }

}
