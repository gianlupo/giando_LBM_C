#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "param.h"

// buffer sizes for each pass
// pass 1 (x-normal): Q * LNY     * LNZ      per face
// pass 2 (y-normal): Q * LNX_H   * LNZ      per face  (LNX_H = LNX+2)
// pass 3 (z-normal): Q * LNX_H   * LNY_H    per face  (LNY_H = LNY+2)
#define XBUF  (Q * LNY   * LNZ  )
#define YBUF  (Q * LNX_H * LNZ  )
#define ZBUF  (Q * LNX_H * LNY_H)

// largest buffer needed across all three passes
#define MAXBUF (ZBUF > YBUF ? (ZBUF > XBUF ? ZBUF : XBUF) : (YBUF > XBUF ? YBUF : XBUF))

static double sbuf_lo[MAXBUF];
static double sbuf_hi[MAXBUF];
static double rbuf_lo[MAXBUF];
static double rbuf_hi[MAXBUF];

// same logic for halo exchange of integer fields
#define XIBUF  (LNY   * LNZ  )
#define YIBUF  (LNX_H * LNZ  )
#define ZIBUF  (LNX_H * LNY_H)
#define MAXIBUF (ZIBUF > YIBUF ? (ZIBUF > XIBUF ? ZIBUF : XIBUF) : (YIBUF > XIBUF ? YIBUF : XIBUF))

static int isbuf_lo[MAXIBUF];
static int isbuf_hi[MAXIBUF];
static int irbuf_lo[MAXIBUF];
static int irbuf_hi[MAXIBUF];

void halo_exchange(double f[Q][LNX_H][LNY_H][LNZ_H])
{
  MPI_Request req[4];
  int n;

  // pass 1: x-normal faces: j=1..LNY, k=1..LNZ
  {
    n = 0;
    for (int q = 0; q < Q; ++q) {
      for (int j = 1; j <= LNY; ++j) {
        for (int k = 1; k <= LNZ; ++k) {
          sbuf_lo[n++] = f[q][1][j][k];
        }
      }
    }

    n = 0;
    for (int q = 0; q < Q; ++q) {
      for (int j = 1; j <= LNY; ++j) {
        for (int k = 1; k <= LNZ; ++k) {
          sbuf_hi[n++] = f[q][LNX][j][k];
        }
      }
    }

    MPI_Isend(sbuf_lo, XBUF, MPI_DOUBLE, neighbours[0], 10, cart_comm, &req[0]);
    MPI_Irecv(rbuf_lo, XBUF, MPI_DOUBLE, neighbours[0], 11, cart_comm, &req[1]);
    MPI_Isend(sbuf_hi, XBUF, MPI_DOUBLE, neighbours[1], 11, cart_comm, &req[2]);
    MPI_Irecv(rbuf_hi, XBUF, MPI_DOUBLE, neighbours[1], 10, cart_comm, &req[3]);
    MPI_Waitall(4, req, MPI_STATUSES_IGNORE);

    n = 0;
    for (int q = 0; q < Q; ++q) {
      for (int j = 1; j <= LNY; ++j) {
        for (int k = 1; k <= LNZ; ++k) {
          f[q][0][j][k] = rbuf_lo[n++];
        }
      }
    }

    n = 0;
    for (int q = 0; q < Q; ++q) {
      for (int j = 1; j <= LNY; ++j) {
        for (int k = 1; k <= LNZ; ++k) {
          f[q][LNX+1][j][k] = rbuf_hi[n++];
        }
      }
    }
  }

  // pass 2: y-normal faces: i=0..LNX+1, k=1..LNZ
  {
    n = 0;
    for (int q = 0; q < Q; ++q) {
      for (int i = 0; i <= LNX+1; ++i) {
        for (int k = 1; k <= LNZ; ++k) {
          sbuf_lo[n++] = f[q][i][1][k];
        }
      }
    }

    n = 0;
    for (int q = 0; q < Q; ++q) {
      for (int i = 0; i <= LNX+1; ++i) {
        for (int k = 1; k <= LNZ; ++k) {
          sbuf_hi[n++] = f[q][i][LNY][k];
        }
      }
    }

    MPI_Isend(sbuf_lo, YBUF, MPI_DOUBLE, neighbours[2], 20, cart_comm, &req[0]);
    MPI_Irecv(rbuf_lo, YBUF, MPI_DOUBLE, neighbours[2], 21, cart_comm, &req[1]);
    MPI_Isend(sbuf_hi, YBUF, MPI_DOUBLE, neighbours[3], 21, cart_comm, &req[2]);
    MPI_Irecv(rbuf_hi, YBUF, MPI_DOUBLE, neighbours[3], 20, cart_comm, &req[3]);
    MPI_Waitall(4, req, MPI_STATUSES_IGNORE);

    n = 0;
    for (int q = 0; q < Q; ++q) {
      for (int i = 0; i <= LNX+1; ++i) {
        for (int k = 1; k <= LNZ; ++k) {
          f[q][i][0][k] = rbuf_lo[n++];
        }
      }
    }

    n = 0;
    for (int q = 0; q < Q; ++q) {
      for (int i = 0; i <= LNX+1; ++i) {
        for (int k = 1; k <= LNZ; ++k) {
          f[q][i][LNY+1][k] = rbuf_hi[n++];
        }
      }
    }
  }

  // pass 3: z-normal faces: i=0..LNX+1, j=0..LNY+1
  {
    n = 0;
    for (int q = 0; q < Q; ++q) {
      for (int i = 0; i <= LNX+1; ++i) {
        for (int j = 0; j <= LNY+1; ++j) {
          sbuf_lo[n++] = f[q][i][j][1];
        }
      }
    }

    n = 0;
    for (int q = 0; q < Q; ++q) {
      for (int i = 0; i <= LNX+1; ++i) {
        for (int j = 0; j <= LNY+1; ++j) {
          sbuf_hi[n++] = f[q][i][j][LNZ];
        }
      }
    }

    MPI_Isend(sbuf_lo, ZBUF, MPI_DOUBLE, neighbours[4], 30, cart_comm, &req[0]);
    MPI_Irecv(rbuf_lo, ZBUF, MPI_DOUBLE, neighbours[4], 31, cart_comm, &req[1]);
    MPI_Isend(sbuf_hi, ZBUF, MPI_DOUBLE, neighbours[5], 31, cart_comm, &req[2]);
    MPI_Irecv(rbuf_hi, ZBUF, MPI_DOUBLE, neighbours[5], 30, cart_comm, &req[3]);
    MPI_Waitall(4, req, MPI_STATUSES_IGNORE);

    n = 0;
    for (int q = 0; q < Q; ++q) {
      for (int i = 0; i <= LNX+1; ++i) {
        for (int j = 0; j <= LNY+1; ++j) {
          f[q][i][j][0] = rbuf_lo[n++];
        }
      }
    }

    n = 0;
    for (int q = 0; q < Q; ++q) {
      for (int i = 0; i <= LNX+1; ++i) {
        for (int j = 0; j <= LNY+1; ++j) {
          f[q][i][j][LNZ+1] = rbuf_hi[n++];
        }
      }
    }
  }

}

void halo_exchange_int(int a[LNX_H][LNY_H][LNZ_H])
{
  MPI_Request req[4];
  int n;

  // pass 1: x-normal faces: j=1..LNY, k=1..LNZ
  {
    n = 0;
    for (int j = 1; j <= LNY; ++j) {
      for (int k = 1; k <= LNZ; ++k) {
        isbuf_lo[n++] = a[1][j][k];
      }
    }
    n = 0;
    for (int j = 1; j <= LNY; ++j) {
      for (int k = 1; k <= LNZ; ++k) {
        isbuf_hi[n++] = a[LNX][j][k];
      }
    }

    MPI_Isend(isbuf_lo, XIBUF, MPI_INT, neighbours[0], 110, cart_comm, &req[0]);
    MPI_Irecv(irbuf_lo, XIBUF, MPI_INT, neighbours[0], 111, cart_comm, &req[1]);
    MPI_Isend(isbuf_hi, XIBUF, MPI_INT, neighbours[1], 111, cart_comm, &req[2]);
    MPI_Irecv(irbuf_hi, XIBUF, MPI_INT, neighbours[1], 110, cart_comm, &req[3]);
    MPI_Waitall(4, req, MPI_STATUSES_IGNORE);

    n = 0;
    for (int j = 1; j <= LNY; ++j) {
      for (int k = 1; k <= LNZ; ++k) {
        a[0][j][k] = irbuf_lo[n++];
      }
    }
    n = 0;
    for (int j = 1; j <= LNY; ++j) {
      for (int k = 1; k <= LNZ; ++k) {
        a[LNX+1][j][k] = irbuf_hi[n++];
      }
    }
  }

  // pass 2: y-normal faces: i=0..LNX+1, k=1..LNZ
  {
    n = 0;
    for (int i = 0; i <= LNX+1; ++i) {
      for (int k = 1; k <= LNZ; ++k) {
        isbuf_lo[n++] = a[i][1][k];
      }
    }
    n = 0;
    for (int i = 0; i <= LNX+1; ++i) {
      for (int k = 1; k <= LNZ; ++k) {
        isbuf_hi[n++] = a[i][LNY][k];
      }
    }

    MPI_Isend(isbuf_lo, YIBUF, MPI_INT, neighbours[2], 120, cart_comm, &req[0]);
    MPI_Irecv(irbuf_lo, YIBUF, MPI_INT, neighbours[2], 121, cart_comm, &req[1]);
    MPI_Isend(isbuf_hi, YIBUF, MPI_INT, neighbours[3], 121, cart_comm, &req[2]);
    MPI_Irecv(irbuf_hi, YIBUF, MPI_INT, neighbours[3], 120, cart_comm, &req[3]);
    MPI_Waitall(4, req, MPI_STATUSES_IGNORE);

    n = 0;
    for (int i = 0; i <= LNX+1; ++i) {
      for (int k = 1; k <= LNZ; ++k) {
        a[i][0][k] = irbuf_lo[n++];
      }
    }
    n = 0;
    for (int i = 0; i <= LNX+1; ++i) {
      for (int k = 1; k <= LNZ; ++k) {
        a[i][LNY+1][k] = irbuf_hi[n++];
      }
    }
  }

  // pass 3: z-normal faces: i=0..LNX+1, j=0..LNY+1
  {
    n = 0;
    for (int i = 0; i <= LNX+1; ++i) {
      for (int j = 0; j <= LNY+1; ++j) {
        isbuf_lo[n++] = a[i][j][1];
      }
    }
    n = 0;
    for (int i = 0; i <= LNX+1; ++i) {
      for (int j = 0; j <= LNY+1; ++j) {
        isbuf_hi[n++] = a[i][j][LNZ];
      }
    }

    MPI_Isend(isbuf_lo, ZIBUF, MPI_INT, neighbours[4], 130, cart_comm, &req[0]);
    MPI_Irecv(irbuf_lo, ZIBUF, MPI_INT, neighbours[4], 131, cart_comm, &req[1]);
    MPI_Isend(isbuf_hi, ZIBUF, MPI_INT, neighbours[5], 131, cart_comm, &req[2]);
    MPI_Irecv(irbuf_hi, ZIBUF, MPI_INT, neighbours[5], 130, cart_comm, &req[3]);
    MPI_Waitall(4, req, MPI_STATUSES_IGNORE);

    n = 0;
    for (int i = 0; i <= LNX+1; ++i) {
      for (int j = 0; j <= LNY+1; ++j) {
        a[i][j][0] = irbuf_lo[n++];
      }
    }
    n = 0;
    for (int i = 0; i <= LNX+1; ++i) {
      for (int j = 0; j <= LNY+1; ++j) {
        a[i][j][LNZ+1] = irbuf_hi[n++];
      }
    }
  }
}
