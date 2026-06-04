# =======================================================================
# Makefile for D3Q19 LBM – MPI parallel
#
# Usage examples:
#
#   make                         # 1×1×1  (serial behaviour, 1 MPI rank)
#   make NPROCX=2 NPROCZ=2       # 2×1×2  decomposition  (4 ranks)
#   make NPROCX=4 NPROCY=1 NPROCZ=4  # 4×1×4  decomposition (16 ranks)
#
#   make run                     # build then mpirun with correct nproc
#   make run NPROCX=2 NPROCZ=2
#
#   make debug                       # debug build  (1×1×1)
#   make debug NPROCX=2 NPROCZ=2     # debug build  (2×1×2)
#   make run-debug                   # debug build + mpirun
#
# =======================================================================

CC      = mpicc
CFLAGS  = -O2 -Wall -Wno-unused-variable -Wno-unused-function

# Debug flags (used when DEBUG=1 or via the 'debug' target)
DEBUG_FLAGS = -g -O0 -DDEBUG

# Topology (override on command line)
NPROCX ?= 1
NPROCY ?= 1
NPROCZ ?= 1
NPROC   = $(shell echo $$(($(NPROCX) * $(NPROCY) * $(NPROCZ))))

CFLAGS += -DNPROCX=$(NPROCX) -DNPROCY=$(NPROCY) -DNPROCZ=$(NPROCZ)

SRCS = src/main.c src/mpimod.c src/lbm.c src/bound.c src/force.c src/init.c src/obstacle.c src/output.c src/param.c src/halo.c
OBJS = $(SRCS:.c=.o)
EXE  = runme
EXE_DBG = runme_debug

# Allow DEBUG=1 on the command line as an alternative to the debug target
ifeq ($(DEBUG),1)
  CFLAGS := $(filter-out -O2,$(CFLAGS)) $(DEBUG_FLAGS)
  EXE    := $(EXE_DBG)
endif

all: $(EXE)

# Optimised build
$(EXE): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lm

%.o: %.c src/utils.h src/param.h src/subroutines.h
	$(CC) $(CFLAGS) -c -o $@ $<

# Debug target: rebuild everything with DEBUG_FLAGS, output to lbm_mpi_debug
# CFLAGS_DBG strips -O2 and appends the debug flags exactly once.
CFLAGS_DBG = $(filter-out -O2,$(CFLAGS)) $(DEBUG_FLAGS)

debug:
	$(CC) $(CFLAGS_DBG) -o $(EXE_DBG) $(SRCS) -lm

$(EXE_DBG): $(SRCS) utils.h param.h subroutines.h
	$(CC) $(CFLAGS_DBG) -o $@ $(SRCS) -lm

run: $(EXE)
	mkdir -p data
	mpirun -n $(NPROC) ./$(EXE)

run-debug: $(EXE_DBG)
	mkdir -p data
	mpirun -n $(NPROC) ./$(EXE_DBG)

clean:
	rm -f $(OBJS) $(EXE) $(EXE_DBG)

.PHONY: all debug run run-debug clean
