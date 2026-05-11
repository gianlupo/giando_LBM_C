#include "param.h"

const double PI = 3.14159265358979323846;

// relaxation time
const double TAU = 1.0;

// lattice velocity
const double U0 = 1.0;

// physical parameters
const double LX_PHYS = 200.0;
const double LY_PHYS = 4.0;
const double LZ_PHYS = 100.0;
const double DX_PHYS = LX_PHYS / NX;
const double DY_PHYS = LY_PHYS / NY;
const double DZ_PHYS = LZ_PHYS / NZ;
const double U0_PHYS = U0;
const double DT_PHYS = DX_PHYS / U0_PHYS;

// D3Q19 weights
const double weights[Q] = {
                           12.0 / 36.0, // 0
                            2.0 / 36.0, // 1
                            2.0 / 36.0, // 2
                            2.0 / 36.0, // 3
                            2.0 / 36.0, // 4
                            2.0 / 36.0, // 5
                            2.0 / 36.0, // 6
                            1.0 / 36.0, // 7
                            1.0 / 36.0, // 8
                            1.0 / 36.0, // 9
                            1.0 / 36.0, // 10
                            1.0 / 36.0, // 11
                            1.0 / 36.0, // 12
                            1.0 / 36.0, // 13
                            1.0 / 36.0, // 14
                            1.0 / 36.0, // 15
                            1.0 / 36.0, // 16
                            1.0 / 36.0, // 17
                            1.0 / 36.0  // 18
                                      };

// D3Q19 lattice velocities
const int c[Q][3] = {
                     { 0,  0,  0}, // 0
                     { 1,  0,  0}, // 1
                     {-1,  0,  0}, // 2
                     { 0,  1,  0}, // 3
                     { 0, -1,  0}, // 4
                     { 0,  0,  1}, // 5
                     { 0,  0, -1}, // 6
                     { 1,  1,  0}, // 7
                     { 1, -1,  0}, // 8
                     { 1,  0,  1}, // 9
                     { 1,  0, -1}, // 10
                     {-1,  1,  0}, // 11
                     {-1, -1,  0}, // 12
                     {-1,  0,  1}, // 13
                     {-1,  0, -1}, // 14
                     { 0,  1,  1}, // 15
                     { 0,  1, -1}, // 16
                     { 0, -1,  1}, // 17
                     { 0, -1, -1}  // 18
                                 };
