#include "param.h"

const double PI = 3.14159265358979323846;

// relaxation time
const double TAU = 6.5;

// lattice velocity
const double U0 = 0.1;

// physical parameters
const double LX_PHYS = NX - 1;
const double LY_PHYS = NY - 1;
const double LZ_PHYS = NZ - 1;
const double DX_PHYS = LX_PHYS / (NX - 1);
const double DY_PHYS = LY_PHYS / (NY - 1);
const double DZ_PHYS = LZ_PHYS / (NZ - 1);
const double U0_PHYS = U0;
const double DT_PHYS = DX_PHYS / U0_PHYS;

// obstacle parameters
const double XC = 0.25 * LX_PHYS;
const double YC = LY_PHYS / 2.0;
const double ZC = (LZ_PHYS / 2.0) + 0.015 * LZ_PHYS;
const double OBST_R = LZ_PHYS / 40.0;

// initial condition
const int ic = 2; // 0 = zero
                  // 1 = uniform streamwise
                  // 2 = poiseuille
                  // 3 = sine wave

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

// D3Q19 opposite Q indices
const int opp[Q] = {
                    0 , // opposite of 0
                    2 , // opposite of 1
                    1 , // opposite of 2
                    4 , // opposite of 3
                    3 , // opposite of 4
                    6 , // opposite of 5
                    5 , // opposite of 6
                    12, // opposite of 7
                    11, // opposite of 8
                    14, // opposite of 9
                    13, // opposite of 10
                    8 , // opposite of 11
                    7 , // opposite of 12
                    10, // opposite of 13
                    9 , // opposite of 14
                    18, // opposite of 15
                    17, // opposite of 16
                    16, // opposite of 17
                    15, // opposite of 18
                      };

// D3Q19 populations streaming from outside the boundary plane
const int qo[6][5] =  {
                       { 1, 7, 8, 9,10}, // x-normal plane, imin
                       { 2,11,12,13,14}, // x-normal plane, imax
                       { 3, 7,11,15,16}, // y-normal plane, jmin
                       { 4, 8,12,17,18}, // y-normal plane, jmax
                       { 5, 9,13,15,17}, // z-normal plane, kmin
                       { 6,10,14,16,18}  // z-normal plane, kmax
                                       };

// D3Q19 populations streaming from inside the boundary plane
const int qi[6][14] = {
                       { 0, 2, 3, 4, 5, 6,11,12,13,14,15,16,17,18}, // x-normal plane, imin
                       { 0, 1, 3, 4, 5, 6, 7, 8, 9,10,15,16,17,18}, // x-normal plane, imax
                       { 0, 1, 2, 4, 5, 6, 8, 9,10,12,13,14,17,18}, // y-normal plane, jmin
                       { 0, 1, 2, 3, 5, 6, 7, 9,10,11,13,14,15,16}, // y-normal plane, jmax
                       { 0, 1, 2, 3, 4, 6, 7, 8,10,11,12,14,16,18}, // z-normal plane, kmin
                       { 0, 1, 2, 3, 4, 5, 7, 8, 9,11,12,13,15,17}  // z-normal plane, kmax
                                                                  };

// boundary conditions
const char bc[6] = {
                    'I', // x-normal plane, imin
                    'O', // x-normal plane, imax
                    'P', // y-normal plane, jmin
                    'P', // y-normal plane, jmax
                    'W', // z-normal plane, kmin
                    'W'  // z-normal plane, kmax
                       };

// boundary plane (outward) normals
const int n[6][3] = {
                     {-1, 0, 0}, // x-normal plane, imin
                     { 1, 0, 0}, // x-normal plane, imax
                     { 0,-1, 0}, // y-normal plane, jmin
                     { 0, 1, 0}, // y-normal plane, jmax
                     { 0, 0,-1}, // z-normal plane, kmin
                     { 0, 0, 1}  // z-normal plane, kmax
                               };
