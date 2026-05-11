#include "param.h"

const double PI = 3.14159265358979323846;
const double TAU = 1.0;
const double U0 = 1.0;

const double LX_PHYS = 200.0;
const double LY_PHYS = 100.0;
const double DX_PHYS = LX_PHYS / NX;
const double DY_PHYS = LY_PHYS / NY;
const double U0_PHYS = 1.0;
const double DT_PHYS = DX_PHYS / U0_PHYS;

const double w[Q] = {
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0
};

const int c[Q][2] = {
    { 0,  0},
    { 1,  0},
    { 0,  1},
    {-1,  0},
    { 0, -1},
    { 1,  1},
    {-1,  1},
    {-1, -1},
    { 1, -1}
};
