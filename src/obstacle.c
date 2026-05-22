#include <math.h>
#include "param.h"

void build_obstacle(int obst[NX][NY][NZ])
{

  for (int k = 0; k < NZ; ++k) {
    for (int j = 0; j < NY; ++j) {
      for (int i = 0; i < NX; ++i) {

        double x = (double)i - XC;
        double y = (double)j - YC;
        double z = (double)k - ZC;

        if ((x*x + z*z) <= (OBST_R*OBST_R) ) {
          obst[i][j][k] = 1;
        } else {
          obst[i][j][k] = 0;
        }
      }
    }
  }

}
