void initialize(int ic,
                double rho[NX][NY][NZ],
                double u[NX][NY][NZ],
                double v[NX][NY][NZ],
                double w[NX][NY][NZ]);

void compute_force(double rho[NX][NY][NZ],
                   double u[NX][NY][NZ],
                   double v[NX][NY][NZ],
                   double w[NX][NY][NZ],
                   double fx[NX][NY][NZ],
                   double fy[NX][NY][NZ],
                   double fz[NX][NY][NZ]);

void equilibrium(double rho[NX][NY][NZ],
                 double u[NX][NY][NZ],
                 double v[NX][NY][NZ],
                 double w[NX][NY][NZ],
                 double feq[Q][NX][NY][NZ]);

void stream(double f[Q][NX][NY][NZ],
            double fstar[Q][NX][NY][NZ]);

void moments(double fstar[Q][NX][NY][NZ],
             double fx[NX][NY][NZ],
             double fy[NX][NY][NZ],
             double fz[NX][NY][NZ],
             double rho[NX][NY][NZ],
             double u[NX][NY][NZ],
             double v[NX][NY][NZ],
             double w[NX][NY][NZ]);

void guo(double fx[NX][NY][NZ],
         double fy[NX][NY][NZ],
         double fz[NX][NY][NZ],
         double u[NX][NY][NZ],
         double v[NX][NY][NZ],
         double w[NX][NY][NZ],
         double fguo[Q][NX][NY][NZ]);

void collide(double fstar[Q][NX][NY][NZ],
             double feq[Q][NX][NY][NZ],
             double fguo[Q][NX][NY][NZ],
             double f[Q][NX][NY][NZ]);

void boundary(double fstar[Q][NX][NY][NZ]);

void triperiodic_bc(double f[Q][NX][NY][NZ],
                    double fstar[Q][NX][NY][NZ]);

void write_field(int iter,
                 double field[NX][NY][NZ],
                 const char* name);

void write_0d(int iter,
              double time,
              double var,
              const char* name);

void write_field_slice(int iter,
                       double slice[NX][NY][NZ],
                       int direction,
                       int position,
                       const char* name);

double max_u(double u[NX][NY][NZ]);

double max_kinetic(double rho[NX][NY][NZ],
                   double u[NX][NY][NZ],
                   double v[NX][NY][NZ],
                   double w[NX][NY][NZ]);

void compute_kinetic(double rho[NX][NY][NZ],
                     double u[NX][NY][NZ],
                     double v[NX][NY][NZ],
                     double w[NX][NY][NZ],
                     double Kin[NX][NY][NZ]);
