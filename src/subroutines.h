void initialize(double rho[NX][NY],
                double u[NX][NY],
                double v[NX][NY]);

void equilibrium(double rho[NX][NY],
                 double u[NX][NY],
                 double v[NX][NY],
                 double feq[Q][NX][NY]);

void stream(double f[Q][NX][NY],
            double fstar[Q][NX][NY]);

void moments(double fstar[Q][NX][NY],
             double rho[NX][NY],
             double u[NX][NY],
             double v[NX][NY]);

void collide(double fstar[Q][NX][NY],
             double feq[Q][NX][NY],
             double f[Q][NX][NY]);

void periodic_bc(double f[Q][NX][NY],
                 double fstar[Q][NX][NY]);

void write_field(int iter,
                 double field[NX][NY],
                 const char* name);

void write_0d(int iter,
              double time,
              double var,
              const char* name);

double max_u(double u[NX][NY]);

double max_kinetic(double rho[NX][NY],
                   double u[NX][NY],
                   double v[NX][NY]);

void compute_kinetic(double rho[NX][NY],
                     double u[NX][NY],
                     double v[NX][NY],
                     double K[NX][NY]);
