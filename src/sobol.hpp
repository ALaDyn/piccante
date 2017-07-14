
int i8_bit_hi1 ( long long int n );
int i8_bit_lo0 ( long long int n );
long long int i8_max ( long long int i1, long long int i2 );
long long int i8_min ( long long int i1, long long int i2 );
void i8_sobol ( int dim_num, long long int *seed, double quasi[ ] );
double *i8_sobol_generate ( int m, int n, int skip );
long long int i8_uniform ( long long int b, long long int c, int *seed );

double r8_abs ( double x );
int r8_nint ( double x );
double r8_uniform_01 ( int *seed );

void r8mat_write ( std::string output_filename, int m, int n, double table[] );

int tau_sobol ( int dim_num );
void timestamp ( );

