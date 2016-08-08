/* test the effect of boundary conditions to the Langevin equation */
#include "../mtrand.h"

double xc = 0.3;
double k = 10;
double dt = 0.001;
double tp = 1;
int boundary = 0;

double Kb = 10;

#define N 1000
double dx = 1.0/N;
double hist[3*N] = {0};


/* get the value */
double gety(double x)
{
  double dx = x - xc;
  return exp(-0.5 * k * dx * dx / tp);
}

/* get the partition function */
double getint(int n)
{
  double q = 0, dx = 1.0 / n, x, y;
  int i;

  for ( i = 0; i <= n; i++ ) {
    x = i * dx;
    y = gety(x);
    if ( i == 0 || i == n ) y *= 0.5;
    q += y;
  }
  q *= dx;

  return q;
}

double getforce(double x)
{
  if ( x <= 0 ) {
    return Kb;
  }
  if ( x >= 1 ) {
    return -Kb;
  }
  return -k * (x - xc);
}

#define TMAX 100000000

int main(int argc, char **argv)
{
  double x = 0.5, x1, f, q, tot;
  int t, i;

  if ( argc > 1 ) {
    boundary = atoi(argv[1]);
  }

  for ( t = 1; t <= TMAX; t++ ) {
    f = getforce(x);
    x1 = x + f * dt + sqrt(2 * tp * dt) * randgaus();
    if ( boundary == 0 ) {
      /* ignore out-of-boundary moves */
      if ( x1 < 0 || x1 > 1 ) {
        x1 = x;
      }
    } else if ( boundary == 1 ) {
      /* reflective boundary conditions */
      if ( x1 < 0 ) {
        x1 = -x1;
      } else if ( x1 > 1 ) {
        x1 = 2 - x1;
      }
    } else if ( boundary == 2 ) {
      /* do nothing, use artifical boundary potential to do the work */
    }
    x = x1;
    if ( x >= -1 && x < 3 ) {
      hist[(int) ((x + 1) / dx)] += 1;
    }
  }

  q = getint(N*10);
  for ( tot = 0, i = 0; i < 3*N; i++ ) tot += hist[i];
  for ( i = -N; i < 2*N; i++ ) {
    if ( hist[i+N] <= 0 ) continue;
    x = (i + 0.5) / N;
    printf("%g %g %g\n", x, hist[i+N]*N/tot, gety(x)/q);
  }

  return 0;
}
