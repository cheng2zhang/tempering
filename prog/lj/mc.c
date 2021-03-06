/* basic Monte Carlo simulation */
#include "lj.h"



int n = 100;
int nequil = 100000;
int nsteps = 1000000;
double rho = 0.5;
double tp = 2.5;
double rcdef = 1e9;
double amp = 0.2; /* Monte Carlo move size */
const char *fnpos = "lj.pos";



int main(void)
{
  int t, acc;
  lj_t *lj;
  double epsm = 0, accsm = 0;

  lj = lj_open(n, rho, rcdef);
  lj_energy(lj);
  for ( t = 1; t <= nequil + nsteps; t++ ) {
    acc = lj_metro(lj, amp, 1/tp);
    if ( t <= nequil ) continue;
    epsm += lj->epot;
    accsm += acc;
  }
  lj_writepos(lj, lj->x, lj->v, fnpos);
  lj_close(lj);
  printf("rho %g, tp %g, ep %g, acc %g%%\n",
      rho, tp, epsm/nsteps/n, 100.*accsm/nsteps);
  return 0;
}

