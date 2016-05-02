/* continuous tempering code */
#include "lj.h"
#include "ljeos.h"
#include "simtemp.h"



int n = 100;
int nequil = 10000; /* number of equilibration steps */
int nsteps = 500000; /* number of simulation steps */
double rho = 0.5;
double rcdef = 2.5;
double dt = 0.002;

enum {
  THERMOSTAT_LANGEVIN,
  THERMOSTAT_VRESCALE,
  THERMOSTAT_COUNT
};
int thermostat_type = THERMOSTAT_VRESCALE;
double thdt = 0.02; /* strength of the thermostat */

int nsttemp = 10; /* number of steps of trying tempering */
double tpmin = 2.5; /* minimal temperature */
double tpdel = 0.02; /* temperature increment */
int ntp = 60; /* number of temperatures */

double betamin;
double betamax;
double betadel;
double ensfac = 2; /* flat-T distribution: 2, dT/T: 1, flat-beta: 0 */
double langdt = 0.0001; /* Langevin step size */


/* thermostat for half MD step */
static void thermostat(lj_t *lj, double temp, double thermdt)
{
  if ( thermostat_type == THERMOSTAT_LANGEVIN ) {
    lj->ekin = lj_langevin(lj, temp, thermdt);
  } else {
    lj->ekin = lj_vrescale(lj, temp, thermdt);
  }
}



/* initialize simulated tempering */
static simtemp_t *init_conttemp(void)
{
  int itp;
  double beta, tp, Fex;
  simtemp_t *st;

  st = simtemp_open(ntp);
  betamin = 1 / ( tpmin + ntp * tpdel );
  betamax = 1 / ( tpmin );
  betadel = (betamax - betamin) / ntp;
  for ( itp = 0; itp < ntp; itp++ ) {
    beta = betamin + (itp + 0.5) * betadel;
    tp = 1 / beta;
    st->beta[itp] = beta;
    st->uref[itp] = n * ljeos3d_get(rho, tp, NULL, &Fex, NULL);
    st->lnz[itp] = -st->beta[itp] * Fex * n;
  }
  return st;
}



/* simulated tempering */
static int simtempmd(void)
{
  lj_t *lj;
  simtemp_t *st;
  int t, itp = 0;
  double beta, tp; /* current temperature */

  lj = lj_open(n, rho, rcdef);

  st = init_conttemp();
  itp = 0; /* initial temperature index */
  beta = (st->beta[itp] + st->beta[itp+1]) / 2;
  tp = 1 / beta;

  for ( t = 1; t <= nequil + nsteps; t++ ) {
    /* standard MD step */
    thermostat(lj, tp, thdt * 0.5);
    lj_vv(lj, dt);
    thermostat(lj, tp, thdt * 0.5);

    /* add a data point */
    simultemp_add(st, itp, lj->epot);

    if ( t % nsttemp == 0 ) {
      /* temperature transition */
      double tp1, beta1, s;
      int i;
     
      tp1 = tp + (lj->epot - st->uref[itp] + ensfac * tp) * langdt
          + sqrt(2 * langdt) * tp * randgaus();
      beta1 = 1 / tp1;
      if ( beta1 > betamin && beta1 < betamax ) {
        /* commit to the new temperature */
        s = sqrt( tp1 / tp );
        for ( i = 0; i < lj->n; i++ )
          vsmul(lj->v[i], s);
        lj->ekin = lj_ekin(lj->v, n);
        beta = beta1;
        tp = tp1;
        itp = (int) ((beta - betamin) / betadel);
        if ( itp < 0 || itp >= ntp ) {
          fprintf(stderr, "beta out of range, beta %g (%g, %g)\n", beta, betamin, betamax);
          exit(1);
        }
      }
    }

    if ( t <= nequil ) continue;
  }
  simtemp_dump(st);
  simtemp_close(st);
  lj_close(lj);
  return 0;
}



int main(int argc, char **argv)
{
  simtempmd();
  return 0;
}

