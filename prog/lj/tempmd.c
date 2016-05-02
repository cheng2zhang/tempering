/* basic molecular dynamics simulation */
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
double tpdel = 0.1; /* temperature increment */
int ntp = 12; /* number of temperatures */
int scale_velocity_after_tempering = 1; /* always set it to 1 */



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
static simtemp_t *init_simtemp(void)
{
  int itp;
  double tp, Fex;
  simtemp_t *st;

  st = simtemp_open(ntp);
  for ( itp = 0; itp < ntp; itp++ ) {
    tp = tpmin + tpdel * itp;
    st->beta[itp] = 1.0 / tp;
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
  int t, itp = 0, jtp;
  double tp; /* current temperature */

  lj = lj_open(n, rho, rcdef);

  st = init_simtemp();
  itp = 0; /* initial temperature index */
  tp = 1 / st->beta[itp];

  for ( t = 1; t <= nequil + nsteps; t++ ) {
    /* standard MD step */
    thermostat(lj, tp, thdt * 0.5);
    lj_vv(lj, dt);
    thermostat(lj, tp, thdt * 0.5);

    /* add a data point */
    simultemp_add(st, itp, lj->epot);

    if ( t % nsttemp == 0 ) {
      /* temperature transition */
      jtp = simtemp_jump(st, itp, lj->epot,
          SIMTEMP_JUMP_GLOBAL);

      /* commit to the new temperature */
      if ( jtp != itp ) {
        if ( scale_velocity_after_tempering ) {
          simtemp_vscale(st, itp, jtp, n * D, (double *) lj->v);
        }
        lj->ekin = lj_ekin(lj->v, n);
        itp = jtp;
        tp = 1 / st->beta[itp];
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
  if ( argc > 0 ) {
    scale_velocity_after_tempering = atoi( argv[1] );
  }
  simtempmd();
  return 0;
}

