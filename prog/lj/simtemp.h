#ifndef SIMTEMP_H__
#define SIMTEMP_H__



#include <math.h>
#include <float.h>
#include "mtrand.h"



/* simulated tempering object */
typedef struct {
  int ntp; /* number of temperatures */
  double *beta;
  double *lnz;
  double *hist;
  double *usum;
  double *uref;
} simtemp_t;



static simtemp_t *simtemp_open(int ntp)
{
  simtemp_t *st;
  int itp;

  xnew(st, 1);
  st->ntp = ntp;
  xnew(st->beta, ntp);
  xnew(st->lnz, ntp);
  xnew(st->hist, ntp);
  xnew(st->usum, ntp);
  xnew(st->uref, ntp);
  for ( itp = 0; itp < ntp; itp++ ) {
    st->hist[itp] = DBL_MIN;
    st->usum[itp] = 0.0;
  }
  return st;
}


static void simtemp_close(simtemp_t *st)
{
  free(st->beta);
  free(st->lnz);
  free(st->hist);
  free(st->usum);
  free(st->uref);
  free(st);
}



static void simultemp_add(simtemp_t *st, int itp, double ep)
{
  st->hist[itp] += 1;
  st->usum[itp] += ep;
}


enum {
  SIMTEMP_JUMP_GLOBAL,
  SIMTEMP_JUMP_LOCAL,
  SIMTEMP_JUMP_COUNT
};


__inline static int simtemp_jump(simtemp_t *st, int itp, double ep, int type)
{
  int jtp, ntp = st->ntp;
  double x, r;

  if ( type == SIMTEMP_JUMP_GLOBAL ) {
    /* choose any temperature other than `itp` */
    jtp = ( itp + 1 + (int) (rand01() * (ntp - 1)) ) % ntp;
  } else {
    /* jump to one of the neighbors */
    if ( rand01() > 0.5 ) {
      jtp = itp + 1;
      if ( jtp >= ntp ) return itp;
    } else {
      jtp = itp - 1;
      if ( jtp < 0 ) return itp;
    }
  }

  /* compute the acceptance probability */
  x = (st->beta[itp] - st->beta[jtp]) * ep
    + st->lnz[itp] - st->lnz[jtp];
  if ( x > 0 ) {
    return jtp;
  } else {
    r = rand01();
    return r < exp(x) ? jtp : itp;
  }
}



__inline static void simtemp_vscale(simtemp_t *st, int itp, int jtp,
    int n, double *v)
{
  double s;
  int i;

  s = sqrt( st->beta[itp] / st->beta[jtp] );
  for ( i = 0; i < n; i++ ) {
    v[i] *= s;
  }
}



/* show statistics */
static void simtemp_dump(simtemp_t *st)
{
  int itp, ntp = st->ntp;
  double u, tp;

  printf("#id\t  temper.\t    uave  \t   uref  \t     hist\n");
  for ( itp = 0; itp < ntp; itp++ ) {
    tp = 1.0 / st->beta[itp];
    u = st->usum[itp] / st->hist[itp];

    printf("%3d\t%9.5f\t%9.3f\t%9.3f\t%9.0f\n",
        itp, tp, u, st->uref[itp], st->hist[itp]);
  }

}

#endif /* SIMTEMP_H__ */

