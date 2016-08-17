/* compute the derivative use in the MC-corrected Langevin equation */
#include <stdio.h>
#include <math.h>

#define kB 0.001987

double r = 1.2;

double energy(double beta, double *dedb, double *lnz)
{
  const double B = 2900, T0 = 4100;
  // const double B = 8100, T0 = 7100;
  double tp = 1/(kB*beta);
  *lnz = -B * beta * (log(tp/T0) + 1);
  *dedb = -B / beta;
  return B * log(tp/T0);
}

double getxp(double de0, double tp, double a, double r, double x,
    double *nbeta_, double *dnb, double *nr_, double *dnr, double *dxp)
{
  double dt = a*a/2, ntp, dtp, ep, nr;
  double beta, de, epave, dedb, lnz;
  double nbeta, nde, nepave, ndedb, nlnz;
  beta = 1./(kB*tp);
  epave = energy(beta, &dedb, &lnz);
  ep = epave + de0;
  de = ep - epave + kB*tp*x;
  dtp = dt*de/kB + r*a*tp;
  ntp = tp + dtp;
  nbeta = 1./(kB*ntp);
  *nbeta_ = nbeta;
  *dnb = -nbeta*nbeta/beta*(r + beta*de*a);
  nepave = energy(nbeta, &ndedb, &nlnz);
  nde = ep - nepave + kB*ntp*x;
  nr = (-dtp - dt*nde/kB) / (a * ntp);
  *nr_ = nr;
  *dnr = -(nbeta*nde+nr/a) + (1/beta/a - (ep-nepave-nbeta*ndedb)*a/2) * (*dnb);
  *dxp = -((nde - 3./nbeta)*(*dnb) + nr * (*dnr));
  return (beta - nbeta)*ep - nlnz + lnz + (x-3)*log(beta/nbeta) + (r*r - nr*nr)/2;
}

int main(void)
{
  double de0 = 27, tp = 300, r = 2.5, x = 1.0;
  double a1 = 0.02, a2 = 0.020001;
  double nbeta1, nbeta2, dnb1, dnb2;
  double nr1, nr2, dnr1, dnr2;
  double xp1, xp2, dxp1, dxp2;
  double dxp, dnb, dnr;
  xp1 = getxp(de0, tp, a1, r, x, &nbeta1, &dnb1, &nr1, &dnr1, &dxp1);
  xp2 = getxp(de0, tp, a2, r, x, &nbeta2, &dnb2, &nr2, &dnr2, &dxp2);
  dnb = (nbeta2 - nbeta1)/(a2 - a1);
  dnr = (nr2 - nr1)/(a2 - a1);
  dxp = (xp2 - xp1)/(a2 - a1);
  printf("beta %g, beta1 %g, beta2 %g, dbeta %g (%g, %g)\n",
      1./(kB*tp), nbeta1, nbeta2, dnb, dnb1, dnb2);
  printf("r %g, nr1 %g, nr2 %g, dnr %g (%g, %g)\n",
      r, nr1, nr2, dnr, dnr1, dnr2);
  printf("xp1 %g; xp2 %g; dxp %g (%g, %g)\n", xp1, xp2, dxp, dxp1, dxp2);
  return 0;
}
