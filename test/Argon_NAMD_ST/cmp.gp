#!/usr/bin/env gnuplot
set terminal pngcairo
set output "cmp.png"

set xtics 0.5
set mxtics 5
set xlabel "Temperature"
set ytics 0.1
set mytics 5
set ylabel "Potential energy per particle"

set key left spacing 1.5

# from id to reduced temperature
eps = 0.238
sig = 3.405
rc = 2.5
kB = 0.0019872041
Tmin = 290
Tmax = 410
nbeta = 40
npart = 100
boxlen = 20.0 # in angstroms

betamin = 1/(kB*410)
betamax = 1/(kB*290)
dbeta = (betamax - betamin)/nbeta;
id2redtemp(x) = 1/(eps*(betamin + dbeta*(x-0.5)));


vol = (boxlen/sig)**3
den = npart/vol
utail = 8*pi/3*den*(1/rc**9 - 1/rc**3)
# potential energy per particle
upp(x) = x/(npart*eps) + utail

load "ljeos/ljeosJZG.gp"

plot [2.4:3.4][:] \
  "<awk '{print FNR,$0}' novscale.txt" u (id2redtemp($1)):(upp($2)) w lp t "No velocity scaling", \
  "<awk '{print FNR,$0}' vscale.txt"   u (id2redtemp($1)):(upp($2)) w lp t "Velocity scaling", \
  U(den, x) t "Equation of state"


unset output
