#!/usr/bin/env gnuplot

# for the new format with the first column being the inverse temperature

set terminal push
set terminal pngcairo
set output "cmp2.png"

set xtics 0.2
set mxtics 4
set xlabel "Reduced temperature"

set ytics 0.1
set mytics 5
set ylabel "Potential energy per particle"

set key left spacing 1.5

# from id to reduced temperature
eps = 0.238
sig = 3.405
rc = 2.5
#kB = 0.0019872041
kB = 0.001987191  # the value used by VMD
Tmin = 290
Tmax = 410
nbeta = 60
npart = 100
boxlen = 20.0 # in angstroms

betamin = 1/(kB*Tmax)
betamax = 1/(kB*Tmin)
dbeta = (betamax - betamin)/nbeta;
invt2redtemp(x) = kB/(eps*x);


vol = (boxlen/sig)**3
den = npart/vol
# tail correction for the potential energy
utail = 8*pi/3*den*(1.0/rc**9 - 1.0/rc**3)
# potential energy per particle
upp(x) = x/(npart*eps) + utail

load "ljeos/ljeosKN.gp"

set key left Left reverse

plot [2.4:3.45][:-2.65] \
  "a.txt"   u (invt2redtemp($1)):(upp($2)) w lp t "Velocity scaling", \
  U(den, x) t "Equation of state, Kolafa and Nezbeda (1994)"


unset output
set terminal pop
