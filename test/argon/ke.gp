#!/usr/bin/env gnuplot
set terminal push
set terminal pngcairo
set output "ke.png"

set xtics 20
set mxtics 4
set xlabel "Kinetic energy"

set ytics 0.02
set mytics 2
set ylabel "Distribution"

set key left spacing 1.5

set key left Left reverse

plot [60:120][:] \
  "ke.dat"       u 1:2  w lp t "Observed", \
  "ke.dat"       u 1:3  w lp t "Reference"


unset output
set terminal pop
