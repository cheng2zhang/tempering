#!/usr/bin/env gnuplot
# Demonstration of the importance of velocity scaling

set encoding cp1250 # make the minus sign longer
set terminal pngcairo enhanced size 640, 480
set output "vscale.png"
set xlabel "Reduced temperature"
set ylabel "Potential energy"
set key left spacing 1.5

plot "good.dat" u 2:3 w lp t "With velocity scaling", \
     "bad.dat"  u 2:3 w lp t "Without velocity scaling", \
     "good.dat" u 2:4 w lp t "Reference"

unset output
reset
