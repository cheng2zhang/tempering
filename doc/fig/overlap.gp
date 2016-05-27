#!/usr/bin/env gnuplot


set encoding cp1250 # make the minus sign longer
set terminal push
# dl 3 make dashed line longer
set terminal postscript eps enhanced dl 3 size 3.7, 3.7 font "Times, 24"
set output "overlap.eps"
set multiplot

T1 = 1.0
T2 = 2.0
s = 2 # half of the number of degrees of freedom

f(x, y, Tx, Ty) = (x/Tx)**(s-1) * exp(-x/Tx) * (y/Ty)**(s-1) * exp(-y/Ty)


set isosamples 100
set view map
set contour
unset surface
set cntrparam levels 1

#set palette maxcolors 2
#set pm3d map interpolate 4,4
#set style fill transparent solid 0.50

unset key
#unset colorbox

set border 10
set bmargin 0

set xrange [0:4.8]
set yrange [0:4.8]
unset xtics
unset ytics
set xlabel "Potential energy" offset 0, 2
set ylabel "Kinetic energy"

set label at 2.5, 2.10 "{/Times-Italic T} = 2"
set label at 0.6, 4.50 "{/Times-Italic T} = 1,velocity scaled"
set label at 0.7, 0.62 "{/Times-Italic T} = 1"

set arrow from 1.0, 4.3 to 1, 4

splot f(x, y, T1, T1), f(x, y, T1, T2), f(x, y, T2, T2)


unset multiplot
unset output
set terminal pop
reset
