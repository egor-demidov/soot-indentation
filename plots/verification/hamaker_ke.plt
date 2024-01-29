#!/usr/bin/env gnuplot

set terminal postscript eps colour enhanced font 'Helvetica,14'
set output "hamaker_ke.eps"

set style line  1 lc rgb '#000000' pt 5 ps 0.2 lt 1 lw 2

set tics scale 0.5
set xlabel 'D, nm'
set ylabel 'E, J'
set size 0.5, 0.5

unset key

plot "05_hamaker.dat" using 1:2 with lines ls 1;