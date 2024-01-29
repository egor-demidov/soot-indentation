#!/usr/bin/env gnuplot

set terminal postscript eps colour enhanced font 'Helvetica,14'
set output "hamaker_parity.eps"

set style line  1 lc rgb '#000000' pt 5 ps 0.2 lt 1 lw 2

set tics scale 0.5
set xlabel 'E, J'
set ylabel 'E, J'
set size square 0.5, 0.5
set xtics 2.5e-21
set yrange [0:7.5e-21]
set xrange [0:7.5e-21]

unset key

plot "05_hamaker.dat" using 2:3 with lines ls 1;