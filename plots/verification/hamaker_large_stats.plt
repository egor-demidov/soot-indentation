#!/usr/bin/env gnuplot

set terminal postscript eps colour enhanced font 'Helvetica,14'
set output "hamaker_large_stats.eps"

set style line  1 lc rgb '#000000' pt 5 ps 0.2 lt 1 lw 2

set tics scale 0.5
set xlabel 'D, nm'
set ylabel 'E, J'
set size 0.5, 0.5
set xrange [0:5e-9]

unset key

plot "06_hamaker.dat" using 1:2 with lines ls 1,\
    "06_hamaker.dat" using 1:4 with lines ls 1,\
    "06_hamaker.dat" using 1:8 with lines ls 1;