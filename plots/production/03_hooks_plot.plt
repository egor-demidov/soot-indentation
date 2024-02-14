#!/usr/bin/env gnuplot

set terminal postscript eps colour enhanced font 'Helvetica,14'
set output "03_hooks_plot.eps"

set style line  1 lc rgb '#000000' pt 5 ps 0.2 lt 1 lw 2
set style line  2 lc rgb '#000000' dashtype 2 pt 5 ps 0.2 lt 1 lw 2
set style line  3 lc rgb '#000000' dashtype 3 pt 5 ps 0.2 lt 1 lw 2

#set tics scale 0.5
#set xlabel 't, s'
#set ylabel 'E, J'
#set size 0.5, 0.5
#set xtics 2e-8
#set yrange[0:2.5e-11]
#set xrange[0:1.6e-6]

unset key

plot "03_hooks.dat" using 1:9 with lines ls 1;
