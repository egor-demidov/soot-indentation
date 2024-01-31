#!/usr/bin/env gnuplot

set terminal postscript eps colour enhanced font 'Helvetica,14'
set output "sintering_plot.eps"

set style line  1 lc rgb '#000000' pt 5 ps 0.2 lt 1 lw 2
set style line  2 lc rgb '#000000' dashtype 2 pt 5 ps 0.2 lt 1 lw 2
set style line  3 lc rgb '#000000' dashtype 3 pt 5 ps 0.2 lt 1 lw 2

set tics scale 0.5
set xlabel 't, s'
set ylabel 'E, J'
set size 0.5, 0.5
set xtics 2e-8
set yrange[-0.1:1.1]

unset key

plot "08_sintering.dat" using 1:2 with lines ls 1,\
    "08_sintering.dat" using 1:3 with lines ls 2,\
    "08_sintering.dat" using 1:7 with lines ls 3;