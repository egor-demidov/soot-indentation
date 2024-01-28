#!/usr/bin/env gnuplot

set encoding iso_8859_1
set term postscript eps colour enhanced 'Helvetica' 20
set output "03_particles_colliding.eps"
set ylabel "{/Helvetica=28 u/u_0}"
set xlabel "{/Helvetica=28 t, us}"
set xr[0:4e-6]
set yr[0.5:1.05]
set key right bottom

plot "03_particles_colliding.dat" using 1:2 with lines title "Kinetic energy",\
     "03_particles_colliding.dat" using 1:3 with lines dashtype 2 linewidth 2 title "Linear momentum"