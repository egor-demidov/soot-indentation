#!/usr/bin/env gnuplot

set encoding iso_8859_1
set term postscript eps colour enhanced 'Helvetica' 20
set output "04_particles_colliding.eps"
set ylabel "{/Helvetica=28 u/u_0}"
set xlabel "{/Helvetica=28 t, us}"
set xr[0:4e-6]
set yr[-0.1:1.05]
set key right bottom

plot "04_particles_colliding.dat" using 1:2 with lines title "Translational kinetic energy",\
     "04_particles_colliding.dat" using 1:3 with lines title "Rotational kinetic energy",\
     "04_particles_colliding.dat" using 1:4 with lines dashtype 2 linewidth 2 title "Linear momentum",\
     "04_particles_colliding.dat" using 1:5 with lines dashtype 1 linewidth 2 title "Angular momentum"