#!/usr/bin/env gnuplot

reset
set terminal postscript eps colour enhanced font 'Helvetica,14'
set output "sintering_multiplot.eps"

# color definitions
set style line  1 lc rgb '#000000' pt 5 ps 0.2 lt 1 lw 2    # blue
set style line  2 lc rgb '#000000' dashtype 2 pt 5 ps 0.2 lt 1 lw 2    # orange
set style line  3 lc rgb '#000000' dashtype 3 pt 5 ps 0.2 lt 1 lw 2    # orange

set macros

set tics scale 0.5
set ytics 0.2
set xtics 2e-8
set xrange [-1.5e-8:11.5e-8]
set yrange [-0.15:1.15]

# MACROS
# x- and ytics for each row resp. column
NOXTICS = "set xtics format ''; unset xlabel"
XTICS = "set xtics rotate by 45 right format '%.0sx10^{%01S}'; set xlabel 't, s'"
NOYTICS = "set format y ''; unset ylabel"
YTICS = "set format y '%.1f'; set ylabel 'u_{norm}'"
KEY = "set key at graph 0.92, 0.75"
NOKEY = "unset key"
# Margins for each row resp. column
TMARGIN = "set tmargin at screen 0.95; set bmargin at screen 0.55"
BMARGIN = "set tmargin at screen 0.55; set bmargin at screen 0.15"
LMARGIN = "set lmargin at screen 0.1; set rmargin at screen 0.525"
RMARGIN = "set lmargin at screen 0.525; set rmargin at screen 0.95"
# Placement of the a,b,c,d labels in the graphs
POS = "at graph 0.92,0.89 font 'Helvetica Bold,12'"

### Start multiplot (2x2 layout)
set multiplot layout 2,2 rowsfirst
# --- GRAPH a
@TMARGIN; @LMARGIN; @NOXTICS; @YTICS; @NOKEY
set label 1 'a' @POS
plot "07_sintering.dat" using 1:2 with lines ls 1,\
    "07_sintering.dat" using 1:3 with lines ls 2,\
    "07_sintering.dat" using 1:7 with lines ls 3;
# --- GRAPH b
@TMARGIN; @RMARGIN; @NOXTICS; @NOYTICS; @KEY
set label 1 'b' @POS
plot "08_sintering.dat" using 1:2 with lines ls 1 title "E_{trs}",\
    "08_sintering.dat" using 1:3 with lines ls 2 title "E_{rot}",\
    "08_sintering.dat" using 1:7 with lines ls 3 title "P";
# --- GRAPH c
@BMARGIN; @LMARGIN; @XTICS; @YTICS; @NOKEY
set label 1 'c' @POS
plot "09_sintering.dat" using 1:2 with lines ls 1,\
    "09_sintering.dat" using 1:3 with lines ls 2,\
    "09_sintering.dat" using 1:7 with lines ls 3;
# --- GRAPH d
@BMARGIN; @RMARGIN; @XTICS; @NOYTICS; @NOKEY
set label 1 'd' @POS
plot "10_sintering.dat" using 1:2 with lines ls 1,\
    "10_sintering.dat" using 1:3 with lines ls 2,\
    "10_sintering.dat" using 1:7 with lines ls 3;
unset multiplot
