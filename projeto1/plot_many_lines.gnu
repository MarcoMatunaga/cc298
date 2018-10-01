#!/usr/bin/gnuplot

# Script gnuplot to iteratively plot residuals of bru3d.
# example also retired from the site below
# http://gnuplot.sourceforge.net/demo_5.0/lines_arrows.html

set term png
set output "residue.png"

set xtics font "Times-Roman, 10"
set ytics font "Times-Roman, 10"

set xlabel "iterations"
set ylabel "residue"
set title "Projeto1"

# set style line 1 \
    linecolor rgb '#0060ad' \
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 1.5
#
# reset linetypes to base dash patterns
#
set for [i=1:5] linetype i dt i

#
# define line styles using explicit rgbcolor names
#
set style line 1 lt 2 lc rgb "red" lw 3
set style line 2 lt 2 lc rgb "orange" lw 2
set style line 3 lt 2 lc rgb "yellow" lw 3
set style line 4 lt 2 lc rgb "green" lw 2

#
set label 1 'set style line 1 lt 2 lc rgb "red" lw 3'    at -0.4, -0.25 tc rgb "red"
set label 2 'set style line 2 lt 2 lc rgb "orange" lw 2' at -0.4, -0.35 tc rgb "orange"
set label 3 'set style line 3 lt 2 lc rgb "yellow" lw 3' at -0.4, -0.45 tc rgb "yellow"
set label 4 'set style line 4 lt 2 lc rgb "green" lw 2'  at -0.4, -0.55 tc rgb "green"
set label 5 'plot ... lt 1 lc 3 ' at -0.4, -0.65 tc lt 3
set label 6 'plot ... lt 3 lc 3 ' at -0.4, -0.75 tc lt 3
set label 7 'plot ... lt 5 lc 3 ' at -0.4, -0.85 tc lt 3
#
set xlabel "You will only see dashed lines if your current terminal setting permits it"
#
show style line
#
# draw some plots
#
plot cos(x)     ls 1 title 'ls 1',   \
     cos(x-.2)  ls 2 title 'ls 2',\
     cos(x-.4)  ls 3 title 'ls 3',\
     cos(x-.6)  ls 4 title 'ls 4', \
     cos(x-.8)  lt 1 lc 3 title 'lt 1 lc 3',  \
     cos(x-1.)  lt 3 lc 3 title 'lt 3 lc 3',  \
     cos(x-1.2) lt 5 lc 3 title 'lt 5 lc 3'
