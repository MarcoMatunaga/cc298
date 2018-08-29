#!/usr/bin/gnuplot

# Script gnuplot to iteratively plot residuals of bru3d.

set term png
set output "residue.png"

set xtics font "Times-Roman, 10"
set ytics font "Times-Roman, 10"

set xlabel "iterations"
set ylabel "residue"
set title "Projeto1"


set style line 1 \
    linecolor rgb '#0060ad' \
    linetype 1 linewidth 2 \
    pointtype 7 pointsize 1.5

plot 'residue.dat' with linespoints linestyle 1
