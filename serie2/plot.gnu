#!/usr/bin/gnuplot

# Script gnuplot to iteratively plot residuals of bru3d.

set term png
set output "output.png"

set xtics font "Times-Roman, 10"
set ytics font "Times-Roman, 10"

set xlabel "x position"
set ylabel "properties"
set title "Série2"

set style line 1 lt rgb "black" 
set style line 2 lt rgb "red" 

plot "diferenças_segundas" u 1:2 with lines title "initial",\
     "diferenças_segundas" u 1:3 with lines title "f2",\
     "diferenças_quartas" u 1:3 with lines title "f4"
