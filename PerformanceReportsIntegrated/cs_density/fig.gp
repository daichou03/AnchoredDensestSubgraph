reset
set terminal postscript eps color "Times-Roman" 18 enhanced
set size 1,1

set style fill solid border 
set style fill pattern
set boxwidth 1

unset mytics

set key off

set xlabel "deg(V)"
set ylabel "Density"

array colors[3]
colors[1] = "#CC0000"
colors[2] = "#00CC00"
colors[3] = "#0000CC"

set xrange [1:41]

set output 'fig.eps'

set encoding utf8

set style data points

plot for [i=1:3] 'fig.txt' using (column(0)):i w linespoints notitle lt rgb colors[i] dashtype 2**(i-1) ps 0.5 pt 1
