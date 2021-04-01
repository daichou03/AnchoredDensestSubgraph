reset
set terminal postscript eps color "Times-Roman" 18 enhanced
set size 1,1

set style fill solid border 
set style fill pattern
set boxwidth 1

set key top horizontal center outside
set key height 1
set key width 5


unset mytics


set xlabel "% Nodes with Degree < vol(R)"
set ylabel "IADS Speed Up"


set output 'fig.eps'

unset key
set encoding utf8

set offset 0,0,0.05,0.25
set style data points
plot for [i=1:9] 'fig.txt' using 4:($2 == i ? $5 : 1/0):(sprintf("%d", $3)) lt -1 ps 2 with labels point pt i offset char 0,1
