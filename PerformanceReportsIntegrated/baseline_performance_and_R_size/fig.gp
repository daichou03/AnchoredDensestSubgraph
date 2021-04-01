reset
set terminal postscript eps color "Times-Roman" 18 enhanced
set size 1,0.65

set style fill solid border 
set style fill pattern
set boxwidth 1

set logscale x
set logscale y

set xtics ('10^{-3}' 0.001, '10^{-2}' 0.01, '10^{-1}' 0.1)
set ytics ('10^{-1}' 0.1, '10^{0}' 1, '10^{1}' 10, '10^{2}' 100)

set xrange [0.00005:0.5]
set yrange [1:100]


set key off

set xlabel "% Nodes Explored"
set ylabel "SLADS Speed Up"



set output 'fig.eps'

set encoding utf8

set offset 0,0,0,0

set style data points

plot 'fig.txt' using 2:3:(sprintf("%d", $1)) every 3::0 notitle lt -1 ps 2 with labels point pt 1 offset char 0,1,\
    'fig.txt' using 2:3:(sprintf("%d", $1)) every 3::1 notitle lt -1 ps 2 with labels point pt 1 offset char 0,-1,\
    'fig.txt' using 2:3:(sprintf("%d", $1)) every 3::2 notitle lt -1 ps 2 with labels point pt 1 offset char 0,1