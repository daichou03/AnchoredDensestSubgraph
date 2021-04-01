reset
set terminal postscript eps color "Times-Roman" 18 enhanced
set size 1,1

set style fill solid border 
set style fill pattern
set boxwidth 1

set logscale x
set logscale y

set xtics ('10^{2}' 100, '10^{3}' 1000, '10^{4}' 10000, '10^{5}' 100000, '10^{6}' 1000000)
set ytics ('10^{-1}' 0.1, '10^{0}' 1, '10^{1}' 10, '10^{2}' 100, '10^{3}' 1000)

set xrange [100:1000000]
set yrange [0.1:1000]


set key off

set xlabel "N / density(G)"
set ylabel "SLADS Speed Up"


set output 'fig.eps'

set encoding utf8

set offset 0,0,0,0

set style data points

plot 'fig.txt' using 2:3:(sprintf("%d", $1)) notitle lt -1 ps 2 with labels point pt 1 offset char 0,1
