reset
set terminal postscript eps color "Times-Roman" 18 enhanced
set size 1,0.65

set style fill solid border 
set style fill pattern
set boxwidth 1

set logscale x
set logscale y

set xtics ('10^{-3}' 0.001, '10^{-2}' 0.01, '10^{-1}' 0.1)
set ytics ('10^{-1}' 0.1, '10^{0}' 1, '10^{1}' 10, '10^{2}' 100, '10^{3}' 1000)

set xrange [0.0005:0.2]
set yrange [1:2000]


set key left bottom
set key box 1
set key title "Number of Edges in the graph"

set xlabel "% of Nodes Explored by LA"
set ylabel "LA Space Reduction Ratio"


set output 'fig.eps'

set encoding utf8

set offset 0,0,0,0

set style data points

plot 'fig_l.txt' using 2:3:(sprintf("%d", $1)) title '>= 5000000' lt -1 with labels point ps 2 pt 5 offset char 0,1.2,\
    'fig_s.txt' using 2:3:(sprintf("%d", $1)) title '< 5000000' lt -1 ps 2 pt 4,\