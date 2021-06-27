reset
set terminal postscript eps color "Times-Roman" 18 enhanced
set size 1,1

set style fill solid border 
set style fill pattern
set boxwidth 1

unset mytics

set key off

set xlabel "d_c"
set ylabel "Speedup VS d_c = 16"

set xrange [1:16]

set output 'fig.eps'

set encoding utf8

set style data points
set label front 'digg' at graph 0.45, graph 0.77
set label front 'trec' at graph 0.45, graph 0.57
set label front 'google' at graph 0.45, graph 0.38
set label front 'youtube' at graph 0.45, graph 0.15

plot for [i=2:25] 'fig.txt' using (2**(column(0))):i w linespoints notitle lt -1 ps 1 pt 1
