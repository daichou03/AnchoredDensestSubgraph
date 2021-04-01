reset
set terminal postscript eps color "Times-Roman" 18 enhanced
set size 1,1

set style fill solid border 
set style fill pattern
set boxwidth 1

unset mytics

set key off

set xlabel "d_c"
set ylabel "Space Reduction Ratio VS d_c = 16"


set output 'fig.eps'

set encoding utf8

set style data points

plot for [i=2:25] 'fig.txt' using (2**(column(0))):i w linespoints notitle lt -1 ps 1 pt 1
