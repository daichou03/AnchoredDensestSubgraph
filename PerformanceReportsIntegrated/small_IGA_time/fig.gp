reset
set terminal postscript eps color "Times-Roman" 18 enhanced
set size 1,1

set style fill solid border 
set style fill pattern
set boxwidth 1

unset mytics

set key right bottom
set key box 1
set key title "Size of R"

set xlabel "% Nodes with Degree >= vol(R)"
set ylabel "IGA Speedup"


set output 'fig.eps'

set encoding utf8

set offset 0,0,0.05,0.25

set style data points

array datasize[5]
datasize[1] = 4
datasize[2] = 6
datasize[3] = 8
datasize[4] = 10
datasize[5] = 12

plot for [i=1:5] 'fig.txt' using 4:($3 == datasize[i] ? $5 : 1/0) title sprintf("%d", datasize[i]) lt -1 ps 2 pt i
