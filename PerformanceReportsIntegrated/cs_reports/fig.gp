reset
set terminal postscript eps color "Times-Roman" 18 enhanced
set size 1,1

set style fill solid border 
set style fill pattern
set boxwidth 1

unset mytics

set key off


array colors[3]
colors[1] = "#CC0000"
colors[2] = "#00CC00"
colors[3] = "#0000CC"

set xrange [1:41]

set encoding utf8

set style data points

set xlabel "deg(V)"

set ylabel "Density"
set output 'density.eps'
plot for [i=1:3] 'density.txt' using (column(0)):i w linespoints notitle lt rgb colors[i] dashtype 2**(i-1) ps 0.5 pt 1

set ylabel "Conductance"
set output 'conductance.eps'
plot for [i=1:3] 'conductance.txt' using (column(0)):i w linespoints notitle lt rgb colors[i] dashtype 2**(i-1) ps 0.5 pt 1

set ylabel "L"
set output 'lscore.eps'
plot for [i=1:3] 'lscore.txt' using (column(0)):i w linespoints notitle lt rgb colors[i] dashtype 2**(i-1) ps 0.5 pt 1

set ylabel "R-Subgraph Density"
set output 'rsdensity.eps'
plot for [i=1:3] 'rsdensity.txt' using (column(0)):i w linespoints notitle lt rgb colors[i] dashtype 2**(i-1) ps 0.5 pt 1

set yrange [0:1]
set ylabel "% R in S"
set output 'rins.eps'
plot for [i=1:3] 'rins.txt' using (column(0)):i w linespoints notitle lt rgb colors[i] dashtype 2**(i-1) ps 0.5 pt 1

set yrange [0:0.5]
set ylabel "% S out of R"
set output 'soutofr.eps'
plot for [i=1:3] 'soutofr.txt' using (column(0)):i w linespoints notitle lt rgb colors[i] dashtype 2**(i-1) ps 0.5 pt 1

set yrange [0:60]
set ylabel "Output Set Size"
set output 'length.eps'
plot for [i=1:3] 'length.txt' using (column(0)):i w linespoints notitle lt rgb colors[i] dashtype 2**(i-1) ps 0.5 pt 1

set yrange [0:1.1]
set ytics ('0' 0.0, '0.2' 0.2, '0.4' 0.4, '0.6' 0.6, '0.8' 0.8, '1.0' 1.0, 'Inf' 1.1)
set ylabel "Local Conductance"
set output 'lconductance.eps'
plot for [i=1:3] 'lconductance.txt' using (column(0)):i w linespoints notitle lt rgb colors[i] dashtype 2**(i-1) ps 0.5 pt 1