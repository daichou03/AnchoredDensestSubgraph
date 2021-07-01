reset
set terminal postscript eps color "Times-Roman" 18 enhanced
set size 1,1

set style fill solid border 
set style fill pattern
set boxwidth 1

unset mytics

set key off


array colors[4]
colors[1] = "#CC0000"
colors[2] = "#00CC00"
colors[3] = "#0000CC"
colors[4] = "#777777"

set xrange [1:41]

set encoding utf8

set style data points

set lmargin 10.5
set bmargin 4.5

set xlabel "d(x)" font ", 40" center offset 0,-0.5
set ylabel font ", 40" offset -2
set xrange [1:41]
set xtics ('1' 1, '5' 5, '10' 10, '15' 15, '20' 20, '25' 25, '30' 30, '35' 35, '' 40, '41+' 41)
set tics font ", 40"

set ylabel "Density"
set output 'density.eps'
plot for [i=1:3] 'density.txt' using (column(0)):i w linespoints notitle lt rgb colors[i] lw 3 dashtype 2**(i-1) ps 0.5 pt 1

set ylabel "L"
set output 'lscore.eps'
plot for [i=1:3] 'lscore.txt' using (column(0)):i w linespoints notitle lt rgb colors[i] lw 3 dashtype 2**(i-1) ps 0.5 pt 1

set ylabel "R-Subgraph Density"
set output 'rsdensity.eps'
plot for [i=1:3] 'rsdensity.txt' using (column(0)):i w linespoints notitle lt rgb colors[i] lw 3 dashtype 2**(i-1) ps 0.5 pt 1

set yrange [0:0.7]
set ylabel "Conductance"
set output 'conductance.eps'
plot for [i=1:3] 'conductance.txt' using (column(0)):i w linespoints notitle lt rgb colors[i] lw 3 dashtype 2**(i-1) ps 0.5 pt 1

set yrange [0:1]
set ylabel "% R in S"
set output 'rins.eps'
plot for [i=1:3] 'rins.txt' using (column(0)):i w linespoints notitle lt rgb colors[i] lw 3 dashtype 2**(i-1) ps 0.5 pt 1

set yrange [0:1]
set ylabel "% S in R"
set output 'sinr.eps'
plot for [i=1:3] 'sinr.txt' using (column(0)):i w linespoints notitle lt rgb colors[i] lw 3 dashtype 2**(i-1) ps 0.5 pt 1

set yrange [0:50]
set ylabel "Output Set Size"
set output 'length.eps'
plot for [i=1:4] 'length.txt' using (column(0)):i w linespoints notitle lt rgb colors[i] lw (i==4 ? 2 : 3) dashtype 2**(i-1) ps 0.5 pt 1

set yrange [0:1.1]
set ytics ('0' 0.0, '0.2' 0.2, '0.4' 0.4, '0.6' 0.6, '0.8' 0.8, '1.0' 1.0, 'Inf' 1.1)
set ylabel "Local Conductance" offset -3
set output 'lconductance.eps'
plot for [i=1:3] 'lconductance.txt' using (column(0)):i w linespoints notitle lt rgb colors[i] lw 3 dashtype 2**(i-1) ps 0.5 pt 1
