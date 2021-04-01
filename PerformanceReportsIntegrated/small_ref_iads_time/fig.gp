reset
set terminal postscript eps color "Times-Roman" 18 enhanced
set size 1,1

set style fill solid border 
set style fill pattern
set boxwidth 1

unset mytics

set key right bottom
set key box 1

set xlabel "% Nodes with Degree < vol(R)"
set ylabel "IADS Speed Up"


set output 'fig.eps'

set encoding utf8

set offset 0,0,0.05,0.25

set style data points

array datanames[6]
datanames[1] = "dblp"
datanames[2] = "github"
datanames[3] = "grqc"
datanames[4] = "hepph"
datanames[5] = "livemocha"
datanames[6] = "youtube"


plot for [i=1:6] 'fig.txt' using 4:($2 == i ? $5 : 1/0):(sprintf("%d", $3)) title datanames[i] lt -1 ps 2 with labels point pt i offset char 0,1
