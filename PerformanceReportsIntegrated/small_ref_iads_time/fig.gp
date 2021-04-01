reset
set terminal postscript eps color "Times-Roman" 18 enhanced
set size 2,0.66

set style data histogram
set style histogram cluster gap 2
set style fill solid border 
set style fill pattern
set boxwidth 1

set key top horizontal center outside
set key height 1
set key width 5


unset mytics


unset xlabel
unset ylabel

set xtics rotate by -30 font "Times-Roman,18"

set yrange [0:1.2]


set output 'fig.eps'

plot 'fig.txt' using 2:xtic(1) title '% Nodes with Degree < vol(R)' lt -1 fs pattern 8,\
	'' using 3 title 'Time(IADS) / Time(ADS)' lt -1 fs pattern 4
