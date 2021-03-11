reset
set terminal postscript eps color "Times-Roman" 18 enhanced
set size 2,0.65

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
set ylabel "Memory Consumption (byte)" offset 0.5 font "Times-Roman,18"

set logscale y
set format y "10^{%L}

set ytics ('10^{0}' 1, '10^{1}' 10, '10^{2}' 100, '10^{3}' 1000, '10^{4}' 10000, '10^{5}' 100000)
set xtics rotate by -30 font "Times-Roman,18"

set yrange [1:500000]


set output 'fig.eps'
plot 	'fig.txt' using 2:xtic(1) title 'AnchoredDS' lt -1 fs pattern 8,\
	'' using 3 title 'ImprovedAnchoredDS' lt -1 fs pattern 4,\
	'' using 4 title 'LocalAnchoredDS' lt -1 fs pattern 2