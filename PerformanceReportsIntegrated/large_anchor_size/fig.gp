reset
set terminal postscript eps color "Times-Roman" 18 enhanced
set size 1,0.65

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
set ylabel "Time Consumption (sec)" offset 0.5 font "Times-Roman,18"

set logscale y
set format y "10^{%L}

set ytics ('10^{-1}' 0.1, '10^{0}' 1, '10^{1}' 10, '10^{2}' 100, 'Inf' 360)

set yrange [0.1:360]


set output 'time.eps'
plot 	'time.txt' using 2:xtic(1) title 'GA' lt -1 fs pattern 8,\
	'' using 3 title 'IGA' lt -1 fs pattern 4,\
	'' using 4 title 'LA' lt -1 fs pattern 2
