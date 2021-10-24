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

set encoding utf8

unset xlabel
set xtics rotate by -30 font "Times-Roman,18"

set yrange [0:10]
set output 'density.eps'
plot 'density.txt' using 2:xtic(1) notitle lt -1 fs pattern 8,\
	'' using 3 notitle lt -1 fs pattern 4,\
	'' using 4 notitle lt -1 fs pattern 2


set yrange [0:10]
set output 'rsdensity.eps'
plot 'rsdensity.txt' using 2:xtic(1) notitle lt -1 fs pattern 8,\
	'' using 3 notitle lt -1 fs pattern 4,\
	'' using 4 notitle lt -1 fs pattern 2


set yrange [0:1]
set output 'conductance.eps'
plot 'conductance.txt' using 2:xtic(1) notitle lt -1 fs pattern 8,\
	'' using 3 notitle lt -1 fs pattern 4,\
	'' using 4 notitle lt -1 fs pattern 2


set yrange [0:1]
set output 'rins.eps'
plot 'rins.txt' using 2:xtic(1) notitle lt -1 fs pattern 8,\
	'' using 3 notitle lt -1 fs pattern 4,\
	'' using 4 notitle lt -1 fs pattern 2


set yrange [0:1]
set output 'sinr.eps'
plot 'sinr.txt' using 2:xtic(1) notitle lt -1 fs pattern 8,\
	'' using 3 notitle lt -1 fs pattern 4,\
	'' using 4 notitle lt -1 fs pattern 2


set yrange [0:1]
set output 'f1score.eps'
plot 'f1score.txt' using 2:xtic(1) notitle lt -1 fs pattern 8,\
	'' using 3 notitle lt -1 fs pattern 4,\
	'' using 4 notitle lt -1 fs pattern 2


set yrange [0:60]
set output 'length.eps'
plot 'length.txt' using 2:xtic(1) notitle lt -1 fs pattern 8,\
	'' using 3 notitle lt -1 fs pattern 4,\
	'' using 4 notitle lt -1 fs pattern 2,\
    '' using 5 notitle lt -1 fs pattern 3


# set yrange [-1:10]
# set ytics ('' -1, '0' 0, '2' 2, '4' 4, '6' 6, '8' 8, '10' 10)
# set ylabel "R-Subgraph Density"
# set output 'rsdensity.eps'
# plot 'rsdensity.txt' using 2:xtic(1) notitle lt -1 fs pattern 8,\
# 	'' using 3 notitle lt -1 fs pattern 4,\
# 	'' using 4 notitle lt -1 fs pattern 2


set yrange [0:1.1]
set ytics ('0' 0.0, '0.2' 0.2, '0.4' 0.4, '0.6' 0.6, '0.8' 0.8, '1.0' 1.0, '' 1.1)
set output 'lconductance.eps'
plot 'lconductance.txt' using 2:xtic(1) notitle lt -1 fs pattern 8,\
	'' using 3 notitle lt -1 fs pattern 4,\
	'' using 4 notitle lt -1 fs pattern 2


set logscale y
set format y "10^{%L}
set yrange [0.001:100]
set ytics ('10^{-3}' 0.001, '10^{-2}' 0.01, '10^{-1}' 0.1, '10^{0}' 1, '10^{1}' 10, '10^{2}' 100)
set output 'time.eps'
plot 'time.txt' using 2:xtic(1) notitle lt -1 fs pattern 8,\
	'' using 3 notitle lt -1 fs pattern 4,\
	'' using 4 notitle lt -1 fs pattern 2