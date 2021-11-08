reset
set terminal postscript eps color "Times-Roman" 20 enhanced
set size 1,0.67

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
set xtics rotate by -30 font "Times-Roman,18" nomirror

set yrange [0:12]
set output 'density.eps'
plot 'density.txt' using 2:xtic(1) notitle lt -1 fs pattern 8,\
	'' using 3 notitle lt -1 fs pattern 4,\
	'' using 4 notitle lt -1 fs pattern 2,\
	'' using 0:($2 > 12 ? 12 : 1/0):(sprintf('%.2f', $2)) with labels font "Times-Roman,12" offset -0.3,0.5 title " "


set yrange [0:12]
set output 'rsdensity.eps'
plot 'rsdensity.txt' using 2:xtic(1) notitle lt -1 fs pattern 8,\
	'' using 3 notitle lt -1 fs pattern 4,\
	'' using 4 notitle lt -1 fs pattern 2,\
	'' using 0:($2 > 12 ? 12 : 1/0):(sprintf('%.2f', $2)) with labels font "Times-Roman,12" offset -0.3,0.5 title " "


set yrange [0:1]
set output 'f1score.eps'
plot 'f1score.txt' using 2:xtic(1) notitle lt -1 fs pattern 8,\
	'' using 3 notitle lt -1 fs pattern 4,\
	'' using 4 notitle lt -1 fs pattern 2,\
	'' using 0:($2 > 1 ? 1 : 1/0):(sprintf('%.2f', $2)) with labels font "Times-Roman,12" offset 0,0.5 title " "


set yrange [0:60]
set output 'length.eps'
plot 'length.txt' using 2:xtic(1) notitle lt -1 fs pattern 8,\
	'' using 3 notitle lt -1 fs pattern 4,\
	'' using 4 notitle lt -1 fs pattern 2,\
    '' using 5 notitle lt -1 fs pattern 3,\
	'' using 0:($3 > 60 ? 60 : 1/0):(sprintf('%.1f', $3)) with labels font "Times-Roman,12" offset 0.05,0.5 title " "


set yrange [0:1.1]
set ytics ('0' 0.0, '0.2' 0.2, '0.4' 0.4, '0.6' 0.6, '0.8' 0.8, '1.0' 1.0, '' 1.1)
set output 'conductance.eps'
plot 'conductance.txt' using 2:xtic(1) notitle lt -1 fs pattern 8,\
	'' using 3 notitle lt -1 fs pattern 4,\
	'' using 4 notitle lt -1 fs pattern 2,\
	'' using 0:($4 > 1.1 ? 1.1 : 1/0):(sprintf('%.3f', $4)) with labels font "Times-Roman,12" offset 0.15,0.5 title " "


set yrange [0:1.1]
set ytics ('0' 0.0, '0.2' 0.2, '0.4' 0.4, '0.6' 0.6, '0.8' 0.8, '1.0' 1.0, '' 1.1)
set output 'lconductance.eps'
plot 'lconductance.txt' using 2:xtic(1) notitle lt -1 fs pattern 8,\
	'' using 3 notitle lt -1 fs pattern 4,\
	'' using 4 notitle lt -1 fs pattern 2,\
	'' using 0:($4 > 1.1 ? 1.1 : 1/0):(sprintf('%.3f', $4)) with labels font "Times-Roman,12" offset 0.15,0.5 title " "


set logscale y
set format y "10^{%L}
set yrange [0.001:2500]
set ytics ('10^{-3}' 0.001, '10^{-2}' 0.01, '10^{-1}' 0.1, '10^{0}' 1, '10^{1}' 10, '10^{2}' 100, '10^{3}' 1000)
set output 'time.eps'
plot 'time.txt' using 2:xtic(1) notitle lt -1 fs pattern 8,\
	'' using 3 notitle lt -1 fs pattern 4,\
	'' using 4 notitle lt -1 fs pattern 2


set logscale y
set format y "10^{%L}
set yrange [0.01:100000]
set ytics ('10^{-2}' 0.01, '10^{-1}' 0.1, '10^{0}' 1, '10^{1}' 10, '10^{2}' 100, '10^{3}' 1000, '10^{4}' 10000, '10^{5}' 100000)
set output 'space.eps'
plot 'space.txt' using 2:xtic(1) notitle lt -1 fs pattern 8,\
	'' using 3 notitle lt -1 fs pattern 4,\
	'' using 4 notitle lt -1 fs pattern 2