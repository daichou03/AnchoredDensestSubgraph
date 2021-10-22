reset
set terminal postscript eps color "Times-Roman" 18 enhanced
set size 2,0.65

set style fill solid 0.25 border -1
set style boxplot outliers pointtype 7
set style data boxplot

unset mytics

unset xlabel
set ylabel "Time Consumption (sec)" offset 0.5 font "Times-Roman,18"

set logscale y
set format y "10^{%L}

set ytics ('10^{-3}' 0.001, '10^{-2}' 0.01, '10^{-1}' 0.1, '10^{0}' 1, '10^{1}' 10, '10^{2}' 100)
set xtics rotate by -30 font "Times-Roman,18"

set yrange [0.001:100]


set xtics ('LA' 1, 'FS' 2, 'MRW' 3)

set output 'fig.eps'
plot for [i=1:3] 'fig.txt' using (i):i notitle
