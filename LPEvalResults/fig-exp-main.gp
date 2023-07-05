reset
set terminal postscript eps color "Times-Roman" 12 enhanced
set size 1,0.3

set datafile separator ","
set style data histogram
set style histogram cluster gap 2
set style fill solid border -1
set boxwidth 0.9

unset xlabel
unset ylabel
unset key
set xtics rotate by -30
set xtics format ""
set xtics nomirror
set ytics nomirror

set xrange [0.4:18.6]
set grid y

y_max = 40
set yrange [0:y_max]
set output 'rhoDS.eps'
plot 'masked-smartL-ρDS' using 2:xtic(1), \
     for [i=3:8] '' using i, \
     for [i=2:8] '' using ($0-0.55+i*0.11):(column(i) > y_max ? y_max : 1/0):(sprintf('%d', column(i))) with labels offset char 0,0.5 font "Times-Roman,8" rotate by -30 notitle

set yrange [0:10]
set output 'rhoADS.eps'
plot 'masked-smartL-ρADS' using 2:xtic(1), \
     for [i=3:8] '' using i

set output 'rhoADSL.eps'
plot 'masked-smartL-ρADSL' using 2:xtic(1), \
     for [i=3:8] '' using i

set output 'rhoADSLS.eps'
plot 'masked-smartL-ρADSLS' using 2:xtic(1), \
     for [i=3:8] '' using i

set yrange [0:10]
set output 'rhoADSF.eps'
plot 'masked-smartL-ρADSF' using 2:xtic(1), \
     for [i=3:8] '' using i

set yrange [0:15]
set output 'rhoADSFS.eps'
plot 'masked-smartL-ρADSFS' using 2:xtic(1), \
     for [i=3:8] '' using i

y_max = 40
set yrange [0:y_max]
set output 'rhoADSI.eps'
plot 'masked-smartL-ρADSI' using 2:xtic(1), \
     for [i=3:8] '' using i, \
     for [i=2:8] '' using ($0-0.55+i*0.11):(column(i) > y_max ? y_max : 1/0):(sprintf('%d', column(i))) with labels offset char 0,0.5 font "Times-Roman,8" rotate by -30 notitle

set output 'rhoADSIS.eps'
plot 'masked-smartL-ρADSIS' using 2:xtic(1), \
     for [i=3:8] '' using i, \
     for [i=2:8] '' using ($0-0.55+i*0.11):(column(i) > y_max ? y_max : 1/0):(sprintf('%d', column(i))) with labels offset char 0,0.5 font "Times-Roman,8" rotate by -30 notitle


set yrange [0:1]
set output 'f1score.eps'
plot 'masked-smartL-f1score' using 2:xtic(1), \
     for [i=3:8] '' using i

y_max = 100
set yrange [0:y_max]
set output 'ssize.eps'
plot 'masked-smartL-ssize' using 2:xtic(1), \
     for [i=3:8] '' using i, \
     for [i=2:8] '' using ($0-0.55+i*0.11):(column(i) > y_max ? y_max : 1/0):(sprintf('%d', column(i))) with labels offset char 0,0.5 font "Times-Roman,8" rotate by -30 notitle

set yrange [*:*]
set output 'iters.eps'
plot 'masked-smartL-iters' using 2:xtic(1), \
     for [i=3:8] '' using i

# density
# cond?

set logscale y
set yrange [10:200000]
set ytics ('10^{1}' 10, '10^{2}' 100, '10^{3}' 1000, '10^{4}' 10000, '10^{5}' 100000)
set output 'lnsize.eps'
plot 'masked-smartL-lnsize' using 2:xtic(1), \
     for [i=3:8] '' using i

set yrange [10:2000000]
set ytics ('10^{1}' 10, '10^{2}' 100, '10^{3}' 1000, '10^{4}' 10000, '10^{5}' 100000, '10^{6}' 1000000)
set output 'lmsize.eps'
plot 'masked-smartL-lmsize' using 2:xtic(1), \
     for [i=3:8] '' using i


set yrange [0.001:300]
set ytics ('10^{-3}' 0.001, '10^{-2}' 0.01, '10^{-1}' 0.1, '10^{0}' 1, '10^{1}' 10, '10^{2}' 100, 'Inf' 300)
set output 'exttime.eps'
plot 'masked-smartL-ext_time' using 2:xtic(1), \
     for [i=3:8] '' using i

set yrange [0.001:300]
set ytics ('10^{-3}' 0.001, '10^{-2}' 0.01, '10^{-1}' 0.1, '10^{0}' 1, '10^{1}' 10, '10^{2}' 100, 'Inf' 300)
set output 'inttime.eps'
plot 'masked-smartL-int_time' using 2:xtic(1), \
     for [i=3:8] '' using i


# Legend only
set key at screen 0.5, screen 0.5 center center horizontal spacing 0.5 font "Times-Roman,12"
set noxtics
set noytics
set noborder

unset title
unset grid

set size 1,1

set output 'legend.eps'
plot NaN title "ADS" with boxes, \
     NaN title "LPA" with boxes, \
     NaN title "LNA" with boxes, \
     NaN title "LIA" with boxes, \
     NaN title "LPA+" with boxes, \
     NaN title "LNA+" with boxes, \
     NaN title "LIA+" with boxes