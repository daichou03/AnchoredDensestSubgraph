reset
set terminal postscript eps color "Times-Roman" 20 enhanced
set size 1,0.67

set datafile separator ","
set style data histogram
set style histogram cluster gap 1
set style fill solid border -1
set boxwidth 0.9
set key outside

set xlabel "Data Name"
set ylabel "Values"
set xtics rotate by -45
set xtics format ""
set xtics nomirror
set ytics nomirror

set grid y

set yrange [*:*]
set output 'rhoDS.eps'
plot 'adjusted-18graph-ρDS' using 2:xtic(1) title "FN100", \
     '' using 3 title "ADSL100C", \
     '' using 4 title "ADSF100C", \
     '' using 5 title "ADSI100C", \
     '' using 6 title "ADSLS100C", \
     '' using 7 title "ADSFS100C", \
     '' using 8 title "ADSIS100C"

set output 'rhoADS.eps'
plot 'adjusted-18graph-ρADS' using 2:xtic(1) title "FN100", \
     '' using 3 title "ADSL100C", \
     '' using 4 title "ADSF100C", \
     '' using 5 title "ADSI100C", \
     '' using 6 title "ADSLS100C", \
     '' using 7 title "ADSFS100C", \
     '' using 8 title "ADSIS100C"

set output 'rhoADSL.eps'
plot 'adjusted-18graph-ρADSL' using 2:xtic(1) title "FN100", \
     '' using 3 title "ADSL100C", \
     '' using 4 title "ADSF100C", \
     '' using 5 title "ADSI100C", \
     '' using 6 title "ADSLS100C", \
     '' using 7 title "ADSFS100C", \
     '' using 8 title "ADSIS100C"

set output 'rhoADSLS.eps'
plot 'adjusted-18graph-ρADSLS' using 2:xtic(1) title "FN100", \
     '' using 3 title "ADSL100C", \
     '' using 4 title "ADSF100C", \
     '' using 5 title "ADSI100C", \
     '' using 6 title "ADSLS100C", \
     '' using 7 title "ADSFS100C", \
     '' using 8 title "ADSIS100C"

set output 'rhoADSF.eps'
plot 'adjusted-18graph-ρADSF' using 2:xtic(1) title "FN100", \
     '' using 3 title "ADSL100C", \
     '' using 4 title "ADSF100C", \
     '' using 5 title "ADSI100C", \
     '' using 6 title "ADSLS100C", \
     '' using 7 title "ADSFS100C", \
     '' using 8 title "ADSIS100C"

set output 'rhoADSFS.eps'
plot 'adjusted-18graph-ρADSFS' using 2:xtic(1) title "FN100", \
     '' using 3 title "ADSL100C", \
     '' using 4 title "ADSF100C", \
     '' using 5 title "ADSI100C", \
     '' using 6 title "ADSLS100C", \
     '' using 7 title "ADSFS100C", \
     '' using 8 title "ADSIS100C"

set output 'rhoADSI.eps'
plot 'adjusted-18graph-ρADSI' using 2:xtic(1) title "FN100", \
     '' using 3 title "ADSL100C", \
     '' using 4 title "ADSF100C", \
     '' using 5 title "ADSI100C", \
     '' using 6 title "ADSLS100C", \
     '' using 7 title "ADSFS100C", \
     '' using 8 title "ADSIS100C"

set output 'rhoADSIS.eps'
plot 'adjusted-18graph-ρADSIS' using 2:xtic(1) title "FN100", \
     '' using 3 title "ADSL100C", \
     '' using 4 title "ADSF100C", \
     '' using 5 title "ADSI100C", \
     '' using 6 title "ADSLS100C", \
     '' using 7 title "ADSFS100C", \
     '' using 8 title "ADSIS100C"


set yrange [0:1]
set output 'f1score.eps'
plot 'adjusted-18graph-f1score' using 2:xtic(1) title "FN100", \
     '' using 3 title "ADSL100C", \
     '' using 4 title "ADSF100C", \
     '' using 5 title "ADSI100C", \
     '' using 6 title "ADSLS100C", \
     '' using 7 title "ADSFS100C", \
     '' using 8 title "ADSIS100C"

set yrange [*:*]
set output 'ssize.eps'
plot 'adjusted-18graph-ssize' using 2:xtic(1) title "FN100", \
     '' using 3 title "ADSL100C", \
     '' using 4 title "ADSF100C", \
     '' using 5 title "ADSI100C", \
     '' using 6 title "ADSLS100C", \
     '' using 7 title "ADSFS100C", \
     '' using 8 title "ADSIS100C"

set yrange [*:*]
set output 'iters.eps'
plot 'adjusted-18graph-iters' using 2:xtic(1) title "FN100", \
     '' using 3 title "ADSL100C", \
     '' using 4 title "ADSF100C", \
     '' using 5 title "ADSI100C", \
     '' using 6 title "ADSLS100C", \
     '' using 7 title "ADSFS100C", \
     '' using 8 title "ADSIS100C"

# density
# cond?

set logscale y
set yrange [10:200000]
set ytics ('10^{1}' 10, '10^{2}' 100, '10^{3}' 1000, '10^{4}' 10000, '10^{5}' 100000)
set output 'lnsize.eps'
plot 'adjusted-18graph-lnsize' using 2:xtic(1) title "FN100", \
     '' using 3 title "ADSL100C", \
     '' using 4 title "ADSF100C", \
     '' using 5 title "ADSI100C", \
     '' using 6 title "ADSLS100C", \
     '' using 7 title "ADSFS100C", \
     '' using 8 title "ADSIS100C"

set yrange [10:5000000]
set ytics ('10^{1}' 10, '10^{2}' 100, '10^{3}' 1000, '10^{4}' 10000, '10^{5}' 100000, '10^{6}' 1000000)
set output 'lmsize.eps'
plot 'adjusted-18graph-lmsize' using 2:xtic(1) title "FN100", \
     '' using 3 title "ADSL100C", \
     '' using 4 title "ADSF100C", \
     '' using 5 title "ADSI100C", \
     '' using 6 title "ADSLS100C", \
     '' using 7 title "ADSFS100C", \
     '' using 8 title "ADSIS100C"


set yrange [0.001:1000]
set ytics ('10^{-3}' 0.001, '10^{-2}' 0.01, '10^{-1}' 0.1, '10^{0}' 1, '10^{1}' 10, '10^{2}' 100, '10^{3}' 1000)
set output 'exttime.eps'
plot 'average-18graph-ext_time' using 2:xtic(1) title "FN100", \
     '' using 3 title "ADSL100C", \
     '' using 4 title "ADSF100C", \
     '' using 5 title "ADSI100C", \
     '' using 6 title "ADSLS100C", \
     '' using 7 title "ADSFS100C", \
     '' using 8 title "ADSIS100C"

set yrange [0.001:1000]
set ytics ('10^{-3}' 0.001, '10^{-2}' 0.01, '10^{-1}' 0.1, '10^{0}' 1, '10^{1}' 10, '10^{2}' 100, '10^{3}' 1000)
set output 'inttime.eps'
plot 'average-18graph-int_time' using 2:xtic(1) title "FN100", \
     '' using 3 title "ADSL100C", \
     '' using 4 title "ADSF100C", \
     '' using 5 title "ADSI100C", \
     '' using 6 title "ADSLS100C", \
     '' using 7 title "ADSFS100C", \
     '' using 8 title "ADSIS100C"