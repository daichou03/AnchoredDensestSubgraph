reset
set terminal postscript eps color "Times-Roman" 18 enhanced
set size 2,0.125

set style data histogram
set style histogram cluster gap 2
set style fill solid border 
set style fill pattern
set boxwidth 1

set key top horizontal center outside
set key height 1
set key width 5

set noborder
set noxtics
set noytics
set notitle
set noxlabel
set noylabel
set xrange[0:1]
set yrange[0:1]

set output 'fig.eps'
plot 	'fig.txt' using 2:xtic(1) title '1' lt -1 fs pattern 8,\
	'' using 3 title '1/2' lt -1 fs pattern 4,\
	'' using 4 title '1/4' lt -1 fs pattern 2,\
	'' using 5 title '1/8' lt -1 fs pattern 3,\
	'' using 6 title '1/16' lt -1 fs pattern 7,\
	'' using 7 title '1/32' lt -1 fs pattern 5
