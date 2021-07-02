reset
set terminal postscript eps color "Times-Roman" 18 enhanced
set size 2,0.1

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
set xrange[1:2]
set yrange[1:2]

set output 'baseline.eps'
plot 	'dummy.txt' using 2:xtic(1) title 'GA' lt -1 fs pattern 8,\
	'' using 3 title 'IGA' lt -1 fs pattern 4,\
	'' using 4 title 'LA' lt -1 fs pattern 2

set output 'halfedge.eps'
plot 	'dummy.txt' using 2:xtic(1) title 'G_0' lt -1 fs pattern 8,\
	'' using 3 title 'G_1' lt -1 fs pattern 4,\
	'' using 4 title 'G_2' lt -1 fs pattern 2,\
	'' using 5 title 'G_3' lt -1 fs pattern 3,\
	'' using 6 title 'G_4' lt -1 fs pattern 7,\
	'' using 7 title 'G_5' lt -1 fs pattern 5

set output 'anchorsize.eps'
plot 	'dummy.txt' using 2:xtic(1) title '8' lt -1 fs pattern 8,\
	'' using 3 title '16' lt -1 fs pattern 4,\
	'' using 4 title '32' lt -1 fs pattern 2,\
	'' using 5 title '64' lt -1 fs pattern 3,\
	'' using 6 title '128' lt -1 fs pattern 5