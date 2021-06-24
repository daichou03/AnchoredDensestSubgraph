reset
set terminal postscript eps color "Times-Roman" 18 enhanced
set size 2,0.125

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

set output 'fig.eps'
plot 	'fig.txt' using 2:xtic(1) w linespoints title 'LA' lt 1 lc rgb "#CC0000" dashtype 1 ps 0.5 pt 1,\
	'' using 3 w linespoints title 'GL' lt 1 lc rgb "#00CC00" dashtype 2 ps 0.5 pt 1,\
	'' using 4 w linespoints title 'FS' lt 1 lc rgb "#0000CC" dashtype 4 ps 0.5 pt 1
