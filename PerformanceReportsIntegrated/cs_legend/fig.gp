reset
set terminal postscript eps color "Times-Roman" 18 enhanced
set size 2,0.1

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
plot 0 w linespoints title 'LA' lt 1 lc rgb "#CC0000" lw 3 dashtype 1 ps 0.5 pt 1,\
	0 w linespoints title 'GL' lt 1 lc rgb "#00CC00" lw 3 dashtype 2 ps 0.5 pt 1,\
	0 w linespoints title 'FS' lt 1 lc rgb "#0000CC" lw 3 dashtype 4 ps 0.5 pt 1,\
	0 w linespoints title 'R' lt 1 lc rgb "#777777" lw 2 dashtype 8 ps 0.5 pt 1
