# In windows, still need to open gnuplot, File > Open, filter type to all files and choose this script to work. Drag and drop etc doesn't do.

folders = system("dir /AD /B")

do for [folder in folders] {
	cd folder
	load "fig.gp"
	cd ".."
}

