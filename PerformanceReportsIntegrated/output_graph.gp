# In windows, still need to open gnuplot, File > Open, filter type to all files and choose this script to work. Drag and drop etc doesn't do.
# 20210328: Double click this script can do too.

folders = system("dir /AD /B")

do for [folder in folders] {
	cd folder
	print folder
	files = system("dir /A /B")
	do for [file in files] {
		if (file eq "fig.gp") {
			load "fig.gp"
		}
	}
	cd ".."
}
