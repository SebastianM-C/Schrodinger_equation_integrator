cd "../Bin"
set terminal gif animate delay 2
set output 'output.gif'
stats 'output.dat' nooutput
set xrange [-20:20]
set yrange [0:0.9]

do for [i=1:int(STATS_blocks)] {
	plot 'output.dat' index (i-1) with lines title "|Psi|^2"
}
