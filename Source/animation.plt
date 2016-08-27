cd "../Bin"
set terminal gif animate delay 1
set output 'output.gif'
stats 'output.dat' nooutput
set xrange [-10:10]
set yrange [0:0.6]

do for [i=1:int(STATS_blocks)] {
	plot 'output.dat' index (i-1) with lines
}
