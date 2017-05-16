set term postscript enhanced color eps
set output "dosloc.eps"
set encoding iso

bint = 2*9.40260
b =  45.4026099999999957
Elo = -8
Ehi = 10
set size 0.7, 0.7
set xrange [0:b+bint]
set yrange [Elo:Ehi]
#set contour
#set title "Density of states"
set xlabel "Distance (\305)"
set ylabel "Energy (eV)"
unset key

set tmargin at screen 0.67
set bmargin at screen 0.10
set lmargin at screen 0.08
set rmargin at screen 0.61

#set tmargin at screen 0.62
#set bmargin at screen 0.10
#set lmargin at screen 0.08
#set rmargin at screen 0.61


set pm3d
set view map
unset surface
#set palette model RGB defined\
# ( 0 '#3366ff', 1 '#99ffcc', 2 '#339900', 3 '#66ff33',\
#   4 '#996633', 5 '#ff9900', 6 '#ffff33' )

load "./moreland.pal"

set arrow from bint, Elo, 0 to bint, Ehi, 0 nohead front lw 1.5
set arrow from b,    Elo, 0 to b,    Ehi, 0 nohead front lw 1.5

splot "dosloc-out" w lines
