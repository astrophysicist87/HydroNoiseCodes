#!/usr/bin/gnuplot

reset

# epslatex
set terminal epslatex standalone color colortext
set output 'plot.tex'

# Line styles
set border linewidth 1.5
set style line 1 linecolor rgb '#0060ad' linetype 1 linewidth 2  # blue
set style line 2 linecolor rgb '#dd181f' linetype 1 linewidth 2  # red
# Legend
set key spacing 2
# Axes label 
set xlabel '$\Delta y$'
set ylabel '$\left< \delta R^2_l(\Delta y) \delta R^2_l(0) \right>/ \left<R^2_l(\Delta y) \right>$  (fm$^2$)'
# Files to plot
f='plot2.out'
# Plot
plot f w lp ls 1 t "", \
	0 lt -1 t ""
