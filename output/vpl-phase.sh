#!/usr/bin/gnuplot
reset

# Terminal eps
set terminal postscript eps size 3.5,2.62 enhanced color \
dashed font 'Sans-Sheriff-Bold,18' lw 2

#set size square
set output 'phase.eps'

# Border style
set border linewidth 1.5

# Line styles
set style line 1 linecolor rgb '#0060ad' linetype 1 linewidth 1.5  # blue
set style line 2 linecolor rgb '#dd181f' linetype 1 linewidth 1.5  # red
set style line 3 linecolor rgb '#006400' linetype 1 linewidth 1.5  # green

# Legend
set key right top #at 6.1,1.3

# Axes label 
set xlabel 'u [m]'
set ylabel 'udot [m/s]' offset 1,0,0

# Axes type
#set logscale y

# Axis ranges
#set xrange[0:25]
#set yrange[1:1.0e-12]

# Axis labels
#set xtics ('-2{/Symbol p}' -2*pi, '-{/Symbol p}' -pi, 0, '{/Symbol p}' pi, \
#'2{/Symbol p}' 2*pi)
#set format y "10^{%.2L}"
#set ytics 1.0,0.01
#set tics scale 0.75
# Functions to plot
#a = 0.9
#f(x) = a * sin(x)
#g(x) = a * cos(x)

# Plot
#plot f(x) title 'sin({/Helvetica-Oblique x})' with lines ls 1, \
#     g(x) notitle with lines ls 2
plot 'abm.dat' u 2:4 with lines ls 2 title '1', \
     'abm.dat' u 3:5 with lines ls 3 title '2'
