#!/usr/bin/gnuplot
reset

# Terminal eps
set terminal postscript eps size 3.5,2.62 enhanced color \
dashed font 'Sans-Sheriff-Bold,18' lw 2

#set size square
set output 'BDFMeshRefinementVsAdjoint.eps'

# Border style
set border linewidth 1.5

# Line styles
set style line 1 linecolor rgb '#0060ad' linetype 1 linewidth 1.5  # blue
set style line 2 linecolor rgb '#dd181f' linetype 2 linewidth 1.5  # red
set style line 3 linecolor rgb '#006400' linetype 3 linewidth 1.5  # green

# Legend
set key right top #at 6.1,1.3

# Axes label 
set xlabel 'Time [s]'
set ylabel 'Adjoint variable {/Symbol l_1}' offset 1,0,0

# Axes type
#set logscale y

# Axis ranges
set xrange[0:10]
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
plot 'bdf0.1.dat' u 1:5 with lines ls 1 title 'h=0.1', \
    'bdf0.01.dat' u 1:5 with lines ls 2 title 'h=0.01', \
    'bdf0.001.dat' u 1:5 with lines ls 3 title 'h=0.001',
