#set term epslatex color
#set terminal postscript eps enhanced 
##set term post color 18
set terminal postscript eps enhanced color 18
#set term epslatex color
set view 60, 40
set yrange[0:1]
set xlabel 'Initial distance factor (D_{I})' font "Times-Roman, 28"
set ylabel "Average success rate" font "Times-Roman, 28"
#set title "Analyzes of Initial Niche Radius"
set output "Tuning_CEC.eps"
plot "cec2016" w linespoints ls 1 lc rgb 'blue' title "CEC 2016", "cec2017" w linespoints ls 1 title "CEC 2017"
