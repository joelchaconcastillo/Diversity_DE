set term post color 18
set view 80, 40
set yrange[0:1]
set xlabel "Niche radius"
set ylabel "Average success rate"
set title "Analyzes of Initial Niche Radius"
set output "Tuning_CEC.eps"
plot "cec2016" w linespoints ls 1 lc rgb 'blue' title "CEC 2016", "cec2017" w linespoints ls 1 title "CEC 2017"
