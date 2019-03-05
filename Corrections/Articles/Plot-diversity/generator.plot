set term post color 18
set view 60, 40
set yrange[0:0.26]
set xlabel "Percentage of evaluations" font "Times-Roman, 28"
set ylabel "Mean distance" font "Times-Roman, 28"
set title "Diversity of the Elite vectors" font "Times-Roman, 28"
set output "Diversity_Elite.eps"
plot "2016/elite_f1.txt" w line ls 1 lc rgb 'blue' lw 4 title "f1 (CEC 2016)", "2016/elite_f30.txt" w line ls 2 lc rgb 'green' lw 4 title "f30 (CEC 2016)", "2017/elite_f1.txt" w line ls 3 lc rgb 'orange' lw 4 title "f1 (CEC 2017)", "2017/elite_f30.txt" w line ls 4 lc rgb 'red' lw 4 title "f30 (CEC 2017)"

set term post color 18
set view 60, 40
set yrange[0:0.26]
set xlabel "Percentage of evaluations" font "Times-Roman, 28"
set ylabel "Mean distance" font "Times-Roman, 28"
set title "Diversity of the Target vectors" font "Times-Roman, 28"
set output "Diversity_Target.eps"
plot "2016/trial_f1.txt" w line ls 1 lc rgb 'blue' lw 4 title "f1 (CEC 2016)", "2016/trial_f30.txt" w line ls 2 lc rgb 'green' lw 4 title "f30 (CEC 2016)", "2017/trial_f1.txt" w line ls 3 lc rgb 'orange' lw 4 title "f1 (CEC 2017)", "2017/trial_f30.txt" w line ls 4 lc rgb 'red' lw 4 title "f30 (CEC 2017)"
