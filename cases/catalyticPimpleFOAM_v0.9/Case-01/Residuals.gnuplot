set logscale y
set title "Residuals"
set ylabel 'Residual'
set xlabel 'Iteration'
plot "< cat log | grep 'Solving for p,' | cut -d' ' -f9 | tr -d ','" title 'p' with lines,\
     "< cat log | grep 'Solving for H2,' | cut -d' ' -f9 | tr -d ','" title 'H2' with lines,\
     "< cat log | grep 'Solving for O2,' | cut -d' ' -f9 | tr -d ','" title 'O2' with lines,\
     "< cat log | grep 'Solving for T,' | cut -d' ' -f9 | tr -d ','" title 'T' with lines
pause 5
reread
