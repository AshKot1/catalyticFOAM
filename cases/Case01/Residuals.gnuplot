set logscale y
set title "Residuals"
set ylabel 'Residual'
set xlabel 'Iteration'
plot "< cat log | grep 'Solving for Ux,' | cut -d' ' -f9 | tr -d ','" title 'Ux' with lines,\
"< cat log | grep 'Solving for Uy,' | cut -d' ' -f9 | tr -d ','" title 'Uy' with lines,\
"< cat log | grep 'Solving for H2,' | cut -d' ' -f9 | tr -d ','" title 'H2' with lines,\
"< cat log | grep 'Solving for O2,' | cut -d' ' -f9 | tr -d ','" title 'O2' with lines,\
"< cat log | grep 'Solving for H2O,' | cut -d' ' -f9 | tr -d ','" title 'H2O' with lines
pause 1
reread
