set title "Shock radius"

set xlabel "log(time) [yr]"
set ylabel "log(r_shock) [pc]"

set key
set autoscale
#set logscale

plot "shock.dat" u 1:2 w lp t "numeric", \
     "shock.dat" u 1:3 w lp lt 2 t "Sedov"

unset logscale
