set multiplot layout 3,2 rowsfirst

set grid

a = "`tail -1 shock.dat | awk '{print $3}'`"

set arrow from first a, graph 0 to first a, graph 1 nohead lt 2

unset xlabel
unset ylabel

#
set title "Density [g/cm^3]"
unset key
set autoscale
plot "sn_b.dat" u 1:2 w lp
#
set title "Pressure [erg s^-1 cm^-3]"
unset key
set autoscale
plot "sn_b.dat" u 1:3 w lp
#
set title "Velocity [km/s]"
unset key
set autoscale
plot "sn_a.dat" u 1:2 w lp
#
set title "Specific energy [cm^2/s^2]"
unset key
set autoscale
plot "sn_b.dat" u 1:($4/$2) w lp
#
set title "Temperature [K]"
unset key
set autoscale
plot "sn_b.dat" u 1:5 w lp
#
load "shock.plt"
#
unset multiplot
unset arrow
