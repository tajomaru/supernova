set multiplot layout 2,2 rowsfirst
#
set title "Density"
unset key
plot "zeus_b.dat" u 1:2 w lp
#
set title "Pressure"
unset key
plot "zeus_b.dat" u 1:3 w lp
#
set title "Velocity"
unset key
plot "zeus_a.dat" u 1:2 w lp
#
set title "Specific energy"
unset key
plot "zeus_b.dat" u 1:($4/$2) w lp
#
unset multiplot
