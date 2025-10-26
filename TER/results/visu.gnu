set xlabel 'time'
set ylabel 'T'
plot 'sensor1' u 1:4 w l, 'sensor2' u 1:4 w l, 'sensor3' u 1:4 w l, 'sensor4' u 1:4 w l, 'sensor5' u 1:4 w l
set term postscript color
set output 'temperature.eps'
rep
set output 'density.eps'
set ylabel 'Density'
plot 'sensor1' u 1:5 w l, 'sensor2' u 1:5 w l, 'sensor3' u 1:5 w l, 'sensor4' u 1:5 w l, 'sensor5' u 1:5 w l
