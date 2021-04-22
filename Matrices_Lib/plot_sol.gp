set log y
set log x
set xlabel 'tamaño de la matriz N'
set ylabel 'tiempo solución t[s]'
set term pdf
set out 'plot_sol.pdf'
plot 'datos.txt' u 1:2 w lp t 'Tiempos para la solución de un sistema lineal'
