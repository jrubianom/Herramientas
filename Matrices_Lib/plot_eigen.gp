set log y
set log x
set xlabel 'tamaño de la matriz N'
set ylabel 'tiempo solución t[s]'
set term pdf
set out 'plot_eigen.pdf'
plot 'datos.txt' u 1:3 w lp t 'Tiempos para hallar eigen valores y vectores'
