all: plot_sol.pdf plot_eigen.pdf

%.x : %.cpp
	g++ -o3 $< -o $@

datos.txt : Tiempo_matrix.x Tiempo_matrix.cpp
	./$< > $@

plot_sol.pdf : plot_sol.gp datos.txt
	gnuplot $<
	xpdf $@

plot_eigen.pdf : plot_eigen.gp datos.txt
	gnuplot $<
	xpdf $@

.PHONY: clean
clean:
	rm -f *.x *.txt *.pdf
