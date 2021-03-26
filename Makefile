all : plot.pdf

sin.x : sin_pade.cpp
	g++ $< -o $@

datos.txt : sin.x sin_pade.cpp
	./$< > $@

plot.pdf : plot.gp datos.txt
	gnuplot $<
	xpdf $@

.PHONY: clean
clean:
	rm -f *.x *.txt
