
all: Histogramas.pdf

main.x: main.cpp Motor.cpp vector_operations.cpp Event.cpp
	g++ -O3 -O0 -fopenmp $^ -o $@

velocities_data.txt: main.x
	./$< $(N) $(n)

Histogramas.pdf: hist.py velocities_data.txt
	python $^
	xpdf -open $@ &

.PHONY: clean
clean:
	rm *.x *.txt
