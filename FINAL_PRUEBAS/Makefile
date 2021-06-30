all: velocities_files

main.x: main.cpp Motor.cpp vector_operations.cpp Event.cpp Particula.cpp
	g++ -O3 -O0 $^ -o $@

velocities_files : main.x
	./$< $(N) $(n) $(v)

.PHONY: clean
clean:
	rm *.x
