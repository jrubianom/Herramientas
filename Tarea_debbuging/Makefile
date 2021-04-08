all: do

log.txt : a.out array_bounds.cpp
	valgrind --tool=memcheck ./$< &>$@
a.out : array_bounds.cpp
	g++ -g $<
do: log.txt
	less -- $<
