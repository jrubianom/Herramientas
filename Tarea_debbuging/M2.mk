all: do

log.txt : a.out ARRAY.cpp
	valgrind --tool=memcheck ./$< &>$@
a.out : ARRAY.cpp
	g++ -g $<
do: log.txt
	less -- $<
