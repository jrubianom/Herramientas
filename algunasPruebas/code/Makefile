#Este makefile se complementa con el run.sh para automatizar la toma de datos

#MODO DE USO:
#Antes de todo, es necesario cargar la libreria gsl con spack (spack load gsl)
#Primero, se debe compilar el programa que se va a medir (e.g., transpuesta_eigen.cpp), de tal forma quel el binario se llame a.out
#Luego, se especifica en el run.sh el numero de filas de la matriz (e.g., 2048), así como el numero de veces que se desea correr el programa (e.g., 10)
#make
#Se crearan varios archivos de texto. El importante será datos_prom-desv.txt, donde estara el promedio de tiempo, mflops y sus desviaciones estandar
#Si se va a medir otro programa, solo sera necesario compilarlo y correr make, no hace falta eliminar los archivos de texto creados.

all: datos_stat.txt


read_file.x: read_file.cpp
	g++ read_file.cpp -o read_file.x -lgsl -lgslcblas


datos_mflops.txt: run.sh a.out read_file.x
	sh run.sh


datos_stat.txt: read_file.cpp datos_mflops.txt read_file.x
	./read_file.x > datos_prom-desv.txt
	echo "Toma de datos para el programa compilado terminada. Revisa datos_prom-desv.txt"
	cat datos_prom-desv.txt

.PHONY: clean_all clean_data


clean_data:
	rm *.txt


clean_all:
	rm *.txt *.x *.out
