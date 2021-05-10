#Este script se complementa con el Makefile para automatizar la toma de datos

#MODO DE USO:
#Antes de todo, es necesario cargar la libreria gsl con spack (spack load gsl)
#Primero, se debe compilar el programa que se va a medir (e.g., transpuesta_eigen.cpp), de tal forma quel el binario se llame a.out
#Luego, se especifica en este archivo run.sh el numero de filas de la matriz (e.g., 2048), así como el numero de veces que se desea correr el programa (e.g., 10)
#make
#Se crearan varios archivos de texto. El importante será datos_prom-desv.txt, donde estara el promedio de tiempo, mflops y sus desviaciones estandar
#Si se va a medir otro programa, solo sera necesario compilarlo y correr make, no hace falta eliminar los archivos de texto creados.


#inputs
N_matriz=4
N_blocking=2
n_iteraciones=1

i=0
> datos_brutos.txt && > datos_linea.txt && > datos_mflops.txt && > datos_prom-desv.txt
while [ $i -lt $n_iteraciones ]
do
	./a.out $N_matriz $N_blocking>> datos_brutos.txt
	((i=i+1))
done

grep -i mflops datos_brutos.txt > datos_linea.txt
cat datos_linea.txt | awk '{print $2, $9}' > datos_mflops.txt
