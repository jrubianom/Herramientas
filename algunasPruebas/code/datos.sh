#Usar este script en caso de que se interrumpa el proceso y solo esté el archivo datos_brutos.txt


#inputs
N_matriz=8192
N_blocking=2
n_iteraciones=100

i=0

grep -i mflops datos_brutos.txt > datos_linea.txt
cat datos_linea.txt | awk '{print $2, $9}' > datos_mflops.txt
