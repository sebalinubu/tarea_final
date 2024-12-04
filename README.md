Guía de Uso: FastaFrags.py

IMPORTANTE:
-modificar DIRECTORIO_BASE = "/mnt/c/Users/sebas/Escritorio/sacc/Tarea_final" en FastaFrags.py
-link con archivos listos para probar el codigo: https://uccl0-my.sharepoint.com/:f:/g/personal/sfdonoso_uc_cl/EsWn35lXLKZBgixXDay4lZMBzzF1Def-jW6jt0jsODZDsQ?e=FJ1Wun

El script FastaFrags.py, ubicado en la carpeta Scripts, ofrece diversas funcionalidades para procesar archivos FASTA y GFF relacionados con transposones. A continuación, se describen sus capacidades y el formato para ejecutarlo correctamente.

El programa cuenta con las siguientes funcionalidades:

-fasta: Cambia el formato de un archivo FASTA para que sea compatible con el programa.
-gff: Modifica el formato de un archivo GFF para que pueda ser utilizado correctamente en el programa.
-gn (Genoma completo - Histograma): Genera un gráfico con el número de transposones por bin y muestra el tamaño del transposón más grande dentro de cada bin.
-gp (Genoma completo - Porcentaje ocupado): Crea un gráfico que muestra el porcentaje de ocupación de transposones por bin, indicando también el tamaño del transposón más grande en cada bin.
-cp (Por cromosoma): Genera un gráfico para cada cromosoma, mostrando el porcentaje de ocupación por bin. Además, genera un archivo GFF que contiene los transposones de cada cromosoma, organizados por bin y tamaño.
Para ejecutar el script, utiliza el siguiente comando en el terminal:
python FastaFrags.py -i <ruta_al_archivo.fasta> -m <modo> -n <nombre_salida> -t <ruta_al_archivo.gff> -b <bins> -p <porcentaje>

Argumentos:

-i: Ruta completa del archivo FASTA de entrada.
-m: Modo de ejecución. Puede ser uno de los siguientes:
fasta: Cambiar formato del archivo FASTA.
gff: Cambiar formato del archivo GFF.
gn: Generar histograma de transposones por bin.
gp: Crear gráfico del porcentaje de ocupación por bin.
cp: Generar gráficos y GFF por cromosoma.
-n: Nombre de la carpeta donde se generarán los resultados.
-t: Ruta completa del archivo GFF con información de los transposones.
-b: Número de bins para los gráficos (opcional; por defecto: 120).
-p: Porcentaje de transposones más grandes a considerar (valor entre 1 y 100, opcional; por defecto: 100).

Ejemplos:

Generar un gráfico de ocupación por bin para todo el genoma (gp):
python FastaFrags.py -i /ruta/genoma.fasta -m gp -n salida_gp -t /ruta/transposones.gff -b (numero) -p (numero) 

Crear gráficos por cromosoma (cp):
python FastaFrags.py -i /ruta/genoma.fasta -m cp -n salida_cp -t /ruta/transposones.gff -b (numero) -p (numero)

Cambiar el formato de un archivo FASTA:
python FastaFrags.py -i /ruta/genoma.fasta -m fasta -n salida_formateada

Cambiar el formato de un archivo GFF:
python FastaFrags.py -i /ruta/transposones.gff -m gff -n salida_formateada

Notas adicionales:

- Los resultados se guardarán en la carpeta Genomas, en una subcarpeta nombrada según el argumento -n. Para los modos gn, gp, y cp, se generarán gráficos y/o archivos GFF que documentan los transposones procesados, todo esto en una carpeta que inicia con la palabra "grafico".

- En la carpeta datos_prueba se encuentran los archivos de Arabidopsis y Saccharomyces cerevisiae para probar el programa. Además, en la carpeta cerevisae_formato_incorrecto se incluyen archivos con formatos incorrectos para probar los modos -fasta y -gff.

- Es importante que los archivos FASTA y GFF estén en formatos compatibles antes de ejecutarlos en los modos gn, gp o cp. Si no lo están, utiliza los modos -fasta y -gff para convertirlos. Los archivos generados en los modos -fasta y -gff se guardarán en la carpeta Genomas con el nombre especificado en -n.