import argparse
import os
import folders
import MapTransposons

# Definir el directorio base
DIRECTORIO_BASE = "/mnt/c/Users/sebas/Escritorio/sacc/Tarea_final"
DIRECTORIO_GENOMAS = os.path.join(DIRECTORIO_BASE, "genomas")

# Configurar argparse para manejar argumentos
parser = argparse.ArgumentParser(description="Procesar archivos FASTA y GFF.")
parser.add_argument("-i", "--entrada", required=True, help="Ruta completa al archivo de entrada.")
parser.add_argument(
    "-m", "--modo", required=True, choices=["gn", "gp", "cp", "fasta", "gff"],
    help="Modo de procesamiento ('gn': histograma genoma, 'gp': porcentaje genoma, 'cp': por cromosoma, 'fasta': cambiar formato FASTA, 'gff': cambiar formato GFF)."
)
parser.add_argument("-n", "--nombre", required=True, help="Nombre del archivo de salida.")
parser.add_argument("-t", "--transposones", help="Archivo con información de los transposones en formato GFF.")
parser.add_argument("-b", "--bins", type=int, default=120, help="Número de divisiones (bins) para el histograma (por defecto: 120).")
parser.add_argument("-p", "--porcentaje", type=int, default=100, help="Porcentaje de transposones más grandes a considerar (1-100).")

# Parsear los argumentos
args = parser.parse_args()

# Validar el archivo de entrada
archivo_entrada = args.entrada
if not os.path.isfile(archivo_entrada):
    print(f"Error: El archivo de entrada '{archivo_entrada}' no existe.")
    exit(1)

# Crear la carpeta de salida dentro de 'genomas/' con el nombre especificado
directorio_salida = os.path.join(DIRECTORIO_GENOMAS, args.nombre)
os.makedirs(directorio_salida, exist_ok=True)

# Determinar el nombre del archivo de salida basado en el modo
archivo_salida = os.path.join(directorio_salida, args.nombre)

# Llamar a la función según el modo especificado
try:
    if args.modo == "gn":  # Histograma del genoma completo
        imagen_salida = os.path.join(directorio_salida, "histograma_transposones.png")
        MapTransposons.mapear_transposones_histograma(
            archivo_gff=args.transposones, 
            imagen_salida=imagen_salida, 
            num_bins=args.bins, 
            porcentaje=args.porcentaje
        )

    elif args.modo == "gp":  # Porcentaje de ocupación en el genoma completo
        MapTransposons.mapear_porcentaje_ocupacion(
            archivo_gff=args.transposones, 
            directorio_salida=directorio_salida, 
            numero_bins=args.bins, 
            porcentaje=args.porcentaje
        )

    elif args.modo == "cp":  # Porcentaje de ocupación por cromosoma
        archivo_tamanos_cromosomas = os.path.join(directorio_salida, "tamanos_cromosomas.sizes")

        # Verificar si el archivo .sizes existe
        if not os.path.exists(archivo_tamanos_cromosomas):
            #print(f"Archivo {archivo_tamanos_cromosomas} no encontrado. Generándolo ahora...")
            folders.procesar_fasta(
                input_fasta=archivo_entrada,
                output_fasta=os.path.join(directorio_salida, args.nombre + "_genoma.fasta"),
                modo="chr"
            )

        directorio_graficos = os.path.join(directorio_salida, "graficos_ocupacion_cromosomas")
        os.makedirs(directorio_graficos, exist_ok=True)

        MapTransposons.mapear_porcentaje_ocupacion_por_cromosoma(
            archivo_gff=args.transposones,
            archivo_tamanos_cromosomas=archivo_tamanos_cromosomas,
            directorio_salida=directorio_graficos,
            numero_bins=args.bins,
            porcentaje=args.porcentaje
        )

    elif args.modo == "fasta":  # Cambiar el formato del archivo FASTA
        folders.cambiar_formato_fasta(
            input_fasta=archivo_entrada, 
            output_fasta=archivo_salida + ".fasta"
        )
        print(f"Formato del archivo FASTA cambiado y guardado en: {archivo_salida}.fasta")

    elif args.modo == "gff":  # Cambiar el formato del archivo GFF
        folders.cambiar_formato_gff(
            input_gff=archivo_entrada, 
            output_gff=archivo_salida + ".gff"
        )
        print(f"Formato del archivo GFF cambiado y guardado en: {archivo_salida}.gff")

    else:
        print("Modo no válido.")

except ValueError as e:
    print(f"Error: {e}")
