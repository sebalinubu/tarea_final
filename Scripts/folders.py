from Bio import SeqIO
import os

def procesar_fasta(input_fasta, output_fasta, modo):
    try:
        secuencias = list(SeqIO.parse(input_fasta, "fasta"))
    except Exception as e:
        raise ValueError(f"Error al leer el archivo FASTA: {e}")

    tamanos = []
    if modo == "scaffold":#esto no lo ocupo en el script, pero si lo elimino el codigo no funciona :()
        secuencias.sort(key=lambda x: len(x.seq), reverse=True)
        for i, record in enumerate(secuencias, start=1):
            record.id = f"scaffold{i}"
            tamanos.append((record.id, len(record.seq)))
            record.description = ""
    elif modo == "chr":
        for i, record in enumerate(secuencias, start=1):
            record.id = f"chr{i}"
            tamanos.append((record.id, len(record.seq)))
            record.description = ""
    else:
        raise ValueError("Modo inválido. Use 'chr' o 'scaffold'.")

    output_dir = os.path.dirname(output_fasta)
    tamanos_file = os.path.join(output_dir, "tamanos_cromosomas.sizes")

    with open(tamanos_file, "w") as file:
        for nombre, tamano in tamanos:
            file.write(f"{nombre}\t{tamano}\n")

    for record in secuencias:
        chrom_dir = os.path.join(output_dir, record.id)
        os.makedirs(chrom_dir, exist_ok=True)
        chrom_fasta = os.path.join(chrom_dir, f"{record.id}.fasta")
        SeqIO.write([record], chrom_fasta, "fasta")

def cambiar_formato_fasta(input_fasta, output_fasta):
    try:
        secuencias = list(SeqIO.parse(input_fasta, "fasta"))
        for i, record in enumerate(secuencias, start=1):
            record.id = f"chr{i}"
            record.description = ""
        SeqIO.write(secuencias, output_fasta, "fasta")
    except Exception as e:
        raise ValueError(f"Error al procesar el archivo FASTA: {e}")

def cambiar_formato_gff(input_gff, output_gff):
    cromosomas = {}
    contador = 1
    try:
        with open(input_gff, "r") as entrada, open(output_gff, "w") as salida:
            for linea in entrada:
                if linea.startswith("##sequence-region"):
                    columnas = linea.strip().split()
                    nombre_original = columnas[1]
                    nuevo_nombre = f"chr{contador}"
                    cromosomas[nombre_original] = nuevo_nombre
                    contador += 1
                    salida.write(linea.replace(nombre_original, nuevo_nombre) + "\n")
                elif linea.startswith("#"):
                    salida.write(linea)
                else:
                    columnas = linea.strip().split("\t")
                    if columnas[0] in cromosomas:
                        columnas[0] = cromosomas[columnas[0]]
                    salida.write("\t".join(columnas) + "\n")
    except Exception as e:
        raise ValueError(f"Error al procesar el archivo GFF: {e}")
