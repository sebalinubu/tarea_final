import matplotlib.pyplot as plt
import numpy as np
import os


def filtrar_transposones_por_porcentaje(transposon_regions, percent):
    """
    Determinado por -p
    Filtra los transposones para obtener solo el porcentaje más grande basado en el tamaño.

    Args:
        transposon_regions (list): Lista de regiones de transposones [(inicio, fin)].
        percent (int): Porcentaje de transposones más grandes a conservar.

    Returns:
        list: Lista filtrada de regiones de transposones.
    """
    if percent == 100:
        return transposon_regions
    transposones_ordenados = sorted(transposon_regions, key=lambda x: (x[1] - x[0] + 1), reverse=True)
    corte = int(len(transposones_ordenados) * (percent / 100))
    return transposones_ordenados[:corte]


import matplotlib.pyplot as plt
import numpy as np
import os


def mapear_transposones_histograma(archivo_gff, imagen_salida, num_bins=120, porcentaje=100):
    """
    -gn
    Mapea gráficamente los transposones como un histograma a partir de un archivo .gff,
    con un número de bins que viene dado por el parametro -b y mostrando la longitud del 
    transposón más grande por bin.

    Args:
        archivo_gff (str): Ruta al archivo GFF que contiene los transposones.
        imagen_salida (str): Ruta para guardar el gráfico generado.
        num_bins (int): Número de bins en el histograma.
        porcentaje (int): Porcentaje de transposones más grandes para los graficos.
    """
    try:
        print(f"Procesando archivo GFF: {archivo_gff}")
        transposones = []
        longitudes_cromosomas = {}

        # Leer transposones del archivo GFF y calcular longitudes de cromosomas, hace lista con informacion importante
        with open(archivo_gff, "r") as archivo:
            for linea in archivo:
                if linea.startswith("#"):
                    continue
                columnas = linea.strip().split("\t")
                if len(columnas) >= 9 and columnas[2] in ["dispersed_repeat"]:
                    cromosoma = columnas[0]
                    inicio = int(columnas[3])
                    fin = int(columnas[4])
                    longitud = fin - inicio + 1
                    transposones.append((cromosoma, inicio, fin, longitud))
                    if cromosoma not in longitudes_cromosomas or fin > longitudes_cromosomas[cromosoma]:
                        longitudes_cromosomas[cromosoma] = fin

        if not transposones:
            print(f"No se encontraron transposones en el archivo {archivo_gff}.")
            return

        print(f"Total de transposones encontrados: {len(transposones)}")

        # asigna a cada cromosoma un desplasamiento inicial
        offsets_cromosomas = {}
        longitud_total_genoma = 0
        for cromosoma, longitud in sorted(longitudes_cromosomas.items()):
            offsets_cromosomas[cromosoma] = longitud_total_genoma
            longitud_total_genoma += longitud

        print(f"Offsets calculados por cromosoma: {offsets_cromosomas}")
        print(f"Tamaño total del genoma: {longitud_total_genoma} bp")

        # Ajustar posiciones al marco global
        posiciones_globales = [
            offsets_cromosomas[cromosoma] + ((inicio + fin) // 2)
            for cromosoma, inicio, fin, _ in transposones
        ]
        longitudes = [longitud for _, _, _, longitud in transposones]

        # Calcular el tamaño de los bins
        tamanio_bin = longitud_total_genoma // num_bins
        bins = np.arange(0, longitud_total_genoma + tamanio_bin, tamanio_bin)
        print(f"Tamaño del bin: {tamanio_bin} bp")

        # Filtrar los transposones por porcentaje
        transposones_ordenados = sorted(zip(posiciones_globales, longitudes), key=lambda x: x[1], reverse=True)
        corte = int(len(transposones_ordenados) * (porcentaje / 100))
        posiciones_filtradas, longitudes_filtradas = zip(*transposones_ordenados[:corte])
        print(f"Transposones considerados tras aplicar el filtro del {porcentaje}%: {len(posiciones_filtradas)}")

        # Calcular el histograma
        conteos, bordes = np.histogram(posiciones_filtradas, bins=bins)

        # Calcular la longitud máxima de los transposones en cada bin
        max_longitudes = [0] * len(conteos)
        for posicion, longitud in zip(posiciones_filtradas, longitudes_filtradas):
            indice_bin = np.searchsorted(bordes, posicion, side='right') - 1
            if indice_bin < len(max_longitudes):
                max_longitudes[indice_bin] = max(max_longitudes[indice_bin], longitud)

        # Generar el gráfico
        fig, ax = plt.subplots(figsize=(15, 5))
        ax.bar(bordes[:-1], conteos, width=tamanio_bin, color="#4CAF50", align="edge", edgecolor="black")  # Verde medio

        # Mostrar el tamaño del bin en el gráfico
        ax.text(0.95, 0.95, f"Tamaño del bin: {tamanio_bin:,} bp",
                transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right',
                bbox=dict(facecolor='white', alpha=0.5))

        # Añadir la longitud del transposón más grande sobre cada barra
        for i in range(len(conteos)):
            if conteos[i] > 0:
                ax.text(bordes[i] + tamanio_bin / 2, conteos[i], f"{max_longitudes[i]}",
                        ha='center', va='bottom', fontsize=8, color='black', rotation=90)

        ax.set_title(f"Distribución de Transposones")# (Tamaño del bin: {tamanio_bin:,} bp)")
        ax.set_xlabel("Posición en el genoma (bp)")
        ax.set_ylabel("Número de transposones")
        ax.grid(axis='y', linestyle="--", alpha=0.7)

        plt.tight_layout()
        plt.savefig(imagen_salida)
        print(f"Histograma de transposones guardado en {imagen_salida}.")
        plt.close()

    except Exception as e:
        print(f"Error al generar el histograma: {e}")


def mapear_porcentaje_ocupacion(archivo_gff, directorio_salida, numero_bins=120, porcentaje=100):
    """"
    -gp
    Genera un gráfico de porcentaje de ocupación en el genoma completo y un archivo GFF
    con los transposones ordenados por tamaño para cada bin.

    Args:
        archivo_gff (str): Ruta al archivo GFF que contiene los transposones.
        directorio_salida (str): Directorio donde se guardará el gráfico y el archivo GFF.
        numero_bins (int): Número total de bins.
        porcentaje (int): Porcentaje de transposones más grandes a considerar.
    """
    try:
        # Leer transposones del archivo GFF
        regiones_transposones = []
        with open(archivo_gff, "r") as archivo:
            for linea in archivo:
                if linea.startswith("#"):
                    continue
                columnas = linea.strip().split("\t")
                if len(columnas) >= 9 and columnas[2] in ["dispersed_repeat"]:
                    inicio = int(columnas[3])
                    fin = int(columnas[4])
                    longitud = fin - inicio + 1
                    regiones_transposones.append((inicio, fin, longitud))

        if not regiones_transposones:
            print(f"No se encontraron transposones en el archivo {archivo_gff}.")
            return

        # Filtrar por porcentaje
        regiones_transposones = filtrar_transposones_por_porcentaje(regiones_transposones, porcentaje)

        # Calcular bins
        posicion_maxima = max(fin for _, fin, _ in regiones_transposones)
        tamanio_bin = posicion_maxima / numero_bins
        bins = [[] for _ in range(numero_bins)]  # Lista para almacenar transposones por bin
        max_longitudes_transposones = [0] * numero_bins  # Longitud del transposón más grande en cada bin

        # Asignar transposones a bins
        for inicio, fin, longitud in regiones_transposones:
            bin_inicio = min(int(inicio / tamanio_bin), numero_bins - 1)
            bin_fin = min(int(fin / tamanio_bin), numero_bins - 1)
            for indice_bin in range(bin_inicio, bin_fin + 1):
                inicio_bin = indice_bin * tamanio_bin
                fin_bin = inicio_bin + tamanio_bin
                superposicion_inicio = max(inicio, inicio_bin)
                superposicion_fin = min(fin, fin_bin)
                if superposicion_inicio < superposicion_fin:
                    bins[indice_bin].append((inicio, fin, longitud))
                    max_longitudes_transposones[indice_bin] = max(max_longitudes_transposones[indice_bin], longitud)

        # Generar archivo GFF con transposones ordenados por bin
        os.makedirs(directorio_salida, exist_ok=True)
        ruta_gff_salida = os.path.join(directorio_salida, "genoma_transposones_por_bins.gff")
        with open(ruta_gff_salida, "w") as archivo_gff_salida:
            archivo_gff_salida.write("# Transposones ordenados por bins y tamaño para el genoma\n")
            for indice_bin, transposones_bin in enumerate(bins):
                transposones_ordenados = sorted(transposones_bin, key=lambda x: x[2], reverse=True)
                for inicio, fin, longitud in transposones_ordenados:
                    archivo_gff_salida.write(f"genoma\tdispersed_repeat\t.\t{int(inicio)}\t{int(fin)}\t.\t.\t.\tlongitud={longitud}\tbin={indice_bin}\n")
        print(f"Archivo GFF por bins guardado en {ruta_gff_salida}")

        # Calcular porcentaje de ocupación por bin
        porcentajes = [(sum(t[2] for t in bin) / tamanio_bin) * 100 for bin in bins]

        # Generar gráfico
        fig, ax = plt.subplots(figsize=(15, 5))
        ax.bar(range(numero_bins), porcentajes, width=1, color="skyblue", edgecolor="black", align="center")

        # Añadir información del transposón más grande en cada bin
        for i, max_longitud in enumerate(max_longitudes_transposones):
            if max_longitud > 0:
                ax.text(i, porcentajes[i] + 2, f"{max_longitud}", ha="center", va="bottom", fontsize=8, color="black")

        ax.set_title("Porcentaje de ocupación de transposones en el genoma completo")
        ax.set_xlabel("Bin")
        ax.set_ylabel("Porcentaje de ocupación")
        ax.set_ylim(0, 100)
        ax.grid(axis='y', linestyle="--", alpha=0.7)

        ax.text(0.95, 0.95, f"Tamaño del bin: {int(tamanio_bin):,} bp",
                transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right',
                bbox=dict(facecolor='white', alpha=0.5))

        imagen_salida = os.path.join(directorio_salida, "ocupacion_genoma.png")
        plt.tight_layout()
        plt.savefig(imagen_salida)
        print(f"Gráfico de porcentaje de ocupación guardado en {imagen_salida}.")
        plt.close()

    except Exception as e:
        print(f"Error al generar el gráfico de ocupación: {e}")



def mapear_porcentaje_ocupacion_por_cromosoma(archivo_gff, archivo_tamanos_cromosomas, directorio_salida, numero_bins=10, porcentaje=100):
    """
    cp
    Genera gráficos del porcentaje de ocupación de transposones por bins para cada cromosoma.
    Además, genera un archivo GFF con los transposones ordenados por tamaño para cada bin por cromosoma e incluye
    la información del transposón más grande en cada bin en el gráfico.

    Args:
        archivo_gff (str): Ruta al archivo GFF que contiene los transposones.
        archivo_tamanos_cromosomas (str): Archivo que contiene los tamaños de los cromosomas.
        directorio_salida (str): Directorio donde se guardarán los gráficos generados y los archivos de transposones por bin.
        numero_bins (int): Número de bins para cada cromosoma.
        porcentaje (int): Porcentaje de transposones más grandes a considerar.
    """
    try:
        # Leer tamaños de los cromosomas
        tamanos_cromosomas = {}
        with open(archivo_tamanos_cromosomas, "r") as archivo:
            for linea in archivo:
                cromosoma, tamano = linea.strip().split("\t")
                tamanos_cromosomas[cromosoma] = int(tamano)

        # Leer transposones del archivo GFF
        transposones_por_cromosoma = {}
        with open(archivo_gff, "r") as archivo:
            for linea in archivo:
                if linea.startswith("#"):
                    continue
                columnas = linea.strip().split("\t")
                if len(columnas) >= 9 and columnas[2] == "dispersed_repeat":
                    cromosoma = columnas[0]
                    inicio = int(columnas[3])
                    fin = int(columnas[4])
                    longitud = fin - inicio + 1
                    if cromosoma not in transposones_por_cromosoma:
                        transposones_por_cromosoma[cromosoma] = []
                    transposones_por_cromosoma[cromosoma].append((inicio, fin, longitud))

        os.makedirs(directorio_salida, exist_ok=True)

        for cromosoma, tamano in tamanos_cromosomas.items():
            if cromosoma not in transposones_por_cromosoma:
                print(f"No se encontraron transposones para {cromosoma}.")
                continue

            transposones = filtrar_transposones_por_porcentaje(transposones_por_cromosoma[cromosoma], porcentaje)

            tamanio_bin = tamano / numero_bins
            bins = [[] for _ in range(numero_bins)]  # Lista para almacenar transposones por bin
            max_longitudes_transposones = [0] * numero_bins  # Longitud del transposón más grande en cada bin

            # Asignar transposones a bins
            for inicio, fin, longitud in transposones:
                inicio_bin = max(0, min(int(inicio / tamanio_bin), numero_bins - 1))  # Limitar índices al rango válido
                fin_bin = max(0, min(int(fin / tamanio_bin), numero_bins - 1))        # Limitar índices al rango válido
                for indice_bin in range(inicio_bin, fin_bin + 1):
                    inicio_bin_actual = indice_bin * tamanio_bin
                    fin_bin_actual = inicio_bin_actual + tamanio_bin
                    superposicion_inicio = max(inicio, inicio_bin_actual)
                    superposicion_fin = min(fin, fin_bin_actual)
                    if superposicion_inicio < superposicion_fin:
                        bins[indice_bin].append((inicio, fin, longitud))
                        max_longitudes_transposones[indice_bin] = max(max_longitudes_transposones[indice_bin], longitud)

            # Guardar transposones por bin en un archivo para este cromosoma
            ruta_salida_bins = os.path.join(directorio_salida, f"{cromosoma}_transposones_por_bins.gff")
            with open(ruta_salida_bins, "w") as archivo_bins:
                archivo_bins.write(f"# Transposones ordenados por bins y tamaño para el cromosoma {cromosoma}\n")
                for indice_bin, transposones_bin in enumerate(bins):
                    transposones_ordenados = sorted(transposones_bin, key=lambda x: x[2], reverse=True)
                    for inicio, fin, longitud in transposones_ordenados:
                        archivo_bins.write(f"{cromosoma}\tdispersed_repeat\t.\t{int(inicio)}\t{int(fin)}\t.\t.\t.\tlongitud={longitud}\tbin={indice_bin}\n")
            print(f"Archivo GFF por bins para {cromosoma} guardado en {ruta_salida_bins}")

            # Generar gráfico de ocupación por bins
            porcentajes = [(sum(t[2] for t in bin) / tamanio_bin) * 100 for bin in bins]
            fig, ax = plt.subplots(figsize=(10, 5))
            ax.bar(range(numero_bins), porcentajes, width=1, color="blue", edgecolor="black", align="center")

            # Añadir el transposón más grande en cada bin al gráfico
            for i, max_longitud in enumerate(max_longitudes_transposones):
                if max_longitud > 0:
                    ax.text(i, porcentajes[i] + 2, f"{max_longitud:,}", ha="center", va="bottom", fontsize=8, color="black")

            ax.set_title(f"Porcentaje de ocupación de transposones - {cromosoma}")
            ax.set_xlabel("Bin")
            ax.set_ylabel("Porcentaje de ocupación")
            ax.set_ylim(0, max(porcentajes) + 10)  # Aumentar el límite del eje Y para evitar solapamiento
            ax.grid(axis='y', linestyle="--", alpha=0.7)

            ax.text(0.95, 0.95, f"Tamaño del bin: {int(tamanio_bin):,} bp",
                    transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right',
                    bbox=dict(facecolor='white', alpha=0.5))

            imagen_salida = os.path.join(directorio_salida, f"{cromosoma}_ocupacion.png")
            plt.tight_layout()
            plt.savefig(imagen_salida)
            print(f"Gráfico para {cromosoma} guardado en {imagen_salida}")
            plt.close()

    except Exception as e:
        print(f"Error al generar gráficos por cromosoma: {e}")
