#!/usr/bin/python3.9
# -*- coding: utf-8 -*-

# Módulos
import re,os,sys,argparse
import colorama as col
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

col.init()  # Inicializamos colorama

# Variables
original_seq = {} # Diccionario cuya clave es el nombre/id de la secuencia y el valor es una tupla con la secuencia y su longitud
trimmed_seq = {}  # Diccionario cuya clave es el nombre/id de la secuencia y el valor es una tupla con la secuencia y su longitud tras recortarla
registros = []  # Lista que va a contener los objetos SeqRecord utilizados para escribir el fichero FASTA
fichero_tsv = "len_seq.tsv"  # Nombre del fichero tsv que se va a crear
fichero_fasta = "trimmed_seq.fasta"  # Nombre del fichero fasta que va contener las secuencias recortadas
fasta_extensions = [".fasta", ".fas", ".fa", ".fna", ".ffn", ".faa", ".mpfa", ".frn" ]  # Lista que contiene las posibles extensiones que puede presentar un fichero fasta

# Funciones
def input_value(x):  # Definimos una función que impida la introducción de números negativos en "-start" y "-end"
    num = int(x)
    if num < 0:
        raise ValueError  # Si el número es negativo el programa devolverá un error
    else:
        return num  # En el caso de que sea positivo devuelve el número


def trimming(secuencia):  # Definimos la función que realizará el recorte de las secuencias
    if results.trim_N:  # Si se usa el argumento '-N' entonces 'trim_N' == True
        recorte = re.sub("^N+","", secuencia)  # Haciendo uso de la expresión regular "la secuencia comienza por una o más N" usamos la función del módulo re para sustituir nuestras coincidencias con nada
    else:
        recorte = secuencia  # Si no encontramos Ns al principio de la secuencia asignamos a la variable recorte la secuencia original
    trim_inicio = recorte[results.trim_start:]  # Hacemos slicing al comienzo de la secuencia (si no se introduce ningún argumento, el valor por defecto es 0)
    if results.trim_end > 0:  # Si el slicing fuera [:-0], que es lo mismo que [:0], cortaríamos la secuencia desde el inicio hasta el index 0 por lo que perderíamos la secuencia (por eso el uso del condicional)
        trim_final = trim_inicio[:-results.trim_end]  # Asignamos a la variable 'trim_final' el resultado de hacer slicing al final de la cadena (recortada o no en 5') usamos el índice negativo para que corte desde atrás la cadena
        res_recorte = trim_final  # Al final asignamos a la variable res_recorte la secuencia ya procesada
    else:
        res_recorte = trim_inicio  # En el caso de que sea 0, la secuencia no se toca en el extremo 3'
    return (Seq(res_recorte), len(res_recorte))  # La función devuelve una tupla con la nueva secuencia (formateada de str a Seq object) y su longitud
        

# Creamos un parseador
parser = argparse.ArgumentParser(
    prog="trimmfa",
    formatter_class=argparse.RawDescriptionHelpFormatter,  # Al usar RawDescription nos da más libertad al aceptar los strings como se escriben y no los formatea
    description="""
 _        _                      __      
| |      (_)                    / _|     
| |_ _ __ _ _ __ ___  _ __ ___ | |_ __ _ 
| __| '__| | '_ ` _ \| '_ ` _ \|  _/ _` |
| |_| |  | | | | | | | | | | | | || (_| |
 \__|_|  |_|_| |_| |_|_| |_| |_|_| \__,_|
___________________________________________                                                                                                                      
 >>Trimmer de ficheros fasta
 ##################################
 El programa usa Biopython para leer el fichero fasta proporcionado y
 elimina las Ns del extremo 5' si se incluye el parámetro.
 Posteriormente, tanto si se eliminan las Ns del principio,
 como si no, elimina los nucleótidos del inicio y fin de la secuencia
 según se pida por consola
""",
    epilog="Python3.9.18; biopython==1.81; colorama==0.4.6"  # Versiones de Python y los módulos utilizados
    )    
# Añadimos los argumentos                               
parser.add_argument("-f",  # Primer argumento, que acepta un fichero fasta
                    action="store",
                    dest="infile",  # Se almacena en la variable/atributo 'infile'
                    help="Fichero FASTA cuyas secuencias van a ser trimmeadas")
parser.add_argument("-start",  # Segundo argumento, acepta un número positivo que indica cuantos nucleótidos se quieren cortar en el extremo 5'
                    action="store",
                    dest="trim_start",  # Es almacenado en la variable/atributo 'trim_start'
                    type=input_value,  # Aquí uso la función que evita la introducción de números negativos
                    default=0,  # Si no se usa el argumento o se usa pero se deja sin valor (ej: -start) el valor por defecto es 0
                    help="Número de nucleótidos a eliminar de 5' (inicio de la secuencia) [Debe ser >= 0]")
parser.add_argument("-end",  # Tercer argumento, acepta un número positivo que indica cuantos nucleótidos se quieren cortar en el extremo 3'
                    action="store",
                    dest="trim_end",  # Se almacena en la variable/atributo 'trim_end'
                    type=input_value, # Lo mismo que en el argumento '-start', con la función evitamos la introducción de números negativos
                    default=0,  # Por defecto el valor es 0
                    help="Número de nucleótidos a eliminar de 3' (final de la secuencia) [Debe ser >= 0]")
parser.add_argument("-N",  # Cuarto argumento (sin valor asociado), si se usa se almacena el valor True en la variable/atributo 'trimm_N'
                    action="store_true",
                    dest="trim_N",
                    help="Elimina las Ns (sólo las que se encuentren en 5')")
parser.add_argument("--version", "-V", action="version", version="%(prog)s 1.0.0")  # Quinto argumento: la versión del programa

# Parseamos los argumentos y almacenamos el resultado en la variable 'result'
results = parser.parse_args()
# A continuación usamos la gestión de errores (try; except; else) para evitar que el programa se 'rompa' o funcione de manera incorrecta
try:
    records = SeqIO.parse(results.infile, "fasta")  # Creamos un iterador (records) con la función parse de SeqIO         
except FileNotFoundError as error:  # Si el fichero introducido no existe te devuelve por pantalla FileNotFoundError
    print(col.Fore.RED+f"FileNotFoundError: Fichero no encontrado {error}"+col.Fore.RESET)
except AttributeError:  # Si no se introduce nada (ej: programa.py o programa.py -N) te imprime la ayuda en pantalla
    parser.print_help()
else:  # Si no hay errores, pasamos a comprobar que la extensión del fichero es la adecuada
    if os.path.splitext(results.infile)[1] in fasta_extensions:  # Con el método split() separamos el nombre del fichero en una lista y extraemos el índice 1 (extensión del fichero) y comprobamos si está en la lista de extensiones
        for record in records:  # Si la extensión se encuentra en la lista de extensiones, iteramos por el iterador
            original_seq[record.id] = (str(record.seq), len(record.seq))  # Convertimos el objeto Seq a string para poder realizar el slicing y los añadimos, junto con la longitud de la secuencia, como una tupla que será el valor del diccionario y la clave será el nombre/id   
    # print(original_seq)
    # Iteramos por el diccionario creado
        for name, seq in original_seq.items():  # seq es una tupla (secuencia, longitud secuencia)
            resultado = trimming(seq[0])  # Llamamos a la funcion trimming que acepta como argumento la secuencia y asignamos el output a la variable resultado. Para acceder a la secuencia buscamos la primera posición de la tupla (index 0)
            trimmed_seq[name] = resultado  # Añadimos al diccionario el nombre/id como clave (name) y como valor la variable resultado, que es una tupla (secuencia trimmeada, longitud secuencia)
            registros.append(SeqRecord(resultado[0], id=name, description=""))  # Añadimos a la lista 'registros' el objeto SeqRecord, que se compone de la secuencia trimmeada, el id/nombre de la secuencia u una descripción que en este caso es un string vacío
            # print(trimmed_seq)
        print(col.Style.BRIGHT+f"Argumentos utilizados: {sys.argv}"+col.Style.RESET_ALL)  # Mostramos en pantalla los argumentos introducidos en el programa
    # Creación fichero tsv
        print("Creando el fichero tsv...")
        with open(fichero_tsv, "w") as f:  # Abrimos/creamos el fichero en modo escritura como f
                for name, len_seq in original_seq.items():  # Iteramos por el diccionario que contiene las secuencias originales
                    f.write(f"{name}\t{len_seq[1]}\t{trimmed_seq[name][1]}\n")  # Como la clave de ambos diccionarios es la misma, podemos extraer del diccionario trimmed_seq el tamaño de la secuencia recortada
        print(col.Fore.GREEN+f"El fichero {fichero_tsv} se ha creado."+col.Fore.RESET)
    
    # Creación fichero fasta con las nuevas secuencias
        print("Creando el fichero fasta...")
        SeqIO.write(registros, fichero_fasta, "fasta")  # Utilizando la función write de SeqIO
        print(col.Fore.GREEN+f"El fichero {fichero_fasta} se ha creado."+col.Fore.RESET)
        print(col.Style.BRIGHT+"Programa finalizado"+col.Style.RESET_ALL)
    else:
        raise ValueError(col.Fore.RED+"Extensión de fichero inválida"+col.Fore.RESET)  # Si la extensión del fichero no se encuentra en nuestra lista de extensiones, se devuelve ValueError en pantalla
