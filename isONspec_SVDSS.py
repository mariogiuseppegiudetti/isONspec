import sys
import os
import time
import threading
import argparse
import shutil
from subprocess import STDOUT, check_output
import shutil

num_threads = '16'

def rename_file(input_file_path):
    # Estrai il nome del file e l'estensione
    file_name, file_extension = os.path.splitext(input_file_path)

    # Mappa le estensioni da rinominare
    extension_mapping = {'.fq': '.fastq'}

    # Verifica se l'estensione è nella mappatura
    if file_extension.lower() in extension_mapping:
        # Costruisci il nuovo percorso del file rinominato
        new_file_path = file_name + extension_mapping[file_extension.lower()]

        # Rinomina il file
        try:
            os.rename(input_file_path, new_file_path)
            return os.path.abspath(new_file_path), os.path.basename(input_file_path)
        except Exception as e:
            return f"Errore durante il cambio di estensione: {str(e)}", None
    elif file_extension.lower() == '.fastq':
        return input_file_path
    else:
        return f"Estensione non supportata: {file_extension}", None


def write_reads_to_file(output_filename, reads):
    with open(output_filename, 'w') as output_file:
        for read in reads:
            output_file.write(read)

def split_reads(input_file):
    print(f'questo è {input_file}')
    # Crea la cartella tmp_reads se non esiste
    output_folder = 'tmp_reads'
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Leggi le linee dal file di input
    with open(input_file, 'r') as input_file_handle:
        lines = input_file_handle.readlines()

    # Divide le linee in gruppi di 4
    for i in range(0, len(lines), 4):
        output_filename = os.path.join(output_folder, f'read_{i//4 + 1}.fastq')
        write_reads_to_file(output_filename, lines[i:i+4])

    # Restituisci il percorso della cartella tmp_reads
    return output_folder


def delete_folder(cartella):
    try:
        shutil.rmtree(cartella)
        print(f"Cartella {cartella} cancellata con successo.")
    except Exception as e:
        print(f"Errore durante la cancellazione della cartella {cartella}: {e}")

def count_row(file_path):
    try:
        with open(file_path, 'r') as file:
            numero_righe = sum(1 for line in file)
        return numero_righe
    except FileNotFoundError:
        return "File non trovato"
    except Exception as e:
        return f"Si è verificato un errore: {e}"

def rimuovi_apice(nome_file):
    return nome_file.replace("'", '')

def calculate_index(representative_path, id_read):
    output_folder = 'tmp_index'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    index_output_path = os.path.join(output_folder, f'{id_read}.fmd')

    # Comando per SVDSS index
    index_command = ['SVDSS', 'index', '--threads', num_threads, '--fastq', representative_path, '--index', index_output_path]
    
    try:
        # Esegui il comando SVDSS index
        print(f'Esecuzione di: {" ".join(index_command)}')
        result_index = check_output(index_command, stderr=STDOUT, timeout=600)  # Timeout di 10 minuti
        #if result_index.returncode != 0:
        #    print(f"Errore nell'esecuzione dell'index SVDSS: {result_index.stderr.decode()}")
        #    return None
        #nuovo_nome = rimuovi_apice(index_output_path)
        #os.rename(index_output_path, nuovo_nome)        
        print(f"Index SVDSS per {id_read}, eseguita con successo.")
        print(f"File {index_output_path}")
        return index_output_path
    #except subprocess.TimeoutExpired:
    #    print("Timeout durante l'esecuzione dell'index SVDSS.")
    #    return None
    #except subprocess.CalledProcessError as e:
    #    print(f"Errore nell'esecuzione dell'index SVDSS: {str(e)}")
    #    return None 
    except Exception as e:
        return f"Errore: {e}"

def extract_info_read(file_path):
    try:
        with open(file_path, 'r') as file:
            # Leggi la prima riga
            prima_riga = file.readline().strip()

            # Elimina il primo carattere
            risultato = prima_riga[1:]

            return risultato
    except FileNotFoundError:
        return "File non trovato"
    except Exception as e:
        return f"Errore: {e}"

def create_search_folder(id_read):
    # Create the main search folder if it doesn't exist
    main_search_folder = 'tmp_search'
    if not os.path.exists(main_search_folder):
        os.makedirs(main_search_folder)

    # Create the subfolder search_{id_read} within the main folder
    subfolder_name = f'search_{id_read}'
    subfolder_path = os.path.join(main_search_folder, subfolder_name)
    if not os.path.exists(subfolder_path):
        os.makedirs(subfolder_path)

    return subfolder_path

def execute_search(index_output_path,read_path,id_read):
    search_output_folder=create_search_folder(id_read)

    # Comando per SVDSS search
    search_command = ['SVDSS', 'search', '--threads', num_threads, '--index', index_output_path, '--fastq', read_path, '--workdir', search_output_folder]

    try:
        # Esegui il comando SVDSS search
        print(f'Esecuzione di: {" ".join(search_command)}')
        result_search = check_output(search_command, stderr=STDOUT,timeout=600)  # Timeout di 10 minuti
        #if result_search.returncode != 0:
        #    print(f"Errore nell'esecuzione della ricerca SVDSS: {result_search}")
        #    return None
        print(f"Search SVDSS per {id_read}, eseguita con successo.")
        print(f'Percorso del risultato: {search_output_folder}/solution_batch_0.sfs')        
        return f'{search_output_folder}/solution_batch_0.sfs'
    #except subprocess.TimeoutExpired:
    #    print("Timeout durante l'esecuzione della ricerca SVDSS.")
    #    return None
    #except subprocess.CalledProcessError as e:
    #    print(f"Errore nell'esecuzione della ricerca SVDSS: {str(e)}")
    #    return None
    except Exception as e:
         return f"Errore: {e}"

def popolate_read_sequences(folder_path):
    # Dizionario per le informazioni ID, lunghezza e nome del file delle read
    read_info = {}

    try:
        # Lista dei file nella cartella
        lista_file = os.listdir(folder_path)

        for file in lista_file:
            # Verifica se il file è di tipo FASTQ o FQ
            if file.lower().endswith(('.fastq', '.fq')):
                file_path = os.path.join(folder_path, file)

                with open(file_path, 'r') as file:
                    # Inizializza una lista vuota per le informazioni di ogni file
                    info_file = []

                    # Leggi le linee del file
                    lines = file.readlines()

                    # Itera sulle linee del file
                    i = 0
                    while i < len(lines):
                        # Verifica se la linea è un ID di read
                        if lines[i].startswith('@'):
                            id_read = lines[i].strip().replace('@', '')
                            i += 1

                            # Leggi la sequenza
                            sequence = lines[i].strip()
                            length_read = len(sequence)

                            # Salva le informazioni nel dizionario read_info
                            read_info[id_read] = {'length_read': length_read, 'file_name': file_path}

                        i += 1

        return read_info

    except Exception as e:
        return f"Errore durante la lettura dei file: {e}"

def process_reads_new(folder_path):

    # Dizionario per i rappresentanti
    representatives = {}

    # Dizionario per i cluster
    clusters = {}

    # Dizionario per le sequenze delle read
    read_sequences = {}

    # Dizionario delle qualità per le sequenze delle read
    quality_lines = {}

    confronto_values = {}

    read_sequences = popolate_read_sequences(folder_path)
    print('Start')

    firstFile=False
    valore_massimo=0;
    valore_soglia =0.005

    for id_read, info in read_sequences.items():
        length_read = info['length_read']
        file_name = info['file_name']        

        if not firstFile:
            index_path=calculate_index(file_name,id_read)
            representatives[id_read] = {'length_read': length_read, 'index_path': index_path}
            firstFile=True
            continue
        else:
            best_match_representative = None
            best_match_distance = float('inf')
            best_num_sfs = 0
            number_max_substrings_obt = (length_read*(length_read+1))/2;

            for id_representative, info_rep in representatives.items():
                length_representative=info_rep['length_read']
                index_representative=info_rep['index_path']
                print(f'index rapp {index_representative}')
                print(f'length_representative {length_representative}')
                file_search_out=execute_search(index_representative,file_name,id_read)
                print(f'file di output dal search {file_search_out}')                
                number_sfs = count_row(file_search_out)
                print(f'number di sfs {number_sfs}')                
                length_check = (length_read / (length_read + length_representative))
                print(f'numero di check lunghezza {length_check}')
                sfs_check = (number_sfs/number_max_substrings_obt)
                print(f'numero di check lunghezza {sfs_check}')
                result_check = 2*length_check*sfs_check
                print(f'result_check soglia da impostare {result_check}')
                valore_confronto = round(result_check, 6)
                if  valore_confronto < best_match_distance:
                    best_match_distance = number_sfs
                    best_match_representative = id_representative
                    best_num_sfs = number_sfs
                    
                    if valore_confronto > valore_massimo:
                        valore_massimo = valore_confronto
                        confronto_values[id_read] = valore_confronto

            if valore_confronto <= valore_soglia:
                if best_match_representative in clusters:
                    clusters[best_match_representative].append(id_read)
                else:
                    clusters[best_match_representative] = [id_read]
            else:
                index_path=calculate_index(file_name,id_read)
                # La differenza di lunghezza è superiore alla soglia, inserisci la read come nuovo rappresentante
                representatives[id_read] = {'length_read': length_read, 'index_path': index_path}

                # Crea un nuovo cluster con la read come rappresentante
                clusters[id_read] = [id_read]
            delete_folder("tmp_search") 
    return representatives, clusters, read_sequences, quality_lines, confronto_values

def main():

    # Utilizza argparse per gestire gli argomenti da riga di comando
    parser = argparse.ArgumentParser(description='Processo di clustering delle reads.')
    parser.add_argument('file_path', metavar='file_path', type=str, help='Il percorso del file FASTQ da processare')
    parser.add_argument('--threads', type=int, default=1, help='Il numero di threads da utilizzare (predefinito: 1)')

    args = parser.parse_args()

    # Ottieni il percorso del file e il numero di thread dall'argomento della riga di comando
    file_path = args.file_path
    numero_thread = args.threads

    # Imposta il numero desiderato di thread
    threading._DummyThread._maxsize = numero_thread
    
    fasta_file_path=rename_file(file_path)
    output_folder_path = split_reads(fasta_file_path)
    representatives, clusters, read_sequences, quality_lines, confronto_values = process_reads_new(output_folder_path)

    print(f'ci sono {len(representatives)} rappresentanti')

    delete_folder("tmp_search")
    delete_folder("tmp_index")

if __name__ == "__main__":
    main()
