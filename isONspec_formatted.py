import sys
import os
import time
import threading
import argparse
import shutil
import subprocess
from colorama import Fore, Style

def calculate_distance(frequencies_rep, info_read):
    freq_A_diff = abs(frequencies_rep['freq_A'] - info_read['freq_A'])
    freq_C_diff = abs(frequencies_rep['freq_C'] - info_read['freq_C'])
    freq_G_diff = abs(frequencies_rep['freq_G'] - info_read['freq_G'])
    freq_T_diff = abs(frequencies_rep['freq_T'] - info_read['freq_T'])

    return freq_A_diff + freq_C_diff + freq_G_diff + freq_T_diff

def create_cluster_folders(clusters, original_names):
    # Create the results_cluster folder if it doesn't exist
    #results_cluster_folder = os.path.join(results_folder, "results_cluster")
    
    results_cluster_folder = "results_cluster"

    if not os.path.exists(results_cluster_folder):
        os.makedirs(results_cluster_folder)

    # Variable di conteggio per il numero dei cluster
    cluster_counter = 1

    # Iterate through clusters
    for cluster_representative, cluster_members in clusters.items():
        cluster_folder = os.path.join(results_cluster_folder, f"cluster_{cluster_counter}")

        # Create a folder for each cluster
        if not os.path.exists(cluster_folder):
            os.makedirs(cluster_folder)

        # Create a file named "all_clustered_reads.fastq" in each cluster folder
        cluster_fastq_path = os.path.join(cluster_folder, "all_clustered_reads.fastq")
        with open(cluster_fastq_path, 'w') as cluster_fastq:
            for read_id in cluster_members:
                # Use the original name if available, otherwise use the modified name
                original_name = original_names.get(read_id, read_id)

                # Write the original name of the read to the file
                cluster_fastq.write(f"{original_name}")

                # Get the file path of the read and copy its content to the cluster file
                read_relative_path = os.path.join("tmp_reads", f"{read_id}.fastq")
                #read_full_path = os.path.join(results_folder, read_relative_path)
                with open(read_relative_path, 'r') as read_file:
                    # Skip the first line (header line) of the read file
                    next(read_file, None)

                    # Copy the remaining content to the cluster file
                    shutil.copyfileobj(read_file, cluster_fastq)
        cluster_counter += 1

    print(f"Cluster folders created in {results_cluster_folder}")



def sum_overlap_sfs(file_path):
    total_sum = 0

    with open(file_path, 'r') as file:
        for line in file:
            # Dividi la riga utilizzando spazi multipli come delimitatori
            values = line.split()

            # Se ci sono abbastanza valori e il secondo valore è un numero, sommalo al totale
            if len(values) > 2:
                try:
                    second_number = int(values[2])
                    total_sum += second_number
                except ValueError:
                    # Gestisci il caso in cui il secondo valore non è un intero
                    print(f"Skipping line due to non-integer second number: {line}")

    return total_sum

def write_comparison_values_to_log(comparison_values, log_file_path):
    with open(log_file_path, 'w') as log_file:
        for read_id, comparison_value in comparison_values.items():
            log_file.write(f"{read_id}: {comparison_value}\n")

def modifica_prima_riga(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    if lines:

        id_read_originale = lines[0]

        # Ottieni il nome del file senza estensione
        nome_file_originale = os.path.splitext(os.path.basename(file_path))[0]

        # Modifica la prima riga con il formato desiderato
        lines[0] = f"@{nome_file_originale}\n"

        # Scrivi le modifiche nel file
        with open(file_path, 'w') as file:
            file.writelines(lines)

    return id_read_originale

def modifica_file_txt_cartella(cartella):
    # Dizionario per i nomi originali delle read
    original_names = {}

    # Verifica che il percorso della cartella esista
    if not os.path.exists(cartella):
        print(f"La cartella {cartella} non esiste.")
        return

    # Itera attraverso tutti i file nella cartella
    for filename in os.listdir(cartella):
        file_path = os.path.join(cartella, filename)

        # Verifica se il file è in formato fastq
        if filename.lower().endswith('.fastq') and os.path.isfile(file_path):
            # Ottieni il nome originale del file senza estensione prima della modifica
            nome_file_originale = os.path.splitext(os.path.basename(file_path))[0]
            print(f"questo è il nome originale {nome_file_originale}")
            # Modifica la prima riga del file
            id_originale = modifica_prima_riga(file_path)
            
            # Ottieni il nuovo nome dell'id dopo la modifica
            nuovo_nome_id = extract_read_info(file_path)

            # Salva il nome originale nel dizionario dei nomi originali
            original_names[nuovo_nome_id] = id_originale
    
    return original_names


def wait_for_file(file_path, timeout=10, polling_interval=1):
    start_time = time.time()

    while not os.path.exists(file_path):
        if time.time() - start_time > timeout:
            raise TimeoutError(f"Timeout while waiting for {file_path} to be created.")
        time.sleep(polling_interval)


def rename_file(input_file_path):
    # Extract the file name and extension
    file_name, file_extension = os.path.splitext(input_file_path)

    # Map the extensions to be renamed
    extension_mapping = {'.fq': '.fastq'}

    # Check if the extension is in the mapping
    if file_extension.lower() in extension_mapping:
        # Construct the new file path
        new_file_path = file_name + extension_mapping[file_extension.lower()]

        # Rename the file
        try:
            os.rename(input_file_path, new_file_path)
            return os.path.abspath(new_file_path), os.path.basename(input_file_path)
        except Exception as e:
            return f"Error while changing extension: {str(e)}", None
    elif file_extension.lower() == '.fastq':
        return input_file_path
    else:
        return f"Unsupported extension: {file_extension}", None

def write_reads_to_file(output_filename, reads):
    with open(output_filename, 'w') as output_file:
        for read in reads:
            output_file.write(read)

def split_reads(input_file):
    print(f'This is {input_file}')
    # Create the tmp_reads folder if it doesn't exist
    output_folder = 'tmp_reads'
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Read lines from the input file
    with open(input_file, 'r') as input_file_handle:
        lines = input_file_handle.readlines()

    # Split lines into groups of 4
    for i in range(0, len(lines), 4):
        output_filename = os.path.join(output_folder, f'read_{i//4 + 1}.fastq')
        write_reads_to_file(output_filename, lines[i:i+4])

    # Return the path of the tmp_reads folder
    return output_folder

def delete_folder(folder):
    try:
        shutil.rmtree(folder)
        print(f"Folder {folder} deleted successfully.")
    except Exception as e:
        print(f"Error while deleting folder {folder}: {e}")

def delete_file(file_path):
    try:
        os.remove(file_path)
        print(f"File {file_path} deleted successfully.")
    except Exception as e:
        print(f"Error while deleting file {file_path}: {e}")

def count_rows(file_path):
    try:
        with open(file_path, 'r') as file:
            row_count = sum(1 for line in file)
        return row_count
    except FileNotFoundError:
        return "File not found"
    except Exception as e:
        return f"An error occurred: {e}"

def remove_apostrophe(file_name):
    return file_name.replace("'", '')

def calculate_index(representative_path, read_id):
    output_folder = 'tmp_index'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    index_output_path = os.path.join(output_folder, f'{read_id}.fmd')

    # Command for SVDSS index
    index_command = ['SVDSS', 'index', '--fastq', representative_path, '--index', index_output_path]

    try:
        # Execute the SVDSS index command
        print(f'Executing: {" ".join(index_command)}')
        while not os.path.exists(index_output_path):
            result_index = subprocess.run(index_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=600)  # 10 minutes timeout

        if result_index.returncode != 0:
            print(f"Error executing SVDSS index: {result_index.stderr.decode()}")
            return None
        new_name = remove_apostrophe(index_output_path)
        os.rename(index_output_path, new_name)        
        print(f"SVDSS index for {read_id} executed successfully.")
        print(f"File {new_name}")
        return new_name
    except subprocess.TimeoutExpired:
        print("Timeout during the execution of SVDSS index.")
        return None
    except subprocess.CalledProcessError as e:
        print(f"Error executing SVDSS index: {str(e)}")
        return None 

def extract_read_info(file_path):
    try:
        with open(file_path, 'r') as file:
            # Read the first line
            first_line = file.readline().strip()
            
            # Remove the first character
            result = first_line[1:]
            
            return result
    except FileNotFoundError:
        return "File not found"
    except Exception as e:
        return f"Error: {e}"

def execute_search(index_output_path, read_path, read_id):
    output_folder = 'tmp_search'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    search_output_folder = os.path.join(output_folder, f'search_{read_id}')

    # Command for SVDSS search
    #search_command = ['SVDSS', 'search', '--index', index_output_path, '--fastq', read_path,'--workdir', search_output_folder]
    search_command = ['SVDSS', 'search', '--index', index_output_path, '--fastq', read_path, '--assemble','--workdir', search_output_folder]
    try:
        # Execute the SVDSS search command
        #print(f'Executing: {" ".join(search_command)}')
        while not os.path.exists(f'{search_output_folder}/solution_batch_0.assembled.sfs'):
            result_search = subprocess.run(search_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=600)  # 10 minutes timeout
        #while not os.path.exists(f'{search_output_folder}/solution_batch_0.sfs'):
        #    result_search = subprocess.run(search_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=600)  # 10 minutes timeout
        
        # if result_search.returncode != 0:
        #     print(f"Error executing SVDSS search: {result_search.stderr.decode()}")
        #     return None
        #print(f"SVDSS search for {read_id} executed successfully.")
        #print(f'Result path: {search_output_folder}/solution_batch_0.sfs')        
        #La Funzione wait_for_file non dovrebbe più servire, aggiunto controllo sul metodo
        #wait_for_file(f'{search_output_folder}/solution_batch_0.assembled.sfs')
        return f'{search_output_folder}/solution_batch_0.assembled.sfs'
        #return f'{search_output_folder}/solution_batch_0.sfs'
    except subprocess.TimeoutExpired:
        print("Timeout during the execution of SVDSS search.")
        return None
    except subprocess.CalledProcessError as e:
        print(f"Error executing SVDSS search: {str(e)}")
        return None

def populate_read_sequences(folder_path):
    # Dictionary for ID, length, file name, and nucleotide frequency information of reads
    read_info = {}

    try:
        # List of files in the folder
        file_list = os.listdir(folder_path)

        for file_name in file_list:
            # Check if the file is of type FASTQ or FQ
            if file_name.lower().endswith(('.fastq', '.fq')):
                file_path = os.path.join(folder_path, file_name)

                with open(file_path, 'r') as file:
                    # Read lines from the file
                    lines = file.readlines()

                    # Iterate over the lines in groups of 4
                    for i in range(0, len(lines), 4):
                        # Extract read ID, sequence, and quality information
                        read_id = lines[i].strip().replace('@', '')
                        sequence = lines[i + 1].strip()
                        length_read = len(sequence)

                        # Calculate nucleotide frequencies
                        freq_A = round(sequence.count('A') / length_read * 100, 2)
                        freq_C = round(sequence.count('C') / length_read * 100, 2)
                        freq_G = round(sequence.count('G') / length_read * 100, 2)
                        freq_T = round(sequence.count('T') / length_read * 100, 2)

                        # Save the information in the read_info dictionary
                        read_info[read_id] = {
                            'length_read': length_read,
                            'file_name': file_path,
                            'freq_A': freq_A,
                            'freq_C': freq_C,
                            'freq_G': freq_G,
                            'freq_T': freq_T
                        }

        # Sort the read_info dictionary by the length of the reads
        read_info = dict(sorted(read_info.items(), key=lambda x: x[1]['length_read'], reverse=True))

        return read_info

    except Exception as e:
        return f"Error while reading files: {e}"

def process_reads(folder_path):
    # Dictionary for representatives
    representatives = {}

    # Dictionary for clusters
    clusters = {}

    # Dictionary for read sequences
    read_sequences = {}

    # Dictionary for quality lines of read sequences
    quality_lines = {}

    comparison_values = {}
    soglia_meno_cinquanta ={}


    read_sequences = populate_read_sequences(folder_path)
    #modifica_file_txt_cartella('tmp_reads')
    print('Start')

    first_file = False
    max_value = 0
    threshold_value = 0.01
    min_percentuale = 100

    for read_id, info in read_sequences.items():
        length_read = info['length_read']
        file_name = info['file_name']

        if not first_file:
            index_path = calculate_index(file_name, read_id)
            # Read the frequencies from the index output
            frequencies = {'freq_A': info['freq_A'], 'freq_C': info['freq_C'], 'freq_G': info['freq_G'], 'freq_T': info['freq_T']}
            representatives[read_id] = {'length_read': length_read, 'index_path': index_path, 'frequencies': frequencies}
            clusters[read_id] = [read_id]
            print(f'La frequenza è A:{info['freq_A']} C:{info['freq_C']} G:{info['freq_G']} T:{info['freq_T']}')
            first_file = True
            continue
        else:
            best_match_representative = None
            best_match_distance = float('inf')
            best_num_sfs = 0
            best_match_percentuale = 0
            #number_max_substrings_obt = (length_read * (length_read + 1)) / 2;
            print(f'La frequenza è A:{info['freq_A']} C:{info['freq_C']} G:{info['freq_G']} T:{info['freq_T']}')


            for rep_id, rep_info in sorted(representatives.items(), key=lambda x: calculate_distance(x[1]['frequencies'], info), reverse=True):
                length_representative = rep_info['length_read']
                index_representative = rep_info['index_path']
                print(f'Rep index: {index_representative}')
                #print(f'Rep length: {length_representative}')
                
                file_search_out = execute_search(index_representative, file_name, read_id)
                #print(f'Search output file: {file_search_out}')                
                
                #number_sfs = count_rows(file_search_out)
                #print(f'Number of sfs: {number_sfs}')                
                
                #length_check = (length_read / (length_read + length_representative))
                #print(f'Length check: {length_check}')
                
                #sfs_check = (number_sfs / number_max_substrings_obt)
                #print(f'Sfs check: {sfs_check}')
               
                #result_check = 2 * length_check * sfs_check
                #print(f'Result check, threshold to set: {result_check}')
                #comparison_value = round(result_check, 6)

                # Altra Logica su OVERLAP SFS
                length_overlap_sfs = sum_overlap_sfs(file_search_out)

                diff_on_total = length_read - length_overlap_sfs
                print(f'la lunghezza della stringa è {length_read}')
                print(f'la differenza tra lunghezza e overlap è {diff_on_total}')
                percentuale_diff_length = round(diff_on_total/length_read*100,2)

                delete_file(file_search_out)
                
                #if comparison_value < best_match_distance:
                #    best_match_distance = comparison_value
                #    best_match_representative = rep_id
                #    best_num_sfs = number_sfs

                if (percentuale_diff_length >= best_match_percentuale):
                    best_match_percentuale = percentuale_diff_length;
                    best_match_representative = rep_id;
                    print(f'{Fore.MAGENTA}Questo è il branch di selezione best {best_match_percentuale} {Style.RESET_ALL}')
                #if comparison_value > max_value:
                    #print(f'Valere attuale {comparison_value}, valore massimo old {max_value}')
                #    max_value = comparison_value
                #    comparison_values[read_id] = comparison_value
                if(percentuale_diff_length<min_percentuale):
                    min_percentuale = percentuale_diff_length

                if(percentuale_diff_length>=10):
                    print(f"{Fore.GREEN}SIMILE! Con percentuale di somiglianza {percentuale_diff_length}{Style.RESET_ALL}")
                    break
                else:
                   print(f"{Fore.RED}NON SIMILE! Con percentuale di somiglianza {percentuale_diff_length}{Style.RESET_ALL}")
                   soglia_meno_cinquanta[read_id] = percentuale_diff_length
                

            #if best_match_distance <= threshold_value:
            if (best_match_percentuale>=10):
                print(f'{Fore.YELLOW}opto per il cluste già esistente, perché {best_match_percentuale}{Style.RESET_ALL}')
                if best_match_representative in clusters:
                    clusters[best_match_representative].append(read_id)
                else:
                    clusters[best_match_representative] = [read_id]
            else:
                print(f'{Fore.BLUE}opto per crea un nuovo rappresentante. {best_match_percentuale}{Style.RESET_ALL}')
                index_path = calculate_index(file_name, read_id)
                # The length difference is above the threshold, insert the read as a new representative
                representatives[read_id] = {'length_read': length_read, 'index_path': index_path, 'frequencies': frequencies}

                # Create a new cluster with the read as the representative
                clusters[read_id] = [read_id]

    return min_percentuale, soglia_meno_cinquanta, max_value, representatives, clusters, read_sequences, quality_lines, comparison_values

def main():

    # Use argparse to handle command line arguments
    parser = argparse.ArgumentParser(description='Read clustering process.')
    parser.add_argument('file_path', metavar='file_path', type=str, help='The path of the FASTQ file to process')
    parser.add_argument('--threads', type=int, default=1, help='The number of threads to use (default: 1)')

    args = parser.parse_args()

    # Get the file path and number of threads from the command line argument
    file_path = args.file_path
    num_threads = args.threads

    # Set the desired number of threads
    threading._DummyThread._maxsize = num_threads
    
    start_time_total = time.time()

    start_time_modifica_file = time.time()
    renamed_file_path = rename_file(file_path)
    output_folder_path = split_reads(renamed_file_path)
    original_names = modifica_file_txt_cartella('tmp_reads')
    end_time_modifica_file = time.time()

    start_time_process_reads = time.time()
    min_percentuale, lista_meno_cinq, max_value, representatives, clusters, read_sequences, quality_lines, comparison_values = process_reads(output_folder_path)
    end_time_process_reads = time.time()

    log_file_path = "comparison_values_log.txt"  # Imposta il percorso del file di log
    write_comparison_values_to_log(comparison_values, log_file_path)

    start_time_create_clusters = time.time()
    #results_folder = os.path.dirname(os.path.abspath(file_path))
    #print(f"Questa è la folder di risultato {results_folder}")
    create_cluster_folders(clusters, original_names)
    end_time_create_clusters = time.time()

    for chiave, valore in clusters.items():
        lunghezza_elemento = len(valore)
        print(f"Lunghezza dell'elemento associato a {chiave}: {lunghezza_elemento}")

    print(f'There are {len(representatives)} representatives')
    print(f'Questo è il valore max di distanza {max_value}')
    delete_folder("tmp_search")
    delete_folder("tmp_index")
    delete_folder("tmp_reads")
    print(f'Ci sono {len(lista_meno_cinq)} reads differenti')
    
    end_time_total = time.time()

    print(f"Tempo per modifica_file_txt_cartella: {end_time_modifica_file - start_time_modifica_file} secondi")
    print(f"Tempo per process_reads: {end_time_process_reads - start_time_process_reads} secondi")
    print(f"Tempo per create_cluster_folders: {end_time_create_clusters - start_time_create_clusters} secondi")
    print(f"Tempo totale di esecuzione: {end_time_total - start_time_total} secondi")

    print(f"Percentuale minima intracluster {min_percentuale}")
if __name__ == "__main__":
    main()
