import sys
import os
import time
import threading
import argparse

# Dizionario per associare il formato del file al carattere di inizio di ogni read
format_to_start_char = {'fastq': '@', 'fq': '@', 'fasta': '>', 'fa': '>'}

def count_reads(file_path):
    # Verifica l'estensione del file per determinare il formato
    file_extension = get_file_extension(file_path)

    # Verifica che l'estensione del file sia supportata
    if file_extension not in format_to_start_char:
        print("Formato del file non supportato.")
        return -1

    start_char = format_to_start_char[file_extension]

    # Contatore per il numero di reads
    reads_count = 0

    # Apri il file e conta le reads
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith(start_char):
                reads_count += 1

    return reads_count

def get_file_extension(file_path):
    return file_path.split('.')[-1].lower()

def read_sequences_from_file(file_path, file_format):
    """
    Legge un file FASTA o FASTQ e restituisce un dizionario con gli identificatori di read come chiavi e le sequenze come valori.
    """
    read_sequences = {}
    quality_lines = {}  # Contenitore per le linee di qualità (solo per FASTQ)

    with open(file_path, 'r') as file:
        current_read_id = None
        current_read_sequence = ""
        current_quality_sequence = ""
        next_is_quality_line = False  # Indica se la linea corrente è una linea di qualità (solo per FASTQ)

        for line in file:
            line = line.strip()

            if file_format.lower() in ['fasta', 'fa']:
                if line.startswith(">"):
                    # Inizio di una nuova sequenza in formato FASTA
                    current_read_id = line[1:]  # Rimuovi il carattere di inizio
                    current_read_sequence = ""
                else:
                    if current_read_id is not None:
                        # Sequenza o linea vuota
                        current_read_sequence += line
                        # Salva la sequenza nel dizionario
                        read_sequences[current_read_id] = current_read_sequence
            elif file_format.lower() in ['fastq', 'fq']:
                if line.startswith("@"):
                    # Inizio di una nuova sequenza in formato FASTQ
                    current_read_id = line[1:]  # Rimuovi il carattere di inizio
                    current_read_sequence = ""
                    current_quality_sequence = ""  # Azzeramento della sequenza di qualità
                    # Controlli per la Quality Line
                    line_count = 0; #Se anche la quality line inizia con +, non da problemi con la terza riga                
                    next_is_quality_line = False
                elif next_is_quality_line and line_count == 1:
                    if current_read_id is not None:
                        # Linea di qualità
                        current_quality_sequence += line
                        # Salva la sequenza di qualità nel contenitore associato all'ID
                        quality_lines[current_read_id] = current_quality_sequence
                        # Debug: stampa il read_id e la sequenza di qualità
                        #print(f"DEBUG: {current_read_id} - Quality: {current_quality_sequence}")
                    line_count = 0
                    next_is_quality_line = False  # Imposta a False per evitare il salvataggio delle linee di qualità
                elif line.startswith("+") and line_count < 1:
                    # linea delimitatore alla qualità
                    line_count += 1
                    next_is_quality_line = True
                else:
                    if current_read_id is not None:
                        # Sequenza o linea vuota
                        current_read_sequence += line
                        # Salva la sequenza nel dizionario
                        read_sequences[current_read_id] = current_read_sequence

                
        # Salva l'ultima sequenza
        if current_read_id is not None:
            read_sequences[current_read_id] = current_read_sequence
            quality_lines[current_read_id] = current_quality_sequence

    return read_sequences, quality_lines

def process_reads(file_path):
    # Verifica l'estensione del file per determinare il formato
    file_extension = get_file_extension(file_path)

    # Verifica che l'estensione del file sia supportata
    if file_extension not in format_to_start_char:
        print("Formato del file non supportato.")
        return

    start_char = format_to_start_char[file_extension]

    # Dizionario per i rappresentanti
    representatives = {}

    # Dizionario per i cluster
    clusters = {}

    # Dizionario per le sequenze delle read
    read_sequences = {}

    # Dizionario delle qualità per le sequenze delle read
    quality_lines = {}

    confronto_values = {}

    # Leggi le sequenze da file
    read_sequences,quality_lines = read_sequences_from_file(file_path, file_extension)

    # Inizializza i rappresentanti con la prima read
    initial_read_id = list(read_sequences.keys())[0]
    representatives[initial_read_id] = read_sequences[initial_read_id]
    clusters[initial_read_id] = [initial_read_id]

    # Ciclo sulle read e confronto con i rappresentanti
    # Definisci la soglia per la differenza di lunghezza
    valore_soglia = 0.004
    contatore = 0
    valore_massimo = 0
    # Ciclo sulle read e confronto con i rappresentanti
    for read_id, read_sequence in read_sequences.items():

        # Se la read è già un rappresentante, passa alla prossima iterazione
        contatore += 1
        print(f"Siamo al {contatore}/{len(read_sequences)} e attualmente i rappresentanti sono {len(representatives)}")

        if read_id in representatives:
            continue
        best_match_representative = None
        best_match_distance = float('inf')
        best_num_sfs = 0

        # Confronto con tutti i rappresentanti
        for representative_id, representative_sequence in representatives.items():
            sfs_read_corr, num_sfs_read_corr = stringhe_specifiche(read_sequence, representative_sequence)
            num_max_sottostringhe_ottenibili = (len(read_sequence)*(len(read_sequence)+1))/2;
            result_check = 2 * (len(read_sequence) / (len(read_sequence) + len(representative_sequence)))*(num_sfs_read_corr/num_max_sottostringhe_ottenibili)
            valore_confronto = round(result_check, 6)
            print(f"Lunghezza read {len(read_sequence)} ")
            print(f"Lunghezza rapp {len(representative_sequence)} ")
            print(f"Num sottostringhe ottenute {num_sfs_read_corr} ")
            print(f"num max sottostringhe {num_max_sottostringhe_ottenibili} ")
            print(f"calcolo {valore_confronto} ")
            # Aggiorna il miglior rappresentante se la distanza è minore
            if  valore_confronto < best_match_distance:
                best_match_distance = num_sfs_read_corr
                best_match_representative = representative_id
                best_num_sfs = num_sfs_read_corr
                #confronto_values[read_id] = valore_confronto
                if valore_confronto > valore_massimo:
                    valore_massimo = valore_confronto
                    confronto_values[read_id] = valore_confronto
            #print(f"Questo è il valore del sfs {valore_confronto} ")
            
        # Verifica la lunghezza e se soddisfa la regola
        if valore_confronto <= valore_soglia:
            # La differenza di lunghezza è inferiore o uguale alla soglia, assegna la read al rappresentante
            # In questo caso, assegna la read al cluster del rappresentante
            if best_match_representative in clusters:
                clusters[best_match_representative].append(read_id)
            else:
                clusters[best_match_representative] = [read_id]
        else:
            # La differenza di lunghezza è superiore alla soglia, inserisci la read come nuovo rappresentante
            representatives[read_id] = read_sequence
    
            # Crea un nuovo cluster con la read come rappresentante
            clusters[read_id] = [read_id]

    return representatives, clusters, read_sequences, quality_lines, file_extension, confronto_values

def save_clusters_to_file(clusters, read_sequences, quality_lines, file_extension):
    for representative_id, cluster_reads in clusters.items():
        folder_path = f"cluster_{representative_id}"
        os.makedirs(folder_path, exist_ok=True)

        # Determina il formato del file iniziale e il carattere di inizio delle read
        original_format = file_extension.lower()
        start_char = format_to_start_char[original_format]

        # Salva le read associate al cluster nel file
        file_path = os.path.join(folder_path, f"cluster_{representative_id}_reads.{original_format}")
        with open(file_path, 'w') as cluster_file:
            for read_id in cluster_reads:
                cluster_file.write(f"{start_char}{read_id}\n{read_sequences[read_id]}\n")

                # Se il formato è FASTQ o FQ, aggiungi le qualità dalla read originale
                if original_format in ['fastq', 'fq']:
                    # Verifica se la chiave esiste in quality_lines prima di accedere
                    quality_line = quality_lines.get(read_id, "")
                    cluster_file.write(f"+\n{quality_line}\n")

def get_quality_line(file_path, read_id):
    """
    Restituisce la linea di qualità associata a un dato identificatore di read da un file FASTQ.
    """
    with open(file_path, 'r') as file:
        in_sequence = False
        for line in file:
            if in_sequence:
                # Se la riga corrente inizia con "+" indica l'inizio delle linee di qualità
                if line.startswith("+"):
                    # Leggi la linea successiva come linea di qualità
                    quality_line = next(file).strip()
                    return quality_line
            if line.startswith(read_id):
                in_sequence = True

    return ""

def main():

    tempo_iniziale = time.time()

    # Utilizza argparse per gestire gli argomenti da riga di comando
    parser = argparse.ArgumentParser(description='Processo di clustering delle reads.')
    parser.add_argument('file_path', metavar='file_path', type=str, help='Il percorso del file da processare')
    parser.add_argument('--threads', type=int, default=1, help='Il numero di threads da utilizzare (predefinito: 1)')

    args = parser.parse_args()

    # Ottieni il percorso del file e il numero di thread dall'argomento della riga di comando
    file_path = args.file_path
    numero_thread = args.threads

    # Imposta il numero desiderato di thread
    threading._DummyThread._maxsize = numero_thread

    reads_count = count_reads(file_path)

    if reads_count == -1:
        return

    reads_count = count_reads(file_path)

    if reads_count == -1:
        return

    representatives, clusters, read_sequences, quality_lines, file_extension, confronto_values = process_reads(file_path)
    
    save_clusters_to_file(clusters, read_sequences, quality_lines, file_extension)

    log_file_path = 'log.txt'
    with open(log_file_path, 'w') as log_file:
        for read_id, confronto_value in confronto_values.items():
            log_file.write(f"{read_id}: {confronto_value}\n")

    print(f"Il numero di reads nel file '{file_path}' è: {len(read_sequences)}")
    print(f"Il numero di cluster è: {len(clusters)}")

    # Stampa il numero di reads nei file creati in output
    print("\nNumero di reads nei file creati in output:")
    for representative_id, cluster_reads in clusters.items():
        file_path = f"cluster_{representative_id}/cluster_{representative_id}_reads.{file_extension}"
        num_reads = len(cluster_reads)
        print(f"{file_path}: {num_reads} reads")

    tempo_finale = time.time()
    tempo_trascorso = tempo_finale - tempo_iniziale
    print(f"Tempo impiegato: {tempo_trascorso} secondi")

def stringhe_specifiche_old(stringa1, stringa2):
    # Inizializziamo una lista vuota per memorizzare le stringhe specifiche
    stringhe_specifiche_list = []

    # Iteriamo attraverso tutte le sottostringhe della prima stringa
    for i in range(len(stringa1)):
        for j in range(i + 1, len(stringa1) + 1):
            sottostringa = stringa1[i:j]

            # Verifichiamo se la sottostringa non è presente nella seconda stringa
            if sottostringa not in stringa2:
                # Aggiungiamo la sottostringa alla lista delle stringhe specifiche
                stringhe_specifiche_list.append(sottostringa)

    return stringhe_specifiche_list, len(stringhe_specifiche_list)

def stringhe_specifiche(stringa1, stringa2):
    # Inizializziamo una lista vuota per memorizzare le stringhe specifiche
    old_St = []
    print("Start: Calcolo Stringhe Spec")
    # Iteriamo attraverso tutte le sottostringhe della prima stringa
    for i in range(len(stringa1)):
        for j in range(i + 1, len(stringa1) + 1):
            sottostringa = stringa1[i:j]

            # Verifichiamo se la sottostringa non è presente nella seconda stringa
            if sottostringa not in stringa2:
                # Verifichiamo se nessuna sottostringa più corta è una sottostringa di questa
                if all(sottostringa not in old_str for old_str in old_St):
                    # Aggiungiamo la sottostringa alla lista delle stringhe specifiche
                    old_St.append(sottostringa)
    # Troviamo St come l'insieme di tutte le stringhe da old_St per cui nessuna sottostringa più corta è una sottostringa
    St = [s for s in old_St if all(s not in other_str for other_str in old_St if len(other_str) < len(s))]

    return St, len(St)

if __name__ == "__main__":
    main()
