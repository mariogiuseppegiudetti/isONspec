# README - isONspec

## Descrizione
Questo script Python è progettato per eseguire il clustering di sequenze di read genomiche. Il clustering viene eseguito utilizzando un algoritmo basato su confronti di stringhe specifiche (SFS) tra le sequenze di read. L'algoritmo confronta ogni read con i rappresentanti dei cluster esistenti e decide se assegnare la read a un cluster esistente o crearne uno nuovo.

## Requisiti
- Python 3.x
- Libreria `argparse` (solitamente inclusa nelle distribuzioni standard di Python)

## Utilizzo

```bash
python nome_script.py <percorso_del_file_da_processare> --threads <numero_di_threads>
```

### Argomenti
- `<percorso_del_file_da_processare>`: Il percorso completo del file contenente le sequenze di read da processare. Il file può essere nel formato FASTA o FASTQ.
- `--threads <numero_di_threads>` (Opzionale): Il numero di threads da utilizzare per il calcolo del clustering (predefinito: 1).

## Output

L'output dello script sarà una serie di directory, ognuna contenente i file relativi a un cluster specifico. Ogni directory avrà il nome "cluster_\<rappresentante>" e conterrà un file di reads associate al cluster (cluster_\<rappresentante>_reads.fasta o cluster_\<rappresentante>_reads.fastq) e, se applicabile, un file di qualità associato (cluster_\<rappresentante>_reads.qual).

Un file di log (`log.txt`) verrà creato nella directory di esecuzione contenente i valori di confronto tra le read e i rappresentanti.

## Parametri

L'algoritmo di clustering utilizza il numero massimo di stringhe specifiche (SFS) e una soglia per la differenza di lunghezza delle sequenze. Questi parametri possono essere adattati nel codice sorgente.

```python
# Nel corpo dello script
valore_soglia = 0.004  # Modificare questo valore per regolare la soglia di differenza di lunghezza
```
