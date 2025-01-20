import csv
from Bio import Entrez, SeqIO
from typing import List, Dict


#---------------------------------------------------------------------------
#  Global variables
#---------------------------------------------------------------------------

# CHANGE: Plik z nazwami taksonów
TAXON_FILE = "Giezgala.txt"

# Rozpatrywane markery molekularne
GENES = ["ABCA1", "ADORA3", "AFF2", "AFF2.2",
        "APP", "ATXN7", "AXIN1", "BCHE",
        "BCOR", "BDNF", "BRCA2", "CFTR",
        "CHRNA1", "CNR1", "CREM", "DACH1",
        "DCTN2", "DMRT1", "EDG1", "ERC2",
        "FAM123B", "FBN1", "FES", "FOXP1",
        "GHR", "KCNMA1", "LRPPRC_169",
        "LRPPRC_171", "LUC7L", "MAPKAP1",
        "MBD5", "NEGR1", "NPAS3", "NPAS3.2",
        "PLCB4", "PNOC", "POLA1", "RAB6IP1",
        "RAG1", "RAG2", "RPGRIP1", "SGMS1",
        "SIM1", "SMCX", "SMCY", "SRY", "TEX2",
        "TTR", "TYR", "USH2A", "UTY", "ZFX",
        "ZFY", "ZIC3"
        ]


CHROMOSOM_MARKER = {gen: None for gen in GENES}
# CHROMOSOM_MARKER = {'ABCA1': 'A', 'ADORA3': None, 'AFF2': 'X', 'AFF2.2': None, 'APP': 'A', 'ATXN7': 'A', 'AXIN1': 'A', 'BCHE': 'A', 'BCOR': 'X', 'BDNF': 'A', 'BRCA2': 'A', 'CFTR': 'A', 'CHRNA1': 'A', 'CNR1': 'A', 'CREM': 'A', 'DACH1': 'A', 'DCTN2': 'A', 'DMRT1': 'A', 'EDG1': None, 'ERC2': 'A', 'FAM123B': None, 'FBN1': 'A', 'FES': 'A', 'FOXP1': 'A', 'GHR': 'A', 'KCNMA1': 'A', 'LRPPRC_169': None, 'LRPPRC_171': None, 'LUC7L': 'A', 'MAPKAP1': 'A', 'MBD5': 'A', 'NEGR1': 'A', 'NPAS3': 'A', 'NPAS3.2': None, 'PLCB4': 'A', 'PNOC': 'A', 'POLA1': 'X', 'RAB6IP1': None, 'RAG1': 'A', 'RAG2': 'A', 'RPGRIP1': 'A', 'SGMS1': 'A', 'SIM1': 'A', 'SMCX': None, 'SMCY': 'Y', 'SRY': 'Y', 'TEX2': 'A', 'TTR': 'A', 'TYR': 'A', 'USH2A': 'A', 'UTY': 'Y', 'ZFX': 'X', 'ZFY': 'Y', 'ZIC3': 'X'}

def read_taxa_from_file(file_path: str) -> List[str]:
    """Funkcja odczytuje taxony z pliku tekstowego.
    
    :param file_path: Ścieżka do pliku, zawierającego nazwy gatunków
    :return: Lista taksonów
    """
    taxons = []
    try:
        with open(file_path, "r") as file:
            for line in file:
                taxon = line.strip()

                if taxon:
                    taxons.append(taxon.lower())

        print(f"W pliku było {len(taxons)} taksonów.")

    except FileNotFoundError:
        print(f"Plik {file_path} nie znaleziony.")
    except Exception as e:
        print(f"ERROR: {e}")
    
    return taxons


def set_gene_location(taxon: str, gene: str) -> bool:
    """ Funkcja sprawdza czy marker jest autosomalny,
    czy sprzężony z odpowiednim chromosomem płci.

    :param organism: Nazwa gatunku
    :param gene: Nazwa genu
    :return: bool czy znaleziono gen dla danego gatunku
    """
    try:
        # Wyszukanie ID genu
        query = f"{gene}[Gene] AND {taxon}[Organism]"
        handle = Entrez.esearch(db="gene", term=query, retmax=1)
        record = Entrez.read(handle)
        handle.close()

        if record["IdList"]:
            gene_id = record["IdList"][0]
            
            # Informacja o chromosomie
            handle = Entrez.efetch(db="gene", id=gene_id, rettype="xml", retmode="text")
            records = Entrez.read(handle)
            handle.close()

            chromosome = records[0]["Entrezgene_source"]["BioSource"]["BioSource_subtype"][0]["SubSource_name"]
            if chromosome == 'Un': return False
            if chromosome not in ('X', 'Y'): chromosome = 'A'

            # Marker Autosomalny, X, Y
            if CHROMOSOM_MARKER[gene] is None:
                CHROMOSOM_MARKER[gene] = chromosome
            elif CHROMOSOM_MARKER[gene] != chromosome :
                print(f"UWAGA: Inny chromosom dla gatunku {taxon}")
                print(f"Zaleca się sprawdzenie genu {gene} w bazie NCBI")
        else:
            return False

    except Exception as e:
        return f"ERROR: {e}"


def fetch_gene_sequence(taxon: str,
                        gene: str,
                        output_file: str = "sequences.fasta",
                        ) -> bool:
    """ Funkcja zapisuje znalezioną sekwencję w pliku FASTA
    
    :param taxon: Nazwa gatunku
    :param gene: Nazwa genu
    :param sequences: Słownik przyporządkowanych genów do gatunku
    :return: bool czy znaleziono sekwencje nukleotydów genu dla danego gatunku
    """
    try:
        # Wyszukaj gen w bazie Nucleotide
        query = f"{gene}[Gene] AND {taxon}[Organism] AND biomol_genomic[PROP] NOT whole genome shotgun"
        handle = Entrez.esearch(db="nucleotide", term=query, retmax=1)
        record = Entrez.read(handle)
        handle.close()

        if record["IdList"]:
            # Zapisanie sekwencji
            seq_id = record["IdList"][0]
            handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
            seq_record = SeqIO.read(handle, "fasta")
            handle.close()

            with open(output_file, "a") as fasta_file:
                SeqIO.write(seq_record, fasta_file, "fasta")

            return True
        else:
            return False

    except Exception as e:
        print(f"ERROR: {e}")
        return False


def fasta_to_csv(output_csv: str,
                 genes_list: List[str],
                 species_list: List[str], 
                 fasta_file: str = "sequences.fasta") -> None:
    """ Funkcja tworzy plik csv zliczeń wystąpień w pliku FASTA. Dla każdego
    genu w liście genes_list policzy, dla ilu gatunków występuję sekwencja.
    W każdej kolumnie (zawierającej także aliasy) będzie 1, jeśli znajduje się
    sekwencja oraz 0 jeśli nie.
    
    :param output_csv: Nazwa tworzonego pliku CSV
    :param genes_list: Lista genów do sprawdzenia
    :param species_list: Lista gatunków do sprawdzenia
    :param fasta_file: plik FASTA z zapisanymi sekwencjami
    """
    # CHANGE: Analysis taxons in file
    species_list.append('ateles belzebuth hybridus')
    species_list.append('carlito syrichta')
    
    # Słownik zliczeń
    data = {gene: {species: 0 for species in species_list} for gene in genes_list}

    # Przetwarzanie pliku FASTA
    for record in SeqIO.parse(fasta_file, "fasta"):
        description = record.description.lower()
        for gene in genes_list:
            if gene.lower() in description:
                for species in species_list:
                    if species.lower() in description:
                        data[gene][species] = 1

    # Tworzenie pliku CSV
    with open(output_csv, mode="w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile)
        
        # Nagłówki
        headers = ["Gene", "Chromosome Type", "Species Count"] + species_list
        writer.writerow(headers)

        # Wiersze danych
        for gene, species_dict in data.items():
            row = [gene]
            chromosome_type = CHROMOSOM_MARKER.get(gene, "NONE")
            row.append(chromosome_type)
            species_count = sum(species_dict.values())
            row.append(species_count)
            row.extend(species_dict.values())
            writer.writerow(row)

    print(f"Plik zapisany jako: {output_csv}")


def select_best_genes(output_fasta_prefix: str = "marker_genes_",
                      genes_file: str = "genes.csv", 
                      fasta_file: str = "sequences.fasta") -> Dict[str, str]:
    """ Funkcja wybiera geny o największej liczbie zliczeń z każdego rodzaju (A, X, Y)
    i zapisuje ich sekwencje do osobnych plików FASTA.
    
    :param genes_file: Plik CSV z danymi genów (2 kolumna ma sumę zliczeń)
    :param fasta_file: Plik FASTA zawierający wszystkie analizowane sekwencje
    :param output_fasta_prefix: Prefiks dla plików wynikowych FASTA (np. "top_gene_").
    :return: Słownik {"rodzaj_chromosomu": "plik.fasta"}.
    """
    # Przechowywanie najlepszych genów
    best_genes = {"A": None, "X": None, "Y": None}
    max_counts = {"A": 0, "X": 0, "Y": 0}

    # Wczytanie pliku CSV
    with open(genes_file, "r", encoding="utf-8") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            gene = row["Gene"]
            chromosome = row["Chromosome Type"]
            count = int(row["Species Count"])
            
            if chromosome in best_genes and count > max_counts[chromosome]:
                best_genes[chromosome] = gene
                max_counts[chromosome] = count

    # Zapisanie sekwencji do plików FASTA
    output_files = {}
    for chromosome, gene in best_genes.items():
        if gene:
            output_file = f"{output_fasta_prefix}{chromosome}.fasta"
            output_files[chromosome] = output_file

            with open(output_file, "w") as outfile:
                for record in SeqIO.parse(fasta_file, "fasta"):
                    if gene in record.description:
                        SeqIO.write(record, outfile, "fasta")

    return output_files


def main():
    # Dostęp do NCBI
    Entrez.email = "KołoNaukowe@Bioinformatyki.UW"
    
    # Wszystkie taxony biorące udział w analizie
    taxons = read_taxa_from_file(TAXON_FILE)
    
    for i in range(len(taxons)):
        for j in range(len(GENES)):
            # fetch_gene_sequence(taxons[i], GENES[j])
            set_gene_location(taxons[i], GENES[j])

    fasta_to_csv("genes.csv", GENES, taxons)

    print(select_best_genes())


if __name__ == "__main__":
    main()
    print(CHROMOSOM_MARKER)
