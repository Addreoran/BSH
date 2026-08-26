def read_tab_hmmer(file):
  pass

def download_fasta(proteins):
  pass


def significant_bsh():
  return {2:"C", 18:"R", 21:"D", 82:"N", 228:"R"}


def select_proteins_with_positions(protein_positions):
  pass

def save_fasta(proteins):
  pass

def calc_mafft(proteins):
  pass

def run_hmmer():
  #?
  pass



#assign proteins and count reads, compare them and calc corr


def run(proteins_files):
  #to jest do szczegółowej analizy BSH. 
  # weryfikacja hmmer czy wynik jest ok
  # 1) wczytanie wyników
  hmmer_proteins_uniprotacc=read_tab_hmmer(proteins_files)
  # 2) pobranie fasta
  # 3) mafft
  # 4) lokalizacja białek po significant_bsh
  # 5) ew. lokalizacja białek z opcjonalnymi aminokwasami
  # 6) policzenie maffta 
  # 7) zapisanie do fasta
  
