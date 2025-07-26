# raw/download_spike_fasta.py
from Bio import SeqIO
from Bio import Entrez

Entrez.email = "avaneesh6404@gmail.com"
handle = Entrez.efetch(db="protein", id="YP_009724390.1", rettype="fasta", retmode="text")
with open("data/raw/spike_protein.fasta", "w") as f:
    f.write(handle.read())
