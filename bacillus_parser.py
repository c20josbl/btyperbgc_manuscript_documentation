
import os
import sys 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pathlib import Path

genbank_file = sys.argv[1]
p = Path(genbank_file)
genbank_filename = p.stem
input = SeqIO.parse(genbank_file, "genbank")
with open("bacillus_antismash.gbk",'a') as output:
        for record in input:
                record.name = "AS." + genbank_filename
                SeqIO.write(record, output, "genbank")

