'''
Armon

Script for querying prosite
'''

from Bio import SeqIO
from Bio.ExPASy import ScanProsite


file = SeqIO.parse("newDNA.txt", "fasta")

output = open("output.txt", 'w')


for sequence in file:
    
   output.write('>' + sequence.id + '\n')
   
   handle = ScanProsite.scan(seq=sequence.seq)
   
   result = ScanProsite.read(handle)
   
   output.write(str(result) + '\n')





