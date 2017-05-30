'''
Armon

Script for querying prosite
'''

from Bio import SeqIO
from Bio.ExPASy import ScanProsite
from Bio import ExPASy
from Bio.ExPASy import Prodoc



records = open("prosite.doc")

keys = {}

for line in records:
    if line.startswith("{PS"):
        mid = line.find(';')
        keys[line[1:mid]] = line[mid+2:len(line)-2]
        

file = SeqIO.parse("newDNA.txt", "fasta")

output = open("output.txt", 'w')

for sequence in file:
    
   output.write('>' + sequence.id + '\n')
   
   print '>' + sequence.id
   
   handle = ScanProsite.scan(seq=sequence.seq)
   
   text = ScanProsite.scan(seq=sequence.seq, output='txt')
   
   print text
   
   with open("prosite_html/"+(sequence.id)[9:]+".html", "w") as outfile:
       print >>outfile, text.read()
   
   result = ScanProsite.read(handle)
   
   for record in result:
       output.write(keys[record['signature_ac']]+'\n')
       
       print keys[record['signature_ac']]





