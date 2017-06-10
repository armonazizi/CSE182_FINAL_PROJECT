'''
Make keyword tsv file from script output
'''

import subprocess
import os

subprocess.call(['python', 'prosite_script.py'])
subprocess.call(['python', 'BLAST_Pfam_Script.py', '-f', 'newDNA.txt', '-d', 'BLAST_tmp.txt'])
subprocess.call(['python', 'BLAST_Pfam_Script.py', '-f', 'newDNA.txt', '-p', 'Pfam_tmp.txt'])


prosite_keys = open("prosite_keywords.txt", 'r')

pfam_keys = open("Pfam_keywords.txt", 'r')

blast_keys = open("BLAST_keywords.txt", 'r')

#fasta = open("DNA.txt")

prosite_arr = [""]*100
pfam_arr = [""]*100
blast_arr = [""]*100
seq_ids = [""]*100

counter=0

for line in prosite_keys:
    if line.startswith('>'):
        seq_ids[counter] = ((line.split(':')[1]).split(' ')[0]).strip()
    else:
        prosite_arr[counter] = line.strip()
        counter += 1

counter = 0

for line in pfam_keys:
    if not line.startswith('>'):
        pfam_arr[counter] = line.strip()
        counter += 1

counter = 0

for line in blast_keys:
    if not line.startswith('>'):
        blast_arr[counter] = line.strip()
        counter += 1
        
        
out_file = open("keywords.tsv", 'w')

out_file.write("SEQUENCE_ID\tBLAST_KEYS\tPFAM_KEYS\tPROSITE_KEYS\n")

for i in range(0, 100):
    out_file.write(seq_ids[i]+'\t'+blast_arr[i]+
                   '\t'+pfam_arr[i]+'\t'+prosite_arr[i]+'\n')

out_file.close()    

os.remove("BLAST_tmp.txt")
os.remove("Pfam_tmp.txt")