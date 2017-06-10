# BLAST_Pfam_Script 
# Query online protein databases

import sys 
import subprocess
from prody import *
from matplotlib.pylab import *
import numpy as np
import json

def runq(query, output, blast):

	in_seq = False # False until loop reaches first header
	header = ''
	seq = ''
	num_lines = 0
	line_count = 0

	with open(query, 'r') as q:	#	Number of lines in query file
		for line in q:
			num_lines += 1
	with open(query, 'r') as q, open(output, 'w') as o:	#open both query and output files
   		for line in q:
   			line_count += 1
   			if line[0] == '>' or line_count == num_lines:	# Header or last line has been reached
   				if line_count == num_lines:
   					seq += line
   				if in_seq == False:	# First header has been reached
   					header = line
   					in_seq = True
   				else:	# non-first header, print previous header and seq to output 
   					if blast:
   						try:
		   					tp = open('temp.fasta', 'w')
		   					tp.truncate()	#	Delete previous information
		   					tp.write(header)
		   					tp.write(seq)
		   					tp.close()
		   					archive = './BLAST_files/' + header[10:18] + '_raw.txt'

		   					blast_out = subprocess.check_output(['blastp', '-query', 'temp.fasta', '-remote', '-db', 'refseq_protein', '-outfmt', "6 qacc sacc evalue stitle", '-max_target_seqs', '10'])
		   					o.write(header)
		   					o.write(str(blast_out, 'utf-8'))

		   					subprocess.call(['blastp', '-query', 'temp.fasta', '-remote', '-db', 'refseq_protein', '-out', archive])	#	Query and save raw results in archive format
		   				except Exception as e:
		   					print(e)
		   					o.write(header)
		   					o.write('Error Occurred\n')
	   				else:
	   					try:
	   						pfam_out = searchPfam(seq) # Query Pfam online
	   						o.write(header)
	   						for x in pfam_out:
	   							o.write(str(x) + ' ')
	   							o.write(json.dumps(pfam_out) + '\n')
	   							archive = header[10:18]
	   							fetchPfamMSA(x, alignment='seed', format='stockholm', outname=archive, folder='.\Pfam_files')

	   					except Exception as e:
	   						print(e)
	   						o.write(header)
	   						o.write('Error Occurred\n')

   					header = line #	Update header
   					seq = '' #	Reset sequence line 
   			else:
   				seq += line	# Append to the sequence
	return 

def getopts(argv):
	"Parse and Import Command Line Arguments"

	opts = {}	#	Empty dictionary to store key-value pairs.
	while argv: #	While there are arguments left to parse...
		if '-' in argv[0]:	#	Found a "-name value" pair.
			opts[argv[0]] = argv[1]  #	Add key and value to the dictionary.
		argv = argv[1:]  #	Reduce the argument list by copying it starting from index 1.
	return opts

def parsePfam(file):

	header = ' '
	keyword = ' '
	with open(file, 'r') as i, open('Pfam_keywords.txt', 'w') as o:
		for line in i:
			if line[0] == '>':
				o.write(line)
			elif line[0:5] == 'Error':
				o.write('null\n')
			else:
				start = line.find('id') + 6
				end = line.find('"', start)
				keyword = line[start:end]
				o.write(keyword + '\n')
	return

def find_nth(haystack, needle, n):
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+len(needle))
        n -= 1
    return start

def parseBLAST(file):
	length = 0
	with open(file, 'r') as x:
		for line in x:
			length += 1

	with open(file, 'r') as o, open('BLAST_keywords.txt', 'w') as i:
		keyword_dict = {}
		count = 0
		for line in o:
			count += 1
			if count == 1:
				i.write(line)
			elif line[0] == '>' or count == length:
				max_occ = 0
				max_keyword = 'null'
				for key in keyword_dict:
					if int(keyword_dict[key]) > max_occ:
						max_keyword = key
						max_occ = int(keyword_dict[key])
				i.write(max_keyword + '\n')
				keyword_dict.clear() # Reset to dictionary for next blastoutput
				i.write(line)
			else:
				starte = find_nth(line, '\t', 2) + 1
				ende = line.find('\t', starte)
				startk = find_nth(line, '\t', 3) + 1
				endk = line.find('\t', startk)
				if 'e' in line[starte:ende] or '0.0' in line[starte:ende]:
					if line[startk:endk] in keyword_dict:
						value = int(keyword_dict[line[startk:endk]]) + 1
						keyword_dict[line[startk:endk]] = str(value)
					else:
						keyword_dict[line[startk:endk]] = '1'

	return

def main():
	arguments = getopts(sys.argv)

	if '-d' in arguments and '-f' in arguments:
		runq(arguments['-f'], arguments['-d'], True)
		parseBLAST(arguments['-d'])
	elif '-p' in arguments and '-f' in arguments:
		runq(arguments['-f'], arguments['-p'], False)
		parsePfam(arguments['-p'])

	return

if __name__ == '__main__':
	main()