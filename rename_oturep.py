import sys
import re

#an empty list to hold the lines, as strings, of the input fasta file
info = []

#iterate through the input file and file up the list
with open(sys.argv[1], 'r') as fasta_in:
    for line in fasta_in:
        info.append(line.replace('\r', '').rstrip('\n'))

#an empty list to hold the OTU fasta seqs and their IDs as tuples
seqs = []

#iterate through the info list, replace the complex seq IDs with just the OTU ID, and fill up the seqs list
for index, item in enumerate(info):
    if index == 0:
        ID = re.search('Otu[0-9]*',item[1:]).group(0)
    else:
        if index % 2 == 0:
            ID = re.search('Otu[0-9]*',item[1:]).group(0)
        else:
            seqs.append((ID, item))

#write the contents of the seqs list to a new file called "clean_repFasta.fasta"
with open("clean_repFasta.fasta", 'w') as fasta_out:
    for item in seqs:
        fasta_out.write('>' + item[0] + '\n' + item[1] + '\n')
