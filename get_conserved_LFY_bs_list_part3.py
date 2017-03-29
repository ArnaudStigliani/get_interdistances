import ast
import re
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
import sys 

with open('list_of_the_LFY_binding_sites_PART2.txt', 'r') as f:
	mylist = ast.literal_eval(f.read())

print("len(mylist) : ",len(mylist))

FastaFile = "LFY_bound_sequences.fas"

sequences = SeqIO.to_dict(SeqIO.parse(FastaFile, "fasta"))

print "  There are %s sequence(s) to analyze"%(len(sequences))

# We will store in these lists all the occurences of each kind of interdistances between motifs found in all sequences.
Inter = [0] * 21

# We look at all the fasta sequences:
for s in sequences:
	seq_id = sequences[s].id
	chrom = re.split(':',seq_id)
	chrom = chrom[0].replace('chr','')
	#print("chrom : ",chrom)
	pos = re.split(':|-',seq_id)
	#print("pos : ",pos)
	begin = int(pos[1])
	#print("begin : ",begin)
	end = int (pos[2])
	#print("end : ",end)
	new = []
	for l in mylist :
		if l[0] == chrom and l[1] >= begin and l[2] <= end :
			new.append(l[1])
	for first in range(0,len(new)-1) :
		firstSubSeq = new[first]	
		for second in range(first+1,len(new)) :
			secondSubSeq = new[second]
			d = ( int(secondSubSeq) ) - (int(firstSubSeq) + 19 )
			if 20 >= d >= 0 :
				Inter[d] += 1
print("Inter : ",Inter)

fig = plt.figure(1,figsize= (18,10))
indexes1 = range(20 + 1)
width = 1
ax1 = fig.add_subplot(1,1,1)
ax1.set_xlabel("base pairs between LFY 19 nucleotides matrix", size = 16)
indexes1 = range(20 + 1)
ax1.axis([0, 20 + 1, 0, 120])
ax1.bar(indexes1, Inter , width , color = 'cornflowerblue')
ax1.set_ylabel("occurences", color = 'cornflowerblue', size = 16)
indexes1 = np.arange(20 + 1)
plt.xticks(indexes1 + width * 0.5 , indexes1)
plt.text(8, 90, sys.argv)
plt.show()