
### ----------------------------
### LFY BS discovery program PART2
### ----------------------------

'''
This program allows to identify LFY binding sites (BS) that are conserved among the 1001 A. thaliana genomes 
(http://1001genomes.org/data/GMI-MPI/releases/v3.1/intersection_snp_short_indel_vcf/).
You need a matrix with frequency values and fasta sequences (bound sequences (i.e. peaks)). Here we use ChIP-Seq sequences for LFY.
This program was written by Adrien Bessy, Arnaud Stigliani and Francois Parcy, and was inspired by Morpheus program written by Eugenio Gomez Minguet
and (4-scores.py and CalculateScoreAboveTh_and_Interdistances.py) programs written by Laura Gregoire.
'''

import interdistances_functions
import ast
import numpy as np
from Bio import SeqIO
import argparse
import os
import glob
from tqdm import tqdm
import operator

parser = argparse.ArgumentParser()                                               

parser.add_argument("--factor", "-fac", type=str, default= "LFY_scores_matrix_19nucl")
parser.add_argument("--threshold", "-th",type=int,default= -23)
args = parser.parse_args()

#python get_conserved_LFY_bs_list_part2.py -fac "LFY_scores_matrix_19nucl" -th -23 

factorTranscription = args.factor
threshold = args.threshold
                    
###################Parameters we can change#################
	
if factorTranscription == "LFY_scores_matrix_19nucl" :
	MatrixFile = "LFY_scores_matrix_19nucl.txt" 
	dependencyFile = "interdependent_bases_matrix_for_LFY.txt"

####################### This is useful to calculate interdependances ###########

codigo = { 'ACC': 5, 'ATG': 14, 'AAG': 2, 'AAA': 0, 'ATC': 13, 'AAC': 1, 'ATA': 12,
           'AGG': 10, 'CCT': 23, 'ACT': 7, 'AGC': 9, 'ACA': 4, 'AGA': 8, 'CAT': 19, 
           'AAT': 3, 'ATT': 15, 'CTG': 30, 'CTA': 28, 'CTC': 29, 'CAC': 17, 'ACG': 6,
           'CAA': 16, 'AGT': 11, 'CCA': 20, 'CCG': 22, 'CCC': 21, 'TAT': 51, 'GGT': 43, 
           'TGT': 59, 'CGA': 24, 'CAG': 18, 'CGC': 25, 'GAT': 35, 'CGG': 26, 'CTT': 31, 
           'TGC': 57, 'GGG': 42, 'TAG': 50, 'GGA': 40, 'TAA': 48, 'GGC': 41, 'TAC': 49,
           'GAG': 34, 'TCG': 54, 'TTA': 60, 'GAC': 33, 'CGT': 27, 'TTT': 63,
           'TCA': 52, 'GCA': 36, 'GTA': 44, 'GCC': 37, 'GTC': 45, 'GCG': 38,
           'GTG': 46, 'TTC': 61, 'GTT': 47, 'GCT': 39, 'TGA': 56, 'TTG': 62,
           'TCC': 53, 'TGG': 58, 'GAA': 32, 'TCT': 55}
           
codigoi = { "A" : "T", "C" : "G", "G" : "C", "T" : "A"}
##################################

# These 3 lines allows to retrieve the matrix from the file
F = open(MatrixFile,"r")
matrix = F.read()
F.close()

# These 3 lines allows to retrieve all the individual frequency values (retrieve only floats) from the matrix and put them in order into a list
import re
num = re.compile(r"([+-]?\d+[.,]\d+)")
Mdata = num.findall(matrix)
matScore, lenMotif = get_score_matrix(Mdata)

'''The following line allows to produce the reversed matrix
if we take the example given before : A T G C
			Position 1:      0.4444  0.155  0.654   0.645
			Position 2:      0.1645  0.1565 0.21614 0.16456
Now, we can notice that scores change between the positions 1 and 2, and between A and T, and between G and C.
So we can calculate the score of the complementary strand with this matRev matrix.
'''
matRev = list(reversed(matScore))

'''The LFY matrix gets 3 * 3 bases that are interdependent, so we have to retrieve the interdependances in a list:
dependency_matrix = [[4, 5, 6], [-2.482, -2.8048, -5.8493, -1.9992, -4.463, -5.8493, -6.0225, -6.0225, ...], 
			[9, 10, 11], [-4.3503, -5.1942, -4.3503, -5.1942, -5.1942, -5.1942, -3.434, -5.1942, -3.9448, ...], 
			[14, 15, 16], [-6.0225, -6.0225, -6.0225, -6.0225, -1.1954, -5.8493, -5.8493, -5.8493, -6.0225, -6.0225, -6.0225, ...]]
'''
dependency_matrix = get_dependency_matrix(dependencyFile)

############### new part###########################
with open('list_of_the_LFY_binding_sites.txt', 'r') as f:
	mylist = ast.literal_eval(f.read())

print("len(mylist) : ",len(mylist))
	
conserved_sequences = []	

for i in mylist : 
	done = False
	if len(i) > 4 :
		#print(i)
		for j in range (4,len(i)) :
			seq = i[j]
			#print("seq : ",seq)
			for c in range(len(seq) - (lenMotif -1)):
				strandPos = seq[c:c+lenMotif].upper()
				#print("strandPos : ",strandPos)
				test = 0
				for nu in strandPos :
					if nu not in ["A","C","G","T"]:
						test = 1
				if test == 1:
					score = "NA"
				else :
					#These lines allows to calculate a score for one sub-sequence
					index = 0
					scoreStrandPos = 0
					scoreStrandNeg = 0 
					while index < lenMotif:
						if strandPos[index] == 'A':
							scoreStrandPos = scoreStrandPos + matScore[index*4]
							scoreStrandNeg = scoreStrandNeg + matRev[index*4]
						elif strandPos[index] == 'C':
							scoreStrandPos = scoreStrandPos + matScore[index*4+1]
							scoreStrandNeg = scoreStrandNeg + matRev[index*4+1]
						elif strandPos[index] == 'G':
							scoreStrandPos = scoreStrandPos + matScore[index*4+2]
							scoreStrandNeg = scoreStrandNeg + matRev[index*4+2]
						elif strandPos[index] == 'T':
							scoreStrandPos = scoreStrandPos + matScore[index*4+3]
							scoreStrandNeg = scoreStrandNeg + matRev[index*4+3]			
						index += 1
					
					# This function allows to add scores that are associated with interdependent positions
					scoreStrandPos, scoreStrandNeg = add_scores_associated_with_interdependent_positions(dependency_matrix,scoreStrandPos,scoreStrandNeg,strandPos)
					if scoreStrandPos < threshold and scoreStrandNeg < threshold:
						done = True
						break
			if done:
				break

			if j == len(i)-1 :
				conserved_sequences.append(i)
					
	else : 
		conserved_sequences.append(i)

print("len(conserved_sequences) : ",len(conserved_sequences))
text_file = open("list_of_the_LFY_binding_sites_PART2.txt", "w")
text_file.write(str(conserved_sequences))
text_file.close()						



