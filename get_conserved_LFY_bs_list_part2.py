
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

import re
from interdistances_functions import *
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
	matrixType = "score"

# These 3 lines allows to retrieve the matrix from the file
F = open(MatrixFile,"r")
matrix = F.read()
F.close()

# These 3 lines allows to retrieve all the individual frequency values (retrieve only floats) from the matrix and put them in order into a list

num = re.compile(r"([+-]?\d+[.,]\d+)")
Mdata = num.findall(matrix)
matScore, lenMotif = get_score_matrix(Mdata,matrixType)

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
dependency_matrix = get_dependency_matrix(dependencyFile,num)

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



