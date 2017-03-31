
### ----------------------------
### LFY BS discovery program
### ----------------------------

'''
This program allows to identify LFY binding sites (BS) that are conserved among the 1001 A. thaliana genomes 
(http://1001genomes.org/data/GMI-MPI/releases/v3.1/intersection_snp_short_indel_vcf/).
You need a matrix with frequency values and fasta sequences (bound sequences (i.e. peaks)). Here we use ChIP-Seq sequences for LFY.
This program was written by Adrien Bessy, Arnaud Stigliani and Francois Parcy, and was inspired by Morpheus program written by Eugenio Gomez Minguet
and (4-scores.py and CalculateScoreAboveTh_and_Interdistances.py) programs written by Laura Gregoire.
'''

import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 
from interdistances_functions import *
import re
import vcf
import pysam
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
parser.add_argument("--pseudoCount", "-pc",type = float, default = 0.001)
args = parser.parse_args()

#python get_conserved_LFY_bs_list_part1.py -fac "ARF2" -th -12

factorTranscription = args.factor
threshold = args.threshold
pseudoCount = args.pseudoCount
                    
###################Parameters we can change#################
	
if factorTranscription == "LFY_scores_matrix_19nucl" :
	FastaFile = "LFY_bound_sequences.fas" 
	MatrixFile = "LFY_scores_matrix_19nucl.txt" 
	dependencyFile = "interdependent_bases_matrix_for_LFY.txt"
	matrixType = "score"
	
if factorTranscription == "ARF2" :
	FastaFile = "../sequences/ARF2_bound_sequences.fas" 
	MatrixFile = "../matrix/ARF2_OMalley_matrixC.txt" 
	matrixType = "freq" 
	
if factorTranscription == "ARF5" :
	FastaFile = "../sequences/ARF5_bound_sequences.fas" 
	MatrixFile = "../matrix/ARF5_OMalley_Cmatrix.txt" 
	matrixType = "freq" 
	
################################################################

def get_LFY_binding_sites(matScore,matRev,FastaFile,threshold,factorTranscription):
	
	# This line allows to retrieve all the sequences from the fasta file
	sequences = SeqIO.to_dict(SeqIO.parse(FastaFile, "fasta"))
	print "  There are %s sequence(s) to analyze"%(len(sequences))
	
	list_of_the_LFY_binding_sites = []
	# We apply a loop on all the fasta sequences:
	for s in sequences:
		
		# We will store in this list all the best scores (see the threshold after) found for subsequences of one sequence
		good_score_positions = [] 
		
		# This line allows to retrieve the DNA sequence
		seq = sequences[s].seq
		seq_id = sequences[s].id
		chrom = re.split(':',seq_id)
		pos = re.split(':|-',seq_id)
		
		# We look at each sub-sequences of the whole sequence. Each sub-sequence has the same length that the matrix length.
		for c in range(len(seq) - (lenMotif -1)):
			strandPos = seq[c:c+lenMotif].upper()
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
				if factorTranscription == "LFY_scores_matrix_19nucl" :
					scoreStrandPos, scoreStrandNeg = add_scores_associated_with_interdependent_positions(get_dependency_matrix(dependencyFile,num),scoreStrandPos,scoreStrandNeg,strandPos)
				
				#These lines allows to retrieve the chromosome and the positions where there is a predicted binding site (score above the threshold fixed by the user) . 
				if scoreStrandPos > threshold or scoreStrandNeg > threshold:
					list_of_the_LFY_binding_sites.append([chrom[0].replace('chr',''),int(pos[1]) + c + 1, int(pos[1]) + c + 1 + 19,str(strandPos[0:19])])

	return(list_of_the_LFY_binding_sites)

########################################### About the main matrix #######################

''' The sens of the matrix is important: The positions are on the vertical sens and the bases are on the horizontal sens as described in the example.
separation between numbers can be spaces, tabulation, comas... DON'T write your floats with comas (0,2). You have to write points (0.2).

                                                                         Example :   A C G T
                                                                  position 1:           0.16456   0.21614       0.1565,0.1645
                                                                  position 2:           0.645; 0.654    0.155 |||||| 0.4444
                                                                        '''
####################################################################################

# These 3 lines allows to retrieve the matrix from the file
F = open(MatrixFile,"r")
matrix = F.read()
F.close()

# These 3 lines allows to retrieve all the individual frequency values (retrieve only floats) from the matrix and put them in order into a list
num = re.compile(r"([+-]?\d+[.,]\d+)")
Mdata = num.findall(matrix)
matScore, lenMotif = get_score_matrix(Mdata,matrixType,pseudoCount)

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

########## get list[chromosome number, first position of a BS, last position of a BS, 'sequence of the BS']:
list_of_binding_sites = get_LFY_binding_sites(matScore,matRev,FastaFile,threshold,factorTranscription)

print("len(list_of_binding_sites) : ",len(list_of_binding_sites))

chr1_list = [s for s in list_of_binding_sites if s[0] == '1' ]
chr1_min = min(i[1] for i in chr1_list)
chr1_max = max(i[2] for i in chr1_list)
chr2_list = [s for s in list_of_binding_sites if s[0] == '2' ]
chr2_min = min(i[1] for i in chr2_list)
chr2_max = max(i[2] for i in chr2_list)
chr3_list = [s for s in list_of_binding_sites if s[0] == '3' ]
chr3_min = min(i[1] for i in chr3_list)
chr3_max = max(i[2] for i in chr3_list)
chr4_list = [s for s in list_of_binding_sites if s[0] == '4' ]
chr4_min = min(i[1] for i in chr4_list)
chr4_max = max(i[2] for i in chr4_list)
chr5_list = [s for s in list_of_binding_sites if s[0] == '5' ]
chr5_min = min(i[1] for i in chr5_list)
chr5_max = max(i[2] for i in chr5_list)

###### STEP 2: Applying a loop on all the vcf files and retrieve the new sequences
print(1)
def add_modified_sequences_from_accessions():
	vcf_reader = vcf.Reader(filename='../../extraPrograms_or_big_files/VCF_files/1001genomes_snp-short-indel_only_ACGTN.vcf.gz')
	for record in vcf_reader.fetch('1',chr1_min,chr1_max):
		for LFY_bs_list in chr1_list:
			if LFY_bs_list[1] <= record.POS < LFY_bs_list[2] :
				LFY_bs_list.append(LFY_bs_list[3][0:record.POS-LFY_bs_list[1]]+str(record.ALT[0])+LFY_bs_list[3][record.POS-LFY_bs_list[1]+1:19])
	print(2)
	for record in vcf_reader.fetch('2',chr2_min,chr2_max):
		for LFY_bs_list in chr2_list:
			if LFY_bs_list[1] <= record.POS < LFY_bs_list[2] :
				LFY_bs_list.append(LFY_bs_list[3][0:record.POS-LFY_bs_list[1]]+str(record.ALT[0])+LFY_bs_list[3][record.POS-LFY_bs_list[1]+1:19])
	print(3)			
	for record in vcf_reader.fetch('3',chr3_min,chr3_max):
		for LFY_bs_list in chr3_list:
			if LFY_bs_list[1] <= record.POS < LFY_bs_list[2] :
				LFY_bs_list.append(LFY_bs_list[3][0:record.POS-LFY_bs_list[1]]+str(record.ALT[0])+LFY_bs_list[3][record.POS-LFY_bs_list[1]+1:19])
	print(4)
	for record in vcf_reader.fetch('4',chr4_min,chr4_max):
		for LFY_bs_list in chr4_list:
			if LFY_bs_list[1] <= record.POS < LFY_bs_list[2] :
				LFY_bs_list.append(LFY_bs_list[3][0:record.POS-LFY_bs_list[1]]+str(record.ALT[0])+LFY_bs_list[3][record.POS-LFY_bs_list[1]+1:19])
	print(5)			
	for record in vcf_reader.fetch('5',chr5_min,chr5_max):
		for LFY_bs_list in chr5_list:
			if LFY_bs_list[1] <= record.POS < LFY_bs_list[2] :
				LFY_bs_list.append(LFY_bs_list[3][0:record.POS-LFY_bs_list[1]]+str(record.ALT[0])+LFY_bs_list[3][record.POS-LFY_bs_list[1]+1:19])
	list_of_the_LFY_binding_sites = chr1_list + chr2_list + chr3_list + chr4_list + chr5_list		
	return(list_of_the_LFY_binding_sites)
	
list_of_binding_sites = add_modified_sequences_from_accessions()

print("len(list_of_binding_sites) : ",len(list_of_binding_sites))

text_file = open("list_of_the_LFY_binding_sites.txt", "w")
text_file.write(str(list_of_ARF2_binding_sites))
text_file.close()

# loop on list_of_the_LFY_binding_sites to retrieve the sequences and write to fasta file.		
# Redo the get_list_LFY_binding_sites function with the new list
# list_of_the_LFY_binding_sites = get_list_LFY_binding_sites(matScore,matRev,FastaFile,dependency_matrix,threshold,factorTranscription)

# You will only have the conserved LFY BS ! 

