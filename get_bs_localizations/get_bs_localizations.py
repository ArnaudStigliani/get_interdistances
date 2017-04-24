### ----------------------------
### TFBS Interdistances program
### ----------------------------

'''
This program allows to localize transcription factor binding sites on sequences such as promoters.
You need a matrix with frequency values and fasta sequences (bound sequences (i.e. peaks)).
This program was written by Adrien Bessy, Arnaud Stigliani and Francois Parcy, and was inspired by Morpheus program written by Eugenio Gomez Minguet
and (4-scores.py and CalculateScoreAboveTh_and_Interdistances.py) programs written by Laura Gregoire.
'''
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 
from interdistances_functions import *
import numpy as np
from Bio import SeqIO
import time
import sys 
from operator import truediv
import operator
from collections import Counter
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker
import matplotlib.patches as mpatches
from matplotlib import pylab
import types
import argparse
import logging
from optparse import OptionParser
from scipy import stats
import PIL
from PIL import ImageFont
from PIL import Image
from PIL import ImageDraw
from PIL import ImageEnhance
#import pygame

parser = argparse.ArgumentParser()                                               

parser.add_argument("--factor", "-fac", type=str, default= "ARF2")
parser.add_argument("--pseudoCount", "-pc",type=float,default= 0.001)
parser.add_argument("--threshold", "-th",nargs='+',type=int,default= -12)
parser.add_argument("--Interdistance_maxValue", "-maxInter", type=int,default= 20)
args = parser.parse_args()

#python get_bs_localizations.py -fac "ARF5" -pc 0.001 -maxInter 20 -th -12
#python get_bs_localizations.py -fac "ARF2" -pc 0.001 -maxInter 20 -th -15

factorTranscription = args.factor
pseudoCount = args.pseudoCount
threshold = args.threshold
Interdistance_maxValue = args.Interdistance_maxValue
                    
###################Parameters we can change#################

if factorTranscription == "ARF2" :
	#FastaFile = "IR7.fas" 
	FastaFile = "test.fas" 
	MatrixFile = "../matrix/ARF2_OMalley_matrixC.txt" 
	matrixType = "freq" 
	
if factorTranscription == "ARF5" :
	FastaFile = "TMO5_promoter.fas" 
	#FastaFile = "test.fas" 
	MatrixFile = "../matrix/ARF5_OMalley_Cmatrix.txt" 
	matrixType = "freq" 

############################################################################### 
	
def interdistance_calcul(DR,ER,IR,good_score_positions) :
	for first in range(0,len(good_score_positions)-1) :
		firstSubSeq = good_score_positions[first]
			
		for second in range(first+1,len(good_score_positions)) :
			secondSubSeq = good_score_positions[second]

			'''
			Here we do different substractions according to we get a Direct Repeat (DR), an Inverted Repeat (IR) or an Everted Repeat (ER).
			Because In the litterature, The interdistance calculation was made between subsequences from our matrix
			and not between the whole sequences from our matrix.
			So according to your specific interdistance calculations you can change these lines.
			'''
			
			if factorTranscription == "ARF2" :
				if firstSubSeq[1] == secondSubSeq[1] :
					d = ( int(secondSubSeq[0]) +2 ) -( int(firstSubSeq[0]) + lenMotif -2)
					if Interdistance_maxValue >= d >= 0 :
						DR = DR + "DR" + str(d) + ", "
				if firstSubSeq[1] == ">" and secondSubSeq[1] == "<" :
					d = ( int(secondSubSeq[0]) +2 ) -( int(firstSubSeq[0]) + lenMotif -2 )
					if Interdistance_maxValue >= d >= 0 :
						ER = ER + "ER" + str(d) + ", "
				if firstSubSeq[1] == "<" and secondSubSeq[1] == ">" :
					d = ( int(secondSubSeq[0]) +2 ) -( int(firstSubSeq[0]) + lenMotif -2 )
					if Interdistance_maxValue >= d >= 0 :
						IR = IR + "IR" + str(d) + ", "

			if factorTranscription == "ARF5" :
				if firstSubSeq[1] == ">" and secondSubSeq[1] == ">" :
					d = (int(secondSubSeq[0])  ) - (int(firstSubSeq[0]) + 6  )
					if Interdistance_maxValue >= d >= 0 :
						DR = DR + "DR" + str(d) + ": " + str([firstSubSeq,secondSubSeq]) + ", "
				if firstSubSeq[1] == "<" and secondSubSeq[1] == "<" :
					d = (int(secondSubSeq[0]) ) - (int(firstSubSeq[0]) + 6 )
					if Interdistance_maxValue >= d >= 0 :
						DR = DR + "DR" + str(d) + ": " + str([firstSubSeq,secondSubSeq]) + ", "
				if firstSubSeq[1] == ">" and secondSubSeq[1] == "<" :
					print(firstSubSeq)
					print(secondSubSeq)
					d = (int(secondSubSeq[0])  ) - (int(firstSubSeq[0]) + 6 )
					print(d)
					if Interdistance_maxValue >= d >= 0 :
						print(2)
						ER = ER + "ER" + str(d) + ": "+ str([firstSubSeq,secondSubSeq]) +", "
				if firstSubSeq[1] == "<" and secondSubSeq[1] == ">" :
					d = (int(secondSubSeq[0])  ) - (int(firstSubSeq[0]) + 6  )
					if Interdistance_maxValue >= d >= 0 :
						IR = IR + "IR" + str(d) + ": " + str([firstSubSeq,secondSubSeq]) + ", "
		
	return(DR,ER,IR)	
	
def get_interdist(matF,matRev,FastaFile,threshold,factorTranscription,Interdistance_maxValue):
	# This line allows to retrieve all the sequences from the fasta file
	sequences = SeqIO.to_dict(SeqIO.parse(FastaFile, "fasta"))
	print(sequences)

	print "  There are %s sequence(s) to analyze"%(len(sequences)) 

	all_data = []
	sequences_ = ""
	explanation = ["Position of the T of the \"TGTCNN\" or of the N of the \"NNGACA\" if the motif is on the complementary strand" , "strand", "score", "binding site"]
	
	# We look at all the fasta sequences:
	for s in sorted(sequences):
		# We will store in this list all the best scores (see the threshold after) found for subsequences of one sequence
		#if type(threshold) is list:
		good_score_positions = [] 
		
		DR = ""
		ER = ""
		IR = ""
		
		# This line allows to retrieve the DNA sequence
		seq = sequences[s].seq
		seq2 = ""
		seq_id = sequences[s].id
		print("seq_id : ",seq_id)
		a = 0
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
				index = 0
				#These lines allows to calculate a score for one sub-sequence
				scoreStrandPos = 0
				scoreStrandNeg = 0 
				while index < lenMotif:
					if strandPos[index] == 'A':
						scoreStrandPos = scoreStrandPos + matF[index*4]
						scoreStrandNeg = scoreStrandNeg + matRev[index*4]
					elif strandPos[index] == 'C':
						scoreStrandPos = scoreStrandPos + matF[index*4+1]
						scoreStrandNeg = scoreStrandNeg + matRev[index*4+1]
					elif strandPos[index] == 'G':
						scoreStrandPos = scoreStrandPos + matF[index*4+2]
						scoreStrandNeg = scoreStrandNeg + matRev[index*4+2]
					elif strandPos[index] == 'T':
						scoreStrandPos = scoreStrandPos + matF[index*4+3]
						scoreStrandNeg = scoreStrandNeg + matRev[index*4+3]			
					index += 1

				#These lines allows to retrieve the position and the strand where there is a predicted binding site. 
				#You can change the threshold.
				#print("scoreStrandPos : ",scoreStrandPos)
				#print("str(strandPos[3:9]) : ",str(strandPos[3:9]))
				#print("scoreStrandNeg : ",scoreStrandNeg)
				#print("str(strandPos[3:9]) : ",str(strandPos[1:7]))
				if scoreStrandPos > threshold:
					if factorTranscription == "ARF2" :
						seq2 = seq2 + seq[a:c+2].upper() + " "
						a = c+2
						good_score_positions.append([c+3,">",round(scoreStrandPos,2),str(strandPos[2:8])])
					if factorTranscription == "ARF5" :
						seq2 = seq2 + seq[a:c+3].upper() + " "
						a = c+3
						print("c+4 : ",c+4)
						good_score_positions.append([c+4,">",round(scoreStrandPos,2),str(strandPos[3:9])])
				if scoreStrandNeg > threshold:
					if factorTranscription == "ARF2" :
						seq2 = seq2 + seq[a:c+2].upper() + " "
						a = c+2
						good_score_positions.append([c+3,"<",round(scoreStrandNeg,2),str(strandPos[2:8])])
					if factorTranscription == "ARF5" :
						seq2 = seq2 + seq[a:c+1].upper() + " "
						a = c+1
						print("c+2 : ",c+2)
						good_score_positions.append([c+2,"<",round(scoreStrandNeg,2),str(strandPos[1:7])])
		seq2 = seq2 + seq[a:].upper()
		#print("good_score_positions : ",good_score_positions)
		DR,ER,IR = interdistance_calcul(DR,ER,IR,good_score_positions)
		if DR == "" :
			DR = "nothing, "
		if ER == "" :
			ER = "nothing, "
		if IR == "" :
			IR = "nothing, "
			
		sequences_ = sequences_ + seq_id + "\n" + seq2 + "\n\n" + str(explanation) + "\nDR : " + str(DR) + "\nER : "+ str(ER) + "\nIR : " + str(IR) + "\n"  + "\n\nAll the sites:\n" + str(good_score_positions) + "\n\n\n"
	return(sequences_)

########################################### About the main matrix #######################

''' The sens of the matrix is important: The positions are on the vertical sens and the bases are on the horizontal sens as described in the example.
separation between numbers can be spaces, tabulation, comas...

                                                                         Example :   A C G T
                                                                  position 1:           0.16456   0.21614       0.1565,0.1645
                                                                  position 2:           0.645; 0.654    0.155 |||||| 0.4444
                                                                                        ...
                                                                        '''
####################################################################################


###############################################################################################

# These 3 lines allows to retrieve the matrix from the file
F = open(MatrixFile,"r")
matrix = F.read().replace("\r","\n") + "\n"
F.close()

# These 3 lines allows to retrieve all the individual frequency values from the matrix and put them in order into a list
import re
num = re.compile(r"([+-]?\d+[.,]\d+)")
Mdata = num.findall(matrix)

matScore, lenMotif = get_score_matrix(Mdata,matrixType,pseudoCount)

# The following line allows to produce the reversed matrix
'''if we take the example given before : A T G C
			Position 1:      0.4444  0.155  0.654   0.645
			Position 2:      0.1645  0.1565 0.21614 0.16456
Now, we can notice that scores change between the positions 1 and 2, and between A and T, and between G and C.
So we can calculate with this reverse matrix, the score of the complementary strand.
'''
matRev = list(reversed(matScore))
	
########## get INTERDISTANCE VALUES for POSITIVE sets:

sequences_ = get_interdist(matScore,matRev,FastaFile,threshold,factorTranscription,Interdistance_maxValue)

output_presentation = "Here the binding sites are present after a white space. For example, with the sequence \n\"AGTAGTCAT TGT CAGATAGAAAGAAAGAG CAGACAAAAGGATCG\"\nBecause a binding site has a length of 6 bp, this means there is a first binding site: \nTGTCAGA, a second one: CAGATA and a last one: CAGACA.\nThe matrix used is: "+MatrixFile + "\n\n\n\n"
		
#text_file = open("localization_of_the_"+factorTranscription+"_binding_sites_on_"+FastaFile+".txt", "w")
#text_file.write(output_presentation+str(sequences_))
#text_file.close()
black = (0,0,0)
white = (255,255,255)

width = 2000
height = 1000
def draw_underlined_text(draw, pos, text, **options):    
	twidth, theight = draw.textsize(text)
	lx, ly = pos[0], pos[1] + theight
	draw.text(pos, text, **options)
	draw.line((lx, ly, lx + twidth, ly), **options)
img = Image.new("RGBA", (width,height),white)
draw = ImageDraw.Draw(img)
w, h = draw.textsize(output_presentation+str(sequences_))
draw_underlined_text(draw, (50, 150), output_presentation, fill=0)
draw_underlined_text(draw, (50, 300), str(sequences_), fill=128)
draw = ImageDraw.Draw(img)
img.show
img.save("result.png")



