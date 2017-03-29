### ----------------------------
### TFBS Interdistances program
### ----------------------------

'''
This program allows to calculate interdistances between transcription factor binding sites.
You need a matrix with frequency values and fasta sequences (bound sequences (i.e. peaks), and unbound sequences).
You can display a plot for several thresholds.
This program was written by Adrien Bessy, Arnaud Stigliani and Francois Parcy, and was inspired by Morpheus program written by Eugenio Gomez Minguet
and (4-scores.py and CalculateScoreAboveTh_and_Interdistances.py) programs written by Laura Gregoire.
'''

import re
import numpy as np
from Bio import SeqIO
import time
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 
from interdistances_functions import *
#sys.path.insert(0, '../')
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

parser = argparse.ArgumentParser()                                               

parser.add_argument("--factor", "-fac", type = str, default= "ARF2")
parser.add_argument("--negative_sets", "-neg", nargs='*')
parser.add_argument("--pseudoCount", "-pc",type = float, default = 0.001)
parser.add_argument("--threshold", "-th",nargs='+',type = int, default= -25 -23- 21)
parser.add_argument("--Interdistance_maxValue", "-maxInter", type = int, default= 20)
parser.add_argument("--histo", "-histo", type = bool, default= False)
parser.add_argument("--sum_threshold", "-sum_threshold", type = bool, default= False)

args = parser.parse_args()

#python get_interdistances.py -fac "ARF2" -pc 0.001 -maxInter 20 -th -20 -histo True -neg ARF2_neg1.fas ARF2_neg2.fas -sum_threshold True
#python get_interdistances.py -fac "ARF5" -pc 0.001 -maxInter 20 -th -15 -histo True -neg
# python get_interdistances.py -fac "LFY_matrix_19nucl" -th -23 -histo True

factorTranscription = args.factor
negative_sets = args.negative_sets
pseudoCount = args.pseudoCount
threshold = args.threshold
Interdistance_maxValue = args.Interdistance_maxValue 
histo = args.histo
sum_threshold = args.sum_threshold
    
if histo == True and len(threshold) > 1 : 
	print("Impossible to display an histogram with several thresholds, you have to change a parameter.")
	sys.exit(0)
	
###################Parameters we can change#################

if factorTranscription == "ARF2" :
	FastaFile = "ARF2_bound_sequences.fas" 
	MatrixFile = "ARF2_OMalley_matrixC.txt" 
	matrixType = "freq" 
	
if factorTranscription == "ARF5" :
	FastaFile = "ARF5_bound_sequences.fas" 
	MatrixFile = "ARF5_OMalley_Cmatrix.txt" 
	matrixType = "freq" 
	
if factorTranscription == "LFY_matrix_19nucl" :
	FastaFile = "../sequences/LFY_bound_sequences.fas"
	MatrixFile = "../matrix/LFY_scores_matrix_19nucl.txt" 
	dependencyFile = "../matrix/interdependent_bases_matrix_for_LFY.txt"
	matrixType = "score" 

###############################################################################       
      
def interdistance_calcul(InterDR,InterER,InterIR,sum_thresold,good_score_positions) :
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
				if sum_threshold == True:
					if int(firstSubSeq[2]) + int(secondSubSeq[2]) > sum_thresold :
						if firstSubSeq[1] == secondSubSeq[1] :
							d = ( int(secondSubSeq[0]) +2 ) -( int(firstSubSeq[0]) + lenMotif -2)
							if Interdistance_maxValue >= d >= 0 :
								InterDR[d] += 1
						if firstSubSeq[1] == ">" and secondSubSeq[1] == "<" :
							d = ( int(secondSubSeq[0]) +2 ) -( int(firstSubSeq[0]) + lenMotif -2 )
							if Interdistance_maxValue >= d >= 0 :
								InterER[d] += 1
						if firstSubSeq[1] == "<" and secondSubSeq[1] == ">" :
							d = ( int(secondSubSeq[0]) +2 ) -( int(firstSubSeq[0]) + lenMotif -2 )
							if Interdistance_maxValue >= d >= 0 :
								InterIR[d] += 1
				else:
					if firstSubSeq[1] == secondSubSeq[1] :
						d = ( int(secondSubSeq[0]) +2 ) -( int(firstSubSeq[0]) + lenMotif -2)
						if Interdistance_maxValue >= d >= 0 :
							InterDR[d] += 1
					if firstSubSeq[1] == ">" and secondSubSeq[1] == "<" :
						d = ( int(secondSubSeq[0]) +2 ) -( int(firstSubSeq[0]) + lenMotif -2 )
						if Interdistance_maxValue >= d >= 0 :
							InterER[d] += 1
					if firstSubSeq[1] == "<" and secondSubSeq[1] == ">" :
						d = ( int(secondSubSeq[0]) +2 ) -( int(firstSubSeq[0]) + lenMotif -2 )
						if Interdistance_maxValue >= d >= 0 :
							InterIR[d] += 1

			if factorTranscription == "ARF5" :
				if sum_threshold == True :
					if int(firstSubSeq[2]) + int(secondSubSeq[2]) > sum_thresold :
						if firstSubSeq[1] == ">" and secondSubSeq[1] == ">" :
							d = ( int(secondSubSeq[0]) +3 ) -( int(firstSubSeq[0]) + lenMotif -1 )
							if Interdistance_maxValue >= d >= 0 :
								InterDR[d] += 1
						if firstSubSeq[1] == "<" and secondSubSeq[1] == "<" :
							d = ( int(secondSubSeq[0]) +1 ) -( int(firstSubSeq[0]) + lenMotif -3)
							if Interdistance_maxValue >= d >= 0 :
								InterDR[d] += 1
						if firstSubSeq[1] == ">" and secondSubSeq[1] == "<" :
							d = ( int(secondSubSeq[0]) +1 ) -( int(firstSubSeq[0]) + lenMotif -1 )
							if Interdistance_maxValue >= d >= 0 :
								InterER[d] += 1
						if firstSubSeq[1] == "<" and secondSubSeq[1] == ">" :
							d = ( int(secondSubSeq[0]) +3 ) -( int(firstSubSeq[0]) + lenMotif -3 )
							if Interdistance_maxValue >= d >= 0 :
								InterIR[d] += 1
				else :
					if firstSubSeq[1] == ">" and secondSubSeq[1] == ">" :
						d = ( int(secondSubSeq[0]) +3 ) -( int(firstSubSeq[0]) + lenMotif -1 )
						if Interdistance_maxValue >= d >= 0 :
							InterDR[d] += 1
					if firstSubSeq[1] == "<" and secondSubSeq[1] == "<" :
						d = ( int(secondSubSeq[0]) +1 ) -( int(firstSubSeq[0]) + lenMotif -3)
						if Interdistance_maxValue >= d >= 0 :
							InterDR[d] += 1
					if firstSubSeq[1] == ">" and secondSubSeq[1] == "<" :
						d = ( int(secondSubSeq[0]) +1 ) -( int(firstSubSeq[0]) + lenMotif -1 )
						if Interdistance_maxValue >= d >= 0 :
							InterER[d] += 1
					if firstSubSeq[1] == "<" and secondSubSeq[1] == ">" :
						d = ( int(secondSubSeq[0]) +3 ) -( int(firstSubSeq[0]) + lenMotif -3 )
						if Interdistance_maxValue >= d >= 0 :
							InterIR[d] += 1
							
			if factorTranscription == "LFY_matrix_19nucl" :
				d = int(secondSubSeq[0]) - (int(firstSubSeq[0])+ lenMotif)
				if Interdistance_maxValue >= d >= 0 :
					InterDR[d] += 1
					InterER[d] += 1
					InterIR[d] += 1
		
	return(InterDR,InterER,InterIR)
	
def get_interdist(matF,matRev,FastaFile,threshold,factorTranscription,Interdistance_maxValue):
	# This line allows to retrieve all the sequences from the fasta file
	sequences = SeqIO.to_dict(SeqIO.parse(FastaFile, "fasta"))

	print "  There are %s sequence(s) to analyze"%(len(sequences))
	
	# We will store in these lists all the occurences of each kind of interdistances between motifs found in all sequences.
	DR = [] 
	ER = [] 
	IR = [] 
	for a in threshold :
		DR.append([0] * (Interdistance_maxValue + 1) )
		ER.append([0] * (Interdistance_maxValue + 1) )
		IR.append([0] * (Interdistance_maxValue + 1) )
	
	score_occurence = 0
	
	# We look at all the fasta sequences:
	for s in sequences:
		# We will store in this list all the best scores (see the threshold after) found for subsequences of one sequence
		#if type(threshold) is list:
		good_score_positions = [] 
		if sum_threshold == False :
			for a in threshold :
				good_score_positions.append( [] )
		bestScore = 0
		positionOfTheBestScore = 0
		# This line allows to retrieve the DNA sequence
		seq = sequences[s].seq
		id=sequences[s].id

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
				n = 0
				#These lines allows to calculate a score for one sub-sequence
				scoreStrandPos = 0
				scoreStrandNeg = 0 
				while n<lenMotif:
					if strandPos[n] == 'A':
						scoreStrandPos = scoreStrandPos + matF[n*4]
						scoreStrandNeg = scoreStrandNeg + matRev[n*4]
					elif strandPos[n] == 'C':
						scoreStrandPos = scoreStrandPos + matF[n*4+1]
						scoreStrandNeg = scoreStrandNeg + matRev[n*4+1]
					elif strandPos[n] == 'G':
						scoreStrandPos = scoreStrandPos + matF[n*4+2]
						scoreStrandNeg = scoreStrandNeg + matRev[n*4+2]
					elif strandPos[n] == 'T':
						scoreStrandPos = scoreStrandPos + matF[n*4+3]
						scoreStrandNeg = scoreStrandNeg + matRev[n*4+3]			
					n += 1
		
				if dependencyFile : 			
					scoreStrandPos, scoreStrandNeg = add_scores_associated_with_interdependent_positions(get_dependency_matrix(dependencyFile,num),scoreStrandPos,scoreStrandNeg,strandPos)

		
				#These lines allows to retrieve the position and the strand where there is a predicted binding site. 
				#You can change the threshold.
				if sum_threshold == True :
					if scoreStrandPos > min(threshold):
						good_score_positions.append([c+1,">",scoreStrandPos])
					if scoreStrandNeg > min(threshold):
						good_score_positions.append([c+1,"<",scoreStrandNeg])
				else :
					for a , b in zip(good_score_positions, threshold) :
						if scoreStrandPos > b:
							score_occurence = score_occurence + 1
							a.append([c+1,">",scoreStrandPos])
						if factorTranscription != "LFY_matrix_19nucl" :
							if scoreStrandNeg > b:
								score_occurence = score_occurence + 1
								a.append([c+1,"<",scoreStrandNeg])
				
		# Once we have stored all the positions, we calculate all the interdistances:
		if sum_threshold == True :
			for interDIR, interEVER, interINVER,sum_thresold in zip(DR,ER,IR,threshold) :
				interdistance_calcul(interDIR,interEVER,interINVER,sum_thresold,good_score_positions)
		else :
			for goodScores, interDIR, interEVER, interINVER in zip(good_score_positions,DR,ER,IR) :
				InterDR,InterER,InterIR = interdistance_calcul(interDIR,interEVER,interINVER,threshold,goodScores)
	print("score_occurence : ",score_occurence)
	return(DR,ER,IR)

########################################### About the main matrix #######################

''' The sens of the matrix is important: The positions are on the vertical sens and the bases are on the horizontal sens as described in the example.
separation between numbers can be spaces, tabulation, comas...

                                                                         Example :   A C G T
                                                                  position 1:           0.16456   0.21614       0.1565,0.1645
                                                                  position 2:           0.645; 0.654    0.155 |||||| 0.4444
                                                                                        ...
                                                                        '''
####################################################################################

################### To capture file names where there are unbound sequences ###################

#FastaFileNnumber = input("\nHow many fasta files with unbound sequences ? ")
#d = {}
#for i in range (1,FastaFileNnumber+1) :
	#d["FastaFileN{0}".format(i)] = raw_input("\nName of fasta file with unbound sequences: ")

###############################################################################################

# These 3 lines allows to retrieve the matrix from the file
F = open(MatrixFile,"r")
matrix = F.read().replace("\r","\n") + "\n"
F.close()

# These 3 lines allows to retrieve all the individual frequency values from the matrix and put them in order into a list

num = re.compile(r"([+-]?\d+[.,]\d+)")
Mdata = num.findall(matrix)

matScore, lenMotif = get_score_matrix(Mdata,matrixType)

# The following line allows to produce the reversed matrix
'''if we take the example given before : A T G C
			Position 1:      0.4444  0.155  0.654   0.645
			Position 2:      0.1645  0.1565 0.21614 0.16456
Now, we can notice that scores change between the positions 1 and 2, and between A and T, and between G and C.
So we can calculate with this reverse matrix, the score of the complementary strand.
'''
matRev = list(reversed(matScore))
	
########## get INTERDISTANCE VALUES for POSITIVE sets:

InterDR, InterER, InterIR = get_interdist(matScore,matRev,FastaFile,threshold,factorTranscription,Interdistance_maxValue)

##### For the negative set:
InterDR_N = []
InterER_N = []
InterIR_N = []
lenThr = 0
listThr = []
for a in threshold :
	InterDR_N.append( [0] * (Interdistance_maxValue + 1) )
	InterER_N.append( [0] * (Interdistance_maxValue + 1) )
	InterIR_N.append( [0] * (Interdistance_maxValue + 1) )	
	listThr.append(lenThr)
	lenThr = lenThr + 1

########## get INTERDISTANCE VALUES for NEGATIVE sets				
		
if negative_sets :
	for fastafileN in negative_sets :
		print("fastafileN : ",fastafileN)
		InterDR_N_temp, InterER_N_temp, InterIR_N_temp = get_interdist(matScore,matRev,fastafileN,threshold,factorTranscription,Interdistance_maxValue)

		for a,b,c,d in zip(InterDR_N_temp,InterER_N_temp,InterIR_N_temp,listThr) :
			InterDR_N[d] = [x + y for x, y in zip(InterDR_N[d], a)]
			InterER_N[d] = [x + y for x, y in zip(InterER_N[d], b)]
			InterIR_N[d] = [x + y for x, y in zip(InterIR_N[d], c)]

	if len(negative_sets) > 0 :
		for a,b,c,d in zip(InterDR_N,InterER_N,InterIR_N,listThr) :
			InterDR_N[d] = [x / float(len(negative_sets)) for x in a]
			InterER_N[d] = [x / float(len(negative_sets)) for x in b]
			InterIR_N[d] = [x / float(len(negative_sets)) for x in c]
		
interdist_sum = []
interdist_sum_N = []
if negative_sets :
	for a,b,c,d,e,f in zip(InterDR,InterER,InterIR,InterDR_N,InterER_N,InterIR_N) :
		interdist_sum.append(sum(a) + sum(b) + sum(c))
		interdist_sum_N.append(sum(d) + sum(e) + sum(f))	
	
relative_DR = []
relative_ER = []
relative_IR = []
relative_DR_neg = []
relative_ER_neg = []
relative_IR_neg = []

if negative_sets :
	for a,b,c,d,e in zip(threshold,InterDR,interdist_sum,InterDR_N,interdist_sum_N) :
		relative_DR.append( [x / float(c) for x in b] )
		if len(negative_sets) > 0 :
			relative_DR_neg.append( [x / float(e) for x in d] )
	for a,b,c,d,e in zip(threshold,InterER,interdist_sum,InterER_N,interdist_sum_N) :
		relative_ER.append( [x / float(c) for x in b] )
		if len(negative_sets) > 0 :
			relative_ER_neg.append( [x / float(e) for x in d] )
	for a,b,c,d,e in zip(threshold,InterIR,interdist_sum,InterIR_N,interdist_sum_N) :
		relative_IR.append( [x / float(c) for x in b] )
		if len(negative_sets) > 0 :
			relative_IR_neg.append( [x / float(e) for x in d] )
		
fig = plt.figure(1,figsize= (18,10))
indexes1 = range(Interdistance_maxValue + 1)
width = 1

if negative_sets and len(negative_sets) > 0 :
	for a,b in zip(relative_DR,relative_DR_neg) :
		ax1 = fig.add_subplot(1,3,1)
		ax1.set_xlabel("base pairs between "+factorTranscription+" direct repeats", size = 16)
		if histo == True :
			indexes1 = range(Interdistance_maxValue + 1)
			ax1.axis([0, Interdistance_maxValue + 1, 0, 5.5])
			ax1.bar(indexes1, map(divide, a, b) , width , color = 'cornflowerblue')
			ax1.set_ylabel(' ', color = 'cornflowerblue', size = 16)
			ax2 = ax1.twinx()
			ax2.axis([0, Interdistance_maxValue + 1, 0, 5.5])
			ax2.bar(indexes1,[x * float(10) for x in a] , width , color = 'b')
			ax2.bar(indexes1,[x * float(10) for x in b] , width, color = 'brown', alpha = 0.85)
			ax2.set_ylabel('D$Rn_+$ frequence / D$Rn_-$ frequence', color = 'cornflowerblue', size = 16)
			ybox1 = TextArea("D$Rn_+$ frequence (X10)", textprops = dict(color = "b", size = 16,rotation = 90,ha = 'left',va = 'bottom'))
			ybox3 = TextArea("D$Rn_-$ frequence (X10), ", textprops = dict(color = "brown", size = 16,rotation = 90,ha = 'left',va = 'bottom'))
			ybox = VPacker(children = [ybox1, ybox3],align = "bottom", pad = 0, sep = 8)
			anchored_ybox = AnchoredOffsetbox(loc = 8, child = ybox, pad = 0., frameon = False, bbox_to_anchor = (-0.08, 0.3), 
					bbox_transform = ax2.transAxes, borderpad = 0.)
			ax2.add_artist(anchored_ybox)
			indexes1 = np.arange(Interdistance_maxValue + 1)
			plt.xticks(indexes1 + width * 0.5 , indexes1)
			plt.text(2, 5, "D$Rn_+$ = DRn number in the bound set")
			plt.text(2, 4, "D$Rn_-$ = DRn number in the unbound set")
			plt.text(2, 3, "DRn frequence = DRn / $\sum_{i=0}^{i=20} DR_i + ER_i + IR_i$")
			plt.text(2, 6, sys.argv)
		else :
			ax1.axis([0, Interdistance_maxValue + 1, 0, 5.5])
			plt.plot(indexes1, map(divide, a, b), lw=2)
			ax1.set_ylabel("$DRn_+$ frequence / $DRn_-$ frequence", size = 16)
		l = plt.axhline(y = 1)
		
	for a, b,c in zip(relative_ER,relative_ER_neg,threshold)  :	
		ax1 = fig.add_subplot(1,3,2)
		ax1.set_xlabel("base pairs between "+factorTranscription+" everted repeats", size = 16)
		if histo == True :
			indexes1 = range(Interdistance_maxValue + 1)
			ax1.axis([0, Interdistance_maxValue + 1, 0, 5.5])
			ax1.bar(indexes1, map(divide, a, b) , width , color = 'cornflowerblue')
			ax1.set_ylabel(' ', color = 'cornflowerblue', size = 16)
			ax2 = ax1.twinx()
			ax2.axis([0, Interdistance_maxValue + 1, 0, 5.5])
			ax2.bar(indexes1,[x * float(10) for x in a] , width , color = 'b')
			ax2.bar(indexes1,[x * float(10) for x in b] , width, color = 'brown', alpha = 0.85)
			ax2.set_ylabel('E$Rn_+$ frequence / D$Rn_-$ frequence', color = 'cornflowerblue', size = 16)
			ybox1 = TextArea("E$Rn_+$ frequence (X10)", textprops = dict(color = "b", size = 16,rotation = 90,ha = 'left',va = 'bottom'))
			ybox3 = TextArea("E$Rn_-$ frequence (X10), ", textprops = dict(color = "brown", size = 16,rotation = 90,ha = 'left',va = 'bottom'))
			ybox = VPacker(children = [ybox1, ybox3],align = "bottom", pad = 0, sep = 8)
			anchored_ybox = AnchoredOffsetbox(loc = 8, child = ybox, pad = 0., frameon = False, bbox_to_anchor = (-0.08, 0.3), 
					bbox_transform = ax2.transAxes, borderpad = 0.)
			ax2.add_artist(anchored_ybox)
			indexes1 = np.arange(Interdistance_maxValue + 1)
			plt.xticks(indexes1 + width * 0.5 , indexes1)
			plt.text(2, 5, "D$Rn_+$ = DRn number in the bound set")
			plt.text(2, 4, "D$Rn_-$ = DRn number in the unbound set")
			plt.text(2, 3, "DRn frequence = DRn / $\sum_{i=0}^{i=20} DR_i + ER_i + IR_i$")
		else :
			ax1.axis([0, Interdistance_maxValue + 1, 0, 5.5])
			plt.plot(indexes1, map(divide, a, b), lw=2, label= "threshold : "+str(c))
			ax1.set_ylabel("$ERn_+$ frequence / $ERn_-$ frequence", size = 16)
		l = plt.axhline(y = 1)
		plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
	for a,b in zip(relative_IR,relative_IR_neg) :	
		ax1 = fig.add_subplot(1,3,3)
		ax1.set_xlabel("base pairs between "+factorTranscription+" inverted repeats", size = 16)
		if histo == True :
			indexes1 = range(Interdistance_maxValue + 1)
			ax1.axis([0, Interdistance_maxValue + 1, 0, 5.5])
			ax1.bar(indexes1, map(divide, a, b) , width , color = 'cornflowerblue')
			ax1.set_ylabel(' ', color = 'cornflowerblue', size = 16)
			ax2 = ax1.twinx()
			ax2.axis([0, Interdistance_maxValue + 1, 0, 5.5])
			ax2.bar(indexes1,[x * float(10) for x in a] , width , color = 'b')
			ax2.bar(indexes1,[x * float(10) for x in b] , width, color = 'brown', alpha = 0.85)
			ax2.set_ylabel('I$Rn_+$ frequence / D$Rn_-$ frequence', color = 'cornflowerblue', size = 16)
			ybox1 = TextArea("I$Rn_+$ frequence (X10)", textprops = dict(color = "b", size = 16,rotation = 90,ha = 'left',va = 'bottom'))
			ybox3 = TextArea("I$Rn_-$ frequence (X10), ", textprops = dict(color = "brown", size = 16,rotation = 90,ha = 'left',va = 'bottom'))
			ybox = VPacker(children = [ybox1, ybox3],align = "bottom", pad = 0, sep = 8)
			anchored_ybox = AnchoredOffsetbox(loc = 8, child = ybox, pad = 0., frameon = False, bbox_to_anchor = (-0.08, 0.3), 
					bbox_transform = ax2.transAxes, borderpad = 0.)
			ax2.add_artist(anchored_ybox)
			indexes1 = np.arange(Interdistance_maxValue + 1)
			plt.xticks(indexes1 + width * 0.5 , indexes1)
			plt.text(2, 5, "D$Rn_+$ = DRn number in the bound set")
			plt.text(2, 4, "D$Rn_-$ = DRn number in the unbound set")
			plt.text(2, 3, "DRn frequence = DRn / $\sum_{i=0}^{i=20} DR_i + ER_i + IR_i$")
		else :
			ax1.axis([0, Interdistance_maxValue + 1, 0, 5.5])
			plt.plot(indexes1, map(divide, a, b), lw=2)
			ax1.set_ylabel("$IRn_+$ frequence / $IRn_-$ frequence", size = 16)
		l = plt.axhline(y = 1)
	plt.show()

else :
	for a in InterDR :
		if factorTranscription == "LFY_matrix_19nucl" :
			ax1 = fig.add_subplot(1,1,1)
		else :	
			ax1 = fig.add_subplot(1,3,1)
		ax1.set_xlabel("base pairs between "+factorTranscription+" direct repeats", size = 16)
		if not negative_sets and histo == True :
			ax1.set_xlabel("base pairs between LFY 19 nucleotides matrix", size = 16)
			indexes1 = range(20 + 1)
			ax1.axis([0, 20 + 1, 0, 350])
			ax1.bar(indexes1, a , width , color = 'cornflowerblue')
			ax1.set_ylabel("occurences", color = 'cornflowerblue', size = 16)
			indexes1 = np.arange(20 + 1)
			plt.xticks(indexes1 + width * 0.5 , indexes1)
			plt.text(8, 230, sys.argv)
		else :
			plt.plot(indexes1, a, lw=2)
			ax1.axis([0, Interdistance_maxValue + 1, 0, 0.05])
			ax1.set_ylabel("$DRn_+$ frequence", size = 16)	
	if factorTranscription != "LFY_matrix_19nucl" :
		for a, c in zip(InterER,threshold)  :
			ax1 = fig.add_subplot(1,3,2)
			ax1.set_xlabel("base pairs between "+factorTranscription+" everted repeats", size = 16)
			if not negative_sets and histo == True :
				indexes1 = range(Interdistance_maxValue + 1)
				ax1.axis([0, Interdistance_maxValue + 1, 0, 350])
				ax1.bar(indexes1, a , width , color = 'cornflowerblue')
				indexes1 = np.arange(Interdistance_maxValue + 1)
				plt.xticks(indexes1 + width * 0.5 , indexes1)
			else :
				plt.plot(indexes1, a, lw=2, label= "threshold : "+str(c))
				ax1.axis([0, Interdistance_maxValue + 1, 0, 0.05])
				ax1.set_ylabel("$ERn_+$ frequence", size = 16)
				plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)	
	if factorTranscription != "LFY_matrix_19nucl" :
		for a in InterIR :	
			ax1 = fig.add_subplot(1,3,3)
			ax1.set_xlabel("base pairs between "+factorTranscription+" inverted repeats", size = 16)
			if not negative_sets and histo == True :
				indexes1 = range(Interdistance_maxValue + 1)
				ax1.axis([0, Interdistance_maxValue + 1, 0, 350])
				ax1.bar(indexes1, a , width , color = 'cornflowerblue')
				indexes1 = np.arange(Interdistance_maxValue + 1)
				plt.xticks(indexes1 + width * 0.5 , indexes1)
			else:
				plt.plot(indexes1, a, lw=2)
				ax1.axis([0, Interdistance_maxValue + 1, 0, 0.05])
				ax1.set_ylabel("$IRn_+$ frequence", size = 16)
	plt.show()

