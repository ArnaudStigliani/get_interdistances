import re

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

###############################################################################

def get_score_matrix(Mdata,matrixType) :
	if matrixType == "freq" :
		## These lines allows to transform the frequency values into scores values
		matScore = []
		lenMotif = 0
		a = 0
		for i in range(0,len(Mdata)/4):
			lenMotif = lenMotif + 1
			if i+1 in dependency_positions1 or dependency_positions2 or dependency_positions3:
				print("i+1 : ",i+1)
				for j in range (0,4):
					matScore.append(0.0)
			else :
				fmax = float(max(Mdata[a],Mdata[a+1],Mdata[a+2],Mdata[a+3])) + pseudoCount
				for j in range (0,4):
					matScore.append(np.log(float(float(Mdata[a+j]) + pseudoCount) /fmax))
			a = a + 4
	else :
		matScore = map(float,Mdata)
		lenMotif = len(matScore)/4
	return (matScore, lenMotif)

def get_dependency_matrix(dependencyFile,num) : 
	G = open(dependencyFile,"r")
	dependency_file_content = G.read().replace("\r","\n") + "\n"
	G.close()
	
	num2 = re.compile(r"(?<![\d.])[0-9]+(?![\d.])")
	position_dependency_matrix = num2.findall(dependency_file_content)
	position_dependency_matrix = map(int, position_dependency_matrix)
	
	dependency_matrix = num.findall(dependency_file_content)
	dependency_matrix = map(float, dependency_matrix)
	
	dependency_matrix_associated_with_position = []
	index1 = 0
	index2 = 3
	index3 = 0
	index4 = 64
	for i in range(0, 3):
		dependency_matrix_associated_with_position.append(position_dependency_matrix[index1:index2])
		dependency_matrix_associated_with_position.append(dependency_matrix[index3:index4])
		index1 = index1 + 3
		index2 = index2 + 3
		index3 = index3 + 64
		index4 = index4 + 64
		
	return(dependency_matrix_associated_with_position)
	
def divide(a, b):
    if b == 0:
        return np.nan
    else: 
        return a/b

def seq_c(site):
        site_i = site[-1::-1]
        site_c = ""
        for x in site_i:
		y = codigoi[x]
                site_c = site_c + y
        return site_c   

def add_scores_associated_with_interdependent_positions(dependency_matrix,scoreStrandPos,scoreStrandNeg,strandPos):
	cStrand = ""
	for lettre in seq_c(strandPos):
		cStrand = lettre + cStrand
	cStrand = cStrand[::-1]
	
	site1 = ""
	Csite1 = ""
	site2 = ""
	Csite2 = ""
	site3 = ""
	Csite3 = ""
	for i in dependency_matrix[0]:
		site1 = site1 + strandPos[i-1]
		Csite1 = Csite1 + cStrand[i-1]
	for i in dependency_matrix[2]:
		site2 = site2 + strandPos[i-1]
		Csite2 = Csite2 + cStrand[i-1]
	for i in dependency_matrix[4]:
		site3 = site3 + strandPos[i-1]
		Csite3 = Csite3 + cStrand[i-1]
	#print("dependency_matrix[1][0] : ",dependency_matrix[1][0])
	scoreStrandPos = scoreStrandPos + dependency_matrix[1][codigo[site1]] + dependency_matrix[3][codigo[site2]] + dependency_matrix[5][codigo[site3]]
	scoreStrandNeg = scoreStrandNeg + dependency_matrix[1][codigo[Csite1]] + dependency_matrix[3][codigo[Csite2]] + dependency_matrix[5][codigo[Csite3]]
	return(scoreStrandPos, scoreStrandNeg)
