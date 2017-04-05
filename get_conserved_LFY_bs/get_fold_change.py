import ast
import csv
with open('list_of_the_LFY_regulated_genes_from_prediction_with_expression_data.txt', 'r') as f:
	mylist = ast.literal_eval(f.read())

dataBase = []	
with open('tab_sep_LFY_regulated_genes_from_database.csv', 'rb') as csvfile:
	dialect = csv.Sniffer().sniff(csvfile.read(1024), delimiters="\t")
	csvfile.seek(0)
	reader = csv.reader(csvfile, dialect)
	#print("len(reader) : ",len(reader))
	for line in reader:
		if line[4] == "FC_35sLFY_LER" :
			print(1)
		else :
			if (float(line[4]) < 1 and float(line[6]) < 0.05) or (float(line[4]) > 1 and float(line[8]) < 0.05) :
				dataBase.append(line[2])
				
print("len(dataBase) : ",len(dataBase))

print("len(mylist) : ",len(mylist))
ChIP_seq = []	
for i in mylist :
	ChIP_seq.append(i[2])
print("len(ChIP_seq) : ",len(ChIP_seq))
	
a = 0
for i in dataBase :
	for j in ChIP_seq :
		if j in i : 
			a = a + 1

print str(a*100 / len(dataBase))+" %"