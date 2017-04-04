import csv
import ast

with open('list_of_the_hypothetical_LFY_regulated_genes_with_non_conserved_BS.txt', 'r') as f:
	mylist = ast.literal_eval(f.read())

print("len(mylist) : ",	len(mylist))
	
expression_data = []	
	
with open('tab_sep_LFY_regulated_genes_from_database.csv', 'rb') as csvfile:
	dialect = csv.Sniffer().sniff(csvfile.read(1024), delimiters="\t")
	csvfile.seek(0)
	reader = csv.reader(csvfile, dialect)
	
	for line in reader:
		for i in mylist :
			if i in line[2] :
				expression_data.append(line)

print("len(expression_data) : ",	len(expression_data))				
				
text_file = open("list_of_the_LFY_regulated_genes_from_prediction_with_expression_data_with_non_conserved_BS.txt", "w")
text_file.write(str(expression_data))
text_file.close()
