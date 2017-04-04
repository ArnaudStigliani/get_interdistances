import ast

with open('list_of_the_LFY_binding_sites.txt', 'r') as f:
	all_sites = ast.literal_eval(f.read())

print("len(all_sites) : ",len(all_sites))
	
with open('list_of_the_conserved_LFY_binding_sites.txt', 'r') as f:
	conserved_sites = ast.literal_eval(f.read())
	
print("len(conserved_sites) : ",len(conserved_sites))
	
non_conserved_sites =  [x for x in all_sites if x not in conserved_sites]

print("len(non_conserved_sites) : ",len(non_conserved_sites))

text_file = open("list_of_the_non_conserved_LFY_regulated_genes.txt", "w")
text_file.write(str(non_conserved_sites))
text_file.close()