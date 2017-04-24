#-*- coding: utf-8 -*-
import os

col_tab='p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}p{0cm}'
width=len(col_tab)/6 -1
#print(width)


with open('promoters_dofs.fasta','r') as f1:
    list_lines=list()
    for line in f1:
        line=line.strip("\n")
        list_lines.append(line)
        #print(list_lines[0])

seq=list_lines[3]

with open('Interdistances_ARF5_DOF34.csv','r') as f1:
    list_sites=list()
    for n,site in enumerate(f1):
        site=site.strip("\n")
        site=site.split("\t")
        if (n!=0):
            list_sites.append(site)
        
orientation=list()
i=0
while (i<len(seq)):
    orientation.append(0)
    i=i+1

orientation2 = list (orientation)
orientation3 = list (orientation)
orientation4 = list (orientation)


for elt in list_sites:
    i=0
    while(i < int(elt[1])+ 12 ):
        if elt[0] == 'DR' and int(elt[3]) >= -9 and int(elt[4]) >= -9 :
            orientation[int(elt[2])-1 + i] +=1
#            orientation[int(elt[2]) + int(elt[1]) + i -1] =  True
        if elt[0] == 'DR rev' and int(elt[3]) >= -9 and int(elt[4]) >= -9 :
            orientation4[int(elt[2])-1 + i]   +=1
#           orientation2[int(elt[2]) + int(elt[1]) + i -1] =  True
        if elt[0] == 'ER' and int(elt[3]) >= -9 and int(elt[4]) >= -9 :
            orientation2[int(elt[2])-1 + i]   +=1
#          orientation2[int(elt[2]) + int(elt[1]) + i -1] =  True
        if elt[0] == 'IR' and int(elt[3]) >= -9 and int(elt[4]) >= -9 :
            orientation3[int(elt[2])-1 + i]   +=1
#         orientation2[int(elt[2]) + int(elt[1]) + i -1] =  True
        i = i+1
        
couleur=list()

i=0
while (i<len(seq)):
    couleur.append(0)
    i=i+1

i=0
while (i < len(couleur)):
    if(orientation[i]):
        couleur[i] += orientation[i]
    if(orientation2[i]):
        couleur[i] += orientation2[i]
    if(orientation3[i]):
        couleur[i] += orientation3[i]
    if(orientation4[i]):
        couleur[i] += orientation4[i]
    i+=1
    



init_tex='\\documentclass[10pt]{article}\n\usepackage[utf8]{inputenc}\n\usepackage{longtable}\n\usepackage[table]{xcolor}\n\usepackage[margin=2cm]{geometry}\n\\begin{document}\n'
legend='\\section*{pDOF34}\n\n\\begin{tabular}{|c|p{2cm}|}\n\hline\nDR &  \\cellcolor{pink}\\\\\n\\hline\nDR rev &  \\cellcolor{purple!30}\\\\\n\\hline\\\nER &  \\cellcolor{blue!25}\\\\\n\\hline\nIR &  \\cellcolor{green!30}\\\\\n\\hline\nSuperimposed &  \\cellcolor{red}\\\\\n\\hline\end{tabular}'
init_tab='\\begin{longtable}{'+col_tab+'}\n'
end_tab='\n\\end{longtable}\n'
end_tex='\\end{document}\n'
#print(init_tab)




i=0
j=0
texte = str()
texte2 = str()
while( i < len(seq)):
    if couleur[i] > 1 :
        texte = texte + '\\cellcolor{red}' 
    elif (orientation[i] != 0 ):
        texte = texte + '\\cellcolor{pink}' 
    elif (orientation2[i] != 0 ):
        texte = texte + '\\cellcolor{blue!25}' 
    elif (orientation3[i] != 0 ):
        texte = texte + '\\cellcolor{green!25}' 
    elif (orientation4[i] != 0 ):
        texte = texte + '\\cellcolor{purple!30}' 
    texte = texte + seq[i]
    texte2 = texte2 + '.' 
    i = i+1
    j = j+1
    if(j == width or i == (len(seq) -1)):
        j=0
        texte = texte + ' & \\\\\n'
        #texte = texte + texte2 + ' & \\\\\n'
        texte2 = str()
    else :
        texte = texte + ' & '
        texte2 = texte2 + ' & '


# i=len(list_sites)-1
# while(i>=0):
#     texte_f = texte[int(list_sites[i][2]):len(seq)-1]
#     texte = texte[1:int(list_sites[i][2])]
#     texte = texte  + ' & ' + list_sites[i][0][0]+' & ' + list_sites[i][0][1]
#     if(len(list_sites[i][0] > 2)):
#         texte = texte  + ' & r &'
#     texte = texte + ' & ' + list_sites[i][1][0]
#     if(len(list_sites[i][1] > 2)):
#         texte = texte  + ' & ' + list_sites[i][1][1]
#     texte = texte + ' & ( & ' + list_sites[i][3][0] + ' & ' + list_sites[i][3][1]
#     if(list_sites[i][3][2] != '.' ):
#         texte = texte  + ' & ' + list_sites[i][3][2]
#     texte = texte + ' & : &' + list_sites[i][4][0] + ' & ' + list_sites[i][4][1]
#     if(list_sites[i][4][2] != '.' ):
#         texte = texte  + ' & ' + list_sites[i][4][2]
#     texte = texte + ' & '
    

    

with open('pDOF34.tex','w') as f1:
    f1.write(init_tex)
    f1.write(legend)
    f1.write(init_tab)
    f1.write(texte)
    f1.write(end_tab)
    f1.write(end_tex)

os.system('pdflatex pDOF34.tex')
