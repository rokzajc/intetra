#!/usr/bin/env python3

#Author: Rok Zajc
print('COLIGO - Oligonucleotide frequency comparison\n---------------------------------------------\nAuthor: Rok Zajc\n\n')

import argparse, os, numpy as np, re, pandas as pd
from Bio import SeqIO,Phylo
from time import time
from programi_args import freq, models, najdi
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from matplotlib import pyplot as plt
import matplotlib

parser = argparse.ArgumentParser(description='The program calculates oligonucleotide (dinucleotide, trinucleotide...) frequency, which is based on chosen statistical scores(ZOM, Z-scores, ROF), of selected DNA sequnces. The frequncies are then compared using pearson R correlation.')
parser.add_argument('-i', help='folder where input fasta files are located (Default=<Current working directory>)', type=str, dest='input', default='' )
parser.add_argument('-o', help='Name of output CSV file', type=str, dest='output', default='Output')
parser.add_argument('-m', help='method of calculating frequencies, multiple methods can be chosen (Default=zscr)', dest='method', type=str, nargs='*', default=['zscr'],choices=['zscr', 'zom','rof'])
parser.add_argument('-n', help='lenght of nucleotide words whose occurrences in the sequence are counted(default=4), multiple can be chosen', type=int, dest='nuc_number', nargs='*', default=[4])
parser.add_argument('--join_multifasta', help='multifasta files will be joined in one sequence before analysed. ', dest='multifasta', action='store_true')
parser.add_argument('-t, --tree', help='construct a tree, select a method', dest='tree', type=str, nargs='*', default=[None],choices=[None, 'upgma','nj'])
args = parser.parse_args()

methods_dic={'zscr': models.z_score,'zom': models.zom_score, 'rof': lambda seznam,nuc: seznam[nuc-1] }

if max(args.nuc_number)>8:
    print('Error: -n must not be bigger than 8')
    quit()

t1=time()
print('Please wait...')

output_mapa=os.path.join(os.getcwd(),args.output)
try:
    os.mkdir(output_mapa)
except FileExistsError:
    pass

freq.naredi_komb(max(args.nuc_number)),freq.naredi_kompl(max(args.nuc_number))
frequncies=[]
names=[]
for file in najdi.najdi_datoteke(os.path.join(os.getcwd(),args.input),'*.f*'):
    name=os.path.basename(file).split('.f')[0]
    names.append(name)
    sekvenca=''
    dolzina=0
    counts=None
    with open(file) as f:
        for oznaka in SeqIO.parse(f, "fasta"):
            if args.multifasta==False:
                sekvenca_list=np.array(re.findall('.',str(oznaka.seq.upper())))
                new_counts=freq.count(sekvenca_list,max(args.nuc_number))
                dolzina+=len(sekvenca_list)
                try:
                    for i,vrednosti in enumerate(new_counts):
                        counts[i]+=vrednosti
                except:
                    counts=new_counts
            else:
                sekvenca+=oznaka.seq.upper()
    if args.multifasta==True:
        sekvenca_list=np.array(re.findall('.',str(sekvenca)))
        frequncies.append([count/len(sekvenca_list) for count in freq.count(sekvenca_list,max(args.nuc_number))])
    else:
        frequncies.append([count/dolzina for count in counts])
    print(f'Counted {najdi.slovar_stevil[max(args.nuc_number)]}nucleotide word occurances for {name}')
    #scores.append([[methods_dic[method](frequncies,nucl) for method in args.method]for nucl in args.nuc_number])
df_list=[]
df3=pd.DataFrame([[''],['']],index=['',''])
for nucl in args.nuc_number:
    for method in args.method:
        if method == 'zscr' and nucl<3:
            continue
        df1=pd.DataFrame(index=['Method:','Counted:',''], data=[[najdi.method_names[method]],[f'{najdi.slovar_stevil[nucl]}nucleotides'],['']])
        df2=pd.DataFrame(np.vstack((names,np.around(np.corrcoef(np.vstack([methods_dic[method](frequence,nucl) for frequence in frequncies])),3))),index=['']+names)
        df_list.append(pd.concat([df1, df2, df3]))
        print(f'Calculated correlations of oligonucleotide frequncies between sequences based on {najdi.slovar_stevil[nucl]}nucleotide {najdi.method_names[method]}.')
        for tree_type in args.tree:
            if tree_type==None:
                continue
            matrix=[]
            [matrix.append(list(vrstica[:i+1]))for i,vrstica in enumerate(1-np.array(df2)[1:].astype(np.float))]
            distMatrix=(DistanceMatrix(list(np.array(df2)[0]), matrix))
            if tree_type=='upgma':
                tree=DistanceTreeConstructor().upgma(distMatrix)
            if tree_type=='nj':
                tree=DistanceTreeConstructor().nj(distMatrix)
            fig = plt.figure(figsize=(13, 5), dpi=100) 
            matplotlib.rc('font', size=12)             
            matplotlib.rc('xtick', labelsize=10)       
            matplotlib.rc('ytick', labelsize=10)       
            #turtle_tree.ladderize()		   
            axes = fig.add_subplot(1, 1, 1)
            Phylo.draw(tree, axes=axes,do_show=False)
            fig.savefig(os.path.join(output_mapa,f'{nucl}{method}_{tree_type}_tree.png'))
            plt.close()

pd.concat(df_list).to_csv(os.path.join(output_mapa,'correlations.csv'),header=False,)

print(f'Time: {round((time()-t1),3)} sec')
print('Finished.')
