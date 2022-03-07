#!/usr/bin/env python3

#Author: Rok Zajc
print('COLIGO - Oligonucleotide frequency comparison\n----------------------------------------------\nAuthor: Rok Zajc\n\n')

import argparse, os, numpy as np, re, pandas as pd
from Bio import SeqIO
from time import time
from programi_args import freq, models, najdi

parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', help='folder where input fasta files are located', type=str, dest='input', default='' )
parser.add_argument('-o', help='Name of output file', type=str, dest='output', default='Output')
parser.add_argument('-m', help='method of calculating frequencies, multiple methods can be chosen (Default=zscr)', dest='method', type=str, nargs='*', default=['zscr'],choices=['zscr', 'zom','rof'])
parser.add_argument('-n', help='lenght of nucleotide words whose occurrences in the sequence are counted(default=4)', type=int, dest='nuc_number', nargs='*', default=[4])
args = parser.parse_args()

methods_dic={'zscr': models.z_score,'zom': models.zom_score, 'rof': lambda seznam,nuc: seznam[nuc-1] }

t1=time()
freq.naredi_komb(max(args.nuc_number)),freq.naredi_kompl(max(args.nuc_number))
frequncies=[]
names=[]
for file in najdi.najdi_datoteke(os.path.join(os.getcwd(),args.input),'*.f*'):
    names.append(os.path.basename(file))
    sekvenca=''
    with open(file) as f:
        for oznaka in SeqIO.parse(f, "fasta"):
            sekvenca+=oznaka.seq.upper()
    sekvenca_list=np.array(re.findall('.',str(sekvenca)))
    frequncies.append([count/len(sekvenca_list) for count in freq.count(sekvenca_list,max(args.nuc_number))])
    print(f'Aqured {najdi.slovar_stevil[max(args.nuc_number)]}nucleotide frequencies for {os.path.basename(file)}.')
    #scores.append([[methods_dic[method](frequncies,nucl) for method in args.method]for nucl in args.nuc_number])
df_list=[]
df3=pd.DataFrame([[''],['']],index=['',''])
for nucl in args.nuc_number:
    for method in args.method:
        df1=pd.DataFrame(index=['Method:','Counted:',''], data=[[najdi.method_names[method]],[f'{najdi.slovar_stevil[nucl]}nucleotides'],['']])
        df2=pd.DataFrame(np.vstack((names,np.around(np.corrcoef(np.vstack([methods_dic[method](frequence,nucl) for frequence in frequncies])),3))),index=['']+names)
        df_list.append(pd.concat([df1, df2, df3]))
pd.concat(df_list).to_csv(os.path.join(os.getcwd(),f'{args.output}.csv'),header=False,)

print(f'Time: {round((time()-t1),3)} sec')
print('Calculated correlations of oligonucleotide frequncies between sequences.')