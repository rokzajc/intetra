#!/usr/bin/env python3

#Author: Rok Zajc

import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import fnmatch
import re
from itertools import chain
seznam_fna=[]
for file in os.listdir(os.getcwd()):
    if fnmatch.fnmatch(file, '*.gbff'):
        seznam_fna.append(file)
for file_name in seznam_fna:
    with open(file_name) as file:
        records=[oznaka for oznaka in SeqIO.parse(file, "genbank")]
    features=list(chain.from_iterable([record.features for record in records]))
    rib_genes=[]
    for feature in features:
        try:
            if len(re.findall('ribosomal protein',feature.qualifiers['product'][0]))>0 and len(re.findall('transferase|protease',feature.qualifiers['product'][0]))==0:
                #print(feature.location)
                #print(feature.qualifiers['product'][0])
                rib_genes.append(feature)
        except:
            pass
    sequence=''
    for rec in records:
        sequence=sequence+rec.seq
    aaa=rib_genes[0].extract(sequence)
    with open(os.path.join(os.getcwd(),f'ribosomal_genes.fna'),'w') as file:
        records2=[SeqRecord(rib_feat.extract(sequence),id=rib_feat.qualifiers['locus_tag'][0],description=rib_feat.qualifiers['product'][0]) for rib_feat in rib_genes]
        SeqIO.write(records2, file, "fasta"), file.close()
    #with open(os.path.join(f'analiza_{file_name[:-3]}','5kb',f'window{najdi.stevilo_nicle(krog)}.fna'),'w') as file:
    #   SeqIO.write(window1, file, "fasta")
    #for rib in rib_genes:

print()