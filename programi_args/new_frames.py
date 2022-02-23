#!/usr/bin/env python3

#Author: Rok Zajc

import os
import numpy as np
from time import time
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re
try:
    from najdi import najdi_datoteke
    from najdi import stevilo_nicle
except ModuleNotFoundError:
    from programi_args.najdi import najdi_datoteke
    from programi_args.najdi import stevilo_nicle

def main(bakterija,output,frame_len,read_by,slide_len,increase_slide,maxsplit,minlen,nucl,blockfasta):
    t1=time()
    #bakterija=re.split('Output_|[.]fna',analiza)[1]
    slide=najdi_datoteke(output,'reads')[0]
    fasta_files=najdi_datoteke(slide,'*.fna')
    sekvence=[str(record.seq) for record in [list(SeqIO.parse(file, "fasta"))[0] for file in fasta_files]]
    sekvence_array=np.array(sekvence)
    frekvence_slide=[]
    mosticki=[]
    #frekvence=[]
    for a in range(1,nucl+1):
        frekvence_slide.append(np.load(os.path.join(slide,f'counts{a}.npy')))
        if a>1:
            mosticki.append(np.load(os.path.join(slide,f'counts_most{a}.npy')))
    for krog in range(int(maxsplit/frame_len)):
        win_len=(krog+1)*frame_len
        if win_len<minlen:
            continue
        slide_total=slide_len+increase_slide*(krog)
        try:
            os.mkdir(os.path.join(output,f'{win_len}b'))
        except FileExistsError:
            pass
        mapa=os.path.join(output,f'{win_len}b')
        stevilo_oken1=len(sekvence_array)if blockfasta!=True and len(sekvence_array)<=len(frekvence_slide[0])else len(frekvence_slide[0])
        stevilo_oken=stevilo_oken1-stevilo_oken1%(win_len/read_by)
        mesta=np.arange(stevilo_oken).reshape(int((stevilo_oken*read_by)/(win_len)),int(win_len/read_by))
        if frame_len!=slide_len:
            seznam=[np.arange(i*int(slide_total/read_by),stevilo_oken-int(win_len/read_by)+i*int(slide_total/read_by)if i!=0 else stevilo_oken) for i in range(int(win_len/slide_total))]
            mesta=np.unique(np.reshape(np.concatenate(seznam),(-1,int(win_len/read_by))),axis=0)
        mesta=mesta.astype(int)        
        if blockfasta!=True:
            try:
                os.mkdir(os.path.join(mapa,'windows'))
            except FileExistsError:
                pass
            nove_sek=sekvence_array[mesta]
            for frame in range(len(nove_sek)):
                window=SeqRecord(Seq(''.join(nove_sek[frame].tolist())),id=f'SW{stevilo_nicle(frame)}',description=f'<{bakterija}>locus={frame*slide_total}..{win_len+frame*slide_total}')
                with open(os.path.join(mapa,'windows',f'window{stevilo_nicle(frame)}.fna'),'w') as file:
                    SeqIO.write(window, file, "fasta")
        for a in range(nucl):
            if a>0:
                frekvence=(np.sum(frekvence_slide[a][mesta],axis=1)+np.sum(mosticki[a-1][mesta[:,:-1]],axis=1))/win_len
            else:
                frekvence=np.sum(frekvence_slide[a][mesta],axis=1)/win_len
            np.save(os.path.join(mapa,f'frequences{a+1}.npy'), frekvence)
        print(f'Aquired sliding windows of lenght {win_len}bp, slid by {slide_total}bp.')
    print(f'Time: {round((time()-t1),3)}sec')
    print(f'Aquired {krog} new lengths of sliding windows.')
if __name__ == "__main__":
    main()
