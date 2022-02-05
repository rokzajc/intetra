#!/usr/bin/env python3

#Author: Rok Zajc

import os
import numpy as np
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from time import time

try:
    import najdi
    import freq
except ModuleNotFoundError:
    import programi_args.najdi as najdi
    import programi_args.freq as freq


def main(input,output,frame_len,nucl,blockfasta):
    t1=time()
    freq.naredi_kompl(nucl)
    freq.naredi_komb(nucl)
    sekvenca=''
    with open(input) as file:
        for oznaka in SeqIO.parse(file, "fasta"):
            sekvenca+=oznaka.seq.upper()

    frames=[]
    sekvenca_list=np.array(re.findall('.',str(sekvenca)))
    len_seq=len(sekvenca_list)

    #Naredi vse potrebne mape
    try:
        os.mkdir(os.path.join(os.getcwd(),output))
        os.mkdir(os.path.join(os.getcwd(),output,'reads'))
    except FileExistsError:
        pass

    #Zanka, ki razdeli sekvenco na 5kb dolge odseke, jih shrani in prešteje frekvence tetranukleotidov
    frekv_list=[np.array([])for i in range(nucl)]
    mosticki_list=[np.array([])for i in range(nucl)]
    frekvence_genom=[np.zeros(4**(i+1))for i in range(nucl)]
    for krog in range(int(len_seq/frame_len)+1):
        mesto=krog*frame_len
        if len_seq-mesto>=frame_len:
            window=sekvenca[mesto:frame_len+mesto]
            if blockfasta!=True:
                window1=SeqRecord(window,id=f'SW{najdi.stevilo_nicle(krog)}',description=f'<{input}>locus={mesto}..{frame_len+mesto}')   
                with open(os.path.join(output,'reads',f'window{najdi.stevilo_nicle(krog)}.fna'),'w') as file:
                    SeqIO.write(window1, file, "fasta")
            okno=sekvenca_list[mesto:frame_len+mesto]
            frekvenca=freq.count(okno,nucl)
            mosticki=freq.count_bridge(sekvenca_list[frame_len+mesto-(nucl-1):frame_len+mesto+(nucl-1)],nucl)
        else:
            okno=sekvenca_list[mesto:]
            frekvenca=freq.count(okno,nucl)
            for a in range(nucl):
                frekvence_genom[a]+=frekvenca[a]
        if len_seq-mesto>=frame_len: 
            for a in range(nucl):
                frekv_list[a]=np.concatenate([frekv_list[a],frekvenca[a]],axis=0)
                frekvence_genom[a]+=frekvenca[a]+(mosticki[a-1] if a>0 else 0)
                if a>0:
                    mosticki_list[a-1]=np.concatenate([mosticki_list[a-1],mosticki[a-1]],axis=0)
    for a in range(nucl):
        frekv_list[a]=frekv_list[a].reshape(-1,4**(1+a))
        np.save(os.path.join(output,'reads',f'counts{a+1}.npy'), frekv_list[a])
        frekvence_genom[a]=frekvence_genom[a]/len_seq
        np.save(os.path.join(output,f'frequences_gen{a+1}.npy'), frekvence_genom[a])
        if a>0:
            mosticki_list[a-1]=mosticki_list[a-1].reshape(-1,4**(1+(a)))
            np.save(os.path.join(output,'reads',f'counts_most{a+1}.npy'), mosticki_list[a-1])

        
    print(f'Time: {round((time()-t1),3)} sec')
    print(f'Aquired {najdi.slovar_stevil[nucl]}nucleotide counts of {krog} fragments, of lenght {frame_len}bp for {input}.')

    
if __name__ == "__main__":
    main()