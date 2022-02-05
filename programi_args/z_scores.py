#!/usr/bin/env python3

#Author: Rok Zajc

import os
import numpy as np
from time import time

try:
    import najdi
except ModuleNotFoundError:
    import programi_args.najdi as najdi
try:
    import models
except ModuleNotFoundError:
    import programi_args.models as models

methods_dic={'zscr': models.z_score,'som': models.som_score,'fom': models.fom_score,'zom': models.zom_score, 'rof': lambda seznam,nuc: seznam[nuc-1] }
def main(output,method,nuc_number):
    t1=time()
    for nucl in nuc_number:
        freq_gen=[]
        for f in najdi.najdi_datoteke(output,'frequences_gen*.npy')[:nucl+1]:
            freq_gen.append(np.load(f))
        for method1 in method:
            obseg_podatkov=0
            z_score_gen=methods_dic[method1](freq_gen,nucl)
            obseg_podatkov+=z_score_gen.size
            np.save(os.path.join(output,f'{method1}_gen{nucl}.npy'), z_score_gen)
            for dir in najdi.najdi_datoteke(output,'*b'):
                freq=[]
                for f in najdi.najdi_datoteke(dir,'frequences*.npy')[:nucl+1]:
                    freq.append(np.load(f))
                z_score_win=methods_dic[method1](freq,nucl)
                obseg_podatkov+=z_score_win.size
                np.save(os.path.join(dir,f'{method1}_win{nucl}.npy'), z_score_win)
            print(f'Time: {round((time()-t1),3)} sec')
            print(f'Calculated {najdi.method_names[method1]} for {obseg_podatkov} {najdi.slovar_stevil[nucl]}nucleotide frequencies of diffrent sequences.')

if __name__ == "__main__":
    main()