#!/usr/bin/env python3

#Author: Rok Zajc

try:
    import models
except ModuleNotFoundError:
    import programi_args.models as models
try:
    import najdi
except ModuleNotFoundError:
    import programi_args.najdi as najdi
import pandas as pd
from time import time
import os
import numpy as np
import matplotlib.pyplot as plt
def main(input,output,frame_len,slide_len,inc_sld,maxlen,minlen,method,nuc_number):
    t1=time()
    for krog in range(1,1+int(maxlen/frame_len)):
        dolzina=frame_len*krog
        if dolzina<minlen:
            continue
        mapa=os.path.join(output,f'{dolzina}b')
        slide_total=slide_len+inc_sld*(krog-1)
        fig=plt.figure(num=1,figsize=(12,5))
        ax1=fig.add_subplot(111)
        for nucl in nuc_number:
            for method1 in method:
                genom_score=np.load(os.path.join(output,f'{method1}_gen{nucl}.npy'))
                win_scores=np.load(os.path.join(mapa,f'{method1}_win{nucl}.npy'))
                coor=models.korelacije2(win_scores,genom_score,nucl)
                loc=np.arange(0,len(coor)*int(slide_total),int(slide_total))
                ax1.plot(loc,coor,linewidth=0.5,label=str(nucl)+method1)
                vsebina=[[dolzina],[slide_total],[najdi.method_names[method1]],[najdi.slovar_stevil[nucl]],[''],[f'Correlation to genome']]+list(coor.reshape(len(coor),1))
                pd.DataFrame(vsebina,index=['Window lenghts','Slid by','Model','Frequencies calculated from','','locus']+list(loc)).to_csv(os.path.join(mapa,f'{method1}_autocorrelation{nucl}.csv',),header=False)
        plt.title(f'Autocorrelation {input}')
        plt.ylabel('Correlation')
        plt.xlabel('Locus')
        plt.ylim(0,1)
        plt.legend()
        plt.savefig(os.path.join(mapa,'autocorrelation_plot.png'))
        #plt.show()
        plt.close()
    print(f'Time: {round((time()-t1),3)} sec')
    print(f'Calculated autocorrelation of {input}.')

if __name__ == "__main__":
    main()