#!/usr/bin/env python3

#Author: Rok Zajc

import os
import numpy as np
from time import time
import matplotlib.pyplot as plt
import pandas as pd

try:
    import najdi
except ModuleNotFoundError:
    import programi_args.najdi as najdi
try:
    import models
except ModuleNotFoundError:
    import programi_args.models as models


def main(bakterija,output,method,nuc_number,frame_len,slide_len,inc_sld,minlen,maxlen):
    okna=[f'SW{najdi.stevilo_nicle(i)}'for i in range(5000)]
    for nucl in nuc_number:
        for method1 in method:
            t1=time()
            dolzine=np.array([])
            mediane=np.array([])
            for krog in range(1,1+int(maxlen/frame_len)):
                dolzina=frame_len*krog
                if dolzina<minlen:
                    continue
                b_mapa=os.path.join(output,f'{dolzina}b')
                slide_total=slide_len+inc_sld*(krog-1)
                z_scores=np.load(os.path.join(b_mapa,f'{method1}_win{nucl}.npy'))
                pearson_r=np.corrcoef(z_scores).flatten()
                dolzine=np.append(dolzine,dolzina)
                try:
                    mediane=np.append(mediane,np.median(pearson_r))
                except:
                    mediane=np.append(mediane,np.nan)
                fig=plt.figure(num=1,figsize=(10,10))
                ax1=fig.add_subplot(111)
                dict={'Sliding window':pearson_r}
                try:
                    ax1.boxplot(dict.values())
                    ax1.set_xticklabels(dict.keys())
                    plt.title(f'Boxplot of {najdi.method_names[method1]} for sliding windows {dolzina}bp long for {bakterija}')
                    plt.ylabel('Correlation')
                    plt.savefig(os.path.join(b_mapa,f'{method1}_boxplot{nucl}.png'))
                    plt.close()
                except:
                    plt.close()
                if len(pearson_r)>len(okna):
                    okna=[f'SW{najdi.stevilo_nicle(i)}'for i in range(len(pearson_r))]
                pearson_r=np.around(pearson_r.reshape(len(z_scores),len(z_scores)),3)
                    #pearson_r=np.vstack((okna[:len(z_scores)],np.around(pearson_r.reshape(len(z_scores),len(z_scores)),3)))
                vsebina=[[dolzina],[slide_total],[najdi.method_names[method1]],[f'{najdi.slovar_stevil[nucl]}nucleotides'],[''],okna[:len(z_scores)]]+list(pearson_r)
                pd.DataFrame(vsebina,index=['Window lenghts','Slid by','Model','Frequencies calculated from','','']+okna[:len(pearson_r)]).to_csv(os.path.join(b_mapa,f'{method1}_correlations{nucl}.csv',),header=False)
            try:
                fig=plt.figure(num=1,figsize=(10,10))
                ax1=fig.add_subplot(111)
                #ax1.scatter(dolzine,mediane,c='red',label='Prekrivajoča okna')
                ax1.scatter(dolzine,mediane)
                plt.title(f'Medians of correlation for {bakterija}. Method: {najdi.method_names[method1]}')
                #plt.legend()
                plt.xlabel('Window lenght')
                plt.ylabel('Median of correlation')
                plt.savefig(os.path.join(output,f'{method1}_graf{nucl}.png'))
                plt.close()
            except:
                plt.close()
            try:
                pd.DataFrame(np.rot90(np.vstack((dolzine,mediane))), columns = ['Window lenght','Median of correlations']).sort_values(by=['Window lenght'],ignore_index=True).to_csv(os.path.join(output,f'{method1}_correlation{nucl}_medians.csv'),index=False)
            except:
                pd.DataFrame(np.rot90(np.vstack((dolzine,mediane))), columns = ['Window lenght','Median of correlations']).sort_values(by=['Window lenght']).to_csv(os.path.join(output,f'{method1}_correlation{nucl}_medians.csv'),index=False)

    print(f'Time: {round((time()-t1),3)} sec')
    print(f'Aquired correlations of {bakterija}')
if __name__ == "__main__":
    main()
