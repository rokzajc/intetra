import os
import numpy as np
from time import time
import matplotlib.pyplot as plt
import pandas as pd

from programi_args.models import z_score
try:
    from najdi import najdi_datoteke,stevilo_nicle
except ModuleNotFoundError:
    from programi_args.najdi import najdi_datoteke,stevilo_nicle
try:
    import models
except ModuleNotFoundError:
    import programi_args.models as models

def out(random_arg):
    exit()

func_dict={'fast':models.korelacije_fast,'all':models.vse_korelacije,'none':out}
def main(bakterija,output,method,coor_method,frame_len,slide_len,inc_sld,minlen,maxlen):
    okna=[f'SW{stevilo_nicle(i)}'for i in range(5000)]
    for method1 in models.method_lst[method]:
        t1=time()
        dolzine=np.array([])
        mediane=np.array([])
        for krog in range(1,1+int(maxlen/frame_len)):
            dolzina=frame_len*krog
            if dolzina<minlen:
                continue
            b_mapa=os.path.join(output,f'{dolzina}b')
            slide_total=slide_len+inc_sld*(krog-1)
            z_scores=np.load(najdi_datoteke(b_mapa,f'{method1}_win.npy')[0])
            pearson_r=func_dict[coor_method](z_scores).reshape(len(z_scores),len(z_scores))
            pearson_r2=models.korelacije_fast2(z_scores)[:len(z_scores),:len(z_scores)]
            razlike=np.around(pearson_r-pearson_r2,decimals=5)
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
                plt.title(f'Boxplot of {models.method_names[method1]} for sliding windows {dolzina}bp long for {bakterija}')
                plt.ylabel('Coorelation')
                plt.savefig(os.path.join(b_mapa,f'{method1}_boxplot.png'))
                plt.close()
            except:
                plt.close()
            pearson_r=np.vstack((okna[:len(z_scores)],np.around(pearson_r.reshape(len(z_scores),len(z_scores)),3)))
            vsebina=[[dolzina],[slide_total],[models.method_names[method1]],['']]+list(pearson_r)
            try:
                pd.DataFrame(vsebina,index=['Window lenghts','Slid by','Model','','']+okna[:len(pearson_r)-1]).to_csv(os.path.join(b_mapa,f'{method1}_coorelations.csv',),header=False)
            except:
                okna=[f'SW{stevilo_nicle(i)}'for i in range(len(pearson_r))]
                pd.DataFrame(vsebina,index=['Window lenghts','Slid by','Model','','']+okna[:len(pearson_r)-1]).to_csv(os.path.join(b_mapa,f'{method1}_coorelations.csv',),header=False)
        try:
            fig=plt.figure(num=1,figsize=(10,10))
            ax1=fig.add_subplot(111)
            #ax1.scatter(dolzine,mediane,c='red',label='Prekrivajoča okna')
            ax1.scatter(dolzine,mediane)
            plt.title(f'Medians of coorelation for {bakterija}. Method: {models.method_names[method1]}')
            #plt.legend()
            plt.xlabel('Window lenght')
            plt.ylabel('Median of coorelation')
            plt.savefig(os.path.join(output,f'{method1}_graf.png'))
            plt.close()
        except:
            plt.close()
        try:
            pd.DataFrame(np.rot90(np.vstack((dolzine,mediane))), columns = ['Window lenght','Median of coorelations']).sort_values(by=['Window lenght'],ignore_index=True).to_csv(os.path.join(output,f'{method1}_coorelation_medians.csv'),index=False)
        except:
            pd.DataFrame(np.rot90(np.vstack((dolzine,mediane))), columns = ['Window lenght','Median of coorelations']).sort_values(by=['Window lenght']).to_csv(os.path.join(output,f'{method1}_coorelation_medians.csv'),index=False)

    print(f'Time: {round((time()-t1),3)} sec')
    print(f'Aquired coorelations of {bakterija}')
if __name__ == "__main__":
    main()
