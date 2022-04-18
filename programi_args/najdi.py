#!/usr/bin/env python3

#Author: Rok Zajc

import os
import fnmatch

def najdi_datoteke(pot,beseda):
    seznam=[]
    for dir in sorted(os.listdir(pot)):
        if fnmatch.fnmatch(dir,beseda):
            seznam.append(os.path.join(pot,dir))
    return seznam

def stevilo_nicle(num):
    dol=len(str(num))
    nicle=''
    for a in range(6-dol):
        nicle+='0'
    return nicle+str(num)

slovar_stevil={1:'mono',2:'di',3:'tri',4:'tetra',5:'penta',6:'heksa', 7:'hepta',8:'octa'}
method_names={'zscr':'z-score',
'fom':'first order Markov chain frequencies',
'som':'second order Markov chain frequencies',
'zom':"zero'th order Markov chain frequencies",
'rof': 'relative oligonucleotide frequencies'}
