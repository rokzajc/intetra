#!/usr/bin/env python3

#Author: Rok Zajc

import numpy as np
import sys
import os
nuc_as_num={'A': 0,'C': 1, 'G': 2, 'T': 3}
def nuc_to_index(zaporedje):
    try:
        stevilka=[nuc_as_num[n] for n in zaporedje]
    except:
        return None
    return np.sum([stevilo*(4**i) for i,stevilo in enumerate(stevilka[::-1])])

#Seznam mest, kjer se nahajajo komplementarni n-nukleotidi
def komplement(velikost):
    arr=np.arange(4**velikost-1,-1,-(4**(velikost-1)))
    for krog in range(velikost-1):
        for a in range(3):
            arr=np.concatenate((arr,arr[-(4**(krog+1)):]-(4**(velikost-2-krog))))
    return arr

kompl=[]
def naredi_kompl(nukl):
    try:
        for i in range(nukl):
            kompl.append(np.load(os.path.join(sys.path[0],'programi_args','data',f'kompl{i}')))
    except FileNotFoundError:
        try:
            os.mkdir(os.path.join(sys.path[0],'programi_args','data'))
        except FileExistsError:
            pass
        for i in range(nukl):
            k=komplement(i+1)
            kompl.append(k)
            np.save(os.path.join(sys.path[0],'programi_args','data',f'kompl{i}'),k)

#Seznam vseh možnih tetranukleotidov
def komb(velikost,arr1=np.array([]),arr=np.array([]),baze=['A','C','G','T']):
    for baza in baze:
        if velikost>1:
            arr2=komb(velikost-1,np.concatenate([arr1,np.array([baza])]),baze=baze)
            arr=np.concatenate((arr,arr2))
        else:
            arr2=np.array([baza])
            arr=np.concatenate((arr,arr1,arr2))
    return arr
vse_komb=[]

def naredi_komb(nukl):
    try:
        for i in range(nukl):
            vse_komb.append(np.load(os.path.join(sys.path[0],'programi_args','data',f'vse_komb{i}')))
    except FileNotFoundError:
        try:
            os.mkdir(os.path.join(sys.path[0],'programi_args','data'))
        except FileExistsError:
            pass
        for i in range(nukl):
            k=komb(i+1).reshape(-1,i+1)
            vse_komb.append(k)
            np.save(os.path.join(sys.path[0],'programi_args','data',f'vse_komb{i}'),k)

def count(sekvenca,nukl=4):
    #sekvenca=np.concatenate((sekvenca,np.repeat('',nukl-1)))
    lenght=len(sekvenca)
    list_arr=[]
    naslednji=None
    ostanek=lenght%nukl
    for a in range(nukl):
        list_arr.append(np.array(sekvenca[a:lenght-ostanek-(nukl if a!=0 else 0)+a]).reshape(-1,nukl))
    for a in range(1,ostanek+1):
        list_arr.append(np.array(sekvenca[lenght-ostanek-nukl+a:lenght-ostanek+a]).reshape(1,nukl))
    if lenght>((4**nukl)*5):
        besede_arr,frekv=np.unique(np.concatenate(list_arr),return_counts=True,axis=0)
        zadnji=[]
        for a in range(nukl-1,0,-1):
            zadnji.append(np.zeros(4**a).astype(int))
            for b in range(nukl-a):
                zadnji[-1][nuc_to_index(sekvenca[-a-b:lenght-b])]+=1
        # Clean up counts of all 'N'
        if np.array_equal(vse_komb[nukl-1],besede_arr)==False:
            zaporedja,lokacije,counts=np.unique(np.concatenate((besede_arr,vse_komb[nukl-1])),return_index=True,return_counts=True,axis=0)
            zamik=0
            len_frekv=len(frekv)
            naslednji=np.array([])
            naslednji_freq=np.array([])
            for lokacija in np.argwhere(counts==1):
                if lokacije[lokacija]>len_frekv-1:
                    frekv=np.insert(frekv,lokacija-zamik,0)
                    zaporedja=np.insert(zaporedja,lokacija-zamik,np.repeat([''],nukl),0)
                else:
                    try:
                        naslednji=np.concatenate((naslednji,zaporedja[lokacija-zamik])).reshape(-1,4)
                    except:
                        naslednji=zaporedja[lokacija-zamik]
                    naslednji_freq=np.concatenate((naslednji_freq,frekv[lokacija-zamik]))
                    zaporedja=np.delete(zaporedja,lokacija-zamik,0)
                    frekv=np.delete(frekv,lokacija-zamik)
                    zamik+=1
    else:
        besede_arr,frekv=np.unique(np.concatenate((vse_komb[nukl-1],np.concatenate(list_arr))),return_counts=True,axis=0)
        frekv-=1
        zadnji=[]
        for a in range(nukl-1,0,-1):
            zadnji.append(np.zeros(4**a).astype(int))
            for b in range(nukl-a):
                zadnji[-1][nuc_to_index(sekvenca[-a-b:lenght-b])]+=1
        if len(frekv)!=4**nukl:
            zaporedja,counts=np.unique(np.concatenate((besede_arr,vse_komb[nukl-1])),return_counts=True,axis=0)
            zamik=0
            naslednji=np.array([])
            naslednji_freq=np.array([])
            for lokacija in np.argwhere(counts==1):
                try:
                    naslednji=np.concatenate((naslednji,zaporedja[lokacija-zamik])).reshape(-1,4)
                except:
                    naslednji=zaporedja[lokacija-zamik]
                naslednji_freq=np.concatenate((naslednji_freq,frekv[lokacija-zamik]))
                zaporedja=np.delete(zaporedja,lokacija-zamik,0)
                frekv=np.delete(frekv,lokacija-zamik)
                zamik+=1
    vse_frekv=[frekv]
    for i in range(nukl-1):
        vse_frekv.append(np.sum(vse_frekv[-1].reshape(-1,4),1))
    for i in range(nukl-1):
        vse_frekv[i+1]+=zadnji[i]
    vse_frekv.reverse()
    if type(naslednji)==list:
        for i in range(1,nukl):
            for i2,kandidat in enumerate(naslednji[:,:i]):
                mesto=nuc_to_index(kandidat)
                if mesto == None:
                    continue
                vse_frekv[i-1][mesto]+=naslednji_freq[i2]
    return [arr + arr[kompl[i]]for i, arr in enumerate(vse_frekv)]

def count_bridge(sekvenca,mode=4):
    frekv_list=[]
    for nukl in range(2,mode+1):
        list_arr=[]
        for a in range(nukl-1):
            list_arr.append(np.array(sekvenca[mode+a-nukl:mode+a]))
        frekv=np.zeros(4**nukl)
        for arr in list_arr:
            mesto=nuc_to_index(arr)
            if mesto!=None:
                frekv[mesto]+=1
        frekv_list.append(frekv+frekv[kompl[nukl-1]])
    return frekv_list
