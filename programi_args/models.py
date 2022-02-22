#!/usr/bin/env python3

#Author: Rok Zajc

import numpy as np
import os
import sys

def z_score(freq,nucl):
    for a in range(nucl):
        if freq[a].ndim==1:
            freq[a]=freq[a].reshape(1,4**(a+1))
    Nxyzw=freq[nucl-1]
    st_vrstic=int(len(Nxyzw))
    Nxyz=np.repeat(freq[nucl-2],4).reshape(st_vrstic,4**(nucl))
    Nyzw=np.tile(freq[nucl-2],4).reshape(st_vrstic,4**(nucl))
    Nyz=np.repeat(np.tile(freq[nucl-3],4),4).reshape(st_vrstic,4**(nucl))
    ex_value=(Nxyz*Nyzw)/Nyz
    var=ex_value*(Nyz-Nxyz)*(Nyz-Nyzw)/Nyz**2
    return np.nan_to_num((Nxyzw-ex_value)/(var**0.5))

def som_score(freq):
    for a in range(4):
        if freq[a].ndim==1:
            freq[a]=freq[a].reshape(1,4**(a+1))
    Nxyzw=freq[3]
    st_vrstic=int(len(Nxyzw))
    Nxyz=np.repeat(freq[2],4).reshape(st_vrstic,256)
    Nyzw=np.tile(freq[2],4).reshape(st_vrstic,256)
    Nyz=np.repeat(np.tile(freq[1],4),4).reshape(st_vrstic,256)
    return(Nxyzw*Nyz)/(Nxyz*Nyzw)

def fom_score(freq):
    for a in range(4):
        if freq[a].ndim==1:
            freq[a]=freq[a].reshape(1,4**(a+1))
    Nxyzw=freq[3]
    st_vrstic=int(len(Nxyzw))
    Ny=np.repeat(np.repeat(np.tile(freq[0],4),4),4).reshape(st_vrstic,256)
    Nz=np.repeat(np.tile(np.tile(freq[0],4),4),4).reshape(st_vrstic,256)
    Nyz=np.repeat(np.tile(freq[1],4),4).reshape(st_vrstic,256)
    Nxy=np.repeat(np.repeat(freq[1],4),4).reshape(st_vrstic,256)
    Nzw=np.tile(np.tile(freq[1],4),4).reshape(st_vrstic,256)
    return(Nxyzw*Ny*Nz)/(Nyz*Nzw*Nxy)

def zom_score(freq,nucl):
    if freq[0].ndim==1:
        freq[0]=freq[0].reshape(1,4)
    if freq[nucl-1].ndim==1:
        freq[nucl-1]=freq[nucl-1].reshape(1,4**nucl)
    Nxyzw=freq[nucl-1]
    try:
        Nukleotidi=np.prod(freq[0][:,mesta[-1]].astype(float),axis=2)
        return Nxyzw/Nukleotidi
    except:
        naredi_mesta(nucl)
        Nukleotidi=np.prod(freq[0][:,mesta[-1]].astype(float),axis=2)
        return Nxyzw/Nukleotidi

baze=np.arange(4)
def komb(velikost,arr1=np.array([]),arr=np.array([])):
    for baza in baze:
        if velikost>1:
            arr2=komb(velikost-1,np.concatenate([arr1,np.array([baza])]))
            arr=np.concatenate((arr,arr2))
        else:
            arr2=np.array([int(baza)])
            arr=np.concatenate((arr,arr1,arr2))
    return arr
mesta=[]
def naredi_mesta(nukl):
    try:
        mesta.append(np.load(os.path.join(sys.path[0],'programi_args','data',f'mesta{nukl}')))
    except FileNotFoundError:
        try:
            os.mkdir(os.path.join(sys.path[0],'programi_args','data'))
        except FileExistsError:
            pass
        mesta.append(komb(nukl).astype(int).reshape(-1,nukl))
        np.save(os.path.join(sys.path[0],'programi_args','data',f'mesta{nukl}'),mesta[-1])
#mesta=komb(4).astype(int).reshape(256,4)

def korelacije(z_score1,z_score2):
    n=256
    z_score2=np.tile(z_score2.reshape(n),len(z_score1)).reshape(len(z_score1),n)
    sum_xy=np.sum((z_score1*z_score2),axis=1)
    sum_x=np.sum(z_score1,axis=1)
    sum_y=np.sum(z_score2,axis=1)
    sum_x2=np.sum(z_score1**2,axis=1)
    sum_y2=np.sum(z_score2**2,axis=1)
    return (n*sum_xy-sum_x*sum_y)/(((n*sum_x2-sum_x**2)**0.5)*((n*sum_y2-sum_y**2)**0.5))

#korelacija med eno in tabelo drugih
def korelacije2(z_score1,z_score2,nucl):
    n=4**nucl
    z_score2=np.tile(z_score2.reshape(n),len(z_score1)).reshape(len(z_score1),n)
    x_avg=np.tile(np.average(z_score1,axis=1),n).reshape(len(z_score1),n)
    y_avg=np.tile(np.average(z_score2,axis=1),n).reshape(len(z_score1),n)
    imen=np.sum((z_score1-x_avg)*(z_score2-y_avg),axis=1)
    stev=(np.sum((z_score1-x_avg)**2,axis=1)*np.sum((z_score2-y_avg)**2,axis=1))**0.5
    return imen/stev
