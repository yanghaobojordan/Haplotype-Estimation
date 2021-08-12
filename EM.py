# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 21:54:02 2019

@author: hyang
"""
from itertools import product
import scipy.stats
import random
import math
import sys
def main():
    f=open(sys.argv[1],"r")
    file=open(sys.argv[2],"w")
    genotype=[]
    for i in f:
        i=i.rstrip()
        genotype.append(i)
    #print ("genotype: @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    #print (genotype)
    g_frequency=dict()
    for i in genotype:
        if i not in g_frequency.keys():
            g_frequency[i]=1
        else:
            g_frequency[i]=g_frequency[i]+1
    #print("genotype frequency: @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    #print(g_frequency)
    total_g_frequency=0
    for t_g_f in g_frequency:
        total_g_frequency+=g_frequency[t_g_f]
    temp_haplotype=list(product(['0','1'],repeat=len(genotype[0])))
    haplotype=[]
    for i in temp_haplotype:
        haplotype.append(''.join(i))
    #print(haplotype)
    h_frequency=dict()
    for i in haplotype:
        temp=random.normalvariate(0.5,0.2401)
        h_frequency[i]=(1/(math.sqrt(2*math.pi*0.2401)))*math.exp((-(temp-0.5)**2)/(2*0.2401))
    #print(h_frequency)
    temp_h_frequency=dict()
    for H in h_frequency:
        temp_h_frequency[H]=h_frequency[H]
    for h in temp_h_frequency:
        total=0
        for g in g_frequency:
            temp=0
            for hh, gg in zip(h,g):
                if hh==gg:
                    temp+=1
                elif hh!=gg and gg=='2':
                    temp+=1
                else:
                    temp+=0
            if temp==len(g):
                total+=1
        if total==0:
            del h_frequency[h]
    #print("haplotype frequency: @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    #print(h_frequency)
    phase=[]
    for g in g_frequency:
        phase.append([g])
    cap=0
    for g in g_frequency:
        for h in h_frequency:
            if len(list(filter(lambda xy: xy[0] == xy[1] or (xy[0]=='2' and xy[0]!=xy[1]), zip(g, h))))==len(g):
                phase[cap].append(h)
                phase[cap].append(h_frequency[h])
        cap+=1
    #print("genotype and haplotype frequency: @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    #for i in phase:
        #print (i)
    #print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    for iteration in range(100):
        new_h_frequency=dict()
        for h in h_frequency:
            h_total=0
        #find each haplotype in every row of phase
            for p in phase:
                if h in p and len(p)>5:
                    h_total+=Get_Probability(h,p)*g_frequency[p[0]]
                elif h in p and len(p)==5:
                    h_total+=g_frequency[p[0]]
                elif h in p and len(p)==3:
                    h_total+=2*g_frequency[p[0]]
            h_total=h_total/(2*total_g_frequency)
            #print(h,h_total)
            new_h_frequency[h]=h_total
        h_frequency=dict()
        for nhf in new_h_frequency:
            h_frequency[nhf]=new_h_frequency[nhf]
        for hf in h_frequency:
            for a in range(len(phase)):
                for b in range(len(phase[a])):
                    if phase[a][b]==hf:
                        phase[a][b+1]=h_frequency[hf]
            
        #print ("new frequency: ")
        #print ("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    
    #Get_Probability('0000000000',['2020002002', '0000000000', 0.6693851933374666, '0000000001', 0.7742632285005067, '0000001000', 0.8077578552793877, '0000001001', 0.6341675881326346, '0010000000', 0.80327444926097, '0010000001', 0.8126004414126986, '0010001000', 0.7082952600977317, '0010001001', 0.7560804415516007, '1000000000', 0.8129378427411442, '1000000001', 0.7474184452957279, '1000001000', 0.7472848651927858, '1000001001', 0.8090754532520193, '1010000000', 0.785104981990187, '1010000001', 0.806682111474882, '1010001000', 0.7634002868996845, '1010001001', 0.7927296285429124])
    #print (h_frequency)
    for h in h_frequency:
        if h_frequency[h]>0:
            file.write(str(h))
            file.write('\n')
            print(h,h_frequency[h])
    f.close()
    file.close()
def Get_Probability(h,p):
    #print ("h is: ", h)
    #print ("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    #print ("p is: ", p)
    #print ("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    new_p=[]
    for i in p:
        if isinstance(i,str):
            new_p.append(i)
    genotype=p[0]
    haplotype=h
    for haplotype2 in new_p:
        f=True
        for a, b,c in zip(haplotype,haplotype2,genotype):
            if c=='2' and a==b:
                f=False
            elif c=='0'and (a!='0' or b!='0'):
                f=False
            elif c=='1' and (a!='1' or b!='1'):
                f=False
        if f:
            compliment=haplotype2
    #print(compliment)
    pair_list=[]
    temp_list=[]
    for i in range(1,len(new_p)):
        if new_p[i] not in temp_list:
            for j in new_p:
                f=True
                for a, b, c in zip(new_p[i],j,genotype):
                    if c=='2' and a==b:
                        f=False
                    elif c=='0'and (a!='0' or b!='0'):
                        f=False
                    elif c=='1' and (a!='1' or b!='1'):
                        f=False
                if f:
                    i_compliment=j
            pair_list.append([new_p[i],p[p.index(new_p[i])+1],i_compliment,p[p.index(i_compliment)+1]])
            temp_list.append(new_p[i])
            temp_list.append(i_compliment)
    #print(pair_list)
    total_frequency=0
    for pair in pair_list:
        total_frequency+=(pair[1]*pair[3])
    #print (p[p.index(haplotype)+1]*p[p.index(compliment)+1]/total_frequency)
    return (p[p.index(haplotype)+1]*p[p.index(compliment)+1]/total_frequency)
    #print ("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
main()
