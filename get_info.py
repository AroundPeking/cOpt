#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##
import re

# find "wordcheck" from "fildsource"
# return the whole desired line and number of advent, e.g. {'etot(Ha): -30.679626963461189\n': 1} 
def checklog(filesource, wordcheck):
    size = {}
    try:
        with open(filesource, "r") as file:
            for line in file:
                x = re.search(wordcheck, line)
                if x:
                    tmp = x.group()
                    size[line] = size.get(line, 0) + 1
    
    except FileNotFoundError as e:
        print(e)
    
    return size
    
#read (element,Nu,Nl,Rcut,Ecut) from STRU, read (Ne) from ORBITAL_RESULTS.txt
def get_info_element():
    count=0
    try:
        file=open('./STRU',"r")
        for i in file:
            x=re.search('NUMERICAL_ORBITAL',i)
            count=count+1
            if x:
                break
        file.close()
    except FileExistsError as e:
            print(e)
    with open('./STRU') as f:
        info=f.readlines()
    orb=info[count]
    element=orb.split('_')[0]
    po=orb.split('_')[2].find("a")
    Rcut=int(orb.split('_')[2][:po])
    pos = orb.split('_')[3].find("R")
    Ecut = int(orb.split('_')[3][:pos])
    position = orb.split('_')[4].split()[0].find(".")
    N_orb = orb.split('_')[4].split()[0][:position]
    l_max = N_orb[-1]
    if l_max == 's':
        Nl=1
        Nu=[int(N_orb[0])]
    elif l_max == 'p':
        Nl=2
        Nu=[int(N_orb[0]), int(N_orb[2])]
    elif l_max == 'd':
        Nl=3
        Nu=[int(N_orb[0]), int(N_orb[2]), int(N_orb[4])]
    elif l_max == 'f':
        Nl=4
        Nu=[int(N_orb[0]), int(N_orb[2]), int(N_orb[4]), int(N_orb[6])]
    else:
        raise ValueError("l_max is not one of spdf")
        
    with open('./ORBITAL_RESULTS.txt','r+') as f:
            flist=f.readlines()
    first_index=flist.index('\tType\tL\tZeta-Orbital\n')
    second_index=flist.index('\tType\tL\tZeta-Orbital\n',first_index+1)
    Ne=second_index-first_index-2
    info_element={element: {'index': 0, 'Nu': Nu, 'Nl': Nl, 'Rcut': Rcut, 'dr': 0.01, 'Ecut': Ecut, 'Ne': Ne}}
    print(f"info_element: {info_element}, reading from STRU and ORBITAL_RESULTS.txt.")
    
    return info_element
          
# receive a list of number 'l'
# return string like "3s3p2d"
def get_orb_str(l_list):
    # 定义轨道角动量量子数列表对应的转换字典
    l_dict = {0: 's', 1: 'p', 2: 'd', 3: 'f', 4: 'g'}
    
    # 定义字符串变量，用于存储转换结果
    l_str = ''
    
    # 遍历数字列表，并进行相应的转换操作
    for i in range(len(l_list)):
        l_str += '{}{}'.format(l_list[i], l_dict[i])  # 添加数字和字母的组合
    return l_str

# receive order
# return name of element
def get_atomic_number(element):
    element_dict = {
        'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8,
        'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15,
        'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22,
        'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29,
        'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36,
        'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43,
        'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50,
        'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57,
        'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64,
        'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71,
        'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78,
        'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Th': 90, 'Pa': 91,
        'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98,
        'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103
    }
    return element_dict.get(element.capitalize())

# check convergence
# receive "filename"
# return "Yeah!" or "No!"
def convergence_test(filesource):
    a=checklog(filesource, wordcheck='!! CONVERGENCE HAS NOT BEEN ACHIEVED !!')
    if a=={}:
        return "Y"
    else:
        return "N"
    
# get DFT etot (Ha)
def get_etot(filesource='./single_Ne.out'):
    a=checklog(filesource, wordcheck='etot')
    for i in range(len(str(a))):
        if(str(a)[i]==':' ):
            for j in range(len(str(a))):
                if(str(a)[j]=='n'):
                    if(j>i):
                        e_pbe=str(a)[i+1:j-1]
    if float(e_pbe):
        pass
    else:
        print("cannot get e_pbe")
    return float(e_pbe)

# get Etot_without_rpa (Ha)
# Etot_without_rpa = etot - exc + exx
def get_Etot_without_rpa(filesource='./single_Ne.out'):
    a=checklog(filesource, wordcheck='Etot_without_rpa')
    for i in range(len(str(a))):
        if(str(a)[i]==':' ):
            for j in range(len(str(a))):
                if(str(a)[j]=='n'):
                    if(j>i):
                        e_pbe_without_RPA=str(a)[i+1:j-1]
    if float(e_pbe_without_RPA):
        pass
    else:
        print("cannot get e_pbe_without_RPA")
    return float(e_pbe_without_RPA)

# get cRPA (Ha)
def get_cRPA(filesource='./LibRPA_single_Ne.out'):
    a=checklog(filesource, wordcheck='Total EcRPA')
    for i in range(len(str(a))):
        if(str(a)[i]==':' ):
            break
    for j in range(len(str(a))):
        if(str(a)[j]=='n'):
            break
    rpa=str(a)[i+1:j-1]
    if float(rpa):
        pass
    else:
        print("cannot get E_cRPA")
    return float(rpa)

# get cRPA without gamma point due to mishandle for Gamma in LibRPA
# just for test
def get_cRPA_without_gamma(filesource='./LibRPA_single_Ne.out'):
    a=checklog(filesource, wordcheck='EcRPA without gamma contributing')
    for i in range(len(str(a))):
        if(str(a)[i]==':' ):
            break
    for j in range(len(str(a))):
        if(str(a)[j]=='\\'):
            break
    rpa=str(a)[i+1:j-1]
    if float(rpa):
        pass
    else:
        print("cannot get e_cRPA_without_gamma")
    return float(rpa)

#-------------------------------------------------
    #condition number
    #with open('./band_out') as f:
    #    eig=f.readlines()
    #Norb=int(eig[2].split()[0])#number of orbital per k point
    #eigenvalues0=[]
    #eigenvalues=[]#final eigenvalues arrangement
    #for i in range(len(eig)):
    #    if(len(eig[i].split())==4):
    #        eigenvalues0.append(eig[i].split()[2])
    #for i in range(len(eigenvalues0)):
    #    eigenvalues.append(float(eigenvalues0[i]))
    #eig_Gamma=eigenvalues[0:Norb]
    #eigenmin=min(eig_Gamma)
    #eigenmax=max(eig_Gamma)
    #Ncondition=abs(eigenmax/eigenmin)
#--------------------------------------------------



