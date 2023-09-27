#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##
import subprocess
import os
import IO
import copy
import get_info
from C_to_orb import create_orb_files

# generate .abfs in "./"
def generate_abfs(info_element, abf_dir, abacus_abf, abf, orb):
    pwd = os.getcwd()
    element = list(info_element.keys())
    fix = []
    for i in range(len(orb)):
        fix.append(0)

    x = IO.read_orb(info_element, fix, mod = orb, file = './ORBITAL_RESULTS.txt')
    os.chdir(abf_dir)
    for i in abf:
        if i:
            orb.append(i)
    IO.write_orb(x, info_element, fix = abf, mod = orb, file = './ORBITAL_RESULTS.txt')
    # !!!!!!!!!!!!!!!!!!!!!!!!!
    # copy python dict variables 
    # shoubld be careful
    info_element_abf = copy.deepcopy(info_element)
    info_element_abf[element[0]]['Nu'] = orb
    info_element_abf[element[0]]['Nl'] = len(orb)
    
    create_orb_files(info_element_abf)
    orb_str = get_info.get_orb_str(info_element_abf[element[0]]['Nu'])
    atom_num = get_info.get_atomic_number(element[0])
    shell = '''
    cp ./ORBITAL_{3}U.dat ./{0}_gga_{4}au_{5}Ry_{6}.orb
    {1} > abf.out
    cp ./{0}.abfs {2}
    '''.format(element[0], abacus_abf, pwd, atom_num, info_element[element[0]]['Rcut'], info_element[element[0]]['Ecut'], orb_str)
    subprocess.run( [shell, "--login"], shell=True, text=True, stdin=subprocess.DEVNULL)
    os.chdir(pwd)
    
    
    
    
