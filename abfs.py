#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##
import subprocess
import os
import IO

# generate .abfs in "./"
def generate_abfs(info_element, abf_dir, abacus_abf, abf, orb):
    pwd = os.getcwd()
    element = list(info_element.keys())
    fix = []
    for i in range(len(orb)):
        fix.append(0)

    x = IO.read_orb(info_element, fix, mod = orb, file = './ORBITAL_RESULTS.txt')
    os.system("cd "+abf_dir)
    for i in abf:
        if abf:
            orb.append(abf)
    IO.write_orb(x, info_element, fix = abf, mod = orb, file = './ORBITAL_RESULTS.txt')
    shell = '''
    {1} > abf.out
    cp ./{0}.abfs {2}
    cd {2}
    '''.format(element[0], abacus_abf, pwd)
    subprocess.run( [shell, "--login"], shell=True, text=True, stdin=subprocess.DEVNULL)
    
    
    
    