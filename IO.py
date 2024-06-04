#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##
import re
import numpy as np

# useless now (C_to_orb.py has been written as a function)
# modify info_element in the C_to_orb.py to generate `.orb` from ORBITAL_RESULTS.txt
def write_info_element(info_element, filesource='./C_to_orb.py'):
    i = 0
    with open(filesource,'r+') as f:
        script=f.readlines()
    for line in script:
        i += 1
        x = re.search("info_element=", line)
        if x:
            break
    script[i-1]="info_element="+str(info_element)+"\n"
    with open(filesource,'w+') as f:
        f.writelines(script)
        
# read `opt.json`
def read_json(filename):
    values = {}
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    for line in lines:
        # remove commentary
        line = re.sub(r'#.*', '', line)
        line = line.strip()
        if not line:
            continue
        # first "="
        match = re.match(r'([^=]+)=\s*(.*)', line)
        if match:
            key = match.group(1).strip()
            value = match.group(2).strip()
            values[key] = value
    
    return values
    ############################################
    # change varibles here
    #fix = [0, 0, 0]
    #mod = [2, 2, 1]
    ## augment abf, e.g. 1f
    #abf = [0, 0, 0, 1]
    #maxiter = 5000
    #opt_method = "local opt" # local opt / global opt(basinhopping) 
    #method = "BFGS" # explicit opt method, # 'Nelder-Mead'
    ##dif_sum = sum((x - y) for x, y in zip(mod, fix))
    ##fre_disp = info_element[element[0]]['Ne'] * dif_sum
    #fre_disp = 10
    #work_dir = "/home/ghj/SIAB/ABACUS-orbitals/SIAB/atom_opt/re-opt/DZP_SZ1f/nn_work"
    #abacus = "mpirun -n 1 -ppn 1 -env OMP_NUM_THREADS=48 /home/ghj/abacus/230820/abacus-develop/build/abacus"
    #abacus_abf = "mpirun -n 1 -ppn 1 -env OMP_NUM_THREADS=48 /home/ghj/abacus/abacus_abfs/abacus-develop/build/abacus"
    #librpa = "mpirun -n 1 -ppn 1 -env OMP_NUM_THREADS=48 /home/ghj/abacus/LibRPA/build/chi0_main.exe 16 0"
    #abf_dir = "/home/ghj/SIAB/ABACUS-orbitals/SIAB/atom_opt/re-opt/DZP_SZ1f/gen_abf"
    ############################################
    
        
# TZDP: fix DZP, modify the extra part (3th s 3th p 2th d)
# receive fixing part and modifying part, read ORBITAL_RESULTS.txt
# return x0
def read_orb(info_element, fix = [2, 2, 1], mod = [3, 3, 2], file = './ORBITAL_RESULTS.txt'):
    element = list(info_element.keys())[0]
    size = [ mod_i - fix_i for mod_i, fix_i in zip(mod, fix)]  
    x0 = []
    l = -1 # angular momentum
    
    with open(file) as f:
        flist=f.readlines()
    # i: modified part, 1 in [1, 1, 1] 
    # j: lines in OBITAL_RESULTS.txt, '\t  Ag \t0\t    1\n'
    # l: angular momentum
    for i in size:
        l += 1
        if i:
            is_append = False
            # read initial ORBITAL_RESULTS.txt
            for j in flist:
                if j == '\t  '+element+' \t'+str(l)+'\t    '+str(fix[l]+1)+'\n':
                    is_append = True
                if is_append:
                    try:
                        n = float(j.strip())
                        x0.append(n)
                    except ValueError:
                        pass
                        
                if j == '\t  '+element+' \t'+str(l+1)+'\t    1\n':
                    break
            
    return np.array(x0)
            
# receive array x
# write ORBITAL_RESULTS.txt
def write_orb(x, info_element, fix = [2, 2, 1], mod = [3, 3, 2], file = './ORBITAL_RESULTS.txt'):
    orb_str = []
    # get all coefficients for certain angular momentum
    lorb_str = []
    num = int(len(x))
    # from 1 dimension array to table   
    for i in range(num):
        orb_str.append('\t   '+str(x[i])+'\n')

    #write coefficients to ORBITAL_RESULTS.txt
    size = [ mod_i - fix_i for mod_i, fix_i in zip(mod, fix)]
    tot_size = 0
    element = list(info_element.keys())[0]
    Ne = info_element[element]['Ne']
    
    l = -1 # angular momentum
    with open(file,'r+') as f:
        flist=f.readlines()
    # i: modified part, 2 in [1, 2, 0] 
    # j: lines in OBITAL_RESULTS.txt, '\t  Ne \t0\t    1\n'
    # l: angular momentum
    # nline: get number of line needed to write
    for i in size:
        l += 1
        tot_size += i
        if i:
            lorb_str = orb_str[(tot_size-i)*Ne:tot_size*Ne]
            while i > 0:
                nline = -1
                for j in flist:
                    nline += 1
                    if j == '\t  '+element+' \t'+str(l)+'\t    '+str(mod[l]-i+1)+'\n':
                        break
                flist[nline+1:nline+1+Ne] = lorb_str[(size[l]-i)*Ne:(size[l]-i+1)*Ne]
                i -= 1
    # check length of x and modified coefficients
    if tot_size*Ne != num:
        raise ValueError("Length of array and modified C_i not match!")
    with open(file,'w+') as f:
        f.writelines(flist)
        
        
def write_iter_header(file, dft):
    if(dft == "rpa_pbe"):
        file = open(file, "a")
        header = "{:<6s} {:<10s} {:<15s} {:<15s} {:<15s} {:<15s}".format("Iter", "Convg", "cRPA(eV)", "E_pbe(eV)", "E_tot(eV)", "change(eV, vs iter0)")
        print(header, file=file)
        file.close()
    elif(dft == "hf"):
        file = open(file, "a")
        header = "{:<6s} {:<10s} {:<15s} {:<15s}".format("Iter", "Convg", "E_hf(eV)", "change(eV, vs iter0)")
        print(header, file=file)
        file.close()

    
def write_iter_rpa_pbe(file, flag, convg, crpa, e_pbe, e_tot, change):
    file = open(file, "a")
    line = "{:<6s} {:<10s} {:<15.8f} {:<15.8f} {:<15.8f} {:<15.8f}".format(str(flag), convg, crpa, e_pbe, e_tot, change)
    print(line, file=file)

    file.close()
    
    
def write_iter_hf(file, flag, convg, obj, obj_change):
    file = open(file, "a")
    # eV
    line = "{:<6s} {:<10s} {:<15.8f} {:<15.8f}".format(str(flag), convg, obj, obj_change)
    print(line, file=file)

    file.close()

def write_best_orb(flag, obj, obj_change, orb_dir):
    import subprocess
    file = open(orb_dir+"/best_orb_info.dat", "a")
    line = "{:<6s} {:<15.8f} {:<15.8f}".format(str(flag), obj, obj_change)
    print(line, file=file)

    file.close()
    sys_run_str = '''
cp ./ORBITAL_RESULTS.txt {0}
cp ./ORBITAL_PLOTU.dat {0}
cp ./*_gga_*au_*Ry_*.orb {0}
'''.format(orb_dir)

    #sys.stdout.flush() 
    subprocess.run( [sys_run_str, "--login"], shell=True, text=True, stdin=subprocess.DEVNULL)
    #sys.stdout.flush() 

def write_best_header(orb_dir):
    file = open(orb_dir+"/best_orb_info.dat", "a")
    header = "{:<6s} {:<15s} {:<15s}".format("Iter", "obj(eV)", "change(eV, vs iter0)")
    print(header, file=file)
    file.close()
