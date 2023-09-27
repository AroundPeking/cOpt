#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##
from scipy.optimize import minimize
import numpy as np
import os
import subprocess
# my code
import get_info
import IO
import abfs
from C_to_orb import create_orb_files


def obj(x0, info_element, new_dir, abf_dir, fix, mod, abf, abacus, abacus_abf, librpa, fre_disp, iter_name):
    global obj_ini
    global flag
    
    element = list(info_element.keys())
    os.chdir("./"+new_dir)
    
    #----------------------------------------------------
    #cp all input files in sub-dir,i.e. "opt_DZP_by_add_"+slurm_id
    #re-write ORBITAL_RESULTS.txt and run C_to_orb.py in sub-sub-dir, i.e. number name 1, 2, ...
    sys_run_str = '''
mkdir {0}
cp ../ORBITAL_RESULTS.txt ./{0}
cp ../{1}_ONCV_PBE-1.0.upf ./{0}
cp ../INPUT ./{0}
cp ../KPT ./{0}
cp ../STRU ./{0}
'''.format(str(flag), element[0])
    #sys.stdout.flush() 
    subprocess.run( [sys_run_str, "--login"], shell=True, text=True, stdin=subprocess.DEVNULL)
    #sys.stdout.flush() 
    
    #sub-sub-dir, i.e. number name 1, 2, ...
    os.chdir("./"+str(flag))
    
    #re-write ORBITAL_RESULTS.txt
    IO.write_orb(x0, info_element, fix, mod, file = './ORBITAL_RESULTS.txt')
    #generate .orb file-------------------------------
    create_orb_files(info_element)
    
    abfs.generate_abfs(info_element, abf_dir, abacus_abf, abf, mod)

    # run abacus and librpa
    orb_str = get_info.get_orb_str(info_element[element[0]]['Nu'])
    atom_num = get_info.get_atomic_number(element[0])
    run_abacus = '''
cp ./ORBITAL_{1}U.dat ./{0}_gga_{2}au_{3}Ry_{4}.orb
{5} >single_{0}.out
{6} > LibRPA_single_{0}.out
'''.format(element[0], atom_num, info_element[element[0]]['Rcut'], info_element[element[0]]['Ecut'], orb_str, abacus, librpa)
    #sys.stdout.flush() 
    subprocess.run( [run_abacus, "--login"], shell=True, text=True, stdin=subprocess.DEVNULL)
    #sys.stdout.flush() 
    
    cRPA=get_info.get_cRPA("LibRPA_single_"+element[0]+".out")
    obj=cRPA
    if(flag==0):
        obj_ini=obj
    obj_change=obj-obj_ini
        
    if flag % fre_disp == 0:
        # all units of energy here are Hartree
        convg=get_info.convergence_test("single_"+element[0]+".out")
        E_withoutRPA=get_info.get_Etot_without_rpa("single_"+element[0]+".out")
        E_pbe=get_info.get_etot("single_"+element[0]+".out")
        E_tot=E_withoutRPA+cRPA
        
    
    #-------------print info to iter.out---------------
    os.chdir("../..")
    
    if flag % fre_disp == 0:
        # all units of energy printed are eV
        IO.write_iter(iter_name, flag, convg, cRPA, E_pbe, E_tot, obj_change)
    
    flag += 1

    return float(obj)


if __name__=="__main__":
    global flag
    
    ############################################
    # change varibles here
    fix = [1, 1, 0]
    mod = [2, 2, 1]
    # augment abf, e.g. 1f
    abf = [0, 0, 0, 1]
    maxiter = 5000
    
    #dif_sum = sum((x - y) for x, y in zip(mod, fix))
    #fre_disp = info_element[element[0]]['Ne'] * dif_sum
    fre_disp = 10
    work_dir = "/home/ghj/SIAB/ABACUS-orbitals/SIAB/atom_opt/re-opt/DZP_SZ1f/work"
    abacus = "mpirun -n 1 -ppn 1 -env OMP_NUM_THREADS=42 /home/ghj/abacus/230820/abacus-develop/build/abacus"
    abacus_abf = "mpirun -n 1 -ppn 1 -env OMP_NUM_THREADS=42 /home/ghj/abacus/abacus_abfs/abacus-develop/build/abacus"
    librpa = "mpirun -n 1 -ppn 1 -env OMP_NUM_THREADS=42 /home/ghj/abacus/LibRPA/build/chi0_main.exe 16 0"
    abf_dir = "/home/ghj/SIAB/ABACUS-orbitals/SIAB/atom_opt/re-opt/DZP_SZ1f/gen_abf"
    ############################################
    
    os.chdir(work_dir)
    flag = 0
    slurm_id = os.environ["SLURM_JOB_ID"]
    new_dir = "opt_"+slurm_id
    os.system("mkdir "+new_dir)
    iter_name = "./iter."+slurm_id+".out"
    
    info_element = get_info.get_info_element()
    element = list(info_element.keys())
    
    #IO.write_info_element(info_element, filesource = './C_to_orb.py')
    
    x0 = IO.read_orb(info_element, fix, mod, file = './ORBITAL_RESULTS.txt')
    IO.write_iter_header(iter_name)
    
    res = minimize(obj, x0, method='BFGS', args=(info_element, new_dir, abf_dir, fix, mod, abf, abacus, abacus_abf, librpa, fre_disp, iter_name), options={'maxiter': maxiter, 'gtol': 1e-6})
    print("Local minimum: x = %s , f(x) = %s" % (res.x, res.fun))
    print("number of iteration for local minization: %d (nit)" %res.nit)
    print('number of total iteration:%d'%flag)
    print(res.message)
    
