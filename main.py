#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##
from typing import ByteString
import numpy as np
import os
import subprocess
# my code
import get_info
import IO
import abfs
from C_to_orb import create_orb_files
import run_dft


def obj(x0, info_element, new_dir, abf_dir, fix, mod, abf, abacus, abacus_abf, librpa, fre_disp, iter_name, init_chg, orb_dir, dft):
    global obj_ini
    global best_obj
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
    if(init_chg=="true"):
        add_chg = '''
        cp ../SPIN*_CHG.cube ./{0}
        '''.format(str(flag))
        sys_run_str += add_chg
    #sys.stdout.flush() 
    subprocess.run( [sys_run_str, "--login"], shell=True, text=True, stdin=subprocess.DEVNULL)
    #sys.stdout.flush() 
    
    #sub-sub-dir, i.e. number name 1, 2, ...
    os.chdir("./"+str(flag))
    
    #re-write ORBITAL_RESULTS.txt
    IO.write_orb(x0, info_element, fix, mod, file = './ORBITAL_RESULTS.txt')
    #generate .orb file-------------------------------
    create_orb_files(info_element)
    
    #abfs.generate_abfs(info_element, abf_dir, abacus_abf, abf, mod)

    orb_str = get_info.get_orb_str(info_element[element[0]]['Nu'])
    atom_num = get_info.get_atomic_number(element[0])
    sys_run_str = '''
cp ./ORBITAL_{1}U.dat ./{0}_gga_{2}au_{3}Ry_{4}.orb
'''.format(element[0], atom_num, info_element[element[0]]['Rcut'], info_element[element[0]]['Ecut'], orb_str)
    subprocess.run( [sys_run_str, "--login"], shell=True, text=True, stdin=subprocess.DEVNULL)

    # run abacus and librpa
    if(dft == "rpa_pbe"):
        run_abacus = run_dft.rpa_pbe(element[0], abacus, librpa)
        run_dft.run(flag, run_abacus)
    elif(dft == "hf"):
        run_abacus = run_dft.hf(element[0], abacus)
        run_dft.run(flag, run_abacus)
    obj = run_dft.get_obj(dft, element) # Ha for rpa, eV for hf
    if(flag==0):
            obj_ini = obj
            best_obj = obj
            obj_change = obj-obj_ini
    
    convg=get_info.convergence_test("single_"+element[0]+".out")
    
    if (flag % fre_disp == 0) or (convg == "N"):
        if (dft == "rpa_pbe"):
            # all units of energy here are Hartree
            # convg=get_info.convergence_test("single_"+element[0]+".out")
            E_withoutRPA=get_info.get_Etot_without_rpa("single_"+element[0]+".out")
            E_pbe=get_info.get_etot("single_"+element[0]+".out")
            E_tot=E_withoutRPA+obj
        
    if (obj < best_obj):
        best_obj = obj
        IO.write_best_orb(flag, obj, obj_change, orb_dir)
    
    #-------------print info to iter.out---------------
    os.chdir("../..")
    
    if (flag % fre_disp == 0) or (convg == "N"):
        if (dft == "rpa_pbe"):
            # all units of energy printed are eV
            IO.write_iter_rpa_pbe(iter_name, flag, convg, obj, E_pbe, E_tot, obj_change)
        elif (dft == "hf"):
            IO.write_iter_hf(iter_name, flag, convg, obj, obj_change)
    
    flag += 1

    return obj*27.2113863


if __name__=="__main__":
    global flag
    
    input_info = IO.read_json('opt.json')
    for key, value in input_info.items():
        print(f"{key}: {value}", flush=True)
        exec(f"{key} = {value}") 
    
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
    
    os.chdir(work_dir)
    flag = 0
    slurm_id = os.environ["SLURM_JOB_ID"]
    new_dir = "opt_"+slurm_id
    orb_dir = work_dir+"/best_orb_"+slurm_id
    os.makedirs(new_dir, exist_ok=False)
    os.makedirs(orb_dir, exist_ok=False)
    IO.write_best_header(orb_dir)
    iter_name = "./iter."+slurm_id+".out"
    
    info_element = get_info.get_info_element()
    element = list(info_element.keys())
    
    #IO.write_info_element(info_element, filesource = './C_to_orb.py')
    
    x0 = IO.read_orb(info_element, fix, mod, file = './ORBITAL_RESULTS.txt')
    IO.write_iter_header(iter_name, dft)
    
    args=(info_element, new_dir, abf_dir, fix, mod, abf, abacus, abacus_abf, librpa, fre_disp, iter_name, init_chg, orb_dir, dft)
    # scipy
    if opt_method == "local opt":
            if method == "fmin":
                from scipy.optimize import fmin
                res = fmin(obj, x0, args=args)
                print("Local minimum: x = %s , f(x) = %s" % (res.xopt, res.fopt))
                print('number of total iteration:%d'%flag)
                print(res)
            else:
                from scipy.optimize import minimize
        # 'Nelder-Mead' 'Powell'
                res = minimize(obj, x0, method=method, tol=1e-7, args=args, options={'maxiter': maxiter})
                print("Local minimum: x = %s , f(x) = %s" % (res.x, res.fun))
                print("number of iteration for local minization: %d (nit)" %res.nit)
                print('number of total iteration:%d'%flag)
                print(res.message)
        
    # basinhopping
    elif opt_method == "global opt":
        from scipy.optimize import basinhopping
        minimizer_kwargs={"method": method, "tol":1e-7, "args":args}
        res = basinhopping(obj, x0, niter=maxiter, T=1.0, minimizer_kwargs=minimizer_kwargs)
        print("Glocal minimum: x = %s , f(x) = %s" % (res.x, res.fun))
        print("number of iteration for local minization: %d (nit)" %res.nit)
        print('number of total iteration:%d'%flag)
        print(res.message)

    # torch_optimizer
    elif opt_method == "nn opt":
        import torch
        from torch_optimizer import SWATS
        from torch.autograd import Variable
        
        xnn = Variable(torch.Tensor(x0), requires_grad=True)
        optimizer = SWATS([xnn], lr = 0.1, eps = 1e-20) 

        lr_decay = 1.0 # useless
        for i in range(maxiter):
            optimizer.zero_grad()  
            loss = obj(xnn.detach().numpy(), info_element, new_dir, abf_dir, fix, mod, abf, abacus, abacus_abf, librpa, fre_disp, iter_name, orb_dir)  
            loss.backward()  
            optimizer.step()  # update x
        
            for param_group in optimizer.param_groups:
                param_group['lr'] *= lr_decay
        
        print("Final result:")
        print(f"x = {x.data}")
        print(f"Minimum value of objective function: {obj(x)}")
        
    else:
        raise ValueError("Please check opt_method")
    
