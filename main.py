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
# import abfs
import run_dft


def obj(x0, info_element, new_dir, fix, mod, abacus, librpa, fre_disp, iter_name, init_chg, orb_dir, dft, dimer_len, pp):
    global obj_ini
    global best_obj
    global flag
    
    element = list(info_element.keys())
    os.chdir("./"+new_dir)
    
    if(dft == "rpa_pbe"):
        obj, convg = run_dft.one_iter_rpa(element[0], flag, init_chg, x0, info_element, fix, mod, abacus, librpa, pp)
    elif(dft == "hf"):
        obj, convg = run_dft.one_iter_hf(element[0], flag, init_chg, x0, info_element, fix, mod, abacus, dimer_len, pp)

    if(flag==0):
        obj_ini = obj
        best_obj = obj
        obj_change = obj-obj_ini
    
    if (flag % fre_disp == 0) or (convg == "N"):
        if (dft == "rpa_pbe"):
            # all units of energy here are eV
            E_withoutRPA=get_info.get_Etot_without_rpa("./"+flag+"/"+"single_"+element[0]+".out")
            E_pbe=get_info.get_etot("./"+flag+"/"+"single_"+element[0]+".out")
            E_tot=E_withoutRPA+obj
        
    if (obj < best_obj):
        best_obj = obj
        IO.write_best_orb(flag, obj, obj_change, orb_dir)
    
    #-------------print info to iter.out---------------
    os.chdir("..")
    
    if (flag % fre_disp == 0) or (convg == "N"):
        if (dft == "rpa_pbe"):
            # all units of energy printed are eV
            IO.write_iter_rpa_pbe(iter_name, flag, convg, obj, E_pbe, E_tot, obj_change)
        elif (dft == "hf"):
            IO.write_iter_hf(iter_name, flag, convg, obj, obj_change)
    
    flag += 1

    return obj


if __name__=="__main__":
    global flag
    
    input_info = IO.read_json('opt.json')
    for key, value in input_info.items():
        print(f"{key}: {value}", flush=True)
        exec(f"{key} = {value}") 
    
    os.chdir(work_dir)
    flag = 0
    slurm_id = os.environ["SLURM_JOB_ID"]
    new_dir = "opt_"+slurm_id
    orb_dir = work_dir+"/best_orb_"+slurm_id
    os.makedirs(new_dir, exist_ok=False)
    os.makedirs(orb_dir, exist_ok=False)
    IO.write_best_header(orb_dir)
    iter_name = "./iter."+slurm_id+".out"
    
    info_element, pp = get_info.get_info_element()
    element = list(info_element.keys())
    
    x0 = IO.read_orb(info_element, fix, mod, file = './ORBITAL_RESULTS.txt')
    IO.write_iter_header(iter_name, dft)
    
    args=(info_element, new_dir, fix, mod, abacus, librpa, fre_disp, iter_name, init_chg, orb_dir, dft, dimer_len, pp)
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
            loss = obj(xnn.detach().numpy(), info_element, new_dir, fix, mod, abacus, librpa, fre_disp, iter_name, orb_dir)  
            loss.backward()  
            optimizer.step()  # update x
        
            for param_group in optimizer.param_groups:
                param_group['lr'] *= lr_decay
        
        print("Final result:")
        print(f"x = {x.data}")
        print(f"Minimum value of objective function: {obj(x)}")
        
    else:
        raise ValueError("Please check opt_method")
    
