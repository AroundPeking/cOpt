#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##
import os
import cOpt.io.read_output as ciro
import cOpt.io.write_output as ciwo
import cOpt.object.cal_exe as coce
import cOpt.driver
def obj(x0, info_element, new_dir, iter_name, orb_dir, pp, coef_init: list, user_setting: dict):
    """
    loss function
    """
    global obj_ini
    global best_obj
    global flag
    flag = cOpt.driver.flag

    element = info_element.keys()
    os.chdir(new_dir)
    
    if(user_setting["dft"] == "rpa_pbe"):
        obj, convg = coce.one_iter_rpa(x0, info_element, pp, flag, coef_init, user_setting)
    elif(user_setting["dft"] == "hf"):
        obj, convg = coce.one_iter_hf(x0, info_element, user_setting["fix"], user_setting["mod"], user_setting["abacus"], user_setting["dimer_len"], pp, flag)

    if(flag == 0):
        obj_ini = obj
        best_obj = obj
    obj_change = obj - obj_ini
    
    if (flag % user_setting["freq_disp"] == 0) or (convg == "N"):
        if (user_setting["dft"] == "rpa_pbe"):
            # all units of energy here are eV
            E_withoutRPA=ciro.get_Etot_without_rpa("./"+str(flag)+"/"+"single_"+element+".out")
            E_pbe=ciro.get_etot("./"+str(flag)+"/"+"single_"+element+".out")
            E_tot=E_withoutRPA+obj
    
    #---------------- best_orb ------------------------
    if (obj < best_obj):
        best_obj = obj
        ciwo.write_best_orb(flag, obj, obj_change, orb_dir, user_setting, info_element)
    
    #-------------print info to iter.out---------------
    os.chdir("..")
    
    if (flag % user_setting["freq_disp"] == 0) or (convg == "N"):
        if (user_setting["dft"] == "rpa_pbe"):
            # all units of energy printed are eV
            ciwo.write_iter_rpa_pbe(iter_name, flag, convg, obj, E_pbe, E_tot, obj_change)
        elif (user_setting["dft"] == "hf"):
            ciwo.write_iter_hf(iter_name, flag, convg, obj, obj_change)
    
    flag += 1
    cOpt.driver.flag = flag
    #print(cOpt.driver.flag)

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
    
