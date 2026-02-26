"""
DESCRIPTION
----
This is a correlation basis optimization based on ABACUS numerical atomic orbital.


"""
global flag 
flag = 0

import argparse
def initialize(command_line: bool = True):
    """initialize the whole workflow of orbital optimization:
    1. specify input script
    2. specify test mode
    return
    """
    welcome = """
=================================================================================================== 

Basis set optimization for correlation calculation

Github repo: https://github.com/AroundPeking/cOpt
           
See reference for more information.
===================================================================================================
    """
    print(welcome, flush = True)
    placeholder_1 = ""
    placeholder_2 = ""
    if command_line:
        parser = argparse.ArgumentParser(description=welcome)
        parser.add_argument(
            "-i", "--input", 
            type=str, 
            default="./opt.json",
            action="store",
            help="specify input script after -i or --input, default value is opt.json")
        parser.add_argument(
            "-t", "--test",
            type=bool, 
            default=False,
            action="store",
            help="test mode, default is False")
        parser.add_argument(
            "-v", "--version",
            type=str,
            default="0.1.0",
            action="store",
            help="specify the version of basis optimization, default is 0.1.0")
        args = parser.parse_args()

        placeholder_1 = args.input
        placeholder_2 = args.test
    else:
        placeholder_1 = "./opt.json"
        placeholder_2 = False

    return placeholder_1, placeholder_2

import cOpt.init.initial as cii
import cOpt.utils.listmanip as culm
import numpy as np


def _build_minimize_options(method: str, user_setting: dict):
    """Build method-specific scipy.optimize.minimize options."""
    maxiter = user_setting["maxiter"]
    options = {"maxiter": maxiter, "disp": True}

    if method == "Powell":
        options["ftol"] = user_setting.get("ftol", 1e-4)
        options["xtol"] = user_setting.get("xtol", 1e-4)
    elif method == "Nelder-Mead":
        options["xatol"] = user_setting.get("xatol", 1e-4)
        options["fatol"] = user_setting.get("fatol", 1e-4)
    elif method in ("BFGS", "CG"):
        options["gtol"] = user_setting.get("gtol", 1e-4)
        options["eps"] = user_setting.get("fd_eps", 1e-3)
    elif method == "L-BFGS-B":
        options["ftol"] = user_setting.get("ftol", 1e-9)
        options["gtol"] = user_setting.get("gtol", 1e-4)
        options["maxcor"] = user_setting.get("maxcor", 20)
        options["eps"] = user_setting.get("fd_eps", 1e-3)

    return options


def _build_bounds(method: str, x0, user_setting: dict):
    """Build optional bounds from opt.json.

    Use `coef_bounds: [lower, upper]` in opt.json to enable.
    """
    if "coef_bounds" not in user_setting:
        return None

    bound_methods = {"Nelder-Mead", "L-BFGS-B", "TNC", "SLSQP", "Powell", "trust-constr"}
    if method not in bound_methods:
        print(f"Warning: method {method} does not use bounds, ignore coef_bounds.", flush=True)
        return None

    lower, upper = user_setting["coef_bounds"]
    return [(float(lower), float(upper)) for _ in x0]


def run(fname: str, 
        test: bool = True):
    """run the whole workflow of optimization:
    1. read input and initial coefficient C_q
    2. call abacus and librpa to calculate cRPA of single atom
    3. call optimizer to optimize coefficients of spherical Bessel functions for forming numerical atomic orbitals
    
    Args:
        `fname`: input filename
        `opt_version`: by specifying different version, will package parse input in different way.  
    
    `dft`: str, only support  
    ```python
    "rpa_pbe" or "hf"
    ```
    `dimer_len`: list of float
    like 
    ```python
    [1.0, 1.2, 1.4, 1.8, 2.3]
    ```
    `fix` and `mod`: list of integer
    like 
    ```python
    [1, 1, 0] and [2, 2, 1] mean optimizing DZP while fixing SZ
    ```
    `init_chg`: bool, whether using init_chg to accelerate calculation of ABACUS. `hf` not suppose this.
    `maxiter`: integer, max No. of single atom calculation.
    `opt_method`: str, "local opt" or "global opt", recommand "local opt"
    `method`: str ,optimizing algorithm  
    ```python
    local: BFGS, CG, Nelder-Mead, ...
    global: basinhopping
    ```
    `freq_disp`: integer, frequency of displaying one-shot ABACUS calculation
    `work_dir`: str, directory of input and output files, default: current directory.
    `abacus` and `librpa`: str, like
    ```python
    abacus="mpirun -n 1 -ppn 1 -env OMP_NUM_THREADS=48 /home/ghj/abacus/240704_lcao_exx/abacus-develop-exxlip/build/abacus"
    librpa="mpirun -n 1 -ppn 1 -env OMP_NUM_THREADS=48 /home/ghj/abacus/gw_librpa_240511/LibRPA-develop/build/chi0_main.exe 16 0"
    ```
    """
    global flag
    flag = 0
    
    # read input, for each term, see above annotation for details
    coef_init, info_element, user_setting, pp, new_dir, orb_dir,\
    iter_name = cii.initialize(fname=fname)

    x0 = np.array(culm.flatten(coef_init))

    from cOpt.object.orbio import read_param
    param = read_param(user_setting["abacus_inputs"]+"/ORBITAL_RESULTS.txt")
    method = user_setting["method"]
    bounds = _build_bounds(method, x0, user_setting)
    options = _build_minimize_options(method, user_setting)
    tol = user_setting.get("tol", 1e-7)
    args = (info_element, new_dir, iter_name, orb_dir, pp, param["coeff"], user_setting)

    import cOpt.object.loss as col
    col.reset_runtime_state()
    # scipy
    if user_setting["opt_method"] == "local opt":
            if user_setting["method"] == "fmin":
                from scipy.optimize import fmin
                res = fmin(col.obj, x0, args=args)
                print("Local minimum: x = %s , f(x) = %s" % (res.xopt, res.fopt))
                print('number of total iteration:%d'%flag)
                print(res)
            else:
                from scipy.optimize import minimize
            # 'Nelder-Mead' 'Powell'
                res = minimize(
                    col.obj,
                    x0,
                    method=method,
                    tol=tol,
                    args=args,
                    options=options,
                    bounds=bounds
                )
                print("Local minimum: x = %s , f(x) = %s" % (res.x, res.fun))
                print("number of iteration for local minization: %d (nit)" %res.nit)
                print('number of total iteration:%d'%flag)
                print(res.message)
        
    # basinhopping
    elif user_setting["opt_method"] == "global opt":
        from scipy.optimize import basinhopping
        minimizer_kwargs={"method": user_setting["method"], "tol":1e-7, "args":args}
        res = basinhopping(col.obj, x0, niter=user_setting["maxiter"], T=1.0, minimizer_kwargs=minimizer_kwargs)
        print("Glocal minimum: x = %s , f(x) = %s" % (res.x, res.fun))
        print("number of iteration for local minization: %d (nit)" %res.nit)
        print('number of total iteration:%d'%flag)
        print(res.message)

        """
    # torch_optimizer
    elif user_setting["opt_method"] == "nn opt":
        import torch
        from torch_optimizer import SWATS
        from torch.autograd import Variable
        
        xnn = Variable(torch.Tensor(x0), requires_grad=True)
        optimizer = SWATS([xnn], lr = 0.1, eps = 1e-20) 

        lr_decay = 1.0 # useless
        for i in range(maxiter):
            optimizer.zero_grad()  
            loss = sol.obj(xnn.detach().numpy(), info_element, new_dir, fix, mod, abacus, librpa, fre_disp, iter_name, orb_dir)  
            loss.backward()  
            optimizer.step()  # update x
        
            for param_group in optimizer.param_groups:
                param_group['lr'] *= lr_decay
        
        print("Final result:")
        print(f"x = {x.data}")
        print(f"Minimum value of objective function: {obj(x)}")
        """

    else:
        raise ValueError("Please check opt_method")

import cOpt.include.citation as cicite
def finalize():
    """finalize the whole workflow of optimization:
    """
    print(cicite.citation(), flush = True)

import time
def main(command_line: bool = True):
    """entry point of the whole workflow of optimization"""
    t_start = time.time()
    
    fname, test = initialize(command_line=command_line)
    t_initialize = time.time()
    run(fname=fname, test=test)
    t_run = time.time()
    finalize()
    t_finalize = time.time()
    t_end = time.time()
    # print time statistics with format %.2f
    print(f"""TIME STATISTICS
---------------
{"initialize":20s} {t_initialize - t_start:10.2f} s
{"run":20s} {t_run - t_initialize:10.2f} s
{"finalize":20s} {t_finalize - t_run:10.2f} s
{"total":20s} {t_end - t_start:10.2f} s
""", flush = True)

if __name__ == '__main__':
    
    main(command_line=True)
