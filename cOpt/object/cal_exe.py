#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##
import subprocess
import time
import os
import cOpt.io.read_output as ciro
import cOpt.io.write_output as ciwo
import cOpt.object.coef2orb as coco
from cOpt.object.orbio import read_param

def sys_rpa_pbe(element, abacus, librpa):
    run_abacus = '''
{1} >single_{0}.out
{2} > LibRPA_single_{0}.out
'''.format(element, abacus, librpa)
    return run_abacus

def sys_hf(element, abacus):
    run_abacus = '''
{1} >single_{0}.out
'''.format(element, abacus)
    return run_abacus

def run_ABACUS(flag, run_abacus, max_attempts = 5):
    # Number of attempts allowed

    # Retry loop
    for attempt in range(1, max_attempts + 1):
        try:
            # Run the command
            subprocess.run([run_abacus, "--login"], shell=True, text=True, stdin=subprocess.DEVNULL, check=True)

            # If the command succeeds, break out of the loop
            break
        except subprocess.CalledProcessError as e:
            # If the command fails, print an error message
            print(f"Iter{flag}: Attempt {attempt} failed with exit code {e.returncode}. Retrying...", flush=True)

            # Add a delay before retrying (optional)
            time.sleep(1)
    else:
        # If all attempts fail, print an error message and handle accordingly
        print(f"All {max_attempts} attempts failed. Exiting.")
        # You can raise an exception, log the failure, or take other appropriate actions here.

def get_obj(dft, element):
    if(dft == "rpa_pbe"):
        obj = ciro.get_cRPA("LibRPA_single_"+element+".out") + ciro.get_Etot_without_rpa("single_"+element+".out")
    elif(dft == "hf"):
        obj = ciro.get_hf("OUT.ABACUS/running_scf.log")

    return obj # eV

# RPA_PBE
def one_iter_rpa(x0, info_element, pp, flag, param, user_setting):

    import cOpt.object.orbio as coo
    element = list(info_element.keys())[0]
    cp_files_rpa(element, user_setting["init_chg"], pp, user_setting["abacus_inputs"], flag)
    #sub-sub-dir, i.e. number name 1, 2, ...
    os.chdir("./"+str(flag))
    #re-write ORBITAL_RESULTS.txt
    coo.rewrite_param("./ORBITAL_RESULTS.txt", x0, param, user_setting["fix"], user_setting["mod"], info_element[element]['Rcut'], element)
    #generate .orb file-------------------------------
    param = read_param("./ORBITAL_RESULTS.txt")
    coco._save_orb([param["coeff"]], element, info_element[element]['Ecut'], info_element[element]['Rcut'], info_element[element]['Nu'])

    orb_str = ciro.get_orb_str(info_element[element]['Nu'])
    sys_run_str = '''
cp ./{0}_{3}/{1}au{2}Ry/{0}_gga_{1}au_{2}Ry_{3}.orb .
'''.format(element, info_element[element]['Rcut'], info_element[element]['Ecut'], orb_str)
    subprocess.run( [sys_run_str, "--login"], shell=True, text=True, stdin=subprocess.DEVNULL)
    abacus_dir = sys_rpa_pbe(element, user_setting["abacus"], user_setting["librpa"])
    run_ABACUS(flag, abacus_dir)
    obj = get_obj("rpa_pbe", element) # eV
    convg = ciro.convergence_test("single_"+element+".out")
    os.chdir("..")

    return obj, convg

# RPA_PBE
def cp_files_rpa(element, init_chg, pp, abacus_inputs, flag):

    sys_run_str = '''
mkdir {0}
cp {3}/ORBITAL_RESULTS.txt ./{0}
cp {3}/{2} ./{0}
cp {3}/INPUT ./{0}
cp {3}/KPT ./{0}
cp {3}/STRU ./{0}
'''.format(str(flag), element, pp, abacus_inputs)
    if(init_chg=="true"):
        add_chg = '''
        cp {1}/SPIN*_CHG.cube ./{0}
        '''.format(str(flag), abacus_inputs)
        sys_run_str += add_chg
    #sys.stdout.flush() 
    subprocess.run( [sys_run_str, "--login"], shell=True, text=True, stdin=subprocess.DEVNULL)
    #sys.stdout.flush() 

# Hartree-Fock
def one_iter_hf(x0, info_element, fix, mod, abacus, dimer_len, pp, flag):

    element = list(info_element.keys())[0]
    os.makedirs(str(flag), exist_ok=False)
    os.chdir("./"+str(flag))
    obj = 0.0
    convg = "Y"
    for ilen in dimer_len:
        cp_files_hf(ilen, pp)
        #sub-sub-sub-dir, i.e. number name 1/2.0, 1/3.0, ...
        os.chdir("./"+str(ilen))
        #re-write ORBITAL_RESULTS.txt
        ciwo.write_orb(x0, info_element, fix, mod, file = './ORBITAL_RESULTS.txt')
        #generate .orb file-------------------------------
        param = read_param("./ORBITAL_RESULTS.txt")
        coco._save_orb([param["coefs"]], element, info_element[element]['Ecut'], info_element[element]['Rcut'], info_element[element]['Nu'])

        orb_str = ciro.get_orb_str(info_element[element]['Nu'])
        sys_run_str = '''
cp ./{0}_{3}/{1}au{2}Ry/{0}_gga_{1}au_{2}Ry_{3}.orb .
'''.format(element, info_element[element]['Rcut'], info_element[element]['Ecut'], orb_str)
        subprocess.run( [sys_run_str, "--login"], shell=True, text=True, stdin=subprocess.DEVNULL)
        abacus_dir = sys_hf(element, abacus)
        run_ABACUS(flag, abacus_dir)
        obj += get_obj("hf", element) # eV
        loc_convg = ciro.convergence_test("single_"+element+".out")
        if(loc_convg == "N"):
            convg = "N"
        os.chdir("..")
    os.chdir("..")

    return obj/len(dimer_len), convg

# Hartree-Fock
def cp_files_hf(ilen, pp, abacus_inputs):
    sys_run_str = '''
mkdir {0}
cp {2}/ORBITAL_RESULTS.txt ./{0}
cp {2}/{1} ./{0}
cp {2}/INPUT ./{0}
cp {2}/KPT ./{0}
cp {2}/STRU ./{0}
sed -i "s/distance/$(echo "scale=5; {0}" | bc)/g" ./{0}/STRU
cd ..
'''.format(ilen, pp, abacus_inputs)
    
    #sys.stdout.flush() 
    subprocess.run( [sys_run_str, "--login"], shell=True, text=True, stdin=subprocess.DEVNULL)
    #sys.stdout.flush() 