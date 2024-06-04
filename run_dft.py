#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##
import subprocess
import time
import get_info
import IO
import os
from C_to_orb import create_orb_files

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
        obj = get_info.get_cRPA("LibRPA_single_"+element+".out")
    elif(dft == "hf"):
        obj = get_info.get_hf("OUT.ABACUS/running_scf.log")

    return obj # eV

# RPA_PBE
def one_iter_rpa(element, flag, init_chg, x0, info_element, fix, mod, abacus, librpa):
    cp_files_rpa(element, flag, init_chg)
    #sub-sub-dir, i.e. number name 1, 2, ...
    os.chdir("./"+str(flag))
    #re-write ORBITAL_RESULTS.txt
    IO.write_orb(x0, info_element, fix, mod, file = './ORBITAL_RESULTS.txt')
    #generate .orb file-------------------------------
    create_orb_files(info_element)

    orb_str = get_info.get_orb_str(info_element[element]['Nu'])
    atom_num = get_info.get_atomic_number(element)
    sys_run_str = '''
cp ./ORBITAL_{1}U.dat ./{0}_gga_{2}au_{3}Ry_{4}.orb
'''.format(element, atom_num, info_element[element]['Rcut'], info_element[element]['Ecut'], orb_str)
    subprocess.run( [sys_run_str, "--login"], shell=True, text=True, stdin=subprocess.DEVNULL)
    abacus_dir = sys_rpa_pbe(element, abacus, librpa)
    run_ABACUS(flag, abacus_dir)
    obj = get_obj("rpa_pbe", element) # eV
    convg = get_info.convergence_test("single_"+element+".out")
    os.chdir("..")

    return obj, convg

# RPA_PBE
def cp_files_rpa(element, flag, init_chg):
    sys_run_str = '''
mkdir {0}
cp ../ORBITAL_RESULTS.txt ./{0}
cp ../{1}_ONCV_PBE-1.0.upf ./{0}
cp ../INPUT ./{0}
cp ../KPT ./{0}
cp ../STRU ./{0}
'''.format(str(flag), element)
    if(init_chg=="true"):
        add_chg = '''
        cp ../SPIN*_CHG.cube ./{0}
        '''.format(str(flag))
        sys_run_str += add_chg
    #sys.stdout.flush() 
    subprocess.run( [sys_run_str, "--login"], shell=True, text=True, stdin=subprocess.DEVNULL)
    #sys.stdout.flush() 

# Hartree-Fock
def one_iter_hf(element, flag, init_chg, x0, info_element, fix, mod, abacus, dimer_len):
    os.makedirs(str(flag), exist_ok=False)
    os.chdir("./"+str(flag))
    obj = 0.0
    convg = "Y"
    for ilen in dimer_len:
        cp_files_hf(element, flag, init_chg, ilen)
        #sub-sub-sub-dir, i.e. number name 1/2.0, 1/3.0, ...
        os.chdir("./"+str(ilen))
        #re-write ORBITAL_RESULTS.txt
        IO.write_orb(x0, info_element, fix, mod, file = './ORBITAL_RESULTS.txt')
        #generate .orb file-------------------------------
        create_orb_files(info_element)

        orb_str = get_info.get_orb_str(info_element[element]['Nu'])
        atom_num = get_info.get_atomic_number(element)
        sys_run_str = '''
cp .    /ORBITAL_{1}U.dat ./{0}_gga_{2}au_{3}Ry_{4}.orb
'''.    format(element, atom_num, info_element[element]['Rcut'], info_element[element]['Ec  ut'], orb_str)
        subprocess.run( [sys_run_str, "--login"], shell=True, text=True, stdin=subprocess.  DEVNULL)
        abacus_dir = sys_hf(element, abacus)
        run_ABACUS(flag, abacus_dir)
        obj += get_obj("hf", element) # eV
        loc_convg = get_info.convergence_test("single_"+element+".out")
        if(loc_convg == "N"):
            convg = "N"
        os.chdir("..")
    os.chdir("..")

    return obj/len(dimer_len), convg

# Hartree-Fock
def cp_files_hf(element, flag, init_chg, ilen):
    sys_run_str = '''
cd {0}
mkdir {2}
cp ../../ORBITAL_RESULTS.txt ./{2}
cp ../../{1}_ONCV_PBE-1.0.upf ./{2}
cp ../../INPUT ./{2}
cp ../../KPT ./{2}
cp ../../STRU ./{2}
sed -i "s/distance/$(echo "scale=5; {2}" | bc)/g" ./{2}/STRU
cd ..
'''.format(str(flag), element, ilen)
    if(init_chg=="true"):
        add_chg = '''
        cp ../SPIN*_CHG.cube ./{0}/{1}
        '''.format(str(flag), ilen)
        sys_run_str += add_chg
    #sys.stdout.flush() 
    subprocess.run( [sys_run_str, "--login"], shell=True, text=True, stdin=subprocess.DEVNULL)
    #sys.stdout.flush() 