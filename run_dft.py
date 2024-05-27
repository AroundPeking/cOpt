#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##
import subprocess
import time
import get_info

def rpa_pbe(element, abacus, librpa):
    run_abacus = '''
{1} >single_{0}.out
{2} > LibRPA_single_{0}.out
'''.format(element, abacus, librpa)
    return run_abacus

def hf(element, abacus):
    run_abacus = '''
{1} >single_{0}.out
'''.format(element, abacus)
    return run_abacus

def run(flag, run_abacus, max_attempts = 5):
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
        cRPA = get_info.get_cRPA("LibRPA_single_"+element[0]+".out")
        obj = cRPA
    elif(dft == "hf"):
        ehf = get_info.get_hf("OUT.ABACUS/running_scf.log")
        obj = ehf
    return obj