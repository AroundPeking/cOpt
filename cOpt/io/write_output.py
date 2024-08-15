#!/usr/bin/env python3
# -*- coding: utf-8 -*-
##


def write_orb(x, info_element, fix = [2, 2, 1], mod = [3, 3, 2], file = './ORBITAL_RESULTS.txt'):
    """
    receive array x
    write ORBITAL_RESULTS.txt 
    """
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
    if tot_size*Ne <= num:
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

def write_best_orb(flag, obj, obj_change, orb_dir, user_setting, info_element):
    import subprocess
    import cOpt.io.read_output as ciro

    file = open(orb_dir+"/best_orb_info.dat", "a")
    line = "{:<6s} {:<15.8f} {:<15.8f}".format(str(flag), obj, obj_change)
    print(line, file=file)
    file.close()
    element = list(info_element.keys())[0]
    orb_str = ciro.get_orb_str(info_element[element]['Nu'])
    if(user_setting["dft"] == "rpa_pbe"):
        sys_run_str = '''
cp ./{1}/{2}_{3}/{4}au{5}Ry/* {0}
'''.format(orb_dir, flag, element, orb_str, info_element[element]['Rcut'], info_element[element]['Ecut'])
    elif(user_setting["dft"] == "hf"):
        sys_run_str = '''
cp ./{1}/{2}/{3}_{4}/{5}au{6}Ry/* {0}
'''.format(orb_dir, flag, user_setting["dimer_len"], element, orb_str, info_element[element]['Rcut'], info_element[element]['Ecut'])

    #sys.stdout.flush() 
    subprocess.run( [sys_run_str, "--login"], shell=True, text=True, stdin=subprocess.DEVNULL)
    #sys.stdout.flush() 

def write_best_header(orb_dir):
    file = open(orb_dir+"/best_orb_info.dat", "a")
    header = "{:<6s} {:<15s} {:<15s}".format("Iter", "obj(eV)", "change(eV, vs iter0)")
    print(header, file=file)
    file.close()

if __name__ == "__main__":
    write_best_orb(3, 0.1, 0.1, "./", "rpa_pbe", [1.0], {'C': {'index': 0, 'Nu': [2, 2, 1], 'Nl': 3, 'Rcut': 8, 'dr': 0.01, 'Ecut': 400, 'Ne': 49}})