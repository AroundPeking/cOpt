def read_json(fname: str = "./opt.json"):
    """read user input `opt.json`
    {
    "dft" : "rpa_pbe",
    "dimer_len" : [1.0, 1.2, 1.4, 1.8, 2.3],
    "fix" : [1, 1, 0],
    "mod" : [3, 3, 2],
    "init_chg" : "false", 
    "maxiter" : 10000,
    "opt_method" : "local opt", 
    "method" : "BFGS", 
    "fre_disp" : 20,
    "work_dir" : "/home/ghj/SIAB/ABACUS-orbitals/SIAB/atom_opt/re-opt/hf/F",
    "abacus_inputs" : "/home/ghj/SIAB/ABACUS-orbitals/SIAB/atom_opt/re-opt/hf/F",
    "abacus" : "mpirun -n 1 -ppn 1 -env OMP_NUM_THREADS=48 /home/ghj/abacus/gw_librpa_240511/abacus-develop/build/abacus",
    "librpa" : "mpirun -n 1 -ppn 1 -env OMP_NUM_THREADS=48 /home/ghj/abacus/LibRPA/build/chi0_main.exe 16 0"
    }
    """
    import json
    with open(fname, "r") as f:
        result = json.load(f)
    
    return result


def read_orb(file = './ORBITAL_RESULTS.txt'):
    """ 
    read ORBITAL_RESULTS.txt
    return all c_q
    """
    x0 = []
    
    with open(file) as f:
        flist = f.readlines()
    for line in flist:
        if line == '    Type   L   Zeta-Orbital\n' \
        or line == '</Coefficient>\n':
            if "tmp" in locals().keys():
                x0.append(tmp)
            tmp = []
        else:
            try: 
                tmp.append(float(line.strip()))
            except:
                continue
    return x0

def get_initx(x_all: list, fix: list, mod: list):
    """get x0
    in x_all, remove fix part, reserve mod part.
    TZDP: fix DZP, modify the extra part (3th s 3th p 2th d)
    """
    x0 = []
    current_index = 0
    for f, m in zip(fix, mod):
        start_index = current_index + f
        end_index = current_index + m
        x0.extend(x_all[start_index:end_index])
        print(x_all[start_index:end_index])
        current_index += m

    return x0
                

if __name__ == "__main__":
    x = read_orb("./ORBITAL_RESULTS.txt")
    print(x)
    #for i in x:
    #    print(len(i))
    #x0 = get_initx(x, [6, 2, 1], [6, 3, 2])
    #print(x0)
