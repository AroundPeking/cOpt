import cOpt.io.read_input as ciri
import cOpt.io.read_output as ciro
import cOpt.io.write_output as ciwo
import os
def initialize(fname: str = "./opt.json"):
    """initialization of numerical atomic orbitals opt task, 
    1. read opt.json file named as fname
    2. check the existence of ABACUS all input files, if input_check is True, go on
    3. unpack the input file, return set of parameters:
        - x0: start of optimization, often chosen as SIAB
        generation(https://github.com/kirk0830/abacus_orbital_generation). 
        - element_info: info of element.
        - user_setting: dict of environment settings, contain ABACUS and LibRPA executable, optimization method...
    """
    fname = fname.strip().replace("\\", "/")
    if fname == "":
        raise ValueError("No filename provided")
    user_setting = ciri.read_json(fname)

    os.chdir(user_setting["work_dir"])
    slurm_id = os.environ["SLURM_JOB_ID"]
    new_dir = user_setting["work_dir"]+"/opting_"+slurm_id
    orb_dir = user_setting["work_dir"]+"/best_orb_"+slurm_id
    os.makedirs(new_dir, exist_ok=False)
    os.makedirs(orb_dir, exist_ok=False)
    ciwo.write_best_header(orb_dir)
    iter_name = user_setting["work_dir"]+"/iter."+slurm_id+".out"
    info_element, pp = ciro.get_info_element(user_setting["abacus_inputs"])
    x_all = ciri.read_orb(file = user_setting["abacus_inputs"]+'/ORBITAL_RESULTS.txt')
    x0 = ciri.get_initx(x_all, user_setting["fix"], user_setting["mod"])
    ciwo.write_iter_header(iter_name, user_setting["dft"])
    
    return x0, info_element, user_setting, pp, new_dir, orb_dir, iter_name



if __name__ == "__main__":
    result = initialize("opt.json")
    #result = initialize("SIAB_INPUT.json")
    print(result)