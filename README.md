# cOpt: basis set optimization for correlation calculation

This is an optimization module for basis set, based on [ABACUS](https://github.com/deepmodeling/abacus-develop "ABACUS git repo").

## Installation

```bash
git clone https://github.com/AroundPeking/cOpt.git
cd cOpt
conda create -n cOpt
conda activate cOpt
pip install .
```

## Usage

You should prepare `opt.json` to set optimization parameters. 

A folder `abacus_inputs` containing all inputs for ABACUS calculation in iterations, e.g. single atom RPA calculation. Besides, `ORBITAL_RESULTS.txt ` is needed to be the starting point. Here `SPIN*_CHG.cube` is recommended to speed up SCF calculation in every iteration.

A folder to store all outputs, i.e. `work_dir`.

Then run

```bash
python cOpt/cOpt/driver.py
```

you can get file `iter.$SLURM_JOB_ID.out`, folder `best_orb_$SLURM_JOB_ID` and `opting_$SLURM_JOB_ID.`

## PARAMETERS in opt.json

`dft`: str, support "rpa_pbe" or "hf" now.

`fix` and `mod`: list, e.g. if you want to optimize DZP while fixing SZ, set `fix=[1, 1, 0]` and `mod=[2, 2, 1]`. Don't recommend to optimize SZ by cRPA due to distorted DFT.

`init_chg:` bool, `true` is recommended to speed up SCF calculation in every iteration. If true, density `SPIN*_CHG.cube` should be put in `abacus_inputs`.

`maxiter`: int, maximum number of iterations.

`opt_method` and `method`: str. If `opt_method="local opt"`, "Powell" (recommended), "CG", "BFGS", ... are optional. If `opt_method="global opt"`, `method=basinhopping`.

`freq_disp`: int, freqency of displaying change of object in `iter.$SLURM_JOB_ID.out`.

`work_dir`: str, outputs directory.

`abacus_inputs`: str, ABACUS inputs directory.

`abacus` and `librpa`: str, ABACUS and LibRPA executable program.
