import numpy as np
import matplotlib.pyplot as plt

def _save_orb(coefs, elem, ecut, rcut: int, nzeta, jY_type: str = "reduced"):
    """
    Plot the orbital and save .orb file
    The coefficients should be in the form of
    [it][l][zeta][q]
    """
    from cOpt.object.orbio import write_nao, write_param
    import os, uuid

    dr = 0.01
    r = np.linspace(0, rcut, int(rcut/dr)+1)

    chi = _build_orb(coefs, rcut, 0.01, jY_type)
    syms = "SPDFGHIKLMNOQRTUVWXYZ".lower()
    nz = nzeta
    suffix = "".join([f"{nz[j]}{syms[j]}" for j in range(len(nz))])
    folder, subfolder = f"{elem}_{suffix}", f"{rcut}au{ecut}Ry"
    os.makedirs(f"{folder}/{subfolder}", exist_ok=True)
    fpng = f"{elem}_gga_{rcut}au_{ecut}Ry_{suffix}.png"
    plot_chi(chi, r, save=fpng)
    os.rename(fpng, os.path.join(f"{folder}/{subfolder}",fpng))
    plt.close()
    forb = fpng.replace(".png", ".orb")
    write_nao(forb, elem, ecut, rcut, len(r), dr, chi)
    fparam = str(uuid.uuid4())
    write_param(fparam, coefs[0], rcut, 0.0, elem)
    os.rename(forb, os.path.join(f"{folder}/{subfolder}", forb))
    os.rename(fparam, os.path.join(f"{folder}/{subfolder}", "ORBITAL_RESULTS.txt"))
    #print(f"orbital saved as {forb}")

def _build_orb(coefs, rcut, dr: float = 0.01, jY_type: str = "reduced"):
    """build real space grid orbital based on the coefficients of the orbitals,
    rcut and grid spacing dr. The coefficients should be in the form of
    [it][l][zeta][q].
    
    Args:
    coefs: list of list of list of list of float
        the coefficients of the orbitals
    rcut: float
        the cutoff radius
    dr: float
        the grid spacing
    jY_type: str
        the type of jY basis, can be "reduced", "nullspace", "svd" or "raw"

    Returns:
    np.ndarray: the real space grid data
    """
    from cOpt.object.radial import build_raw, build_reduced, coeff_normalized2raw
    import numpy as np

    r = np.linspace(0, rcut, int(rcut/dr)+1) # hard code dr to be 0.01? no...
    if jY_type in ["reduced", "nullspace", "svd"]:
        chi = build_reduced(coefs[0], rcut, r, True)
    else:
        coefs = coeff_normalized2raw(coefs, rcut)
        chi = build_raw(coefs[0], rcut, r, 0.0, True, True)
    return chi

def plot_chi(chi, r, save=None):
    lmax = len(chi)-1
    nzeta = [len(chi_l) for chi_l in chi]

    fig, ax = plt.subplots(nrows=1, ncols=lmax+1, figsize=((lmax+1)*7, 6), layout='tight', squeeze=False)
    
    for l, chi_l in enumerate(chi):
        for zeta, chi_lz in enumerate(chi_l):
            # adjust the sign so that the largest value is positive
            if chi_lz[np.argmax(np.abs(chi_lz))] < 0:
                chi_lz *= -1
    
            ax[0, l].plot(r, chi_lz, label='$\zeta=%i$'%(zeta))
    
        ax[0, l].legend(fontsize=16)
        ax[0, l].axhline(0, color='black', linestyle=':')
        ax[0, l].set_title('$l=%d$'%(l), fontsize=20)
        ax[0, l].set_xlim([0, r[-1]])

    if save is not None:
        plt.savefig(save)

if __name__ == "__main__":
    from cOpt.object.orbio import read_param
    a = read_param("./ORBITAL_RESULTS.txt")
    print(a)
    _save_orb([a["coeff"]],a["elem"],100,a["rcut"],[6,3,2])
