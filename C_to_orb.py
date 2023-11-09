#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

import sys
import os
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'opt_orb_pytorch_dpsi'))

import orbital
import torch
import numpy as np
from scipy.special import spherical_jn
from scipy.optimize import fsolve
from util import *
import functools

#########derive E#####################################
def find_eigenvalue(Nl,Ne):
	""" E[il,ie] """
	E = np.zeros((Nl,Ne+Nl+1))
	for ie in range(1,Ne+Nl+1):
		E[0,ie] = ie*np.pi
	for il in range(1,Nl):
		jl = functools.partial(spherical_jn,il)
		for ie in range(1,Ne+Nl+1-il):
			E[il,ie] = fsolve( jl, (E[il-1,ie]+E[il-1,ie+1])/2 )
	return E[:,1:Ne+1]
def set_E(info_element):
	""" E[it][il,ie] """
	eigenvalue = { it:find_eigenvalue(info_element[it]['Nl'],info_element[it]['Ne']) for it in info_element }
	E = dict()
	for it in info_element:
		E[it] = torch.from_numpy(( eigenvalue[it]/info_element[it]['Rcut'] ).astype("float64"))
	return E

#########derive C#####################################
def random_C_init(info_element):
	""" C[it][il][ie,iu]	<jY|\phi> """
	C = dict()
	for it in info_element.keys():
		C[it] = ND_list(info_element[it]['Nl'])
		for il in range(info_element[it]['Nl']):
			C[it][il] = torch.tensor(np.random.uniform(-1,1, (info_element[it]['Ne'], info_element[it]['Nu'][il])), dtype=torch.float64, requires_grad=True)
	return C
	
def read_C_init(file_name,info_element):
	""" C[it][il][ie,iu]	<jY|\phi> """
	C = random_C_init(info_element)

	with open(file_name,"r") as file:
	
		for line in file:
			if line.strip() == "<Coefficient>":	
				line=None
				break
		ignore_line(file,1)
	
		C_read_index = set()
		while True:
			line = file.readline().strip()
			if line.startswith("Type"):
				it,il,iu = file.readline().split();	
				il = int(il)
				iu = int(iu)-1
				C_read_index.add((it,il,iu))
				line = file.readline().split()
				for ie in range(info_element[it]['Ne']):
					if not line:	line = file.readline().split()
					C[it][il].data[ie,iu] = float(line.pop(0))
			elif line.startswith("</Coefficient>"):
				break;
			else:
				raise IOError("unknown line in read_C_init "+file_name+"\n"+line)
	return C, C_read_index


periodtable = {   'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7,
                  'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13,
               'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19,
               'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25,
               'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31,
               'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37,
               'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43,
               'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49,
               'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55,
               'Ba': 56, #'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61,
                    ## 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67,
                    ## 'Er': 68, 'Tm': 69, 'Yb': 70, 
                    ## 'Lu': 71, 
               'Hf': 72, 'Ta': 73,
               'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79,
               'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 
                    ## 'Po': 84, #'At': 85,
                    ## 'Rn': 86, #'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91,
                    ## 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97,
                    ## 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103,
                    ## 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108,
                    ## 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112, 'Uut': 113,
                    ## 'Fl': 114, 'Uup': 115, 'Lv': 116, 'Uus': 117, 'Uuo': 118
               } 

def print_orbital(orb,info_element):
	""" orb[it][il][iu][r] """
	for it,orb_t in orb.items():
		#with open("orb_{0}.dat".format(it),"w") as file:
		with open("ORBITAL_{0}U.dat".format( periodtable[it] ),"w") as file:
			print_orbital_head(file,info_element,it)
			for il,orb_tl in enumerate(orb_t):
				for iu,orb_tlu in enumerate(orb_tl):
					print("""                Type                   L                   N""",file=file)
					print("""                   0                   {0}                   {1}""".format(il,iu),file=file)
					for ir,orb_tlur in enumerate(orb_tlu):
						print( '%.14e'%orb_tlur, end="  ",file=file)
						if ir%4==3:	print(file=file)
					print(file=file)
					
					
def plot_orbital(orb,Rcut,dr):
	for it,orb_t in orb.items():
		#with open("orb_{0}_plot.dat".format(it),"w") as file:
		with open("ORBITAL_PLOTU.dat", "w") as file:
			Nr = int(Rcut[it]/dr[it])+1
			for ir in range(Nr):
				print( '%10.6f'%(ir*dr[it]),end="  ",file=file)
				for il,orb_tl in enumerate(orb_t):
					for orb_tlu in orb_tl:
						print( '%18.14f'%orb_tlu[ir],end="  ",file=file)
				print(file=file)
				
				
#Use of attribution of info_element has been modified by Huanjing Gong
def print_orbital_head(file,info_element,it):
	print( "---------------------------------------------------------------------------", file=file )
	print( "Element                     {0}".format(it), file=file )
	print( "Energy Cutoff(Ry)           {0}".format(info_element[it]['Ecut']), file=file )
	print( "Radius Cutoff(a.u.)         {0}".format(info_element[it]['Rcut']), file=file )
	print( "Lmax                        {0}".format(info_element[it]['Nl']-1), file=file )
	l_name = ["S","P","D"]+list(map(chr,range(ord('F'),ord('Z')+1)))
	for il,iu in enumerate(info_element[it]['Nu']):
		print( "Number of {0}orbital-->       {1}".format(l_name[il],iu), file=file )
	print( "---------------------------------------------------------------------------", file=file )
	print( "SUMMARY  END", file=file )
	print( file=file )
	print( "Mesh                        {0}".format(int(info_element[it]['Rcut']/info_element[it]['dr'])+1), file=file )
	print( "dr                          {0}".format(info_element[it]['dr']), file=file )

def generate_orbital(info_element,C,E):
	""" C[it][il][ie,iu] """
	""" orb[it][il][iu][r] = \suml_{ie} C[it][il][ie,iu] * jn(il,ie*r) """
	orb = dict()
	for it in info_element:
		Nr = int(info_element[it]['Rcut']/info_element[it]['dr'])+1
		orb[it] = ND_list(info_element[it]['Nl'])
		for il in range(info_element[it]['Nl']):
			orb[it][il] = ND_list(info_element[it]['Nu'][il])
			for iu in range(info_element[it]['Nu'][il]):
				orb[it][il][iu] = np.zeros(Nr)
				for ir in range(Nr):
					r = ir * info_element[it]['dr']
					for ie in range(info_element[it]['Ne']):
						orb[it][il][iu][ir] += C[it][il][ie,iu].item() * spherical_jn(il,E[it][il,ie].item()*r)
	return orb


# read ./ORBITAL_RESULTS.txt (Coefficents) and write ./ORBITAL_6U.dat (.orb) and ./ORBITAL_PLOTU.dat
def create_orb_files(info_element):
    #info_element={'Ne': {'index': 0, 'Nu': [2, 2, 1], 'Nl': 3, 'Rcut': 6, 'dr': 0.01, 'Ecut': 100, 'Ne': 19}}
    
    E=set_E(info_element)
    C, C_read_index = read_C_init( './ORBITAL_RESULTS.txt', info_element )
    orbital.normalize(generate_orbital(info_element,C,E),{it:info_element[it]['dr'] for it in info_element},C, flag_norm_C=True)
    orb = generate_orbital(info_element,C,E)
    
    if True:#info_opt.cal_smooth:
        orbital.smooth_orbital(
            orb,
            {it:info_element[it]['Rcut'] for it in info_element}, {it:info_element[it]['dr'] for it in info_element},
            0.1)
    orbital.orth(
        orb,
        {it:info_element[it]['dr'] for it in info_element})
    
    print_orbital(orb,info_element)
    plot_orbital(
        orb,
        {it:info_element[it]['Rcut'] for it in info_element},
        {it:info_element[it]['dr'] for it in info_element})
    
