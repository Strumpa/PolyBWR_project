# PyGan script to post-treat results for HOM_UOX_Gd15* calculations
# Author: R. Guasch
# Date: 2024-11-12
# Purpose: Post-treat results for HOM_UOX_Gd155/157 calculations :
#         - Extract absorption rates and cross sections from COMPO_MESHES multi-compo 
#         - Extract absorption rates and cross sections from Serpent2 reference
#
# In COMP_MESHES, the results are tabulated according to : 
#           - the energy mesh used, SHEM281/295/315 identified in EDIHOM_xxx directories
#           - the self-shielding method used, identified through the 'SSH' parameters : A, R, P, S for Autosecol, RSE, PT and SUBG respectively

# In Serpent2, the results for absorption rates are obtained from detector results for MT=101 and MT=102 reactions of Gd157/U238
# The cross sections are obtained from microdepletion output files
# Results are obtained for PyNjoy2016 and oldlib librairies, only on the SHEM295 energy mesh for now.
#

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import lifo
import lcm
import cle2000
from MeshHandler import energyMeshHandler as ENEMESH

#os.chdir("Gd157_rates_XS_proc")
path = os.getcwd()
print(path)
print(os.listdir())
name_compo='COMPO_MESHES'
os.chdir("Gd157_rates_XS_proc")
print(os.listdir())
pyCOMPO=lcm.new('LCM_INP',name_compo,impx=0)
meshes = ["SHEM281","SHEM295","SHEM315"]
SSH_methods = {"A":"Autoecol","R":"RSE","P":"PT","S":"SUBG"}
for mesh in meshes:
    if mesh == "SHEM281":
        DIR = "EDIHOM_281"
    elif mesh == "SHEM295":
        DIR = "EDIHOM_295"
    elif mesh == "SHEM315":
        DIR = "EDIHOM_315"
    energyMESH = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][0]["ENERGY"]
    print(f"The {mesh} mesh has {len(energyMESH)} energy bounds, so {len(energyMESH)-1} energy groups")
    MESH_obj = ENEMESH(mesh,energyMESH,1.0E+07,"eV")

    SSH_methods_par = pyCOMPO[DIR]["GLOBAL"]['pval00000001']
    #for i in range(len(SSH_methods_par)):
        #print(f"SSH method {i} = {SSH_methods[SSH_methods_par[i]]}")
    isotope = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][1]['ALIAS'][0:5]
    print(f"isotope = {isotope}")

    N_U8 = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESDENS'][0]
    N_Gd157 = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESDENS'][1]

    XS_NGAMMA_U8 = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][0]['NG'] 
    XS_NGAMMA_Gd157 = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][1]['NG'] 

    XS_NABS_Gd157 = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][1]['NA'] 

    PHI_U8 = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][0]['NWT0'] 
    PHI_Gd157 = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][1]['NWT0'] 

    diff_NWT0 = np.array(PHI_U8) - np.array(PHI_Gd157) # Same so we good to go!!
    #print(f"diff_NWT0 = {diff_NWT0}")

    diff_NG_NA = np.array(XS_NGAMMA_Gd157) - np.array(XS_NABS_Gd157)
    #print(f"diff_NG_NA = {diff_NG_NA}") 

    # Try to find the energy mesh in the compo : that would be really useful
    # Will need to use the energyMesh handler class to create Serpent2 ene cards for cross section generation
    
    

# Reconstruct reaction rates
