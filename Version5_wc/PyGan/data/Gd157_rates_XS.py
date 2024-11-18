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
from PT_Gd157 import postTreatment_rates_XS_D5 as PT_D5
from PT_Gd157 import postTreatment_rates_XS_S2 as PT_S2

#os.chdir("Gd157_rates_XS_proc")
path = os.getcwd()
name_compo='COMPO_MESHES'
os.chdir("Gd157_rates_XS_proc")
pyCOMPO=lcm.new('LCM_INP',name_compo,impx=0)
os.chdir(path)
# Creation of results directory
path=os.getcwd()
save_path = f"Gd157_Rates_and_XS_results_PyGan"
a=os.path.exists(save_path)
if a==False:
	os.mkdir(save_path)
case_name = "HOM_UOX_Gd157"

meshes = ["SHEM281","SHEM295","SHEM315"]
SSH_methods = {"A":"Autosecol","R":"RSE","P":"PT","S":"SUBG"}
Gd157_ngamma = {"SHEM281":{"Autosecol":{"XS":[], "rates":[]}, "RSE":{"XS":[], "rates":[]}, "PT":{"XS":[], "rates":[]}, "SUBG":{"XS":[], "rates":[]}},
                "SHEM295":{"Autosecol":{"XS":[], "rates":[]}, "RSE":{"XS":[], "rates":[]}, "PT":{"XS":[], "rates":[]}, "SUBG":{"XS":[], "rates":[]}},
                "SHEM315":{"Autosecol":{"XS":[], "rates":[]}, "RSE":{"XS":[], "rates":[]}, "PT":{"XS":[], "rates":[]}, "SUBG":{"XS":[], "rates":[]}}} # build a dictionary of dictionaries to store the rates for each mesh and SSH method
Gd157_abs =  {"SHEM281":{"Autosecol":{"XS":[], "rates":[]}, "RSE":{"XS":[], "rates":[]}, "PT":{"XS":[], "rates":[]}, "SUBG":{"XS":[], "rates":[]}},
                "SHEM295":{"Autosecol":{"XS":[], "rates":[]}, "RSE":{"XS":[], "rates":[]}, "PT":{"XS":[], "rates":[]}, "SUBG":{"XS":[], "rates":[]}},
                "SHEM315":{"Autosecol":{"XS":[], "rates":[]}, "RSE":{"XS":[], "rates":[]}, "PT":{"XS":[], "rates":[]}, "SUBG":{"XS":[], "rates":[]}}}
U8_ngamma = {"SHEM281":{"Autosecol":{"XS":[], "rates":[]}, "RSE":{"XS":[], "rates":[]}, "PT":{"XS":[], "rates":[]}, "SUBG":{"XS":[], "rates":[]}},
                "SHEM295":{"Autosecol":{"XS":[], "rates":[]}, "RSE":{"XS":[], "rates":[]}, "PT":{"XS":[], "rates":[]}, "SUBG":{"XS":[], "rates":[]}},
                "SHEM315":{"Autosecol":{"XS":[], "rates":[]}, "RSE":{"XS":[], "rates":[]}, "PT":{"XS":[], "rates":[]}, "SUBG":{"XS":[], "rates":[]}}}
mesh_objects = {"SHEM281":None,"SHEM295":None,"SHEM315":None}

fluxes = {"SHEM281":{"Autosecol":[], "RSE":[], "PT":[], "SUBG":[]},"SHEM295":{"Autosecol":[], "RSE":[], "PT":[], "SUBG":[]},"SHEM315":{"Autosecol":[], "RSE":[], "PT":[], "SUBG":[]}}

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
    MESH_obj.printnfgCard()
    mesh_objects[mesh] = MESH_obj

    SSH_methods_par = pyCOMPO[DIR]["GLOBAL"]['pval00000001'].strip().split("           ")
    print(SSH_methods_par)
    print(pyCOMPO[DIR]["GLOBAL"]['pval00000001'])
    for i in range(len(SSH_methods_par)):
        print(i)
        print(f"SSH method {i} = {SSH_methods[SSH_methods_par[i]]}")
        isotope = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESLIST'][1]['ALIAS'][0:5]
        print(f"isotope = {isotope}")
        print(f"number of calculations = {len(pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'])}")

        N_U8 = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESDENS'][0]
        N_Gd157 = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESDENS'][1]

        XS_NGAMMA_U8 = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESLIST'][0]['NG'] 
        XS_NGAMMA_Gd157 = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESLIST'][1]['NG'] 


        PHI = pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESLIST'][0]['NWT0']  

        # Try to find the energy mesh in the compo : that would be really useful
        # Will need to use the energyMesh handler class to create Serpent2 ene cards for cross section generation
        
        

        # Reconstruct reaction rates
        NGAMMA_Gd157_rates = np.array(XS_NGAMMA_Gd157)*np.array(PHI)*N_Gd157
        NGAMMA_U8_rates = np.array(XS_NGAMMA_U8)*np.array(PHI)*N_U8

        # Store the results in the dictionaries

        Gd157_ngamma[mesh][SSH_methods[SSH_methods_par[i]]]["rates"] = NGAMMA_Gd157_rates
        U8_ngamma[mesh][SSH_methods[SSH_methods_par[i]]]["rates"] = NGAMMA_U8_rates

        Gd157_ngamma[mesh][SSH_methods[SSH_methods_par[i]]]["XS"] = XS_NGAMMA_Gd157
        U8_ngamma[mesh][SSH_methods[SSH_methods_par[i]]]["XS"] = XS_NGAMMA_U8

        fluxes[mesh][SSH_methods[SSH_methods_par[i]]] = PHI

# Visualize the ditionaries contents
print(f"Gd157_ngamma = {Gd157_ngamma}")
print(f"U8_ngamma = {U8_ngamma}")

# Call post treatment class

DRAGON5_case = PT_D5(case_name, mesh_objects, SSH_methods.values(), save_path)
DRAGON5_case.set_reaction_data("Gd157_ngamma",Gd157_ngamma)
DRAGON5_case.set_reaction_data("U8_ngamma",U8_ngamma)
DRAGON5_case.set_fluxes(fluxes)

# Plot the results

for mesh in meshes:
    DRAGON5_case.plot_fluxes(mesh)
    DRAGON5_case.plot_XS(mesh,"Gd157_ngamma")
    DRAGON5_case.plot_XS(mesh,"U8_ngamma")
    DRAGON5_case.plot_reaction_rates(mesh,"Gd157_ngamma")
    DRAGON5_case.plot_reaction_rates(mesh,"U8_ngamma")
    if mesh == "SHEM281":
        grmin = 1
        grmax = 281
    elif mesh == "SHEM295":
        grmin = 1
        grmax = 295
    elif mesh == "SHEM315":
        grmin = 1
        grmax = 315

    DRAGON5_case.compute_relative_differences_XS("Gd157_ngamma",mesh, "Autosecol", ["RSE","PT","SUBG"])
    DRAGON5_case.compute_relative_differences_XS("U8_ngamma",mesh, "Autosecol", ["RSE","PT","SUBG"])
    DRAGON5_case.compute_relative_differences_Rates("Gd157_ngamma",mesh, "Autosecol", ["RSE","PT","SUBG"])
    DRAGON5_case.compute_relative_differences_Rates("U8_ngamma",mesh, "Autosecol", ["RSE","PT","SUBG"])

    DRAGON5_case.plot_histogram_relative_differences_XS("Gd157_ngamma",mesh, ["RSE","PT","SUBG"],grmin,grmax)
    DRAGON5_case.plot_histogram_relative_differences_XS("U8_ngamma",mesh, ["RSE","PT","SUBG"],grmin,grmax)
    DRAGON5_case.plot_histogram_relative_differences_Rates("Gd157_ngamma",mesh, ["RSE","PT","SUBG"],grmin,grmax)
    DRAGON5_case.plot_histogram_relative_differences_Rates("U8_ngamma",mesh, ["RSE","PT","SUBG"],grmin,grmax)

    DRAGON5_case.plot_histogram_relative_differences_XS("Gd157_ngamma",mesh, ["RSE"],grmin,grmax)
    DRAGON5_case.plot_histogram_relative_differences_XS("U8_ngamma",mesh, ["RSE"],grmin,grmax)
    DRAGON5_case.plot_histogram_relative_differences_Rates("Gd157_ngamma",mesh, ["RSE"],grmin,grmax)
    DRAGON5_case.plot_histogram_relative_differences_Rates("U8_ngamma",mesh, ["RSE"],grmin,grmax)

    DRAGON5_case.plot_histogram_relative_differences_XS("Gd157_ngamma",mesh, ["PT"],grmin,grmax)
    DRAGON5_case.plot_histogram_relative_differences_XS("U8_ngamma",mesh, ["PT"],grmin,grmax)
    DRAGON5_case.plot_histogram_relative_differences_Rates("Gd157_ngamma",mesh, ["PT"],grmin,grmax)
    DRAGON5_case.plot_histogram_relative_differences_Rates("U8_ngamma",mesh, ["PT"],grmin,grmax)

    DRAGON5_case.plot_histogram_relative_differences_XS("Gd157_ngamma",mesh, ["SUBG"],grmin,grmax)
    DRAGON5_case.plot_histogram_relative_differences_XS("U8_ngamma",mesh, ["SUBG"],grmin,grmax)
    DRAGON5_case.plot_histogram_relative_differences_Rates("Gd157_ngamma",mesh, ["SUBG"],grmin,grmax)
    DRAGON5_case.plot_histogram_relative_differences_Rates("U8_ngamma",mesh, ["SUBG"],grmin,grmax)


SERPENT2_case = PT_S2(case_name, mesh_objects, ["PyNjoy2016", "oldlib"], range(0,70))
SERPENT2_case.parse_S2_outputs()


print("Post-treatment done")

