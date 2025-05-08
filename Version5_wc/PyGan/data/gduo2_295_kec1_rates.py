# PyGan procedure to post treat gduo2_295_kec1 procedures output
# 1) Post treat BU=0 case : gduo2_295_kec1.x2m Dragon case
# 2) Post treat evolving BU case : gduo2_295_kec1_evol.py PyGan case
# Purpose : (1) extend validation of RSE/RSE CORR, vs AUTO: and vs Serpent2 [Jeff3.1.1]
#          (2) extend validation of RSE/RSE CORR, vs AUTO: and vs Serpent2 [ENDFb8r1]
#          (3) extend validation to evolving BU case --> gduo2_295_kec1_evol.py


import re
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import lifo
import lcm
import cle2000
from MeshHandler import energyMeshHandler as ENEMESH
from PT_D5S2 import postTreatment_rates_XS_D5 as PT_D5


#
# 1) Post treat BU=0 case : gduo2_295_kec1.x2m Dragon case

# Read Dragon output COMPO
path = os.getcwd()
os.chdir('gduo2_295_kec1_rates_proc')
name_compo='COMPO_gduo2_kec1'
pyCOMPO=lcm.new('LCM_INP',name_compo,impx=0)
os.chdir(path)


save_path_D5vsAUTO='gduo2_295_kec1_rates_results/D5vsAUTO'
if not os.path.exists(save_path_D5vsAUTO):
    os.makedirs(save_path_D5vsAUTO)
save_path_D5vsS2='gduo2_295_kec1_rates_results/D5vsS2'
if not os.path.exists(save_path_D5vsS2):
    os.makedirs(save_path_D5vsS2)

Cpo_keys = ["ENDFb8r1_295", "J311_295"]
ssh_keys = ["RSE", "RSE_C_87", "RSE_C_857", "AUTO", "PT", "PT_C_87", "PT_C_857"]

# Parse COMPO : 
# extract absorption rates for Gd157, Gd155 and U238 for RSE, RSE_CORR and AUTO

mesh_objects = {}
Gd157_ngamma = {}
Gd155_ngamma = {}
U238_ngamma =  {}
fluxes = {}

for eval in Cpo_keys:
    Gd157_ngamma[eval] = {}
    Gd155_ngamma[eval] = {}
    U238_ngamma[eval] = {}
    fluxes[eval] = {}
    energyMESH = pyCOMPO[eval]['MIXTURES'][0]['CALCULATIONS'][0]["ENERGY"]
    MESH_obj = ENEMESH("SHEM295",energyMESH,1.0E+07,"eV")
    mesh_objects["SHEM295"] = MESH_obj
    SSH_methods_par = pyCOMPO[eval]["GLOBAL"]['pval00000001'].strip().split()
    print(SSH_methods_par)
    for i in range(len(SSH_methods_par)):
        print(SSH_methods_par[i])
        U238_ngamma[eval][SSH_methods_par[i]] = {"rates":[], "XS":[]}
        Gd155_ngamma[eval][SSH_methods_par[i]] = {"rates":[], "XS":[]}
        Gd157_ngamma[eval][SSH_methods_par[i]] = {"rates":[], "XS":[]}
        fluxes[eval][SSH_methods_par[i]] = []
        print(pyCOMPO[eval]['MIXTURES'][0])
        isotope0 = pyCOMPO[eval]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESLIST'][0]['ALIAS'][0:5]
        isotope1 = pyCOMPO[eval]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESLIST'][1]['ALIAS'][0:5]
        isotope2 = pyCOMPO[eval]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESLIST'][2]['ALIAS'][0:5]

        print(f"Isotopes are {isotope0}, {isotope1}, {isotope2}")

        N_U238 = pyCOMPO[eval]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESDENS'][0]
        N_Gd155 = pyCOMPO[eval]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESDENS'][1]
        N_Gd157 = pyCOMPO[eval]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESDENS'][2]

        XS_NGAMMA_U238 = pyCOMPO[eval]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESLIST'][0]['NG'] 
        XS_NGAMMA_Gd155 = pyCOMPO[eval]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESLIST'][1]['NG'] 
        XS_NGAMMA_Gd157 = pyCOMPO[eval]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESLIST'][2]['NG']

        PHI = pyCOMPO[eval]['MIXTURES'][0]['CALCULATIONS'][i]['ISOTOPESLIST'][0]['NWT0']  # Flux

        # Reconstruct reaction rates
        NGAMMA_U238_rates = np.array(XS_NGAMMA_U238)*np.array(PHI)*N_U238
        NGAMMA_Gd155_rates = np.array(XS_NGAMMA_Gd155)*np.array(PHI)*N_Gd155
        NGAMMA_Gd157_rates = np.array(XS_NGAMMA_Gd157)*np.array(PHI)*N_Gd157

        # Store the results in the dictionaries
        print(U238_ngamma[eval].keys())
        U238_ngamma[eval][SSH_methods_par[i]]["rates"] = NGAMMA_U238_rates
        Gd155_ngamma[eval][SSH_methods_par[i]]["rates"] = NGAMMA_Gd155_rates
        Gd157_ngamma[eval][SSH_methods_par[i]]["rates"] = NGAMMA_Gd157_rates

        U238_ngamma[eval][SSH_methods_par[i]]["XS"] = XS_NGAMMA_U238
        Gd155_ngamma[eval][SSH_methods_par[i]]["XS"] = XS_NGAMMA_Gd155
        Gd157_ngamma[eval][SSH_methods_par[i]]["XS"] = XS_NGAMMA_Gd157

        fluxes[eval][SSH_methods_par[i]] = PHI
# Sanity check :
for eval in Cpo_keys:
    for ssh in ssh_keys:
        print(f"U238_ngamma[{eval}][{ssh}]['rates'] = {U238_ngamma[eval][ssh]['rates']}")
        print(f"Gd155_ngamma[{eval}][{ssh}]['rates'] = {Gd155_ngamma[eval][ssh]['rates']}")
        print(f"Gd157_ngamma[{eval}][{ssh}]['rates'] = {Gd157_ngamma[eval][ssh]['rates']}")
        print(f"fluxes[{eval}][{ssh}] = {fluxes[eval][ssh]}")
        print(f"U238_ngamma[{eval}][{ssh}]['XS'] = {U238_ngamma[eval][ssh]['XS']}")
        print(f"Gd155_ngamma[{eval}][{ssh}]['XS'] = {Gd155_ngamma[eval][ssh]['XS']}")
        print(f"Gd157_ngamma[{eval}][{ssh}]['XS'] = {Gd157_ngamma[eval][ssh]['XS']}")

 # Call post treatment class
    
DRAGON5_BU0 = PT_D5("gduo2_295_kec1", "EVAL", mesh_objects, ["U238_ngamma", "Gd155_ngamma", "Gd157_ngamma"], Cpo_keys, ssh_keys, save_path_D5vsAUTO)
DRAGON5_BU0.set_reaction_data("U238_ngamma",U238_ngamma)
DRAGON5_BU0.set_reaction_data("Gd155_ngamma",Gd155_ngamma)
DRAGON5_BU0.set_reaction_data("Gd157_ngamma",Gd157_ngamma)
DRAGON5_BU0.set_fluxes(fluxes)

# Plot results for Jeff3.1.1 and ENDFb-VIII-1 for RSE, RSE_CORR and AUTO self-shielding methods

for eval in Cpo_keys:
    # Plot fluxes
    DRAGON5_BU0.plot_fluxes(eval, "SHEM295")
    # Plot cross sections
    DRAGON5_BU0.plot_XS(eval, "SHEM295", "U238_ngamma")
    DRAGON5_BU0.plot_XS(eval, "SHEM295", "Gd155_ngamma")
    DRAGON5_BU0.plot_XS(eval, "SHEM295", "Gd157_ngamma")
    # Plot reaction rates
    DRAGON5_BU0.plot_reaction_rates(eval, "SHEM295", "U238_ngamma")
    DRAGON5_BU0.plot_reaction_rates(eval, "SHEM295", "Gd155_ngamma")
    DRAGON5_BU0.plot_reaction_rates(eval, "SHEM295", "Gd157_ngamma")
    # Compute relative differences on cross sections and reaction rates between AUTO: and RSE and RSE_CORR
    # U238 
    DRAGON5_BU0.compute_relative_differences_XS("U238_ngamma", eval, "AUTO", ["RSE", "RSE_C_87", "RSE_C_857", "PT", "PT_C_87", "PT_C_857"])
    DRAGON5_BU0.compute_relative_differences_Rates("U238_ngamma", eval, "AUTO", ["RSE", "RSE_C_87", "RSE_C_857", "PT", "PT_C_87", "PT_C_857"])
    # Gd155
    DRAGON5_BU0.compute_relative_differences_XS("Gd155_ngamma", eval, "AUTO", ["RSE", "RSE_C_87", "RSE_C_857", "PT", "PT_C_87", "PT_C_857"])
    DRAGON5_BU0.compute_relative_differences_Rates("Gd155_ngamma", eval, "AUTO", ["RSE", "RSE_C_87", "RSE_C_857", "PT", "PT_C_87", "PT_C_857"])
    # Gd157
    DRAGON5_BU0.compute_relative_differences_XS("Gd157_ngamma", eval, "AUTO", ["RSE", "RSE_C_87", "RSE_C_857", "PT", "PT_C_87", "PT_C_857"])
    DRAGON5_BU0.compute_relative_differences_Rates("Gd157_ngamma", eval, "AUTO", ["RSE", "RSE_C_87", "RSE_C_857", "PT", "PT_C_87", "PT_C_857"])
    
# Zoom on self-shielded region : RSE and RSE_CORR with CORR for U238 and Gd157
# U238, Jeff3.1.1
DRAGON5_BU0.plot_zoom_XS_and_errors_D5_AUTO("U238_ngamma", "J311_295", "SHEM295", ["RSE", "RSE_C_87"], grmin=70, grmax=206)
DRAGON5_BU0.plot_zoom_rates_and_errors_D5_AUTO("U238_ngamma", "J311_295", "SHEM295", ["RSE", "RSE_C_87"], grmin=70, grmax=206)
# Gd155, Jeff3.1.1
DRAGON5_BU0.plot_zoom_XS_and_errors_D5_AUTO("Gd155_ngamma", "J311_295", "SHEM295", ["RSE", "RSE_C_87"], grmin=70, grmax=206)
DRAGON5_BU0.plot_zoom_rates_and_errors_D5_AUTO("Gd155_ngamma", "J311_295", "SHEM295", ["RSE", "RSE_C_87"], grmin=70, grmax=206)
# Gd157, Jeff3.1.1
DRAGON5_BU0.plot_zoom_XS_and_errors_D5_AUTO("Gd157_ngamma", "J311_295", "SHEM295", ["RSE", "RSE_C_87"], grmin=70, grmax=206)
DRAGON5_BU0.plot_zoom_rates_and_errors_D5_AUTO("Gd157_ngamma", "J311_295", "SHEM295", ["RSE", "RSE_C_87"], grmin=70, grmax=206)

# Zoom on self-shielded region : RSE and RSE_CORR with CORR for U238 and Gd157
# U238, ENDFb-VIII-1
DRAGON5_BU0.plot_zoom_XS_and_errors_D5_AUTO("U238_ngamma", "ENDFb8r1_295", "SHEM295", ["RSE", "RSE_C_87"], grmin=70, grmax=206)
DRAGON5_BU0.plot_zoom_rates_and_errors_D5_AUTO("U238_ngamma", "ENDFb8r1_295", "SHEM295", ["RSE", "RSE_C_87"], grmin=70, grmax=206)
# Gd155, ENDFb-VIII-1
DRAGON5_BU0.plot_zoom_XS_and_errors_D5_AUTO("Gd155_ngamma", "ENDFb8r1_295", "SHEM295", ["RSE", "RSE_C_87"], grmin=70, grmax=206)
DRAGON5_BU0.plot_zoom_rates_and_errors_D5_AUTO("Gd155_ngamma", "ENDFb8r1_295", "SHEM295", ["RSE", "RSE_C_87"], grmin=70, grmax=206)
# Gd157, ENDFb-VIII-1
DRAGON5_BU0.plot_zoom_XS_and_errors_D5_AUTO("Gd157_ngamma", "ENDFb8r1_295", "SHEM295", ["RSE", "RSE_C_87"], grmin=70, grmax=206)
DRAGON5_BU0.plot_zoom_rates_and_errors_D5_AUTO("Gd157_ngamma", "ENDFb8r1_295", "SHEM295", ["RSE", "RSE_C_87"], grmin=70, grmax=206)

# Zoom on self-shielded region : RSE and RSE_CORR with CORR for U238, Gd155 and Gd157
# U238, ENDFb-VIII-1
DRAGON5_BU0.plot_zoom_XS_and_errors_D5_AUTO("U238_ngamma", "ENDFb8r1_295", "SHEM295", ["RSE", "RSE_C_857"], grmin=70, grmax=206)
DRAGON5_BU0.plot_zoom_rates_and_errors_D5_AUTO("U238_ngamma", "ENDFb8r1_295", "SHEM295", ["RSE", "RSE_C_857"], grmin=70, grmax=206)
# Gd155, ENDFb-VIII-1
DRAGON5_BU0.plot_zoom_XS_and_errors_D5_AUTO("Gd155_ngamma", "ENDFb8r1_295", "SHEM295", ["RSE", "RSE_C_857"], grmin=70, grmax=206)
DRAGON5_BU0.plot_zoom_rates_and_errors_D5_AUTO("Gd155_ngamma", "ENDFb8r1_295", "SHEM295", ["RSE", "RSE_C_857"], grmin=70, grmax=206)
# Gd157, ENDFb-VIII-1
DRAGON5_BU0.plot_zoom_XS_and_errors_D5_AUTO("Gd157_ngamma", "ENDFb8r1_295", "SHEM295", ["RSE", "RSE_C_857"], grmin=70, grmax=206)
DRAGON5_BU0.plot_zoom_rates_and_errors_D5_AUTO("Gd157_ngamma", "ENDFb8r1_295", "SHEM295", ["RSE", "RSE_C_857"], grmin=70, grmax=206)

# Zoom on self-shielded region : PT and PT_CORR with CORR for U238 and Gd157
# U238, ENDFb-VIII-1
DRAGON5_BU0.plot_zoom_XS_and_errors_D5_AUTO("U238_ngamma", "ENDFb8r1_295", "SHEM295", ["PT", "PT_C_87"], grmin=70, grmax=206)
DRAGON5_BU0.plot_zoom_rates_and_errors_D5_AUTO("U238_ngamma", "ENDFb8r1_295", "SHEM295", ["PT", "PT_C_87"], grmin=70, grmax=206)
# Gd155, ENDFb-VIII-1
DRAGON5_BU0.plot_zoom_XS_and_errors_D5_AUTO("Gd155_ngamma", "ENDFb8r1_295", "SHEM295", ["PT", "PT_C_87"], grmin=70, grmax=206)
DRAGON5_BU0.plot_zoom_rates_and_errors_D5_AUTO("Gd155_ngamma", "ENDFb8r1_295", "SHEM295", ["PT", "PT_C_87"], grmin=70, grmax=206)
# Gd157, ENDFb-VIII-1
DRAGON5_BU0.plot_zoom_XS_and_errors_D5_AUTO("Gd157_ngamma", "ENDFb8r1_295", "SHEM295", ["PT", "PT_C_87"], grmin=70, grmax=206)
DRAGON5_BU0.plot_zoom_rates_and_errors_D5_AUTO("Gd157_ngamma", "ENDFb8r1_295", "SHEM295", ["PT", "PT_C_87"], grmin=70, grmax=206)

# Zoom on self-shielded region : PT and PT_CORR with CORR for U238, Gd155 and Gd157
# U238, ENDFb-VIII-1
DRAGON5_BU0.plot_zoom_XS_and_errors_D5_AUTO("U238_ngamma", "ENDFb8r1_295", "SHEM295", ["PT", "PT_C_857"], grmin=70, grmax=206)
DRAGON5_BU0.plot_zoom_rates_and_errors_D5_AUTO("U238_ngamma", "ENDFb8r1_295", "SHEM295", ["PT", "PT_C_857"], grmin=70, grmax=206)
# Gd155, ENDFb-VIII-1
DRAGON5_BU0.plot_zoom_XS_and_errors_D5_AUTO("Gd155_ngamma", "ENDFb8r1_295", "SHEM295", ["PT", "PT_C_857"], grmin=70, grmax=206)
DRAGON5_BU0.plot_zoom_rates_and_errors_D5_AUTO("Gd155_ngamma", "ENDFb8r1_295", "SHEM295", ["PT", "PT_C_857"], grmin=70, grmax=206)
# Gd157, ENDFb-VIII-1
DRAGON5_BU0.plot_zoom_XS_and_errors_D5_AUTO("Gd157_ngamma", "ENDFb8r1_295", "SHEM295", ["PT", "PT_C_857"], grmin=70, grmax=206)
DRAGON5_BU0.plot_zoom_rates_and_errors_D5_AUTO("Gd157_ngamma", "ENDFb8r1_295", "SHEM295", ["PT", "PT_C_857"], grmin=70, grmax=206)

print("Post treatment for BU=0 case completed")