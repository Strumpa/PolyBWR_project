# PyGan script for burnup evolution of homogeneous fuel cell with evolving Gd157
# Cell : AT10 45Gd homogenized
# Constant POWER normalization of rates over time stages in EVO: -- > keyword POWER 38.6 [W/gU] 
# Author: R. Guasch
# Date: 2025-01-30
# Hypothesis to test : Difference (D5-S2) in N_Gd157 is due to differente energy deposition models between DRAGON5 and SERPENT2
# 

import sys
import numpy as np
import os, shutil
import lifo
import lcm
import cle2000
import matplotlib.pyplot as plt

# --- DRAGON5 procedures
from assertS import *
# POST-PROCESSING class
from postproc_cst_pow_evol import DRAGON_case as D5_case
from postproc_cst_pow_evol import Serpent2_case as S2_case
from postproc_cst_pow_evol import multiD5S2_comparisons as multiD5S2
from postproc_cst_pow_evol import D5multiS2_comparisons as D5multiS2 
#from MULTI_SERP_POSTPROC import *

# --- HOMOGENEOUS Gd157 microlib definition
from MixGd157 import *
from LIB_NG0 import *
# --- HOMOGENEOUS Gd157 burnup evolution
from HOM_Gd157_cpow import *

# --- OTHERS
from getLists import *

########################################################################################################################################################################################
# run options :
exec_D5 = True # True : run DRAGON5 calculations, False : skip DRAGON5 calculations and open compo object with lcm

#
# Library = name of the draglib used to build LIB: object
dlib_name = "endfb8r1_295" # or 'nJ311_295' : new JEFF3.1.1 library with more Gd157 dilutions and finer autolib data 

S2_lib = "endfb8r1_pynjoy2012_kerma"#,"oldlib"]
# Serpent2 energy deposition mode : set edep 0 = "Constant energy deposition per fission" --> "PyNjoy2016" and "oldlib" / default setting
# 											 1 = "Local energy deposition based on ENDF MT458 data"  --> "PyNjoy2016" and "oldlib" with _Ecapt_ --> "more" precise
# 											 2 = "Local photon energy deposition" --> "PyNjoy2016" only, using KERMA coefficients, not available in draglib

ssh_module = "USS" #"USS", "AUTO"
ssh_method = "PT" #"PT", "RSE", "SUBG" all supported for USS: but AUTO: only takes SUBG
# for RSE method : test with eps_RSE = 1.0E-1, 5.0E-2, 1.0E-2, 5.0E-3, 1.0E-3, 1.0E-4
correlation = "CORR" # "CORR", "noCORR"


saturation_option = "" # "NSAT", "SAT", "" for short-lived nuclide saturation or not
depl_solution = "RUNG" # "RUNG" : 5th order Runge-Kutta (Cash-Karp) or "KAPS": 4th order Kaps-Rentrop
#
burnup_points = 'BOC_fine_autop5' #'Gd_VBOC_fine2_autop9' #'Gd_BOC_fine', 'Gd_autop4', 'Gd_autop3', 'Gd_VBOC_fine', 'Gd_BOC_fine'
# suffixe = suffixe added to name_geom for creation of figures, MULTICOMPO and BU vector
tracked_nuclides = ["U235","U238","Pu239", 
                    "Gd157","Gd158", "Xe135","Sm149"]

#
# Creation of results directory
path=os.getcwd()

save_dir_D5 = f"{path}/HOM_Gd157_Cst_pow_evol_results/D5"
save_dir_S2 = f"{path}/HOM_Gd157_Cst_pow_evol_results/S2"
save_dir_comparison = f"{path}/HOM_Gd157_Cst_pow_evol_results/Comparison"
if not os.path.exists(save_dir_D5):
    os.makedirs(save_dir_D5)
if not os.path.exists(save_dir_S2):
    os.makedirs(save_dir_S2)
if not os.path.exists(save_dir_comparison):
    os.makedirs(save_dir_comparison)

# Recovering ListBU ListAUTOP ListCOMPO
[ListeBU,ListeAUTOP,ListeCOMPO]=getLists(burnup_points)
BU_lists = {"BU": ListeBU, "AUTOP": ListeAUTOP, "COMPO": ListeCOMPO}
# Create Steplist for BU - SELFSHIELDING - COMPO save 
StepList = lcm.new('LCM','burnup_steps')
StepList['ListBU']    = np.array(ListeBU, dtype='f')
StepList['ListAutop'] = np.array(ListeAUTOP, dtype='f')
StepList['ListCompo'] = np.array(ListeCOMPO, dtype='f')
StepList.close() # close without erasing

# --- Call to DRAGON5 CLE-2000 procedures :

# Nominal case : depletion chain given by draglib, error on Gd157 evolution noticed.
# Hypothesis : energy deposition model difference between DRAGON5 and SERPENT2 is the cause of the difference in Gd157 in constant power evolution
#  				--> normalization of reaction rates over time stages in EVO: with POWER 38.6 [W/gU]
# --- DRAGON5 microlib generation
pyLIB = MixGd157(dlib_name, ssh_method, correlation) # Creation of the microlib, default D5 energy deposition mode
#
# names for exportation
if saturation_option == "SAT":
	SAT = "_SAT"
elif saturation_option == "SATOFF":
	SAT = "_SATOFF"
elif saturation_option == "NSAT":
	SAT = "_NSAT"
else:
	SAT = ""
print(f"State of the calculation : {dlib_name} {ssh_module} {ssh_method} {correlation} {saturation_option} {depl_solution}")
compo_name = f"_COMPO_HOM_Gd157_{dlib_name}_{ssh_module}_{ssh_method}_{correlation}{SAT}_{depl_solution}"
# run DRAGON5 calculation with BU evolution
pyCOMPO = HOM_Gd157_cpow(f"COMPO",pyLIB,StepList,compo_name,ssh_module,saturation_option,depl_solution)
#
# --- Post-processing of DRAGON5 results
D5case = D5_case(pyCOMPO, dlib_name, burnup_points, ssh_module, ssh_method, correlation, saturation_option, depl_solution, tracked_nuclides, BU_lists, save_dir_D5)
D5case.plot_keffs()
for iso in tracked_nuclides:
	D5case.plot_Ni(iso)


# --- Hypothesis to test : Difference (D5-S2) in N_Gd157 is due to differente energy deposition models between DRAGON5 and SERPENT2
# 		--> Try to reproduce energy deposition mode 0 of SERPENT2 with DRAGON5,
# 			--> test LIB_NG0 : generate microlib for RSE+CORR (hypothesised as best option), but modified depletion chain : Q values = 0 for (n,gamma) reactions
		
# --- Call to DRAGON5 CLE-2000 procedures :

# --- DRAGON5 microlib generation
pyLIB_NG0 = LIB_NG0(dlib_name) # Creation of the microlib, default D5 energy deposition mode
#
# names for exportation
if saturation_option == "SAT":
	SAT = "_SAT"
elif saturation_option == "SATOFF":
	SAT = "_SATOFF"
elif saturation_option == "NSAT":
	SAT = "_NSAT"
else:
	SAT = ""
print(f"State of the calculation NG0 : {dlib_name} {ssh_module} {ssh_method} {correlation} {saturation_option} {depl_solution}")
compo_name = f"_COMPO_HOM_Gd157_{dlib_name}_NG0_{ssh_module}_{ssh_method}_{correlation}{SAT}_{depl_solution}"
# run DRAGON5 calculation with BU evolution
pyCOMPO_NG0 = HOM_Gd157_cpow(f"COMPO",pyLIB_NG0,StepList,compo_name,ssh_module,saturation_option,depl_solution)
#
# --- Post-processing of DRAGON5 results
D5case_NG0 = D5_case(pyCOMPO_NG0, dlib_name, burnup_points, ssh_module, ssh_method, correlation, saturation_option, depl_solution, tracked_nuclides, BU_lists, save_dir_D5)
D5case_NG0.plot_keffs()
for iso in tracked_nuclides:
	D5case_NG0.plot_Ni(iso)
