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
import matplotlib

# --- DRAGON5 procedures
from assertS import *
# POST-PROCESSING class
#from POSTPROC_hom import *
#from MULTI_SERP_POSTPROC import *

# --- HOMOGENEOUS Gd157

from HOM_Gd157_cpow import *

# --- OTHERS
from getLists import *

########################################################################################################################################################################################
#
# Library = name of the draglib used to build LIB: object
Libraries = ['J311_295', 'nJ311_295'] # or 'nJ311_295' : new JEFF3.1.1 library with more Gd157 dilutions and finer autolib data 

S2_libs = ["PyNjoy2016"]#,"oldlib"]
# Serpent2 energy deposition mode : set edep 0 = "Constant energy deposition per fission" --> "PyNjoy2016" and "oldlib" / default setting
# 											 1 = "Local energy deposition based on ENDF MT458 data"  --> "PyNjoy2016" and "oldlib" with _Ecapt_ --> "more" precise
# 											 2 = "Local photon energy deposition" --> "PyNjoy2016" only, using KERMA coefficients, not available in draglib
set_edep_mode_to_treat = {"PyNjoy2016":[0,1,2],"oldlib":[0,1]}

ssh_modules = ["USS","AUTO"] #"USS", "AUTO"
ssh_methods = ["RSE","PT","SUBG"] #"PT", "RSE", "SUBG" all supported for USS: but AUTO: only takes SUBG
# for RSE method : test with eps_RSE = 1.0E-1, 5.0E-2, 1.0E-2, 5.0E-3, 1.0E-3, 1.0E-4
correlations = ["CORR","noCORR"] # "CORR", "noCORR"

sats = ["SAT","NSAT",""] # "NSAT", "SAT", "" for short-lived nuclide saturation or not
depl_solutions = ["KAPS","RUNG"] # "RUNG" : 5th order Runge-Kutta (Cash-Karp) or "KAPS": 4th order Kaps-Rentrop

#
# burnup_points = string to use predefine burnup points matching with SERPENT2 results

burnup_points = 'BOC_fine_autop5' #'Gd_VBOC_fine2_autop9' #'Gd_BOC_fine', 'Gd_autop4', 'Gd_autop3', 'Gd_VBOC_fine', 'Gd_BOC_fine'
# suffixe = suffixe added to name_geom for creation of figures, MULTICOMPO and BU vector
#
# Creation of results directory
path=os.getcwd()

a=os.path.exists(f"HOM_Gd157_Cst_pow_evol")
if a==False:
	os.mkdir(f"HOM_Gd157_Cst_pow_evol")
#
# Recovering ListBU ListAUTOP ListCOMPO
[ListeBU,ListeAUTOP,ListeCOMPO]=getLists(burnup_points)
#
# Create Steplist for BU - SELFSHIELDING - COMPO save 
StepList = lcm.new('LCM','burnup_steps')
StepList['ListBU']    = np.array(ListeBU, dtype='f')
StepList['ListAutop'] = np.array(ListeAUTOP, dtype='f')
StepList['ListCompo'] = np.array(ListeCOMPO, dtype='f')
StepList.close() # close without erasing

pyCOMPOS_dict = {}
# --- Call to DRAGON5 CLE-2000 procedures :

for lib in Libraries:
	for ssh_module in ssh_modules:
		if ssh_module == "AUTO":
			ssh_methods = ["SUBG"]
			correlations = ["noCORR"]
		else:
			ssh_methods = ["RSE","PT","SUBG"]
			correlations = ["CORR","noCORR"]
		for ssh_method in ssh_methods:
			for correlation in correlations:
				for sat in sats:
					for depl_sol in depl_solutions:
						#
						# names for exportation
						if sat == "SAT":
							SAT = "_SAT"
						elif sat == "SATOFF":
							SAT = "_SATOFF"
						elif sat == "NSAT":
							SAT = "_NSAT"
						else:
							SAT = ""
						print(f"State of the calculation : {lib} {ssh_module} {ssh_method} {correlation} {sat} {depl_sol}")
						suffixe = f"{lib}_{ssh_module}_{ssh_method}_{correlation}{SAT}_{depl_sol}"
						compo_name = f"_COMPO_HOM_Gd157_{lib}_{ssh_module}_{ssh_method}_{correlation}{SAT}_{depl_sol}"
						#
						pyCOMPOS_dict[suffixe] = HOM_Gd157_cpow(f"COMPO",StepList,compo_name,ssh_module,ssh_method,sat,depl_sol)
							


