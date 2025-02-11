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
from postproc_cst_pow_evol import D5S2_comparisons as D5S2_comparisons
#from MULTI_SERP_POSTPROC import *

# --- HOMOGENEOUS Gd157
from MixGd157 import *
from HOM_Gd157_cpow import *

# --- OTHERS
from getLists import *

########################################################################################################################################################################################
#
# Library = name of the draglib used to build LIB: object
Libraries = ['J311_295'] # or 'nJ311_295' : new JEFF3.1.1 library with more Gd157 dilutions and finer autolib data 

S2_lib = "PyNjoy2016"#,"oldlib"]
# Serpent2 energy deposition mode : set edep 0 = "Constant energy deposition per fission" --> "PyNjoy2016" and "oldlib" / default setting
# 											 1 = "Local energy deposition based on ENDF MT458 data"  --> "PyNjoy2016" and "oldlib" with _Ecapt_ --> "more" precise
# 											 2 = "Local photon energy deposition" --> "PyNjoy2016" only, using KERMA coefficients, not available in draglib
set_edep_mode_to_treat = {"PyNjoy2016":[0,1,2],"oldlib":[0,1]}

ssh_modules = ["USS"] #"USS", "AUTO"
ssh_methods = ["RSE"] #"PT", "RSE", "SUBG" all supported for USS: but AUTO: only takes SUBG
# for RSE method : test with eps_RSE = 1.0E-1, 5.0E-2, 1.0E-2, 5.0E-3, 1.0E-3, 1.0E-4
correlations = ["CORR"] # "CORR", "noCORR"

sats = [""] # "NSAT", "SAT", "" for short-lived nuclide saturation or not
depl_solutions = ["RUNG"] # "RUNG" : 5th order Runge-Kutta (Cash-Karp) or "KAPS": 4th order Kaps-Rentrop
#
burnup_points = 'BOC_fine_autop5' #'Gd_VBOC_fine2_autop9' #'Gd_BOC_fine', 'Gd_autop4', 'Gd_autop3', 'Gd_VBOC_fine', 'Gd_BOC_fine'
# suffixe = suffixe added to name_geom for creation of figures, MULTICOMPO and BU vector
tracked_nuclides = ["U235","U238","Pu239","Pu240","Pu241","Pu242", 
                    "Gd154","Gd155","Gd156","Gd157","Gd158","Gd160",
                    "Xe135","Sm149","Sm154","Eu155","Eu156","Eu157"]

D5_case_options_for_multi_S2_comp = {"Library": "J311_295", "ssh_module": "USS", "ssh_method": "RSE", "correlation": "CORR", "sat": "", "depl_sol": "RUNG"}
#
# Creation of results directory
path=os.getcwd()

save_dir_D5 = f"{path}/HOM_Gd157_Cst_pow_evol/D5"
save_dir_S2 = f"{path}/HOM_Gd157_Cst_pow_evol/S2"
save_dir_comparison = f"{path}/HOM_Gd157_Cst_pow_evol/Comparison"
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

D5_cases = []
# --- Call to DRAGON5 CLE-2000 procedures :

for dlib_name in Libraries:
	for ssh_module in ssh_modules:
		if ssh_module == "AUTO":
			ssh_methods = ["SUBG"]
			correlations = ["noCORR"]
		for ssh_method in ssh_methods:
			for correlation in correlations:
				pyLIB = MixGd157(dlib_name, ssh_method, correlation) # Creation of the microlib, default D5 energy deposition mode
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
						print(f"State of the calculation : {dlib_name} {ssh_module} {ssh_method} {correlation} {sat} {depl_sol}")
						compo_name = f"_COMPO_HOM_Gd157_{dlib_name}_{ssh_module}_{ssh_method}_{correlation}{SAT}_{depl_sol}"
						pyCOMPO = HOM_Gd157_cpow(f"COMPO",pyLIB,StepList,compo_name,ssh_module,sat,depl_sol)
						#
						# --- Post-processing of DRAGON5 results
						D5case = D5_case(pyCOMPO, dlib_name, burnup_points, ssh_module, ssh_method, correlation, sat, depl_sol, tracked_nuclides, BU_lists, save_dir_D5)
						D5case.plot_keffs()
						for iso in tracked_nuclides:
							D5case.plot_Ni(iso)
						D5_cases.append(D5case)

# --- Load SERPENT2 results
S2_cases = []

S2_edep0 = S2_case("HOM_UOX_Gd157", S2_lib, edep_id=0, isEcaptSet=False, pcc_id=1, tracked_nuclides = tracked_nuclides, specific_power = 38.6, save_dir = save_dir_S2)
S2_edep1 = S2_case("HOM_UOX_Gd157", S2_lib, edep_id=1, isEcaptSet=True, pcc_id=1, tracked_nuclides = tracked_nuclides, specific_power = 38.6, save_dir = save_dir_S2)
S2_edep2 = S2_case("HOM_UOX_Gd157", S2_lib, edep_id=2, isEcaptSet=False, pcc_id=1, tracked_nuclides = tracked_nuclides, specific_power = 38.6, save_dir = save_dir_S2)

"""
exemple de comment acceder aux attributs des cas S2
for case in [S2_edep0, S2_edep1, S2_edep2]:
	case.plot_keff()
	for iso in tracked_nuclides:
		case.plot_concentrations(iso)
	S2_cases.append(case)
"""

# --- Comparison of DRAGON5 and SERPENT2 results
Comparison_D5_S2_edep0 = D5S2_comparisons("HOM_Gd157_D5_cst_pow_evol_vs_S2_edep_0", D5_cases, S2_edep0, tracked_nuclides, save_dir_comparison)
Comparison_D5_S2_edep1 = D5S2_comparisons("HOM_Gd157_D5_cst_pow_evol_vs_S2_edep_1", D5_cases, S2_edep1, tracked_nuclides, save_dir_comparison)
Comparison_D5_S2_edep2 = D5S2_comparisons("HOM_Gd157_D5_cst_pow_evol_vs_S2_edep_2", D5_cases, S2_edep2, tracked_nuclides, save_dir_comparison)

for comp in [Comparison_D5_S2_edep0, Comparison_D5_S2_edep1, Comparison_D5_S2_edep2]:
	comp.compare_keffs()
	comp.compare_Ni()
	comp.plot_delta_Keff()
	comp.plot_delta_Ni()

# --- Make custom plots : 1 Dragon5 vs all energy deposition modes
# Pick out a D5 case from the list based on the D5_case_options_for_multi_S2_comp
D5_case_to_compare = None
for case in D5_cases:
	if case.draglib_name == D5_case_options_for_multi_S2_comp["Library"] and case.ssh_module == D5_case_options_for_multi_S2_comp["ssh_module"] and case.ssh_method == D5_case_options_for_multi_S2_comp["ssh_method"] and case.correlation == D5_case_options_for_multi_S2_comp["correlation"] and case.sat == D5_case_options_for_multi_S2_comp["sat"] and case.depl_sol == D5_case_options_for_multi_S2_comp["depl_sol"]:
		D5_case_to_compare = case
		break
if D5_case_to_compare is None:
	raise ValueError(f"No D5 case found with the following options : {D5_case_options_for_multi_S2_comp}")

# Initialize dictionnaries for delta values : key = energy deposition mode, values = delta values
delta_keffs_D5_vs_edepmodes = {}
delta_iso_D5_vs_edepmodes = {}
iso_list = ["Gd157", "Gd158", "U238", "U235", "Pu239"]


if D5_case_to_compare is not None:
	for comp in [Comparison_D5_S2_edep0, Comparison_D5_S2_edep1, Comparison_D5_S2_edep2]:
		comparison_identifier = f"{D5_case_to_compare.draglib_name}_{D5_case_to_compare.ssh_module}_{D5_case_to_compare.ssh_method}_{D5_case_to_compare.correlation}_to_S2_edep{comp.S2_case.edep_id}"
		delta_keffs_D5_vs_edepmodes[comp.S2_case.edep_id] = comp.delta_keffs[comparison_identifier]
		delta_iso_D5_vs_edepmodes[comp.S2_case.edep_id] = {}
		for iso in iso_list:
			delta_iso_D5_vs_edepmodes[comp.S2_case.edep_id][iso] = comp.delta_Niso[iso]

	# Plot delta keffs
	fig, ax = plt.subplots()
	ax.set_title(f"$\\Delta$ keffs : Dragon5 vs Serpent2")
	ax.set_xlabel("Burnup [MWd/kgU]")
	ax.set_ylabel("$\\Delta$  keff (pcm)")
	for edep_id, delta_keffs in delta_keffs_D5_vs_edepmodes.items():
		ax.plot(D5_case_to_compare.COMPO_BU_steps, delta_keffs, label=f"D5- Serpent2 energy deposition mode {edep_id}")
	ax.axhline(y = 200.0, color = 'r', linestyle = '-') 
	ax.axhline(y = -200.0, color = 'r', linestyle = '-')  
	ax.grid()
	ax.legend()
	fig.savefig(f"{save_dir_comparison}/delta_keffs_D5_vs_edepmodes.png")
	plt.close(fig)

	
	for iso in iso_list:
		fig, ax = plt.subplots()
		ax.set_title(f"$\\Delta$ {iso} : Dragon5 vs Serpent2")
		ax.set_xlabel("Burnup [MWd/kgU]")
		ax.set_ylabel(f"$\\Delta$  {iso} (%)")
		for edep_id, delta_iso in delta_iso_D5_vs_edepmodes.items():
			ax.plot(D5_case_to_compare.COMPO_BU_steps, delta_iso[iso][f"{D5_case_to_compare.draglib_name}_{D5_case_to_compare.ssh_module}_{D5_case_to_compare.ssh_method}_{D5_case_to_compare.correlation}_to_S2_edep{edep_id}"], label=f"D5- Serpent2 energy deposition mode {edep_id}")
		ax.axhline(y = 2.0, color = 'r', linestyle = '-')
		ax.axhline(y = -2.0, color = 'r', linestyle = '-')
		ax.grid()
		ax.legend()
		
		fig.savefig(f"{save_dir_comparison}/delta_{iso}_D5_vs_edepmodes.png")
		plt.close(fig)




