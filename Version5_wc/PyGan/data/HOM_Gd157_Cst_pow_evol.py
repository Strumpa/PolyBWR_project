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
from RSEC_NG0 import *
from RSEC_qfiss import *
from RSEC_S2Q import *
# --- HOMOGENEOUS Gd157 burnup evolution
from HOM_Gd157_cpow import *

# --- OTHERS
from getLists import *

########################################################################################################################################################################################
# run options :
exec_D5 = True # True : run DRAGON5 calculations, False : skip DRAGON5 calculations and open compo object with lcm

#
# Library = name of the draglib used to build LIB: object
Libraries = ['J311_295'] # or 'nJ311_295' : new JEFF3.1.1 library with more Gd157 dilutions and finer autolib data 

S2_lib = "PyNjoy2016"#,"oldlib"]
# Serpent2 energy deposition mode : set edep 0 = "Constant energy deposition per fission" --> "PyNjoy2016" and "oldlib" / default setting
# 											 1 = "Local energy deposition based on ENDF MT458 data"  --> "PyNjoy2016" and "oldlib" with _Ecapt_ --> "more" precise
# 											 2 = "Local photon energy deposition" --> "PyNjoy2016" only, using KERMA coefficients, not available in draglib

ssh_modules = ["USS"] #"USS", "AUTO"
ssh_methods = ["RSE"] #"PT", "RSE", "SUBG" all supported for USS: but AUTO: only takes SUBG
# for RSE method : test with eps_RSE = 1.0E-1, 5.0E-2, 1.0E-2, 5.0E-3, 1.0E-3, 1.0E-4
correlations = ["CORR"] # "CORR", "noCORR"

module_to_possible_methods = {"USS":["RSE","SUBG","PT"],"AUTO":["SUBG"]}

saturation_options = [""] # "NSAT", "SAT", "" for short-lived nuclide saturation or not
depl_solutions = ["RUNG"] # "RUNG" : 5th order Runge-Kutta (Cash-Karp) or "KAPS": 4th order Kaps-Rentrop
#
burnup_points = 'BOC_fine_autop5' #'Gd_VBOC_fine2_autop9' #'Gd_BOC_fine', 'Gd_autop4', 'Gd_autop3', 'Gd_VBOC_fine', 'Gd_BOC_fine'
# suffixe = suffixe added to name_geom for creation of figures, MULTICOMPO and BU vector
tracked_nuclides = ["U235","U238","Pu239", 
                    "Gd157","Gd158", "Xe135","Sm149",
					"Sm154", "Eu155", "Eu156", "Eu157"]

D5_case_options_for_multi_S2_comp = {"Library": "J311_295", "ssh_module": "USS", "ssh_method": "RSE", "correlation": "CORR", "sat": "", "depl_sol": "RUNG"}
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

D5_cases = []
# --- Call to DRAGON5 CLE-2000 procedures :

# Nominal case : depletion chain given by draglib, error on Gd157 evolution noticed.
# Hypothesis : energy deposition model difference between DRAGON5 and SERPENT2 is the cause of the difference in Gd157 in constant power evolution
#  				--> normalization of reaction rates over time stages in EVO: with POWER 38.6 [W/gU]
if exec_D5:
	for dlib_name in Libraries:
		for ssh_module in ssh_modules:
			ssh_methods = module_to_possible_methods[ssh_module]
			for ssh_method in ssh_methods:
				for correlation in correlations:
					pyLIB = MixGd157(dlib_name, ssh_method, correlation) # Creation of the microlib, default D5 energy deposition mode
					for sat in saturation_options:
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
else:
	
	name_compo = 'COMPO_HOM_Gd157_J311_295_USS_RSE_CORR_RUNG'
	os.listdir(path)
	pyCOMPO = lcm.new('LCM_INP', name_compo, impx=0)
	print(save_dir_D5)
	D5_case = D5_case(pyCOMPO, "J311_295", burnup_points, "USS", "RSE", "CORR", "", "RUNG", tracked_nuclides, BU_lists, save_dir_D5)
	D5_cases = [D5_case]
# --- Load SERPENT2 results
S2_cases = []

S2_edep0 = S2_case("HOM_UOX_Gd157", S2_lib, edep_id=0, isEcaptSet=False, pcc_id=1, tracked_nuclides = tracked_nuclides, specific_power = 38.6, save_dir = save_dir_S2)
#S2_edep1 = S2_case("HOM_UOX_Gd157", S2_lib, edep_id=1, isEcaptSet=False, pcc_id=1, tracked_nuclides = tracked_nuclides, specific_power = 38.6, save_dir = save_dir_S2)
S2_edep1_Ecapt = S2_case("HOM_UOX_Gd157", S2_lib, edep_id=1, isEcaptSet=True, pcc_id=1, tracked_nuclides = tracked_nuclides, specific_power = 38.6, save_dir = save_dir_S2)
S2_edep2 = S2_case("HOM_UOX_Gd157", S2_lib, edep_id=2, isEcaptSet=False, pcc_id=1, tracked_nuclides = tracked_nuclides, specific_power = 38.6, save_dir = save_dir_S2)


# --- Comparison of DRAGON5 and SERPENT2 results : each comparison = several D5 cases vs 1 Serpent2 case (energy deposition mode)
Comparison_D5_S2_edep0 = multiD5S2("HOM_Gd157_D5_vs_S2_edep_0", D5_cases, S2_edep0, tracked_nuclides, save_dir_comparison)
#Comparison_D5_S2_edep1 = multiD5S2("HOM_Gd157_D5_vs_S2_edep_1", D5_cases, S2_edep1, tracked_nuclides, save_dir_comparison)
Comparison_D5_S2_edep1_Ecapt = multiD5S2("HOM_Gd157_D5_vs_S2_edep_1_Ecapt", D5_cases, S2_edep1_Ecapt, tracked_nuclides, save_dir_comparison)
Comparison_D5_S2_edep2 = multiD5S2("HOM_Gd157_D5_vs_S2_edep_2", D5_cases, S2_edep2, tracked_nuclides, save_dir_comparison)

comparisons_D5_to_S2 = [Comparison_D5_S2_edep0, Comparison_D5_S2_edep1_Ecapt, Comparison_D5_S2_edep2]
for comp in comparisons_D5_to_S2:
	comp.compare_keffs(BU_interp=True)
	comp.compare_Ni(BU_interp=True, abs=False, rel=True)
	comp.compare_keffs(BU_interp=False)
	comp.compare_Ni(BU_interp=False, abs=False, rel=True)


# --- 1 Dragon5 vs all energy deposition modes
# Pick out a D5 case from the list based on the D5_case_options_for_multi_S2_comp
D5_case_to_compare = None
for case in D5_cases:
	if case.draglib_name == D5_case_options_for_multi_S2_comp["Library"] and case.ssh_module == D5_case_options_for_multi_S2_comp["ssh_module"] and case.ssh_method == D5_case_options_for_multi_S2_comp["ssh_method"] and case.correlation == D5_case_options_for_multi_S2_comp["correlation"] and case.sat == D5_case_options_for_multi_S2_comp["sat"] and case.depl_sol == D5_case_options_for_multi_S2_comp["depl_sol"]:
		D5_case_to_compare = case
		break
if D5_case_to_compare is None:
	raise ValueError(f"No D5 case found with the following options : {D5_case_options_for_multi_S2_comp}")

#D5multiS2_case = D5multiS2("HOM_Gd157_D5_cst_pow_evol_vs_S2_all_edep", D5_case_to_compare, [S2_edep0, S2_edep1, S2_edep1_Ecapt, S2_edep2], tracked_nuclides, save_dir_comparison)
D5multiS2_case = D5multiS2("HOM_Gd157_D5_vs_S2_all_edep", D5_case_to_compare, [S2_edep0, S2_edep1_Ecapt, S2_edep2], tracked_nuclides, save_dir_comparison)
D5multiS2_case.compare_keffs(BU_interp=True)
D5multiS2_case.compare_Ni(BU_interp=True, abs=False, rel=True)
D5multiS2_case.compare_keffs(BU_interp=False)
D5multiS2_case.compare_Ni(BU_interp=True, abs=False, rel=True)


# --- Hypothesis to test : Difference (D5-S2) in N_Gd157 is due to differente energy deposition models between DRAGON5 and SERPENT2
# 		--> Try to reproduce energy deposition mode 0 of SERPENT2 with DRAGON5,
# 			--> test RSEC_NG0 : generate microlib for RSE+CORR (hypothesised as best option), but modified depletion chain : Q values = 0 for (n,gamma) reactions
#           --> test RSEC_qfiss : generate microlib for RSE+CORR (hypothesised as best option), but modified depletion chain : fission Q values edited to match S2 edepmode 0 and Q values = 0 for (n,gamma) reactions
		
# --- Call to DRAGON5 CLE-2000 procedures :

NG0_cases = []
QFISS_cases = []
S2Q_fiss_cases = []
for draglib in Libraries:
	# Generate microlib with RSEC_NG0, Q_NG = 0
	pyLIB_NG0 = RSEC_NG0(draglib)
	# Generate microlib with RSEC_qfiss = Q_NG = 0 and Qfiss = Qi/Q_U235 * 202.27 MeV
	pyLIB_qfiss = RSEC_qfiss(draglib)
	# Generate microlib with RSEC_S2Q = Q_NG = 0 and Qfiss = Q-values from Serpent2_oldlib.out
	pyLIB_S2Q = RSEC_S2Q(draglib)
	# loop over BU evolution options
	for sat in saturation_options:
		for sol in depl_solutions:
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
			print(f"State of the calculation : {draglib} {sat} {sol}")
			compo_name_NG0 = f"_COMPO_HOM_Gd157_{draglib}_RSEC_NG0{SAT}_{sol}"
			compo_name_qfiss = f"_COMPO_HOM_Gd157_{draglib}_RSEC_qfiss{SAT}_{sol}"
			compo_name_S2Q = f"_COMPO_HOM_Gd157_{draglib}_RSEC_S2Q{SAT}_{sol}"
			pyCOMPO_NG0 = HOM_Gd157_cpow(f"COMPO",pyLIB_NG0,StepList,compo_name_NG0,"USS",sat,sol)
			pyCOMPO_qfiss = HOM_Gd157_cpow(f"COMPO",pyLIB_qfiss,StepList,compo_name_qfiss,"USS",sat,sol)
			pyCOMPO_S2Q = HOM_Gd157_cpow(f"COMPO",pyLIB_S2Q,StepList,compo_name_S2Q,"USS",sat,sol)
			#
			# --- Post-processing of DRAGON5 results
			dlib_NG0 = f"{draglib}_RSEC_NG0"
			dlib_qfiss = f"{draglib}_RSEC_qfiss"
			dlib_S2Q = f"{draglib}_RSEC_S2Q"
			D5case_NG0 = D5_case(pyCOMPO_NG0, dlib_NG0, burnup_points, "USS", "RSE", "CORR", sat, sol, tracked_nuclides, BU_lists, save_dir_D5)
			D5case_qfiss = D5_case(pyCOMPO_qfiss, dlib_qfiss, burnup_points, "USS", "RSE", "CORR", sat, sol, tracked_nuclides, BU_lists, save_dir_D5)
			D5case_S2Q = D5_case(pyCOMPO_S2Q, dlib_S2Q, burnup_points, "USS", "RSE", "CORR", sat, sol, tracked_nuclides, BU_lists, save_dir_D5)
			#D5case_NG0.plot_keffs()
			#D5case_qfiss.plot_keffs()
			#for iso in tracked_nuclides:
			#	D5case_NG0.plot_Ni(iso)
			#	D5case_qfiss.plot_Ni(iso)
			NG0_cases.append(D5case_NG0)
			QFISS_cases.append(D5case_qfiss)
			S2Q_fiss_cases.append(D5case_S2Q)

# --- Create a comparison between the RSEC_NG0 case and all SERPENT2 edep cases
# Pick out a D5 case from the list based on the D5_case_options_for_multi_S2_comp
D5_case_to_compare = None
for case in NG0_cases:
	if case.draglib_name == dlib_NG0 and case.ssh_module == D5_case_options_for_multi_S2_comp["ssh_module"] and case.ssh_method == D5_case_options_for_multi_S2_comp["ssh_method"] and case.correlation == D5_case_options_for_multi_S2_comp["correlation"] and case.sat == D5_case_options_for_multi_S2_comp["sat"] and case.depl_sol == D5_case_options_for_multi_S2_comp["depl_sol"]:
		D5_case_to_compare = case
		break
if D5_case_to_compare is None:
	raise ValueError(f"No D5 case found with the following options : {D5_case_options_for_multi_S2_comp}")

D5_NG0_to_S2_edepmodes = D5multiS2("HOM_Gd157_D5_RSEC_NG0_vs_S2_all_edep", D5_case_to_compare, [S2_edep0, S2_edep1_Ecapt, S2_edep2], tracked_nuclides, save_dir_comparison)
D5_NG0_to_S2_edepmodes.compare_keffs()
D5_NG0_to_S2_edepmodes.compare_Ni()
D5_NG0_to_S2_edepmodes.plot_delta_Keff()
D5_NG0_to_S2_edepmodes.plot_delta_Ni()

# --- Create a comparison between the RSEC_qfiss case and all SERPENT2 edep cases
# Pick out a D5 case from the list based on the D5_case_options_for_multi_S2_comp
D5_case_to_compare = None
for case in QFISS_cases:
	if case.draglib_name == dlib_qfiss and case.ssh_module == D5_case_options_for_multi_S2_comp["ssh_module"] and case.ssh_method == D5_case_options_for_multi_S2_comp["ssh_method"] and case.correlation == D5_case_options_for_multi_S2_comp["correlation"] and case.sat == D5_case_options_for_multi_S2_comp["sat"] and case.depl_sol == D5_case_options_for_multi_S2_comp["depl_sol"]:
		D5_case_to_compare = case
		break
D5_qfiss_to_S2_edepmodes = D5multiS2("HOM_Gd157_D5_RSEC_qfiss_vs_S2_all_edep", D5_case_to_compare, [S2_edep0, S2_edep1_Ecapt, S2_edep2], tracked_nuclides, save_dir_comparison)
D5_qfiss_to_S2_edepmodes.compare_keffs(BU_interp=True)
D5_qfiss_to_S2_edepmodes.compare_Ni(BU_interp=True, abs=False, rel=True)
D5_qfiss_to_S2_edepmodes.compare_keffs(BU_interp=False)
D5_qfiss_to_S2_edepmodes.compare_Ni(BU_interp=False, abs=False, rel=True)

D5_case_to_compare = None
for case in S2Q_fiss_cases:
	if case.draglib_name == dlib_S2Q and case.ssh_module == D5_case_options_for_multi_S2_comp["ssh_module"] and case.ssh_method == D5_case_options_for_multi_S2_comp["ssh_method"] and case.correlation == D5_case_options_for_multi_S2_comp["correlation"] and case.sat == D5_case_options_for_multi_S2_comp["sat"] and case.depl_sol == D5_case_options_for_multi_S2_comp["depl_sol"]:
		D5_case_to_compare = case
		break

D5_S2Q_to_S2_edepmodes = D5multiS2("HOM_Gd157_D5_RSEC_S2Q_vs_S2_all_edep", D5_case_to_compare, [S2_edep0, S2_edep1_Ecapt, S2_edep2], tracked_nuclides, save_dir_comparison)
D5_S2Q_to_S2_edepmodes.compare_keffs(BU_interp=True)
D5_S2Q_to_S2_edepmodes.compare_Ni(BU_interp=True, abs=False, rel=True)
D5_S2Q_to_S2_edepmodes.compare_keffs(BU_interp=False)
D5_S2Q_to_S2_edepmodes.compare_Ni(BU_interp=False, abs=False, rel=True)

print(f"For burnup point {burnup_points} :")
print(f"BU lists['COMPO'] : {BU_lists['COMPO']}") # should be equal to number of burnup points in evolution (excluding t0)
print(f"BU lists['AUTOP'] : {BU_lists['AUTOP']}")
print(f"BU lists['BU'] : {BU_lists['BU']}") # can be more that number of COMPO points. It is the number of burnup points computed in DRAGON5 (excluding t0)