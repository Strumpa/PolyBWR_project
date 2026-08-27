# PyGan script for burnup evolution of homogeneous fuel cell with evolving Gd157
# Homogenized cell with UOX+(optional Gd) + clad + coolant
# Constant POWER normalization of rates over time stages in EVO: -- > keyword POWER 25.6 [W/gU] 
# Author: R. Guasch
# Date: 2026-05-14
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
from postproc_hfactors import plot_microlib_info, plot_self_shielded_microlib_info
#from MULTI_SERP_POSTPROC import *

# --- HOMOGENEOUS Gd157 microlib definition
from Mix import *
# --- HOMOGENEOUS Gd157 burnup evolution
from HOM_cpow import HOM_cpow_CELI, HOM_cpow_EXTR, flux_295g

# --- OTHERS
from getLists import *

########################################################################################################################################################################################
# run options :
case_name = "HOM_UOX_Gd157" # "HOM_UOX_Gd157"
#
# Library = name of the draglib used to build LIB: object
dlib_name = "e81d1mt18v2" #"endfb8r1_295" #"e81d1mt18v2" #"e81d1mt458" # "endfb8r1_295" # "eb81dec1" or "eb81dec0" or 'nJ311_295' : new JEFF3.1.1 library with more Gd157 dilutions and finer autolib data 

ssh_module = "USS" #"USS", "AUTO"
ssh_method = "PT" #"PT", "RSE", "SUBG" all supported for USS: but AUTO: only takes SUBG
# for RSE method : test with eps_RSE = 1.0E-1, 5.0E-2, 1.0E-2, 5.0E-3, 1.0E-3, 1.0E-4
correlation = "CORR" # "CORR", "noCORR"
time_intergator = "EXTR"

edpmode_flag = "EDP0"
glob_opt = "EDP0" # "GLOB", "NOGL", "EDP0"
#
burnup_points = "Gd_BOC" # "2_steps" #"BOC_fine_autop5" # "VBOC" or "VBOC_finerBU", "BOC_fine_autop5", 'Gd_VBOC_fine2_autop9' #'Gd_BOC_fine', 'Gd_autop4', 'Gd_autop3', 'Gd_VBOC_fine', 'Gd_BOC_fine'
tracked_nuclides = ["U235","U238","Pu239", 
					"Gd157","Gd158", "Xe135","Sm149"]

#
# Creation of results directory
path=os.getcwd()

save_dir_D5 = f"{path}/HOM_TESTS_Cst_pow_evol_results/{case_name}/D5"
save_dir_S2 = f"{path}/HOM_TESTS_Cst_pow_evol_results/{case_name}/S2"
save_dir_comparison = f"{path}/HOM_TESTS_Cst_pow_evol_results/{case_name}/Comparison"
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
#  				--> normalization of reaction rates over time stages in EVO: with POWER 25.6 [W/gU]
# --- DRAGON5 microlib generation
pyLIB = Mix(dlib_name, ssh_method, case_name, edpmode_flag) # Creation of the microlib, default D5 energy deposition mode
#
plot_microlib_info(pyLIB, f"{case_name}_{dlib_name}_{edpmode_flag}", save_dir_D5)
compo_name_295 = f"_CPO_{case_name}_{dlib_name}_{ssh_module}_{ssh_method}_{correlation}_{burnup_points}_{glob_opt}_295"
pyCOMPO_295, pyLIB2 = flux_295g(f"COMPO295",pyLIB,compo_name_295) 
#
plot_self_shielded_microlib_info(pyLIB2, f"{case_name}_{dlib_name}_{edpmode_flag}", save_dir_D5)

# names for exportation
if True:
	print(f"State of the calculation : {dlib_name} {ssh_module} {ssh_method} {correlation} {glob_opt}")
	#compo_name_295 = f"_CPO_{case_name}_{dlib_name}_{ssh_module}_{ssh_method}_{correlation}_{burnup_points}_{glob_opt}_295"
	#pyCOMPO_295 = flux_295g(f"COMPO295",pyLIB,compo_name_295) 

	compo_name = f"_CPO_{case_name}_{dlib_name}_{edpmode_flag}_{ssh_module}_{ssh_method}_{correlation}_{burnup_points}_{glob_opt}_{time_intergator}"
	# run DRAGON5 calculation with BU evolution
	if time_intergator == "CELI":
		pyCOMPO = HOM_cpow_CELI(f"COMPO",pyLIB,StepList,compo_name,glob_opt)
	elif time_intergator == "EXTR":
		pyCOMPO = HOM_cpow_EXTR(f"COMPO",pyLIB,StepList,compo_name,glob_opt)
	#
	D5_NOM_case = D5_case(pyCOMPO = pyCOMPO,
					dlib_name = dlib_name,
					bu_points = burnup_points,
					ssh_opt = "PT",
					correlation = correlation,
					sat = "NODI",
					depl_sol = "RUNG",
					tracked_nuclides = tracked_nuclides,
					BU_lists = BU_lists,
					save_dir = save_dir_D5)
	D5_NOM_case.set_BUscheme("predictor", "LE")

	evaluation = "endfb8r1"
	origin_of_S2_data = "pynjoy2012_kerma"
	S2_ref_case = S2_case(case_name = case_name,
							lib_name = f"{evaluation}_{origin_of_S2_data}",
							edep_id = 0,
							areQfissSet = False,
							isEcaptSet = False,
							pcc_id = 2,
							specific_power = 25.6,
							tracked_nuclides = tracked_nuclides,
							save_dir = save_dir_S2)

	if time_intergator == "CELI":
		S2_ref_case.read_S2_outputs(f"{case_name}_test_endfb8r1_edep0_pcc1_mc")
		pccid = 1
	elif time_intergator == "EXTR":
		S2_ref_case.read_S2_outputs(f"{case_name}_test_endfb8r1_edep0_pcc2_mc")
		pccid = 2

	comparisonD5_EDP0_S2 = multiD5S2(f"D5 {ssh_method}-{correlation} {edpmode_flag} - S2 : {glob_opt} vs edep0 pcc {pccid}", [D5_NOM_case], S2_ref_case, burnup_points, tracked_nuclides, f"{save_dir_comparison}/{dlib_name}/{glob_opt}")

	comparisonD5_EDP0_S2.compare_keffs()
	comparisonD5_EDP0_S2.plot_delta_keff()

	comparisonD5_EDP0_S2.compare_Ni()
	comparisonD5_EDP0_S2.plot_delta_Ni()