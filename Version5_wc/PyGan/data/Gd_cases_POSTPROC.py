## Gd_cases_POSTPROC.py
# Purpose: Post-processing script for Gd cases :
#          - HOM_UOX_Gd157
#          - AT10_45Gd
#          - gduo2_295_kec1
# Aim : study the impact of the energy deposition mode on depleting system with Gd(155/157) 

# Importing libraries
import re
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import lifo
import lcm
import cle2000


from postproc_cst_pow_evol import DRAGON_case as D5_case
from postproc_cst_pow_evol import Serpent2_case as S2_case
from postproc_cst_pow_evol import multiD5S2_comparisons as multiD5S2
from postproc_cst_pow_evol import D5multiS2_comparisons as D5multiS2 
from getLists import getLists


# evaluation and origin of S2 data

evaluation = "endfb8r1"
origin_of_S2_data = "pynjoy2012_kerma"

# isotopes to be tracked 
tracked_nuclides = ["U235","U238","Pu239","Pu240","Pu241","Pu242","Am241","Xe135","Sm149","Gd155","Gd157"]

# build the path to the results directory
path=os.getcwd()
save_dir_HOM_Gd157_VBOC = f"{path}/Gd_cases_POSTPROC/HOM_Gd157_VBOC"
save_dir_HOM_UOX_Gd157 = f"{path}/Gd_cases_POSTPROC/HOM_UOX_Gd157"
print(f"save_dir_HOM_UOX_Gd157 : {save_dir_HOM_UOX_Gd157}")
save_dir_AT10_45Gd = f"{path}/Gd_cases_POSTPROC/AT10_45Gd"
save_dir_gduo2_295_kec1 = f"{path}/Gd_cases_POSTPROC/gduo2_295_kec1"
# create the directories if they do not exist
if not os.path.exists(save_dir_HOM_Gd157_VBOC):
    os.makedirs(save_dir_HOM_Gd157_VBOC)
if not os.path.exists(save_dir_HOM_UOX_Gd157):
    os.makedirs(save_dir_HOM_UOX_Gd157)
if not os.path.exists(save_dir_AT10_45Gd):
    os.makedirs(save_dir_AT10_45Gd)
if not os.path.exists(save_dir_gduo2_295_kec1):
    os.makedirs(save_dir_gduo2_295_kec1)


### Begin post treatment of DRAGON and SERPENT2 results
# Case 0 : HOM_Gd157_VBOC_OMC (already done in HOM_Gd157_VBOC.py from DRAGON5 results (not called from PyGan)) <-- not a priority
# Case 1 : HOM_Gd157_VBOC, focus on this
# Case 2 : HOM_UOX_Gd157, focus on this
# Case 3 : AT10_45Gd, focus on this
# Case 4 : gduo2_295_kec1

# Case 1 : HOM_Gd157_VBOC
    """
    _COMPO_HOM_Gd157_VBOC_endfb8r1_295_PCC0
    _COMPO_HOM_Gd157_VBOC_endfb8r1_295_NG0_PCC0
    """

# Load COMPOS from DRAGON_RESULTS
    
path_to_DRAGON_results = f"{os.getcwd()}/DRAGON_COMPOS_path/"
cwd_path = os.getcwd()
os.chdir(path_to_DRAGON_results)
# Load the data
name_CPO_HOM_Gd157_VBOC_NG0 = f"COMPO_HOM_Gd157_VBOC_{evaluation}_295_NG0"
name_CPO_HOM_Gd157_VBOC_NG0_PCC0 = f"COMPO_HOM_Gd157_VBOC_{evaluation}_295_NG0_PCC0"
name_CPO_HOM_Gd157_VBOC_NG0_DIRA = f"COMPO_HOM_Gd157_VBOC_{evaluation}_295_NG0_DIRA"
name_CPO_HOM_Gd157_VBOC_NG0_PCC0_DIRA = f"COMPO_HOM_Gd157_VBOC_{evaluation}_295_NG0_PCC0_DIRA"
name_CPO_HOM_Gd157_VBOC_NG0_PCC0_SAT = f"COMPO_HOM_Gd157_VBOC_{evaluation}_295_NG0_PCC0_SAT"
name_CPO_HOM_Gd157_VBOC_NG0_SAT = f"COMPO_HOM_Gd157_VBOC_{evaluation}_295_NG0_SAT"

CPO_HOM_Gd157_VBOC_NG0 = lcm.new('LCM_INP', name_CPO_HOM_Gd157_VBOC_NG0, impx=0)
CPO_HOM_Gd157_VBOC_NG0_PCC0 = lcm.new('LCM_INP', name_CPO_HOM_Gd157_VBOC_NG0_PCC0, impx=0)
CPO_HOM_Gd157_VBOC_NG0_DIRA = lcm.new('LCM_INP', name_CPO_HOM_Gd157_VBOC_NG0_DIRA, impx=0)
CPO_HOM_Gd157_VBOC_NG0_PCC0_DIRA = lcm.new('LCM_INP', name_CPO_HOM_Gd157_VBOC_NG0_PCC0_DIRA, impx=0)
CPO_HOM_Gd157_VBOC_NG0_PCC0_SAT = lcm.new('LCM_INP', name_CPO_HOM_Gd157_VBOC_NG0_PCC0_SAT, impx=0)
CPO_HOM_Gd157_VBOC_NG0_SAT = lcm.new('LCM_INP', name_CPO_HOM_Gd157_VBOC_NG0_SAT, impx=0)

os.chdir(cwd_path)    

## Load S2 results with set fission Q-values, edepmode 0, pcc 0

# reminder of S2_case object contrustructor (case_name, lib_name, edep_id, areQfissSet, isEcaptSet, pcc_id, specific_power, tracked_nuclides, save_dir)
S2_edep0_setQfiss_pcc0 = S2_case(case_name = "HOM_Gd157_VBOC",
                                lib_name = "endfb8r1_pynjoy2012_kerma",
                                edep_id = 0, areQfissSet = True, isEcaptSet = False,
                                pcc_id = 0, specific_power = 38.6, tracked_nuclides = tracked_nuclides, save_dir = save_dir_HOM_Gd157_VBOC)
S2_edep0_setQfiss_pcc1 = S2_case(case_name = "HOM_Gd157_VBOC",
                                lib_name = "endfb8r1_pynjoy2012_kerma",
                                edep_id = 0, areQfissSet = True, isEcaptSet = False,
                                pcc_id = 1, specific_power = 38.6, tracked_nuclides = tracked_nuclides, save_dir = save_dir_HOM_Gd157_VBOC)
S2_edep0_setQfiss_pcc2 = S2_case(case_name = "HOM_Gd157_VBOC",
                                lib_name = "endfb8r1_pynjoy2012_kerma",
                                edep_id = 0, areQfissSet = True, isEcaptSet = False,
                                pcc_id = 2, specific_power = 38.6, tracked_nuclides = tracked_nuclides, save_dir = save_dir_HOM_Gd157_VBOC)

# Create D5 cases with and without (n,gamma) energy deposition
# reminder of D5_case object contrustructor (pyCOMPO, dlib_name, bu_points, ssh_opt, correlation, sat, depl_sol, tracked_nuclides, BU_lists, save_dir)


D5_HOM_Gd157_VBOC_NG0 = D5_case(pyCOMPO = CPO_HOM_Gd157_VBOC_NG0,
                            dlib_name = "endfb8r1_295_NG0",
                            bu_points = "VBOC",
                            ssh_opt = "PT",
                            correlation = "CORR",
                            sat = "",
                            depl_sol = "RUNG",
                            tracked_nuclides = tracked_nuclides,
                            BU_lists = getLists("VBOC"),
                            save_dir = save_dir_HOM_Gd157_VBOC)

D5_HOM_Gd157_VBOC_NG0_PCC0 = D5_case(pyCOMPO = CPO_HOM_Gd157_VBOC_NG0_PCC0,
                            dlib_name = "endfb8r1_295_NG0",
                            bu_points = "VBOC",
                            ssh_opt = "PT",
                            correlation = "CORR",
                            sat = "",
                            depl_sol = "RUNG",
                            tracked_nuclides = tracked_nuclides,
                            BU_lists = getLists("VBOC"),
                            save_dir = save_dir_HOM_Gd157_VBOC)

D5_HOM_Gd157_VBOC_NG0_DIRA = D5_case(pyCOMPO = CPO_HOM_Gd157_VBOC_NG0_DIRA,
                            dlib_name = "endfb8r1_295_NG0",
                            bu_points = "VBOC",
                            ssh_opt = "PT",
                            correlation = "CORR",
                            sat = "DIRA",
                            depl_sol = "RUNG",
                            tracked_nuclides = tracked_nuclides,
                            BU_lists = getLists("VBOC"),
                            save_dir = save_dir_HOM_Gd157_VBOC)

D5_HOM_Gd157_VBOC_NG0_PCC0_DIRA = D5_case(pyCOMPO = CPO_HOM_Gd157_VBOC_NG0_PCC0_DIRA,
                            dlib_name = "endfb8r1_295_NG0",
                            bu_points = "VBOC",
                            ssh_opt = "PT",
                            correlation = "CORR",
                            sat = "DIRA",
                            depl_sol = "RUNG",
                            tracked_nuclides = tracked_nuclides,
                            BU_lists = getLists("VBOC"),
                            save_dir = save_dir_HOM_Gd157_VBOC)

D5_HOM_Gd157_VBOC_NG0_PCC0_SAT = D5_case(pyCOMPO = CPO_HOM_Gd157_VBOC_NG0_PCC0_SAT,
                            dlib_name = "endfb8r1_295_NG0",
                            bu_points = "VBOC",
                            ssh_opt = "PT",
                            correlation = "CORR",
                            sat = "SAT",
                            depl_sol = "RUNG",
                            tracked_nuclides = tracked_nuclides,
                            BU_lists = getLists("VBOC"),
                            save_dir = save_dir_HOM_Gd157_VBOC)
D5_HOM_Gd157_VBOC_NG0_SAT = D5_case(pyCOMPO = CPO_HOM_Gd157_VBOC_NG0_SAT,
                            dlib_name = "endfb8r1_295_NG0",
                            bu_points = "VBOC",
                            ssh_opt = "PT",
                            correlation = "CORR",
                            sat = "SAT",
                            depl_sol = "RUNG",
                            tracked_nuclides = tracked_nuclides,
                            BU_lists = getLists("VBOC"),
                            save_dir = save_dir_HOM_Gd157_VBOC)

# Compare the results of DRAGON5 and SERPENT2
BU_points_plot = getLists("VBOC")["BU"]
# compare NG0 NOEX + NODI with S2 edepmode 0 pcc0
delta_keff_NG0_NOEX_HOM_Gd157_VBOC_pcc0 = (D5_HOM_Gd157_VBOC_NG0_PCC0.keff - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm <-- This one works, use as a reference / building bloc for other cases
# compare NG0 NOEX + DIRA with S2 edepmode 0 pcc0
delta_keff_NG0_NOEX_DIRA_HOM_Gd157_VBOC_pcc0 = (D5_HOM_Gd157_VBOC_NG0_PCC0_DIRA.keff - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm

# compare NG0 EXTR + NODI with S2 edepmode 0 pcc2
delta_keff_EXTR_NG0_HOM_Gd157_VBOC_pcc2 = (D5_HOM_Gd157_VBOC_NG0.keff - S2_edep0_setQfiss_pcc2.keff) * 1e5 # pcm
# compare NG0 EXTR + DIRA with S2 edepmode 0 pcc2
delta_keff_EXTR_NG0_DIRA_HOM_Gd157_VBOC_pcc2 = (D5_HOM_Gd157_VBOC_NG0_DIRA.keff - S2_edep0_setQfiss_pcc2.keff) * 1e5 # pcm

# compare NOEX + SAT with S2 edepmode 0 pcc0
delta_keff_NG0_NOEX_SAT_HOM_Gd157_VBOC_pcc0 = (D5_HOM_Gd157_VBOC_NG0_PCC0_SAT.keff - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm
# compare EXTR + SAT with S2 edepmode 0 pcc2
delta_keff_EXTR_SAT_HOM_Gd157_VBOC_pcc2 = (D5_HOM_Gd157_VBOC_NG0_SAT.keff - S2_edep0_setQfiss_pcc2.keff) * 1e5 # pcm

plt.figure(figsize=(10, 6))
plt.plot(BU_points_plot, delta_keff_NG0_NOEX_HOM_Gd157_VBOC_pcc0, label="D5 NG0 NOEX - S2 edepmode 0 - pcc 0", color='red', linestyle='--', marker='x')
plt.plot(BU_points_plot, delta_keff_NG0_NOEX_DIRA_HOM_Gd157_VBOC_pcc0, label="D5 NG0 NOEX DIRA - S2 edepmode 0 - pcc 0", color='blue', linestyle='--', marker='x')
plt.plot(BU_points_plot, delta_keff_EXTR_NG0_HOM_Gd157_VBOC_pcc2, label="D5 NG0 EXTR - S2 edepmode 0 - pcc 2", color='brown', linestyle='--', marker='x')
plt.plot(BU_points_plot, delta_keff_EXTR_NG0_DIRA_HOM_Gd157_VBOC_pcc2, label="D5 NG0 EXTR DIRA - S2 edepmode 0 - pcc 2", color='green', linestyle='--', marker='x')
plt.plot(BU_points_plot, delta_keff_NG0_NOEX_SAT_HOM_Gd157_VBOC_pcc0, label="D5 NG0 NOEX SAT - S2 edepmode 0 - pcc 0", color='orange', linestyle='--', marker='x')
plt.plot(BU_points_plot, delta_keff_EXTR_SAT_HOM_Gd157_VBOC_pcc2, label="D5 NG0 EXTR SAT - S2 edepmode 0 - pcc 2", color='purple', linestyle='--', marker='x')
plt.plot(BU_points_plot, 300*np.ones_like(BU_points_plot),  color='red', linestyle='--')
plt.plot(BU_points_plot, -300*np.ones_like(BU_points_plot),  color='red', linestyle='--')
plt.xlabel("Burnup (MWd/tU)")
plt.ylabel("Delta keff (pcm)")
plt.title("Delta keff between D5 and S2 for HOM_Gd157_VBOC")
plt.legend()
plt.grid()
plt.savefig(f"{save_dir_HOM_Gd157_VBOC}/delta_keff_HOM_Gd157_VBOC.png")
plt.close()

## Compare results on the isotopic inventory

for iso in tracked_nuclides:
    delta_Niso_NOEX = [(D5_HOM_Gd157_VBOC_NG0_PCC0.DRAGON_ISOTOPESDENS[iso][idx] - S2_edep0_setQfiss_pcc0.Ni[iso][idx]) * 100 / S2_edep0_setQfiss_pcc0.Ni[iso][idx]
                    if S2_edep0_setQfiss_pcc0.Ni[iso][idx] != 0 else 0
                    for idx in range(len(S2_edep0_setQfiss_pcc0.Ni[iso]))]
    delta_Niso_EXTR = [(D5_HOM_Gd157_VBOC_NG0.DRAGON_ISOTOPESDENS[iso][idx] - S2_edep0_setQfiss_pcc2.Ni[iso][idx]) * 100 / S2_edep0_setQfiss_pcc2.Ni[iso][idx]
                    if S2_edep0_setQfiss_pcc2.Ni[iso][idx] != 0 else 0
                    for idx in range(len(S2_edep0_setQfiss_pcc2.Ni[iso]))]
    plt.figure(figsize=(10, 6))
    plt.plot(BU_points_plot, delta_Niso_NOEX, label=f"D5 NG0 NOEX - S2 edepmode 0 - pcc 0 - {iso}", color='blue', linestyle='--', marker='x')
    plt.plot(BU_points_plot, delta_Niso_EXTR, label=f"D5 NG0 EXTR - S2 edepmode 0 - pcc 2 - {iso}", color='green', linestyle='--', marker='x')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel(f"Delta N {iso} (%)")
    plt.title(f"Delta N {iso} between D5 and S2 for HOM_Gd157_VBOC")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_HOM_Gd157_VBOC}/delta_Niso_{iso}_HOM_Gd157_VBOC.png")
    plt.close()



# Case 2 : HOM_UOX_Gd157

# Load COMPOS from PYGAN_RESULTS
name_CPO_HOM_UOX_Gd157 = f"COMPO_HOM_Gd157_{evaluation}_295_USS_PT_CORR_RUNG"
name_CPO_HOM_UOX_Gd157_NG0 = f"COMPO_HOM_Gd157_{evaluation}_295_NG0_USS_PT_CORR_RUNG"

path_to_PYGAN_results = f"{os.getcwd()}/PYGAN_COMPOS_path/HOM_Gd157_Cst_pow_evol_results/"
cwd_path = os.getcwd()
os.chdir(path_to_PYGAN_results)
# Load the data
CPO_HOM_UOX_Gd157 = lcm.new('LCM_INP', name_CPO_HOM_UOX_Gd157, impx=0)
CPO_HOM_UOX_Gd157_NG0 = lcm.new('LCM_INP', name_CPO_HOM_UOX_Gd157_NG0, impx=0)
os.chdir(cwd_path)

# Load S2 results with set fission Q-values, edepmode 0, pcc 0, 1 and 2
# reminder of S2_case object contrustructor (case_name, lib_name, edep_id, areQfissSet, isEcaptSet, pcc_id, specific_power, tracked_nuclides, save_dir)
S2_edep0_setQfiss_pcc0 = S2_case(case_name = "HOM_UOX_Gd157", 
                                lib_name = "endfb8r1_pynjoy2012_kerma", 
                                edep_id = 0, areQfissSet = True, isEcaptSet = False, 
                                pcc_id = 0, specific_power = 38.6, tracked_nuclides = tracked_nuclides, save_dir = save_dir_HOM_UOX_Gd157)

S2_edep0_setQfiss_pcc2 = S2_case(case_name = "HOM_UOX_Gd157",
                                lib_name = "endfb8r1_pynjoy2012_kerma",
                                edep_id = 0, areQfissSet = True, isEcaptSet = False,
                                pcc_id = 2, specific_power = 38.6, tracked_nuclides = tracked_nuclides, save_dir = save_dir_HOM_UOX_Gd157)
# Create DRAGON5 cases with and without (n,gamma) energy deposition 
# reminder of D5_case object contrustructor (pyCOMPO, dlib_name, bu_points, ssh_opt, correlation, sat, depl_sol, tracked_nuclides, BU_lists, save_dir)
D5_HOM_UOX_Gd157 = D5_case(pyCOMPO = CPO_HOM_UOX_Gd157,
                            dlib_name = "endfb8r1_295",
                            bu_points = "BOC_fine_autop5",
                            ssh_opt = "RSE",
                            correlation = "CORR",
                            sat = "",
                            depl_sol = "RUNG",
                            tracked_nuclides = tracked_nuclides,
                            BU_lists = getLists("BOC_fine_autop5"),
                            save_dir = save_dir_HOM_UOX_Gd157)

D5_HOM_UOX_Gd157_NG0 = D5_case(pyCOMPO = CPO_HOM_UOX_Gd157_NG0,
                            dlib_name = "endfb8r1_295_NG0",
                            bu_points = "BOC_fine_autop5",
                            ssh_opt = "RSE",
                            correlation = "CORR",
                            sat = "",
                            depl_sol = "RUNG",
                            tracked_nuclides = tracked_nuclides,
                            BU_lists = getLists("BOC_fine_autop5"),
                            save_dir = save_dir_HOM_UOX_Gd157)

# Compare the results of DRAGON5 and SERPENT2
# could use D5multiS2 object to compare D5 and S2 results
# gonna do it manually for now

delta_keff_HOM_UOX_Gd157_edep0_set_qfiss_pcc0 = (D5_HOM_UOX_Gd157.keff - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm
#delta_keff_HOM_UOX_Gd157_edep0_set_qfiss_pcc1 = (D5_HOM_UOX_Gd157.keff - S2_edep0_setQfiss_pcc1.keff) * 1e5 # pcm
delta_keff_HOM_UOX_Gd157_edep0_set_qfiss_pcc2 = (D5_HOM_UOX_Gd157.keff - S2_edep0_setQfiss_pcc2.keff) * 1e5 # pcm
delta_keff_HOM_UOX_Gd157_NG0_edep0_set_qfiss_pcc0 = (D5_HOM_UOX_Gd157_NG0.keff - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm
#delta_keff_HOM_UOX_Gd157_NG0_edep0_set_qfiss_pcc1 = (D5_HOM_UOX_Gd157_NG0.keff - S2_edep0_setQfiss_pcc1.keff) * 1e5 # pcm
delta_keff_HOM_UOX_Gd157_NG0_edep0_set_qfiss_pcc2 = (D5_HOM_UOX_Gd157_NG0.keff - S2_edep0_setQfiss_pcc2.keff) * 1e5 # pcm

# Plot the results
print(f"{save_dir_HOM_UOX_Gd157}/delta_keff_HOM_UOX_Gd157.png")
plt.figure(figsize=(10, 6))
plt.plot(D5_HOM_UOX_Gd157.BU, delta_keff_HOM_UOX_Gd157_edep0_set_qfiss_pcc0, label="D5 HOM_UOX_Gd157 - S2 edepmode 0 - pcc 0", color='blue', linestyle='--', marker='x')
#plt.plot(D5_HOM_UOX_Gd157.DRAGON_BU, delta_keff_HOM_UOX_Gd157_edep0_set_qfiss_pcc1, label="D5 HOM_UOX_Gd157 - S2 edepmode 0 - pcc 1", color='orange')
plt.plot(D5_HOM_UOX_Gd157.BU, delta_keff_HOM_UOX_Gd157_edep0_set_qfiss_pcc2, label="D5 HOM_UOX_Gd157 - S2 edepmode 0 - pcc 2", color='green', linestyle='--', marker='x')
plt.plot(D5_HOM_UOX_Gd157_NG0.BU, delta_keff_HOM_UOX_Gd157_NG0_edep0_set_qfiss_pcc0, label="D5 HOM_UOX_Gd157 NG0 - S2 edepmode 0 - pcc 0", color='red', linestyle='--', marker='x')
#plt.plot(D5_HOM_UOX_Gd157_NG0.DRAGON_BU, delta_keff_HOM_UOX_Gd157_NG0_edep0_set_qfiss_pcc1, label="D5 HOM_UOX_Gd157 NG0 - S2 edepmode 0 - pcc 1", color='purple')
plt.plot(D5_HOM_UOX_Gd157_NG0.BU, delta_keff_HOM_UOX_Gd157_NG0_edep0_set_qfiss_pcc2, label="D5 HOM_UOX_Gd157 NG0 - S2 edepmode 0 - pcc 2", color='green', linestyle='--', marker='x')
plt.xlabel("Burnup (MWd/tU)")
plt.ylabel("Delta keff (pcm)")
plt.axhline(y=300, color='red', linestyle='--')
plt.axhline(y=-300, color='red', linestyle='--')
plt.title("Delta keff between D5 and S2 for HOM_UOX_Gd157")
plt.legend()
plt.grid()
plt.savefig(f"{save_dir_HOM_UOX_Gd157}/delta_keff_HOM_UOX_Gd157.png")
plt.close()



