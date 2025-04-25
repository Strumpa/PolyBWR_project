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
import lcm


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
post_treat_case0 = False
# Case 1 : HOM_Gd157_VBOC, focus on this
post_treat_case1 = False
# Case 2 : HOM_UOX_Gd157, focus on this
post_treat_case2 = False
# Case 3 : AT10_45Gd, focus on this
post_treat_case3 = False
post_treat_case3_1 = True # pcc scheme for AT10_45Gd
# Case 4 : gduo2_295_kec1
post_treat_case4 = False


if post_treat_case1:
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
    name_CPO_HOM_Gd157_VBOC_NG0_EXTR = f"COMPO_HOM_Gd157_VBOC_{evaluation}_295_NG0"
    name_CPO_HOM_Gd157_VBOC_NG0_NOEX = f"COMPO_HOM_Gd157_VBOC_{evaluation}_295_NG0_PCC0"
    name_CPO_HOM_Gd157_VBOC_NG0_EXTR_DIRA = f"COMPO_HOM_Gd157_VBOC_{evaluation}_295_NG0_DIRA"
    name_CPO_HOM_Gd157_VBOC_NG0_NOEX_DIRA = f"COMPO_HOM_Gd157_VBOC_{evaluation}_295_NG0_PCC0_DIRA"
    name_CPO_HOM_Gd157_VBOC_NG0_NOEX_SAT = f"COMPO_HOM_Gd157_VBOC_{evaluation}_295_NG0_PCC0_SAT"
    name_CPO_HOM_Gd157_VBOC_NG0_EXTR_SAT = f"COMPO_HOM_Gd157_VBOC_{evaluation}_295_NG0_SAT"

    CPO_HOM_Gd157_VBOC_NG0_EXTR_NODI = lcm.new('LCM_INP', name_CPO_HOM_Gd157_VBOC_NG0_EXTR, impx=0)
    CPO_HOM_Gd157_VBOC_NG0_NOEX_NODI = lcm.new('LCM_INP', name_CPO_HOM_Gd157_VBOC_NG0_NOEX, impx=0)
    CPO_HOM_Gd157_VBOC_NG0_EXTR_DIRA = lcm.new('LCM_INP', name_CPO_HOM_Gd157_VBOC_NG0_EXTR_DIRA, impx=0)
    CPO_HOM_Gd157_VBOC_NG0_NOEX_DIRA = lcm.new('LCM_INP', name_CPO_HOM_Gd157_VBOC_NG0_NOEX_DIRA, impx=0)
    CPO_HOM_Gd157_VBOC_NG0_NOEX_SAT = lcm.new('LCM_INP', name_CPO_HOM_Gd157_VBOC_NG0_NOEX_SAT, impx=0)
    CPO_HOM_Gd157_VBOC_NG0_EXTR_SAT = lcm.new('LCM_INP', name_CPO_HOM_Gd157_VBOC_NG0_EXTR_SAT, impx=0)

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


    D5_HOM_Gd157_VBOC_NG0_EXTR_NODI = D5_case(pyCOMPO = CPO_HOM_Gd157_VBOC_NG0_EXTR_NODI,
                                dlib_name = "endfb8r1_295_NG0",
                                bu_points = "VBOC",
                                ssh_opt = "PT",
                                correlation = "CORR",
                                sat = "",
                                depl_sol = "RUNG",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("VBOC"),
                                save_dir = save_dir_HOM_Gd157_VBOC)

    D5_HOM_Gd157_VBOC_NG0_NOEX_NODI = D5_case(pyCOMPO = CPO_HOM_Gd157_VBOC_NG0_NOEX_NODI,
                                dlib_name = "endfb8r1_295_NG0",
                                bu_points = "VBOC",
                                ssh_opt = "PT",
                                correlation = "CORR",
                                sat = "",
                                depl_sol = "RUNG",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("VBOC"),
                                save_dir = save_dir_HOM_Gd157_VBOC)

    D5_HOM_Gd157_VBOC_NG0_EXTR_DIRA = D5_case(pyCOMPO = CPO_HOM_Gd157_VBOC_NG0_EXTR_DIRA,
                                dlib_name = "endfb8r1_295_NG0",
                                bu_points = "VBOC",
                                ssh_opt = "PT",
                                correlation = "CORR",
                                sat = "DIRA",
                                depl_sol = "RUNG",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("VBOC"),
                                save_dir = save_dir_HOM_Gd157_VBOC)

    D5_HOM_Gd157_VBOC_NG0_NOEX_DIRA = D5_case(pyCOMPO = CPO_HOM_Gd157_VBOC_NG0_NOEX_DIRA,
                                dlib_name = "endfb8r1_295_NG0",
                                bu_points = "VBOC",
                                ssh_opt = "PT",
                                correlation = "CORR",
                                sat = "DIRA",
                                depl_sol = "RUNG",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("VBOC"),
                                save_dir = save_dir_HOM_Gd157_VBOC)

    D5_HOM_Gd157_VBOC_NG0_NOEX_SAT = D5_case(pyCOMPO = CPO_HOM_Gd157_VBOC_NG0_NOEX_SAT,
                                dlib_name = "endfb8r1_295_NG0",
                                bu_points = "VBOC",
                                ssh_opt = "PT",
                                correlation = "CORR",
                                sat = "SAT",
                                depl_sol = "RUNG",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("VBOC"),
                                save_dir = save_dir_HOM_Gd157_VBOC)

    D5_HOM_Gd157_VBOC_NG0_EXTR_SAT = D5_case(pyCOMPO = CPO_HOM_Gd157_VBOC_NG0_EXTR_SAT,
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
    BU_points_plot = getLists("VBOC")["COMPO"]
    # compare NG0 NOEX + NODI with S2 edepmode 0 pcc0
    delta_keff_NG0_NOEX_NODI_HOM_Gd157_VBOC_pcc0 = (D5_HOM_Gd157_VBOC_NG0_NOEX_NODI.keff - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm <-- This one works, use as a reference / building bloc for other cases
    # compare NG0 NOEX + DIRA with S2 edepmode 0 pcc0
    delta_keff_NG0_NOEX_DIRA_HOM_Gd157_VBOC_pcc0 = (D5_HOM_Gd157_VBOC_NG0_NOEX_DIRA.keff - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm

    # compare NG0 + NOEX + NODI with S2 edepmode 0 pcc1
    delta_keff_NG0_NOEX_NODI_HOM_Gd157_VBOC_pcc1 = (D5_HOM_Gd157_VBOC_NG0_NOEX_NODI.keff - S2_edep0_setQfiss_pcc1.keff) * 1e5 # pcm
    # compare NG0 + NOEX + DIRA with S2 edepmode 0 pcc1
    delta_keff_NG0_NOEX_DIRA_HOM_Gd157_VBOC_pcc1 = (D5_HOM_Gd157_VBOC_NG0_NOEX_DIRA.keff - S2_edep0_setQfiss_pcc1.keff) * 1e5 # pcm
    # compare NG0 EXTR + NODI with S2 edepmode 0 pcc1
    delta_keff_NG0_EXTR_NODI_HOM_Gd157_VBOC_pcc1 = (D5_HOM_Gd157_VBOC_NG0_EXTR_NODI.keff - S2_edep0_setQfiss_pcc1.keff) * 1e5 # pcm
    # compare NG0 EXTR + DIRA with S2 edepmode 0 pcc1
    delta_keff_NG0_EXTR_DIRA_HOM_Gd157_VBOC_pcc1 = (D5_HOM_Gd157_VBOC_NG0_EXTR_DIRA.keff - S2_edep0_setQfiss_pcc1.keff) * 1e5 # pcm


    # compare NG0 EXTR + NODI with S2 edepmode 0 pcc2
    delta_keff_NG0_EXTR_NODI_HOM_Gd157_VBOC_pcc2 = (D5_HOM_Gd157_VBOC_NG0_EXTR_NODI.keff - S2_edep0_setQfiss_pcc2.keff) * 1e5 # pcm
    # compare NG0 EXTR + DIRA with S2 edepmode 0 pcc2
    delta_keff_NG0_EXTR_DIRA_HOM_Gd157_VBOC_pcc2 = (D5_HOM_Gd157_VBOC_NG0_EXTR_DIRA.keff - S2_edep0_setQfiss_pcc2.keff) * 1e5 # pcm

    # compare NOEX + SAT with S2 edepmode 0 pcc0
    delta_keff_NG0_NOEX_SAT_HOM_Gd157_VBOC_pcc0 = (D5_HOM_Gd157_VBOC_NG0_NOEX_SAT.keff - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm
    # compare EXTR + SAT with S2 edepmode 0 pcc2
    delta_keff_NG0_EXTR_SAT_HOM_Gd157_VBOC_pcc2 = (D5_HOM_Gd157_VBOC_NG0_EXTR_SAT.keff - S2_edep0_setQfiss_pcc2.keff) * 1e5 # pcm

    plt.figure(figsize=(10, 6))
    plt.plot(BU_points_plot, delta_keff_NG0_NOEX_NODI_HOM_Gd157_VBOC_pcc0, label="D5 NG0 NOEX NODI - S2 edepmode 0 - pcc 0", color='blue', linestyle='--', marker='x')
    #plt.plot(BU_points_plot, delta_keff_NG0_NOEX_DIRA_HOM_Gd157_VBOC_pcc0, label="D5 NG0 NOEX DIRA - S2 edepmode 0 - pcc 0", color='red', linestyle='--', marker='x')
    plt.plot(BU_points_plot, delta_keff_NG0_NOEX_NODI_HOM_Gd157_VBOC_pcc1, label="D5 NG0 NOEX NODI - S2 edepmode 0 - pcc 1", color='red', linestyle='--', marker='x')
    plt.plot(BU_points_plot, delta_keff_NG0_EXTR_NODI_HOM_Gd157_VBOC_pcc2, label="D5 NG0 EXTR NODI - S2 edepmode 0 - pcc 2", color='green', linestyle='--', marker='D')
    ##plt.plot(BU_points_plot, delta_keff_NG0_EXTR_DIRA_HOM_Gd157_VBOC_pcc2, label="D5 NG0 EXTR DIRA - S2 edepmode 0 - pcc 2", color='orange', linestyle='--', marker='D')
    #plt.plot(BU_points_plot, delta_keff_NG0_NOEX_SAT_HOM_Gd157_VBOC_pcc0, label="D5 NG0 NOEX SAT - S2 edepmode 0 - pcc 0", color='purple', linestyle='--', marker='x')
    #plt.plot(BU_points_plot, delta_keff_NG0_EXTR_SAT_HOM_Gd157_VBOC_pcc2, label="D5 NG0 EXTR SAT - S2 edepmode 0 - pcc 2", color='brown', linestyle='--', marker='D')
    plt.plot(BU_points_plot, 300*np.ones_like(BU_points_plot),  color='red', linestyle='--')
    plt.plot(BU_points_plot, -300*np.ones_like(BU_points_plot),  color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("Delta keff (pcm)")
    plt.title("Delta keff between D5 and S2 for HOM_Gd157_VBOC")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_HOM_Gd157_VBOC}/delta_keff_HOM_Gd157_VBOC_vs_S2_pcc0_pcc1_pcc2.png")
    plt.close()

    ## Compare results on the isotopic inventory

    for iso in tracked_nuclides:
        delta_Niso_NOEX_NODI = [(D5_HOM_Gd157_VBOC_NG0_NOEX_NODI.DRAGON_ISOTOPESDENS[iso][idx] - S2_edep0_setQfiss_pcc0.Ni[iso][idx]) * 100 / S2_edep0_setQfiss_pcc0.Ni[iso][idx]
                        if S2_edep0_setQfiss_pcc0.Ni[iso][idx] != 0 else 0
                        for idx in range(len(S2_edep0_setQfiss_pcc0.Ni[iso]))]
        delta_Niso_NOEX_NODI_pcc1 = [(D5_HOM_Gd157_VBOC_NG0_NOEX_NODI.DRAGON_ISOTOPESDENS[iso][idx] - S2_edep0_setQfiss_pcc1.Ni[iso][idx]) * 100 / S2_edep0_setQfiss_pcc2.Ni[iso][idx]
                        if S2_edep0_setQfiss_pcc1.Ni[iso][idx] != 0 else 0
                        for idx in range(len(S2_edep0_setQfiss_pcc1.Ni[iso]))]
        delta_Niso_EXTR_NODI_pcc1 = [(D5_HOM_Gd157_VBOC_NG0_EXTR_NODI.DRAGON_ISOTOPESDENS[iso][idx] - S2_edep0_setQfiss_pcc1.Ni[iso][idx]) * 100 / S2_edep0_setQfiss_pcc2.Ni[iso][idx]
                        if S2_edep0_setQfiss_pcc1.Ni[iso][idx] != 0 else 0
                        for idx in range(len(S2_edep0_setQfiss_pcc1.Ni[iso]))]
        delta_Niso_EXTR_NODI = [(D5_HOM_Gd157_VBOC_NG0_EXTR_NODI.DRAGON_ISOTOPESDENS[iso][idx] - S2_edep0_setQfiss_pcc2.Ni[iso][idx]) * 100 / S2_edep0_setQfiss_pcc2.Ni[iso][idx]
                        if S2_edep0_setQfiss_pcc2.Ni[iso][idx] != 0 else 0
                        for idx in range(len(S2_edep0_setQfiss_pcc2.Ni[iso]))]
        plt.figure(figsize=(10, 6))
        plt.plot(BU_points_plot, delta_Niso_NOEX_NODI, label=f"D5 NG0 NOEX NODI - S2 edepmode 0 - pcc 0 - {iso}", color='blue', linestyle='--', marker='x')
        plt.plot(BU_points_plot, delta_Niso_NOEX_NODI_pcc1, label=f"D5 NG0 NOEX NODI - S2 edepmode 0 - pcc 1 - {iso}", color='red', linestyle='--', marker='x')
        plt.plot(BU_points_plot, delta_Niso_EXTR_NODI_pcc1, label=f"D5 NG0 EXTR NODI - S2 edepmode 0 - pcc 1 - {iso}", color='orange', linestyle='--', marker='D')
        plt.plot(BU_points_plot, delta_Niso_EXTR_NODI, label=f"D5 NG0 EXTR NODI - S2 edepmode 0 - pcc 2 - {iso}", color='green', linestyle='--', marker='D')
        plt.axhline(y=2.0, color='red', linestyle='--')
        plt.axhline(y=-2.0, color='red', linestyle='--')
        plt.xlabel("Burnup (MWd/tU)")
        plt.ylabel(f"Delta N {iso} (%)")
        plt.title(f"Delta N {iso} between D5 and S2 for HOM_Gd157_VBOC")
        plt.legend()
        plt.grid()
        plt.savefig(f"{save_dir_HOM_Gd157_VBOC}/delta_Niso_{iso}_HOM_Gd157_VBOC.png")
        plt.close()


    ### Compare D5 NG0 NOEX with D5 NG0 EXTR and S2 edepmode 0 pcc0 with pcc2
    # compare keffs
    delta_keff_D5_NG0_NODI_EXTR_vs_NOEX = (D5_HOM_Gd157_VBOC_NG0_EXTR_NODI.keff - D5_HOM_Gd157_VBOC_NG0_NOEX_NODI.keff) * 1e5 # pcm
    delta_keff_D5_NG0_DIRA_EXTR_vs_NOEX = (D5_HOM_Gd157_VBOC_NG0_EXTR_DIRA.keff - D5_HOM_Gd157_VBOC_NG0_NOEX_DIRA.keff) * 1e5 # pcm
    delta_keff_S2_pcc2_vs_pcc0 = (S2_edep0_setQfiss_pcc2.keff - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm

    plt.figure(figsize=(10, 6))
    plt.plot(BU_points_plot, delta_keff_D5_NG0_NODI_EXTR_vs_NOEX, label="D5 NG0 EXTR vs NOEX NODI", color='blue', linestyle='--', marker='x')
    plt.plot(BU_points_plot, delta_keff_D5_NG0_DIRA_EXTR_vs_NOEX, label="D5 NG0 EXTR vs NOEX DIRA", color='green', linestyle='--', marker='D')
    plt.plot(BU_points_plot, delta_keff_S2_pcc2_vs_pcc0, label="S2 edepmode 0 pcc2 vs pcc0", color='red', linestyle='--', marker='o')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("Delta keff (pcm)")
    plt.axhline(y=300, color='red', linestyle='--')
    plt.axhline(y=-300, color='red', linestyle='--')
    plt.title("Delta keff (EXTR-NOEX) D5 and (pcc2-pcc0) S2")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_HOM_Gd157_VBOC}/delta_keff_D5_NG0_NODI_DIRA_EXTR_vs_NOEX_delta_keff_S2_pcc2_vs_pcc0.png")
    plt.close()

    # compare NGd157 for EXTR vs NOEX:
    delta_NGd157_D5_NG0_NODI_EXTR_vs_NOEX = [(D5_HOM_Gd157_VBOC_NG0_EXTR_NODI.DRAGON_ISOTOPESDENS["Gd157"][idx] - D5_HOM_Gd157_VBOC_NG0_NOEX_NODI.DRAGON_ISOTOPESDENS["Gd157"][idx]) * 100 / D5_HOM_Gd157_VBOC_NG0_NOEX_NODI.DRAGON_ISOTOPESDENS["Gd157"][idx]
                        if D5_HOM_Gd157_VBOC_NG0_NOEX_NODI.DRAGON_ISOTOPESDENS["Gd157"][idx] != 0 else 0
                        for idx in range(len(D5_HOM_Gd157_VBOC_NG0_NOEX_NODI.DRAGON_ISOTOPESDENS["Gd157"]))]
    delta_NGd157_D5_NG0_DIRA_EXTR_vs_NOEX = [(D5_HOM_Gd157_VBOC_NG0_EXTR_DIRA.DRAGON_ISOTOPESDENS["Gd157"][idx] - D5_HOM_Gd157_VBOC_NG0_NOEX_DIRA.DRAGON_ISOTOPESDENS["Gd157"][idx]) * 100 / D5_HOM_Gd157_VBOC_NG0_NOEX_DIRA.DRAGON_ISOTOPESDENS["Gd157"][idx]
                        if D5_HOM_Gd157_VBOC_NG0_NOEX_DIRA.DRAGON_ISOTOPESDENS["Gd157"][idx] != 0 else 0
                        for idx in range(len(D5_HOM_Gd157_VBOC_NG0_NOEX_DIRA.DRAGON_ISOTOPESDENS["Gd157"]))]
    delta_NGd157_S2_pcc2_vs_pcc0 = [(S2_edep0_setQfiss_pcc2.Ni["Gd157"][idx] - S2_edep0_setQfiss_pcc0.Ni["Gd157"][idx]) * 100 / S2_edep0_setQfiss_pcc0.Ni["Gd157"][idx]
                        if S2_edep0_setQfiss_pcc0.Ni["Gd157"][idx] != 0 else 0
                        for idx in range(len(S2_edep0_setQfiss_pcc0.Ni["Gd157"]))]
    plt.figure(figsize=(10, 6))
    plt.plot(BU_points_plot, delta_NGd157_D5_NG0_NODI_EXTR_vs_NOEX, label="D5 NG0 EXTR vs NOEX NODI", color='blue', linestyle='--', marker='x')
    plt.plot(BU_points_plot, delta_NGd157_D5_NG0_DIRA_EXTR_vs_NOEX, label="D5 NG0 EXTR vs NOEX DIRA", color='green', linestyle='--', marker='D')
    plt.plot(BU_points_plot, delta_NGd157_S2_pcc2_vs_pcc0, label="S2 edepmode 0 pcc2 vs pcc0", color='red', linestyle='--', marker='o')
    plt.axhline(y=2.0, color='red', linestyle='--')
    plt.axhline(y=-2.0, color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("Delta N Gd157 (%)")
    plt.title("Delta N Gd157 (EXTR-NOEX) D5 and (pcc2-pcc0) S2")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_HOM_Gd157_VBOC}/delta_NGd157_D5_NG0_NODI_DIRA_EXTR_vs_NOEX_delta_NGd157_S2_pcc2_vs_pcc0.png")
    plt.close()


    #### Idea : Add BU points to BU lists in HOM_Gd15_Cst_pow_evol_proc/getLists.py and investigate if :
    #   - switching from BU_points = VBOC to BU_points = VBOC_finerBU affects the results stored in the COMPO.
    #   - refining BU list gives similar results for NG0 NOEX cases in comparison with S2 edepmode 0 pcc 0
    #   - refining BU list gives similar or better results for NG0 EXTR cases in comparison with S2 edepmode 0 pcc 2 

    # format of COMPO name "_COMPO_HOM_Gd157_{dlib_name}_{ssh_module}_{ssh_method}_{correlation}_{burnup_points}_{depl_solution}_{rates_extrapolation}{SAT}_{DIRA_opt}",
    # stored in PYGAN_RESULTS/HOM_Gd157_Cst_pow_evol_results/
    # _COMPO_HOM_Gd157_endfb8r1_295_NG0_USS_PT_CORR_VBOC_RUNG_NOEX_DIRA

    # Case 1.2 : HOM_Gd157_VBOC through PyGan : sanity check --> checked
    # Load COMPOS from PYGAN_RESULTS
    name_CPO_HOM_Gd157_VBOC_NG0_NOEX_DIRA = f"CPO_endfb8r1_295_NG0_USS_PT_CORR_VBOC_RUNG_NOEX_DIRA_550K"
    name_CPO_HOM_Gd157_VBOC_NG0_EXTR_DIRA = f"CPO_endfb8r1_295_NG0_USS_PT_CORR_VBOC_RUNG_EXTR_DIRA_550K"

    path_to_PYGAN_results = f"{os.getcwd()}/PYGAN_COMPOS_path/HOM_Gd157_Cst_pow_evol_results/"
    cwd_path = os.getcwd()
    os.chdir(path_to_PYGAN_results)
    # Load the data
    CPO_HOM_Gd157_VBOC_NG0_NOEX_DIRA = lcm.new('LCM_INP', name_CPO_HOM_Gd157_VBOC_NG0_NOEX_DIRA, impx=0)
    CPO_HOM_Gd157_VBOC_NG0_EXTR_DIRA = lcm.new('LCM_INP', name_CPO_HOM_Gd157_VBOC_NG0_EXTR_DIRA, impx=0)

    os.chdir(cwd_path)

    # Create D5 cases without (n,gamma) energy deposition
    Case_1_2_VBOC_NG0_NOEX_DIRA = D5_case(pyCOMPO = CPO_HOM_Gd157_VBOC_NG0_NOEX_DIRA,
                                            dlib_name = "endfb8r1_295_NG0",
                                            bu_points = "VBOC",
                                            ssh_opt = "PT",
                                            correlation = "CORR",
                                            sat = "",
                                            depl_sol = "RUNG",
                                            tracked_nuclides = tracked_nuclides,
                                            BU_lists = getLists("VBOC"),
                                            save_dir = save_dir_HOM_Gd157_VBOC)
    Case_1_2_VBOC_NG0_EXTR_DIRA = D5_case(pyCOMPO = CPO_HOM_Gd157_VBOC_NG0_EXTR_DIRA,
                                            dlib_name = "endfb8r1_295_NG0",
                                            bu_points = "VBOC",
                                            ssh_opt = "PT",
                                            correlation = "CORR",
                                            sat = "",
                                            depl_sol = "RUNG",
                                            tracked_nuclides = tracked_nuclides,
                                            BU_lists = getLists("VBOC"),
                                            save_dir = save_dir_HOM_Gd157_VBOC)

        

    # Compare the results of DRAGON5 and SERPENT2
    # compare NG0 NOEX + DIRA with S2 edepmode 0 pcc0
    delta_keff_CASE_1_2_VBOC_NG0_NOEX_DIRA = (Case_1_2_VBOC_NG0_NOEX_DIRA.keff - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm
    delta_keff_CASE_1_2_VBOC_NG0_EXTR = (Case_1_2_VBOC_NG0_EXTR_DIRA.keff - S2_edep0_setQfiss_pcc2.keff) * 1e5 # pcm

    plt.figure(figsize=(10, 6))
    plt.plot(BU_points_plot, delta_keff_CASE_1_2_VBOC_NG0_NOEX_DIRA, label="D5 NG0 NOEX DIRA - S2 edepmode 0 - pcc 0", color='blue', linestyle='--', marker='x')
    plt.plot(BU_points_plot, delta_keff_CASE_1_2_VBOC_NG0_EXTR, label="D5 NG0 EXTR - S2 edepmode 0 - pcc 2", color='green', linestyle='--', marker='D')
    plt.plot(BU_points_plot, 300*np.ones_like(BU_points_plot),  color='red', linestyle='--')
    plt.plot(BU_points_plot, -300*np.ones_like(BU_points_plot),  color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("Delta keff (pcm)")
    plt.title("Delta keff between D5 and S2 for HOM_Gd157_VBOC (PyGan)")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_HOM_Gd157_VBOC}/delta_keff_PyGan_CASE_1_2_VBOC_NG0_NOEX_DIRA.png")

    # Case 1.3 : compare with finer BU list, on same COMPO BU points : hypothesis, (1) this should reduce error because of more points
    # Original hypothesis (1) proved wrong by comparing with D5 with VBOC_finerBU vs S2 edepmode 0 pcc 0 with VBOC list.
    # Hypothesis (2) : S2 results on VBOC are not "temporaly converged" : study S2 results with VBOC_finerBU and VBOC_finer2 
    name_CPO_VBOC_finerBU_NG0_NOEX_DIRA = f"CPO_endfb8r1_295_NG0_USS_PT_CORR_VBOC_finerBU_RUNG_NOEX_DIRA_550K"
    name_CPO_VBOC_finerBU_NG0_EXTR_DIRA = f"CPO_endfb8r1_295_NG0_USS_PT_CORR_VBOC_finerBU_RUNG_EXTR_DIRA_550K"

    name_CPO_VBOC_finerBU2_NG0_NOEX_DIRA = f"CPO_endfb8r1_295_NG0_USS_PT_CORR_VBOC_finerBU2_RUNG_NOEX_DIRA_550K"
    name_CPO_VBOC_finerBU2_NG0_EXTR_DIRA = f"CPO_endfb8r1_295_NG0_USS_PT_CORR_VBOC_finerBU2_RUNG_EXTR_DIRA_550K"
    # Case 1.4 : refine the BU list to VBOC_finer2
    name_CPO_VBOC_finer2_NG0_NOEX_DIRA = f"CPO_endfb8r1_295_NG0_USS_PT_CORR_VBOC_finer2_RUNG_NOEX_DIRA_550K"
    name_CPO_VBOC_finer2_NG0_EXTR_DIRA = f"CPO_endfb8r1_295_NG0_USS_PT_CORR_VBOC_finer2_RUNG_EXTR_DIRA_550K"

    os.chdir(path_to_PYGAN_results)
    CPO_VBOC_finerBU_NG0_NOEX_DIRA = lcm.new('LCM_INP', name_CPO_VBOC_finerBU_NG0_NOEX_DIRA, impx=0)
    CPO_VBOC_finerBU_NG0_EXTR_DIRA = lcm.new('LCM_INP', name_CPO_VBOC_finerBU_NG0_EXTR_DIRA, impx=0)
    CPO_VBOC_finerBU2_NG0_NOEX_DIRA = lcm.new('LCM_INP', name_CPO_VBOC_finerBU2_NG0_NOEX_DIRA, impx=0)
    CPO_VBOC_finerBU2_NG0_EXTR_DIRA = lcm.new('LCM_INP', name_CPO_VBOC_finerBU2_NG0_EXTR_DIRA, impx=0)
    CPO_VBOC_finer2_NG0_NOEX_DIRA = lcm.new('LCM_INP', name_CPO_VBOC_finer2_NG0_NOEX_DIRA, impx=0)
    CPO_VBOC_finer2_NG0_EXTR_DIRA = lcm.new('LCM_INP', name_CPO_VBOC_finer2_NG0_EXTR_DIRA, impx=0)
    os.chdir(cwd_path)

    # Create D5 cases without (n,gamma) energy deposition for case 1.3
    # Case 1.3.1 :
    Case_1_3_VBOC_finerBU_NG0_NOEX_DIRA = D5_case(pyCOMPO = CPO_VBOC_finerBU_NG0_NOEX_DIRA,
                                            dlib_name = "endfb8r1_295_NG0",
                                            bu_points = "VBOC_finerBU",
                                            ssh_opt = "PT",
                                            correlation = "CORR",
                                            sat = "",
                                            depl_sol = "RUNG",
                                            tracked_nuclides = tracked_nuclides,
                                            BU_lists = getLists("VBOC_finerBU"),
                                            save_dir = save_dir_HOM_Gd157_VBOC)
    
    Case_1_3_VBOC_finerBU_NG0_EXTR_DIRA = D5_case(pyCOMPO = CPO_VBOC_finerBU_NG0_EXTR_DIRA,
                                            dlib_name = "endfb8r1_295_NG0",
                                            bu_points = "VBOC_finerBU",
                                            ssh_opt = "PT",
                                            correlation = "CORR",
                                            sat = "",
                                            depl_sol = "RUNG",
                                            tracked_nuclides = tracked_nuclides,
                                            BU_lists = getLists("VBOC_finerBU"),
                                            save_dir = save_dir_HOM_Gd157_VBOC)
    
    # Case 1.3.2 :
    Case_1_32_VBOC_finerBU2_NG0_NOEX_DIRA = D5_case(pyCOMPO = CPO_VBOC_finerBU2_NG0_NOEX_DIRA,
                                            dlib_name = "endfb8r1_295_NG0",
                                            bu_points = "VBOC_finerBU2",
                                            ssh_opt = "PT",
                                            correlation = "CORR",
                                            sat = "",
                                            depl_sol = "RUNG",
                                            tracked_nuclides = tracked_nuclides,
                                            BU_lists = getLists("VBOC_finerBU2"),
                                            save_dir = save_dir_HOM_Gd157_VBOC)
    
    Case_1_32_VBOC_finerBU2_NG0_EXTR_DIRA = D5_case(pyCOMPO = CPO_VBOC_finerBU2_NG0_EXTR_DIRA,
                                            dlib_name = "endfb8r1_295_NG0",
                                            bu_points = "VBOC_finerBU2",
                                            ssh_opt = "PT",
                                            correlation = "CORR",
                                            sat = "",
                                            depl_sol = "RUNG",
                                            tracked_nuclides = tracked_nuclides,
                                            BU_lists = getLists("VBOC_finerBU2"),
                                            save_dir = save_dir_HOM_Gd157_VBOC)
    

    # Case 1.4 : refine the BU list to VBOC_finer2
    Case_1_4_VBOC_finer2_NG0_NOEX_DIRA = D5_case(pyCOMPO = CPO_VBOC_finer2_NG0_NOEX_DIRA,
                                            dlib_name = "endfb8r1_295_NG0",
                                            bu_points = "VBOC_finer2",
                                            ssh_opt = "PT",
                                            correlation = "CORR",
                                            sat = "",
                                            depl_sol = "RUNG",
                                            tracked_nuclides = tracked_nuclides,
                                            BU_lists = getLists("VBOC_finer2"),
                                            save_dir = save_dir_HOM_Gd157_VBOC)

    Case_1_4_VBOC_finer2_NG0_EXTR_DIRA = D5_case(pyCOMPO = CPO_VBOC_finer2_NG0_EXTR_DIRA,
                                            dlib_name = "endfb8r1_295_NG0",
                                            bu_points = "VBOC_finer2",
                                            ssh_opt = "PT",
                                            correlation = "CORR",
                                            sat = "",
                                            depl_sol = "RUNG",
                                            tracked_nuclides = tracked_nuclides,
                                            BU_lists = getLists("VBOC_finer2"),
                                            save_dir = save_dir_HOM_Gd157_VBOC)
    # Create S2 Post Treatment objects for VBOC_finer2 BU list
    # 

    # Case 1.3 : S2 results with VBOC_finerBU list
    S2_edep0_setQfiss_pcc0_VBOC_finerBU = S2_case(case_name = "HOM_Gd157_VBOC_finerBU",
                                    lib_name = "endfb8r1_pynjoy2012_kerma",
                                    edep_id = 0, areQfissSet = True, isEcaptSet = False,
                                    pcc_id = 0, specific_power = 38.6, tracked_nuclides = tracked_nuclides, save_dir = save_dir_HOM_Gd157_VBOC)
    S2_edep0_setQfiss_pcc1_VBOC_finerBU = S2_case(case_name = "HOM_Gd157_VBOC_finerBU",
                                    lib_name = "endfb8r1_pynjoy2012_kerma",
                                    edep_id = 0, areQfissSet = True, isEcaptSet = False,
                                    pcc_id = 1, specific_power = 38.6, tracked_nuclides = tracked_nuclides, save_dir = save_dir_HOM_Gd157_VBOC)
    S2_edep0_setQfiss_pcc2_VBOC_finerBU = S2_case(case_name = "HOM_Gd157_VBOC_finerBU",
                                    lib_name = "endfb8r1_pynjoy2012_kerma",
                                    edep_id = 0, areQfissSet = True, isEcaptSet = False,
                                    pcc_id = 2, specific_power = 38.6, tracked_nuclides = tracked_nuclides, save_dir = save_dir_HOM_Gd157_VBOC)

    # Case 1.3.2 : S2 results with VBOC_finerBU list
    S2_edep0_setQfiss_pcc0_VBOC_finerBU2 = S2_case(case_name = "HOM_Gd157_VBOC_finerBU2",
                                    lib_name = "endfb8r1_pynjoy2012_kerma",
                                    edep_id = 0, areQfissSet = True, isEcaptSet = False,
                                    pcc_id = 0, specific_power = 38.6, tracked_nuclides = tracked_nuclides, save_dir = save_dir_HOM_Gd157_VBOC)
    S2_edep0_setQfiss_pcc1_VBOC_finerBU2 = S2_case(case_name = "HOM_Gd157_VBOC_finerBU2",
                                    lib_name = "endfb8r1_pynjoy2012_kerma",
                                    edep_id = 0, areQfissSet = True, isEcaptSet = False,
                                    pcc_id = 1, specific_power = 38.6, tracked_nuclides = tracked_nuclides, save_dir = save_dir_HOM_Gd157_VBOC)
    S2_edep0_setQfiss_pcc2_VBOC_finerBU2 = S2_case(case_name = "HOM_Gd157_VBOC_finerBU2",
                                    lib_name = "endfb8r1_pynjoy2012_kerma",
                                    edep_id = 0, areQfissSet = True, isEcaptSet = False,
                                    pcc_id = 2, specific_power = 38.6, tracked_nuclides = tracked_nuclides, save_dir = save_dir_HOM_Gd157_VBOC)

    BU_points_plot_finerBU = getLists("VBOC_finerBU")["BU"]
    BU_points_plot_finerBU2 = getLists("VBOC_finerBU2")["BU"]
    # Case 1.4 : S2 results with VBOC_finer2 list
    S2_edep0_setQfiss_pcc0_VBOC_finer2 = S2_case(case_name = "HOM_Gd157_VBOC_finer2",
                                    lib_name = "endfb8r1_pynjoy2012_kerma",
                                    edep_id = 0, areQfissSet = True, isEcaptSet = False,
                                    pcc_id = 0, specific_power = 38.6, tracked_nuclides = tracked_nuclides, save_dir = save_dir_HOM_Gd157_VBOC)
    S2_edep0_setQfiss_pcc1_VBOC_finer2 = S2_case(case_name = "HOM_Gd157_VBOC_finer2",
                                    lib_name = "endfb8r1_pynjoy2012_kerma",
                                    edep_id = 0, areQfissSet = True, isEcaptSet = False,
                                    pcc_id = 1, specific_power = 38.6, tracked_nuclides = tracked_nuclides, save_dir = save_dir_HOM_Gd157_VBOC)
    S2_edep0_setQfiss_pcc2_VBOC_finer2 = S2_case(case_name = "HOM_Gd157_VBOC_finer2",
                                    lib_name = "endfb8r1_pynjoy2012_kerma",
                                    edep_id = 0, areQfissSet = True, isEcaptSet = False,
                                    pcc_id = 2, specific_power = 38.6, tracked_nuclides = tracked_nuclides, save_dir = save_dir_HOM_Gd157_VBOC)
    BU_points_plot_finer2 = getLists("VBOC_finer2")["BU"]
    # Compare the results of DRAGON5 and SERPENT2
    # compare NG0 NOEX + DIRA with S2 edepmode 0 pcc0 (original VBOC)
    delta_keff_CASE_1_3_VBOC_finerBU_NG0_NOEX_DIRA = (Case_1_3_VBOC_finerBU_NG0_NOEX_DIRA.keff - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm
    delta_keff_CASE_1_3_VBOC_finerBU_NG0_EXTR = (Case_1_3_VBOC_finerBU_NG0_EXTR_DIRA.keff - S2_edep0_setQfiss_pcc2.keff) * 1e5 # pcm
    delta_keff_CASE_1_32_VBOC_finerBU2_NG0_NOEX_DIRA = (Case_1_32_VBOC_finerBU2_NG0_NOEX_DIRA.keff - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm
    delta_keff_CASE_1_32_VBOC_finerBU2_NG0_EXTR = (Case_1_32_VBOC_finerBU2_NG0_EXTR_DIRA.keff - S2_edep0_setQfiss_pcc2.keff) * 1e5 # pcm
    delta_keff_CASE_1_4_VBOC_finer2_NG0_NOEX_DIRA = (Case_1_4_VBOC_finer2_NG0_NOEX_DIRA.keff - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm
    delta_keff_CASE_1_4_VBOC_finer2_NG0_EXTR = (Case_1_4_VBOC_finer2_NG0_EXTR_DIRA.keff - S2_edep0_setQfiss_pcc2.keff) * 1e5 # pcm

    plt.figure(figsize=(10, 6))
    plt.plot(BU_points_plot, delta_keff_CASE_1_3_VBOC_finerBU_NG0_NOEX_DIRA, label="D5 NOEX (finerBU) - S2 - pcc 0 (VBOC)", color='blue', linestyle='--', marker='x')
    plt.plot(BU_points_plot, delta_keff_CASE_1_3_VBOC_finerBU_NG0_EXTR, label="D5 EXTR - (finerBU) - S2 - pcc 2 (VBOC)", color='green', linestyle='--', marker='D')
    plt.plot(BU_points_plot, delta_keff_CASE_1_32_VBOC_finerBU2_NG0_NOEX_DIRA, label="D5 NOEX (finerBU2) - S2 - pcc 0 (VBOC)", color='purple', linestyle='--', marker='o')
    plt.plot(BU_points_plot, delta_keff_CASE_1_32_VBOC_finerBU2_NG0_EXTR, label="D5 EXTR (finerBU2) - S2 - pcc 2 (VBOC)", color='brown', linestyle='--', marker='o')
    plt.plot(BU_points_plot, delta_keff_CASE_1_4_VBOC_finer2_NG0_NOEX_DIRA, label="D5 NOEX (finer2) - S2 - pcc 0 (VBOC)", color='red', linestyle='--', marker='x')
    plt.plot(BU_points_plot, delta_keff_CASE_1_4_VBOC_finer2_NG0_EXTR, label="D5 EXTR (finer2) - S2 - pcc 2 (VBOC)", color='orange', linestyle='--', marker='D')
    plt.plot(BU_points_plot, 300*np.ones_like(BU_points_plot),  color='red', linestyle='--')
    plt.plot(BU_points_plot, -300*np.ones_like(BU_points_plot),  color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("$\\Delta$ keff (pcm)")
    plt.title("$\\Delta$ keff between D5 (NG0) and S2 edep 0 for HOM_Gd157, refining D5 BU lists")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_HOM_Gd157_VBOC}/delta_keff_PyGan_CASES_1_3_1_4_VBOC_finerBU_finer2_NG0_vs_edepmode0.png")
    plt.close()


    # Compare D5 EXTR to D5 NOEX
    delta_keff_D5_NG0_EXTR_vs_NOEX_Case12 = (Case_1_2_VBOC_NG0_EXTR_DIRA.keff - Case_1_2_VBOC_NG0_NOEX_DIRA.keff) * 1e5 # pcm
    delta_keff_D5_NG0_EXTR_vs_NOEX_Case13 = (Case_1_3_VBOC_finerBU_NG0_EXTR_DIRA.keff - Case_1_3_VBOC_finerBU_NG0_NOEX_DIRA.keff) * 1e5 # pcm
    delta_keff_D5_NG0_EXTR_vs_NOEX_Case14 = (Case_1_4_VBOC_finer2_NG0_EXTR_DIRA.keff - Case_1_4_VBOC_finer2_NG0_NOEX_DIRA.keff) * 1e5 # pcm
    delta_keff_S2_pcc2_vs_pcc0 = (S2_edep0_setQfiss_pcc2.keff - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm
    delta_keff_S2_pcc2_vs_pcc0_finerBU = (S2_edep0_setQfiss_pcc2_VBOC_finerBU.keff - S2_edep0_setQfiss_pcc0_VBOC_finerBU.keff) * 1e5 # pcm
    delta_keff_S2_pcc2_vs_pcc0_finer2 = (S2_edep0_setQfiss_pcc2_VBOC_finer2.keff - S2_edep0_setQfiss_pcc0_VBOC_finer2.keff) * 1e5 # pcm 
    plt.figure(figsize=(10, 6))
    plt.plot(BU_points_plot, delta_keff_D5_NG0_EXTR_vs_NOEX_Case12, label="D5 NG0 EXTR vs NOEX : VBOC", color='blue', linestyle='--', marker='x')
    plt.plot(BU_points_plot, delta_keff_D5_NG0_EXTR_vs_NOEX_Case13, label="D5 NG0 EXTR vs NOEX : VBOC_finerBU", color='green', linestyle='--', marker='D')
    plt.plot(BU_points_plot, delta_keff_D5_NG0_EXTR_vs_NOEX_Case14, label="D5 NG0 EXTR vs NOEX : VBOC_finest", color='orange', linestyle='--', marker='D')
    plt.plot(BU_points_plot, delta_keff_S2_pcc2_vs_pcc0, label="S2 edepmode 0 pcc2 vs pcc0 : VBOC", color='red', linestyle='--', marker='o')
    plt.plot(BU_points_plot_finerBU, delta_keff_S2_pcc2_vs_pcc0_finerBU, label="S2 edepmode 0 pcc2 vs pcc0 : VBOC_finerBU", color='purple', linestyle='--', marker='o')
    plt.plot(BU_points_plot_finer2, delta_keff_S2_pcc2_vs_pcc0_finer2, label="S2 edepmode 0 pcc2 vs pcc0 : VBOC_finest", color='brown', linestyle='--', marker='o')
    plt.axhline(y=300, color='red', linestyle='--')
    plt.axhline(y=-300, color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("Delta keff (pcm)")
    plt.title("Delta keff between (EXTR-NOEX) D5 and (pcc2-pcc0) S2")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_HOM_Gd157_VBOC}/delta_keff_D5_NG0_EXTR_vs_NOEX_Case12_Case13.png")
    plt.close()


    ### Study S2 convergence between VBOC, VBOC_finerBU and VBOC_finer2
    ## 08/04/2025 : Adding VBOC_finerBU2 to the study, corresponds to VBOC_finerBU, with delta BU = 100 MWd/t on the [4000,8000] MWd/t range
    # ie VBOC_finerBU2 = VBOC_finerBU with VBOC_finer2 points on [4000,8000] MWd/t interval. 
    # assume that VBOC_finer2 is the most converged : 112 BU points vs 58 and 40

    keff_interp_VBOC = np.interp(BU_points_plot_finer2, BU_points_plot, S2_edep0_setQfiss_pcc0.keff) # interpolate S2 edepmode 0 pcc0 keffs to VBOC_finer2
    keff_interp_VBOC_finerBU = np.interp(BU_points_plot_finer2, BU_points_plot_finerBU, S2_edep0_setQfiss_pcc0_VBOC_finerBU.keff) # interpolate S2 edepmode 0 pcc0 keffs to VBOC_finer2
    keff_interp_VBOC_finerBU2 = np.interp(BU_points_plot_finer2, BU_points_plot_finerBU2, S2_edep0_setQfiss_pcc0_VBOC_finerBU2.keff) # interpolate S2 edepmode 0 pcc0 keffs to VBOC_finer2
    keff_VBOC_finer2 = S2_edep0_setQfiss_pcc0_VBOC_finer2.keff # S2 edepmode 0 pcc0 keffs on VBOC_finer2
    delta_keff_VBOC_finer2_VBOC = (keff_interp_VBOC - keff_VBOC_finer2) * 1e5 # pcm
    delta_keff_VBOC_finer2_VBOC_finerBU = (keff_interp_VBOC_finerBU - keff_VBOC_finer2) * 1e5 # pcm
    delta_keff_VBOC_finer2_VBOC_finerBU2 = (keff_interp_VBOC_finerBU2 - keff_VBOC_finer2) * 1e5 # pcm

    

    plt.figure(figsize=(10, 6))
    plt.plot(BU_points_plot_finer2, delta_keff_VBOC_finer2_VBOC, label="S2 edepmode 0 pcc0 VBOC - VBOC_finest", color='blue', linestyle='--', marker='x')
    plt.plot(BU_points_plot_finer2, delta_keff_VBOC_finer2_VBOC_finerBU, label="S2 edepmode 0 pcc0 VBOC_finerBU - VBOC_finest", color='green', linestyle='--', marker='D')
    plt.plot(BU_points_plot_finer2, delta_keff_VBOC_finer2_VBOC_finerBU2, label="S2 edepmode 0 pcc0 VBOC_finerBU2 - VBOC_finest", color='purple', linestyle='--', marker='o')
    plt.axhline(y=300, color='red', linestyle='--')
    plt.axhline(y=-300, color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("Delta keff (pcm)")
    plt.title("Delta keff between S2 edepmode 0 pcc0 VBOC and VBOC_finer2")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_HOM_Gd157_VBOC}/delta_keff_S2_edepmode0_pcc0_VBOC_VBOC_finerBU_vs_VBOC_finer2.png")
    plt.close()

    # Try interpolating S2 edepmode 0 pcc0 VBOC_finerBU, VBOC_finerBU2 and VBOC_finer2 to VBOC
    # pcc0
    keff_interp_VBOC_finerBU_pcc0 = np.interp(BU_points_plot, BU_points_plot_finerBU, S2_edep0_setQfiss_pcc0_VBOC_finerBU.keff) # interpolate S2 edepmode 0 pcc0 keffs to VBOC
    keff_interp_VBOC_finerBU2_pcc0 = np.interp(BU_points_plot, BU_points_plot_finerBU2, S2_edep0_setQfiss_pcc0_VBOC_finerBU2.keff) # interpolate S2 edepmode 0 pcc0 keffs to VBOC
    keff_interp_VBOC_finer2_pcc0 = np.interp(BU_points_plot, BU_points_plot_finer2, S2_edep0_setQfiss_pcc0_VBOC_finer2.keff) # interpolate S2 edepmode 0 pcc0 keffs to VBOC
    keff_VBOC_pcc0 = S2_edep0_setQfiss_pcc0.keff # S2 edepmode 0 pcc0 keffs on VBOC
    # pcc1
    keff_interp_VBOC_finer2_pcc1 = np.interp(BU_points_plot, BU_points_plot_finer2, S2_edep0_setQfiss_pcc1_VBOC_finer2.keff) # interpolate S2 edepmode 0 pcc1 keffs to VBOC
    keff_interp_VBOC_finerBU_pcc1 = np.interp(BU_points_plot, BU_points_plot_finerBU, S2_edep0_setQfiss_pcc1_VBOC_finerBU.keff) # interpolate S2 edepmode 0 pcc1 keffs to VBOC
    keff_interp_VBOC_finerBU2_pcc1 = np.interp(BU_points_plot, BU_points_plot_finerBU2, S2_edep0_setQfiss_pcc1_VBOC_finerBU2.keff) # interpolate S2 edepmode 0 pcc1 keffs to VBOC
    keff_VBOC_pcc1 = S2_edep0_setQfiss_pcc1.keff # S2 edepmode 0 pcc1 keffs on VBOC
    # pcc2
    keff_interp_VBOC_finerBU_pcc2 = np.interp(BU_points_plot, BU_points_plot_finerBU, S2_edep0_setQfiss_pcc2_VBOC_finerBU.keff) # interpolate S2 edepmode 0 pcc2 keffs to VBOC
    keff_interp_VBOC_finerBU2_pcc2 = np.interp(BU_points_plot, BU_points_plot_finerBU2, S2_edep0_setQfiss_pcc2_VBOC_finerBU2.keff) # interpolate S2 edepmode 0 pcc2 keffs to VBOC
    keff_interp_VBOC_finer2_pcc2 = np.interp(BU_points_plot, BU_points_plot_finer2, S2_edep0_setQfiss_pcc2_VBOC_finer2.keff) # interpolate S2 edepmode 0 pcc2 keffs to VBOC
    keff_VBOC_pcc2 = S2_edep0_setQfiss_pcc2.keff # S2 edepmode 0 pcc2 keffs on VBOC

    # pcc0
    delta_keff_VBOC_VBOC_finer2_pcc0 = (keff_VBOC_pcc0 - keff_interp_VBOC_finer2_pcc0) * 1e5 # pcm
    delta_keff_VBOC_finerBU_VBOC_finer2_pcc0 = (keff_interp_VBOC_finerBU_pcc0 - keff_interp_VBOC_finer2_pcc0) * 1e5 # pcm
    delta_keff_VBOC_finerBU2_VBOC_finer2_pcc0 = (keff_interp_VBOC_finerBU2_pcc0 - keff_interp_VBOC_finer2_pcc0) * 1e5 # pcm
    # pcc1
    delta_keff_VBOC_VBOC_finer2_pcc1 = (keff_VBOC_pcc0 - keff_interp_VBOC_finer2_pcc1) * 1e5 # pcm
    delta_keff_VBOC_finerBU_VBOC_finer2_pcc1 = (keff_interp_VBOC_finerBU_pcc1 - keff_interp_VBOC_finer2_pcc1) * 1e5 # pcm
    delta_keff_VBOC_finerBU2_VBOC_finer2_pcc1 = (keff_interp_VBOC_finerBU2_pcc1 - keff_interp_VBOC_finer2_pcc1) * 1e5 # pcm
    # pcc2
    delta_keff_VBOC_VBOC_finer2_pcc2 = (keff_VBOC_pcc2 - keff_interp_VBOC_finer2_pcc2) * 1e5 # pcm
    delta_keff_VBOC_finerBU_VBOC_finer2_pcc2 = (keff_interp_VBOC_finerBU_pcc2 - keff_interp_VBOC_finer2_pcc2) * 1e5 # pcm
    delta_keff_VBOC_finerBU2_VBOC_finer2_pcc2 = (keff_interp_VBOC_finerBU2_pcc2 - keff_interp_VBOC_finer2_pcc2) * 1e5 # pcm

    # pcc0
    plt.figure(figsize=(10, 6))
    plt.plot(BU_points_plot, delta_keff_VBOC_VBOC_finer2_pcc0, label="S2 VBOC - VBOC_finest (pcc0)", color='blue', linestyle='--', marker='x')
    plt.plot(BU_points_plot, delta_keff_VBOC_finerBU_VBOC_finer2_pcc0, label="S2 VBOC_finerBU - VBOC_finest (pcc0)", color='red', linestyle='--', marker='x')
    plt.plot(BU_points_plot, delta_keff_VBOC_finerBU2_VBOC_finer2_pcc0, label="S2 VBOC_finerBU2 - VBOC_finest (pcc0)", color='green', linestyle='--', marker='D')
    plt.axhline(y=300, color='red', linestyle='--')
    plt.axhline(y=-300, color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("$\\Delta$ keff (pcm)")
    plt.title("$\\Delta$ keff between S2 pcc 0 VBOC, finerBU, finerBU2 vs VBOC_finer2")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_HOM_Gd157_VBOC}/delta_keff_S2_VBOC_VBOC_finerBU_finerBU2_vs_VBOC_finer2_pcc0.png")
    plt.close()


    # pcc1
    plt.figure(figsize=(10, 6))
    plt.plot(BU_points_plot, delta_keff_VBOC_VBOC_finer2_pcc1, label="S2 VBOC - VBOC_finest (pcc1)", color='green', linestyle='--', marker='D')
    plt.plot(BU_points_plot, delta_keff_VBOC_finerBU_VBOC_finer2_pcc1, label="S2 VBOC_finerBU - VBOC_finest (pcc1)", color='purple', linestyle='--', marker='D')
    plt.plot(BU_points_plot, delta_keff_VBOC_finerBU2_VBOC_finer2_pcc1, label="S2 VBOC_finerBU2 - VBOC_finest (pcc1)", color='orange', linestyle='--', marker='D')
    print(f"delta keff VBOC_finerBU - VBOC_finer2 pcc1 : {delta_keff_VBOC_finerBU_VBOC_finer2_pcc1}")
    plt.axhline(y=300, color='red', linestyle='--')
    plt.axhline(y=-300, color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("$\\Delta$ keff (pcm)")
    plt.title("$\\Delta$ keff between S2 pcc 1 VBOC, finerBU, finerBU2 vs VBOC_finer2")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_HOM_Gd157_VBOC}/delta_keff_S2_VBOC_VBOC_finerBU_finerBU2_vs_VBOC_finer2_pcc1.png")
    plt.close()


    #pcc2
    plt.figure(figsize=(10, 6))
    plt.plot(BU_points_plot, delta_keff_VBOC_VBOC_finer2_pcc2, label="S2 VBOC - VBOC_finest (pcc2)", color='orange', linestyle='--', marker='D')
    plt.plot(BU_points_plot, delta_keff_VBOC_finerBU_VBOC_finer2_pcc2, label="S2 VBOC_finerBU - VBOC_finest (pcc2)", color='brown', linestyle='--', marker='D')
    plt.plot(BU_points_plot, delta_keff_VBOC_finerBU2_VBOC_finer2_pcc2, label="S2 VBOC_finerBU2 - VBOC_finest (pcc2)", color='purple', linestyle='--', marker='D')
    plt.axhline(y=300, color='red', linestyle='--')
    plt.axhline(y=-300, color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("$\\Delta$ keff (pcm)")
    plt.title("$\\Delta$ keff between S2 pcc 2 VBOC, finerBU, finerBU2 vs VBOC_finer2")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_HOM_Gd157_VBOC}/delta_keff_S2_VBOC_VBOC_finerBU_finerBU2_vs_VBOC_finer2_pcc2.png")
    plt.close()

    # Check S2 convergence for VBOC_finer2 : compare to pcc0 and assess variation when changing pcc

    # Compare S2 edepmode 0 pcc0 and pcc1
    delta_keff_S2_pcc1_vs_pcc0 = (S2_edep0_setQfiss_pcc1_VBOC_finer2.keff - S2_edep0_setQfiss_pcc0_VBOC_finer2.keff) * 1e5 # pcm
    delta_keff_S2_pcc2_vs_pcc0 = (S2_edep0_setQfiss_pcc2_VBOC_finer2.keff - S2_edep0_setQfiss_pcc0_VBOC_finer2.keff) * 1e5 # pcm

    plt.figure(figsize=(10, 6))
    plt.plot(BU_points_plot_finer2, delta_keff_S2_pcc1_vs_pcc0, label="S2 edepmode 0 pcc1 - pcc0", color='blue', linestyle='--', marker='x')
    plt.plot(BU_points_plot_finer2, delta_keff_S2_pcc2_vs_pcc0, label="S2 edepmode 0 pcc2 - pcc0", color='green', linestyle='--', marker='D')
    plt.axhline(y=300, color='red', linestyle='--')
    plt.axhline(y=-300, color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("$\\Delta$ keff (pcm)")
    plt.title("$\\Delta$ keff between S2 edepmode 0 pcc1 and pcc0, VBOC_finer2")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_HOM_Gd157_VBOC}/delta_keff_S2_pcc1_pcc2_vs_pcc0_VBOC_finer2.png")
    plt.close()

    # Check S2 convergence for VBOC_finerBU : compare to pcc0 and assess variation when changing pcc

    # Compare S2 edepmode 0 pcc0 and pcc1
    delta_keff_S2_pcc1_vs_pcc0 = (S2_edep0_setQfiss_pcc1_VBOC_finerBU.keff - S2_edep0_setQfiss_pcc0_VBOC_finerBU.keff) * 1e5 # pcm
    delta_keff_S2_pcc2_vs_pcc0 = (S2_edep0_setQfiss_pcc2_VBOC_finerBU.keff - S2_edep0_setQfiss_pcc0_VBOC_finerBU.keff) * 1e5 # pcm

    plt.figure(figsize=(10, 6))
    plt.plot(BU_points_plot_finerBU, delta_keff_S2_pcc1_vs_pcc0, label="S2 edepmode 0 pcc1 - pcc0", color='blue', linestyle='--', marker='x')
    plt.plot(BU_points_plot_finerBU, delta_keff_S2_pcc2_vs_pcc0, label="S2 edepmode 0 pcc2 - pcc0", color='green', linestyle='--', marker='D')
    plt.axhline(y=300, color='red', linestyle='--')
    plt.axhline(y=-300, color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("$\\Delta$ keff (pcm)")
    plt.title("$\\Delta$ keff between S2 edepmode 0 pcc1 and pcc0, VBOC_finerBU")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_HOM_Gd157_VBOC}/delta_keff_S2_pcc1_pcc2_vs_pcc0_VBOC_finerBU.png")
    plt.close()

    # Check S2 convergence for VBOC : compare to pcc0 and assess variation when changing pcc

    # Compare S2 edepmode 0 pcc0 and pcc1
    delta_keff_S2_pcc1_vs_pcc0 = (S2_edep0_setQfiss_pcc1.keff - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm
    delta_keff_S2_pcc2_vs_pcc0 = (S2_edep0_setQfiss_pcc2.keff - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm

    plt.figure(figsize=(10, 6))
    plt.plot(BU_points_plot, delta_keff_S2_pcc1_vs_pcc0, label="S2 edepmode 0 pcc1 - pcc0", color='blue', linestyle='--', marker='x')
    plt.plot(BU_points_plot, delta_keff_S2_pcc2_vs_pcc0, label="S2 edepmode 0 pcc2 - pcc0", color='green', linestyle='--', marker='D')
    plt.axhline(y=300, color='red', linestyle='--')
    plt.axhline(y=-300, color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("$\\Delta$ keff (pcm)")
    plt.title("$\\Delta$ keff between S2 edepmode 0 pcc1 and pcc0, VBOC")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_HOM_Gd157_VBOC}/delta_keff_S2_pcc1_pcc2_vs_pcc0_VBOC.png")
    plt.close()

    ### Compare D5 to S2 now --> Assess convergence in D5 - S2 comparison
    # because getLists("VBOC_finer2")["COMPO"] is the same as VBOC, (ie DRAGON calculations are performed on VBOC_finer2 but only saved if BU_pt in VBOC)
    # --> need to interpolate S2 results to ListCOMPO
    keff_S2_pcc0_VBOC_finer2_interp = np.interp(BU_points_plot, BU_points_plot_finer2, S2_edep0_setQfiss_pcc0_VBOC_finer2.keff) # interpolate S2 edepmode 0 pcc0 keffs to VBOC
    keff_S2_pcc1_VBOC_finer2_interp = np.interp(BU_points_plot, BU_points_plot_finer2, S2_edep0_setQfiss_pcc1_VBOC_finer2.keff) # interpolate S2 edepmode 0 pcc1 keffs to VBOC
    keff_S2_pcc2_VBOC_finer2_interp = np.interp(BU_points_plot, BU_points_plot_finer2, S2_edep0_setQfiss_pcc2_VBOC_finer2.keff) # interpolate S2 edepmode 0 pcc2 keffs to VBOC
    # Compare D5 to S2 for VBOC_finer2
    # compare NOEX results with S2 edepmode 0 pcc0
    delta_keff_D5_NOEX_VBOC_finer2_pcc0 = (Case_1_4_VBOC_finer2_NG0_NOEX_DIRA.keff - keff_S2_pcc0_VBOC_finer2_interp) * 1e5 # pcm
    # compare EXTR results with S2 edepmode 0 pcc2
    delta_keff_D5_EXTR_VBOC_finer2_pcc2 = (Case_1_4_VBOC_finer2_NG0_EXTR_DIRA.keff - keff_S2_pcc2_VBOC_finer2_interp) * 1e5 # pcm
    # compare NOEX results with S2 edepmode 0 pcc1
    delta_keff_D5_NOEX_VBOC_finer2_pcc1 = (Case_1_4_VBOC_finer2_NG0_NOEX_DIRA.keff - keff_S2_pcc1_VBOC_finer2_interp) * 1e5 # pcm
    # compare EXTR results with S2 edepmode 0 pcc1
    delta_keff_D5_EXTR_VBOC_finer2_pcc1 = (Case_1_4_VBOC_finer2_NG0_EXTR_DIRA.keff - keff_S2_pcc1_VBOC_finer2_interp) * 1e5 # pcm

    plt.figure(figsize=(10, 6))
    plt.plot(BU_points_plot, delta_keff_D5_NOEX_VBOC_finer2_pcc0, label="D5 NOEX - S2 edepmode 0 pcc0, VBOC_finer2", color='blue', linestyle='--', marker='x')
    plt.plot(BU_points_plot, delta_keff_D5_EXTR_VBOC_finer2_pcc2, label="D5 EXTR - S2 edepmode 0 pcc2, VBOC_finer2", color='green', linestyle='--', marker='D')
    plt.plot(BU_points_plot, delta_keff_D5_NOEX_VBOC_finer2_pcc1, label="D5 NOEX - S2 edepmode 0 pcc1, VBOC_finer2", color='red', linestyle='--', marker='x')
    plt.plot(BU_points_plot, delta_keff_D5_EXTR_VBOC_finer2_pcc1, label="D5 EXTR - S2 edepmode 0 pcc1, VBOC_finer2", color='purple', linestyle='--', marker='D')
    plt.axhline(y=300, color='red', linestyle='--')
    plt.axhline(y=-300, color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("$\\Delta$ keff (pcm)")
    plt.title("$\\Delta$ keff between D5 and S2, VBOC_finer2")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_HOM_Gd157_VBOC}/delta_keff_D5_NOEX_EXTR_VBOC_finer2_vs_S2_edepmode0_pcc0_pcc1_pcc2.png")
    plt.close()


    # Case 1.5 : refine the BU list to VBOC_finer2 and save all points to CPO
    name_CPO_VBOC_finer2_NG0_NOEX_DIRA = f"CPO_endfb8r1_295_NG0_USS_PT_CORR_VBOC_finer2_CPO_RUNG_NOEX_DIRA_550K"
    name_CPO_VBOC_finer2_NG0_EXTR_DIRA = f"CPO_endfb8r1_295_NG0_USS_PT_CORR_VBOC_finer2_CPO_RUNG_EXTR_DIRA_550K"

    os.chdir(path_to_PYGAN_results)
    CPO_VBOC_finer2_NG0_NOEX_DIRA_test = lcm.new('LCM_INP', name_CPO_VBOC_finer2_NG0_NOEX_DIRA, impx=0)
    CPO_VBOC_finer2_NG0_EXTR_DIRA_test = lcm.new('LCM_INP', name_CPO_VBOC_finer2_NG0_EXTR_DIRA, impx=0)
    os.chdir(cwd_path)

    # Case 1.5 : refine the BU list to VBOC_finer2 and save all points to CPO
    Case_1_5_VBOC_finer2_NG0_NOEX_DIRA = D5_case(pyCOMPO = CPO_VBOC_finer2_NG0_NOEX_DIRA_test,
                                            dlib_name = "endfb8r1_295_NG0",
                                            bu_points = "VBOC_finer2_CPO",
                                            ssh_opt = "PT",
                                            correlation = "CORR",
                                            sat = "",
                                            depl_sol = "RUNG",
                                            tracked_nuclides = tracked_nuclides,
                                            BU_lists = getLists("VBOC_finer2_CPO"),
                                            save_dir = save_dir_HOM_Gd157_VBOC)

    Case_1_5_VBOC_finer2_NG0_EXTR_DIRA = D5_case(pyCOMPO = CPO_VBOC_finer2_NG0_EXTR_DIRA_test,
                                            dlib_name = "endfb8r1_295_NG0",
                                            bu_points = "VBOC_finer2_CPO",
                                            ssh_opt = "PT",
                                            correlation = "CORR",
                                            sat = "",
                                            depl_sol = "RUNG",
                                            tracked_nuclides = tracked_nuclides,
                                            BU_lists = getLists("VBOC_finer2_CPO"),
                                            save_dir = save_dir_HOM_Gd157_VBOC)
    
    delta_keff_NOEX_finer2_CPO_pcc0 = (Case_1_5_VBOC_finer2_NG0_NOEX_DIRA.keff - S2_edep0_setQfiss_pcc0_VBOC_finer2.keff) * 1e5 # pcm
    delta_keff_EXTR_finer2_CPO_pcc2 = (Case_1_5_VBOC_finer2_NG0_EXTR_DIRA.keff - S2_edep0_setQfiss_pcc2_VBOC_finer2.keff) * 1e5 # pcm
    delta_keff_NOEX_finer2_CPO_pcc1 = (Case_1_5_VBOC_finer2_NG0_NOEX_DIRA.keff - S2_edep0_setQfiss_pcc1_VBOC_finer2.keff) * 1e5 # pcm
    delta_keff_EXTR_finer2_CPO_pcc1 = (Case_1_5_VBOC_finer2_NG0_EXTR_DIRA.keff - S2_edep0_setQfiss_pcc1_VBOC_finer2.keff) * 1e5 # pcm

    plt.figure(figsize=(10, 6))
    plt.plot(BU_points_plot_finer2, delta_keff_NOEX_finer2_CPO_pcc0, label="D5 NOEX - S2 edepmode 0 pcc0, VBOC_finer2", color='blue', linestyle='--', marker='x')
    plt.plot(BU_points_plot_finer2, delta_keff_EXTR_finer2_CPO_pcc2, label="D5 EXTR - S2 edepmode 0 pcc2, VBOC_finer2", color='green', linestyle='--', marker='D')
    plt.plot(BU_points_plot_finer2, delta_keff_NOEX_finer2_CPO_pcc1, label="D5 NOEX - S2 edepmode 0 pcc1, VBOC_finer2", color='red', linestyle='--', marker='x')
    plt.plot(BU_points_plot_finer2, delta_keff_EXTR_finer2_CPO_pcc1, label="D5 EXTR - S2 edepmode 0 pcc1, VBOC_finer2", color='purple', linestyle='--', marker='D')
    plt.axhline(y=300, color='red', linestyle='--')
    plt.axhline(y=-300, color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("$\\Delta$ keff (pcm)")
    plt.title("$\\Delta$ keff between D5 and S2, VBOC_finer2")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_HOM_Gd157_VBOC}/delta_keff_D5_NOEX_EXTR_VBOC_finer2_CPO_vs_S2_edepmode0_pcc0_pcc1_pcc2.png")

    ### Similarly to Serpent2 results : analyse convergence of solutions on coarser BU lists to VBOC_finer2

    # NOEX cases : 
    # VBOC vs VBOC_finer2
    delta_keff_D5_NOEX_VBOC_vs_VBOC_finer2 = (Case_1_2_VBOC_NG0_NOEX_DIRA.keff - Case_1_4_VBOC_finer2_NG0_NOEX_DIRA.keff) * 1e5 # pcm
    # VBOC_finerBU vs VBOC_finer2
    delta_keff_D5_NOEX_VBOC_finerBU_vs_VBOC_finer2 = (Case_1_3_VBOC_finerBU_NG0_NOEX_DIRA.keff - Case_1_4_VBOC_finer2_NG0_NOEX_DIRA.keff) * 1e5 # pcm
    # VBOC_finerBU2 vs VBOC_finer2
    delta_keff_D5_NOEX_VBOC_finerBU2_vs_VBOC_finer2 = (Case_1_32_VBOC_finerBU2_NG0_NOEX_DIRA.keff - Case_1_4_VBOC_finer2_NG0_NOEX_DIRA.keff) * 1e5 # pcm
    # EXTR cases :
    # VBOC vs VBOC_finer2
    delta_keff_D5_EXTR_VBOC_vs_VBOC_finer2 = (Case_1_2_VBOC_NG0_EXTR_DIRA.keff - Case_1_4_VBOC_finer2_NG0_EXTR_DIRA.keff) * 1e5 # pcm
    # VBOC_finerBU vs VBOC_finer2
    delta_keff_D5_EXTR_VBOC_finerBU_vs_VBOC_finer2 = (Case_1_3_VBOC_finerBU_NG0_EXTR_DIRA.keff - Case_1_4_VBOC_finer2_NG0_EXTR_DIRA.keff) * 1e5 # pcm
    # VBOC_finerBU2 vs VBOC_finer2
    delta_keff_D5_EXTR_VBOC_finerBU2_vs_VBOC_finer2 = (Case_1_32_VBOC_finerBU2_NG0_EXTR_DIRA.keff - Case_1_4_VBOC_finer2_NG0_EXTR_DIRA.keff) * 1e5 # pcm
    # Plot the results
    plt.figure(figsize=(10, 6))
    plt.plot(BU_points_plot, delta_keff_D5_NOEX_VBOC_vs_VBOC_finer2, label="D5 NOEX VBOC - VBOC_finer2", color='blue', linestyle='--', marker='x')
    plt.plot(BU_points_plot, delta_keff_D5_NOEX_VBOC_finerBU_vs_VBOC_finer2, label="D5 NOEX VBOC_finerBU - VBOC_finer2", color='red', linestyle='--', marker='x')
    plt.plot(BU_points_plot, delta_keff_D5_NOEX_VBOC_finerBU2_vs_VBOC_finer2, label="D5 NOEX VBOC_finerBU2 - VBOC_finer2", color='green', linestyle='--', marker='D')
    plt.plot(BU_points_plot, delta_keff_D5_EXTR_VBOC_vs_VBOC_finer2, label="D5 EXTR VBOC - VBOC_finer2", color='purple', linestyle='--', marker='D')
    plt.plot(BU_points_plot, delta_keff_D5_EXTR_VBOC_finerBU_vs_VBOC_finer2, label="D5 EXTR VBOC_finerBU - VBOC_finer2", color='orange', linestyle='--', marker='D')
    plt.plot(BU_points_plot, delta_keff_D5_EXTR_VBOC_finerBU2_vs_VBOC_finer2, label="D5 EXTR VBOC_finerBU2 - VBOC_finer2", color='brown', linestyle='--', marker='D')
    plt.axhline(y=300, color='red', linestyle='--')
    plt.axhline(y=-300, color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("$\\Delta$ keff (pcm)")
    plt.title("$\\Delta$ keff between D5 BU lists compared to VBOC_finer2, EXTR and NOEX")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_HOM_Gd157_VBOC}/delta_keff_D5_NOEX_EXTR_VBOC_finerBU_finerBU2_vs_VBOC_finer2.png")
    plt.close()
    

if post_treat_case2:
    # Case 2 : HOM_UOX_Gd157

    # Load COMPOS from PYGAN_RESULTS
    name_CPO_HOM_UOX_Gd157 = "CPO_endfb8r1_295_USS_PT_CORR_BOC_fine_autop5_RUNG_EXTR_DIRA_750K"
    name_CPO_HOM_UOX_Gd157_NG0 ="CPO_endfb8r1_295_NG0_USS_PT_CORR_BOC_fine_autop5_RUNG_EXTR_DIRA_750K"

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

        # Plot the results
    print(f"{save_dir_HOM_UOX_Gd157}/delta_keff_HOM_UOX_Gd157.png")
    plt.figure(figsize=(10, 6))
    plt.plot(D5_HOM_UOX_Gd157.BU, delta_keff_HOM_UOX_Gd157_edep0_set_qfiss_pcc0, label="D5 HOM_UOX_Gd157 - S2 edepmode 0 - pcc 0", color='blue', linestyle='--', marker='x')
    #plt.plot(D5_HOM_UOX_Gd157.DRAGON_BU, delta_keff_HOM_UOX_Gd157_edep0_set_qfiss_pcc1, label="D5 HOM_UOX_Gd157 - S2 edepmode 0 - pcc 1", color='orange')
    #plt.plot(D5_HOM_UOX_Gd157.BU, delta_keff_HOM_UOX_Gd157_edep0_set_qfiss_pcc2, label="D5 HOM_UOX_Gd157 - S2 edepmode 0 - pcc 2", color='green', linestyle='--', marker='x')
    plt.plot(D5_HOM_UOX_Gd157_NG0.BU, delta_keff_HOM_UOX_Gd157_NG0_edep0_set_qfiss_pcc0, label="D5 HOM_UOX_Gd157 NG0 - S2 edepmode 0 - pcc 0", color='red', linestyle='--', marker='x')
    #plt.plot(D5_HOM_UOX_Gd157_NG0.DRAGON_BU, delta_keff_HOM_UOX_Gd157_NG0_edep0_set_qfiss_pcc1, label="D5 HOM_UOX_Gd157 NG0 - S2 edepmode 0 - pcc 1", color='purple')
    #plt.plot(D5_HOM_UOX_Gd157_NG0.BU, delta_keff_HOM_UOX_Gd157_NG0_edep0_set_qfiss_pcc2, label="D5 HOM_UOX_Gd157 NG0 - S2 edepmode 0 - pcc 2", color='green', linestyle='--', marker='x')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("Delta keff (pcm)")
    plt.axhline(y=300, color='red', linestyle='--')
    plt.axhline(y=-300, color='red', linestyle='--')
    plt.title("Delta keff between D5 and S2 pcc0 for HOM_UOX_Gd157")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_HOM_UOX_Gd157}/delta_keff_HOM_UOX_Gd157_PCC0.png")
    plt.close()


if post_treat_case3:
    print("Post-treatment of case 3 : AT10_45Gd_Cst_pow_evol")

    name_CPO_45Gd_Gdautop3_EXTR = "CPO_endfb8r1_295_NG0_RSE_SALT_Gd_autop3_RUNG_NODI_EXTR_GLOB"
    name_CPO_45Gd_Gdautop3_NOEX = "CPO_endfb8r1_295_NG0_RSE_SALT_Gd_autop3_RUNG_NODI_NOEX_GLOB"

    name_CPO_45Gd_Gd_EXTR = "CPO_endfb8r1_295_NG0_RSE_SALT_Gd_RUNG_NODI_EXTR_GLOB"
    name_CPO_45Gd_Gd_NOEX = "CPO_endfb8r1_295_NG0_RSE_SALT_Gd_RUNG_NODI_NOEX_GLOB"
    name_CPO_45Gd_Gd_EXTR_DIRA = "CPO_endfb8r1_295_NG0_RSE_SALT_Gd_RUNG_DIRA_EXTR_GLOB"
    name_CPO_45Gd_Gd_NOEX_DIRA = "CPO_endfb8r1_295_NG0_RSE_SALT_Gd_RUNG_DIRA_NOEX_GLOB"

    name_CPO_45Gd_Gd2_autop6_EXTR = "CPO_endfb8r1_295_NG0_RSE_SALT_Gd2_autop6_RUNG_NODI_EXTR_GLOB"
    name_CPO_45Gd_Gd2_autop6_NOEX = "CPO_endfb8r1_295_NG0_RSE_SALT_Gd2_autop6_RUNG_NODI_NOEX_GLOB"
    name_CPO_45Gd_Gd2_autop6_EXTR_DIRA = "CPO_endfb8r1_295_NG0_RSE_SALT_Gd2_autop6_RUNG_DIRA_EXTR_GLOB"
    name_CPO_45Gd_Gd2_autop6_NOEX_DIRA = "CPO_endfb8r1_295_NG0_RSE_SALT_Gd2_autop6_RUNG_DIRA_NOEX_GLOB"

    "_CPO_{draglib_name}_NG0_{ssh_option}_{tracking_option}_{burnup_points}_{solver_option}_{saturation_option}_{rates_extr}_{glob_opt}"
    name_CPO_Gd3_autop6_NOEX_NODI_GLOB = "CPO_endfb8r1_295_NG0_RSE_SALT_Gd3_autop6_KAPS_NODI_NOEX_GLOB"
    name_CPO_Gd3_autop6_EXTR_NODI_GLOB = "CPO_endfb8r1_295_NG0_RSE_SALT_Gd3_autop6_KAPS_NODI_EXTR_GLOB"
    name_CPO_Gd3_autop6_NOEX_NODI_NOGL = "CPO_endfb8r1_295_NG0_RSE_SALT_Gd3_autop6_KAPS_NODI_NOEX_NOGL"
    name_CPO_Gd3_autop6_EXTR_NODI_NOGL = "CPO_endfb8r1_295_NG0_RSE_SALT_Gd3_autop6_KAPS_NODI_EXTR_NOGL"
        
    name_CPO_Gd4_autop6_NOEX_NODI_GLOB = "CPO_endfb8r1_295_NG0_RSE_SALT_Gd4_autop6_KAPS_NODI_NOEX_GLOB"
    name_CPO_Gd4_autop6_EXTR_NODI_GLOB = "CPO_endfb8r1_295_NG0_RSE_SALT_Gd4_autop6_KAPS_NODI_EXTR_GLOB"
    name_CPO_Gd4_autop6_NOEX_NODI_NOGL = "CPO_endfb8r1_295_NG0_RSE_SALT_Gd4_autop6_KAPS_NODI_NOEX_NOGL"
    name_CPO_Gd4_autop6_EXTR_NODI_NOGL = "CPO_endfb8r1_295_NG0_RSE_SALT_Gd4_autop6_KAPS_NODI_EXTR_NOGL"

    path_to_PYGAN_results = f"{os.getcwd()}/PYGAN_COMPOS_path/AT10_45Gd_Cst_pow_evol_results/"
    cwd_path = os.getcwd()
    os.chdir(path_to_PYGAN_results)
    # Load the data
    CPO_45Gd_Gdautop3_EXTR = lcm.new('LCM_INP', name_CPO_45Gd_Gdautop3_EXTR, impx=0)
    CPO_45Gd_Gdautop3_NOEX = lcm.new('LCM_INP', name_CPO_45Gd_Gdautop3_NOEX, impx=0)

    CPO_45Gd_Gd_EXTR = lcm.new('LCM_INP', name_CPO_45Gd_Gd_EXTR, impx=0)
    CPO_45Gd_Gd_NOEX = lcm.new('LCM_INP', name_CPO_45Gd_Gd_NOEX, impx=0)
    CPO_45Gd_Gd_EXTR_DIRA = lcm.new('LCM_INP', name_CPO_45Gd_Gd_EXTR_DIRA, impx=0)
    CPO_45Gd_Gd_NOEX_DIRA = lcm.new('LCM_INP', name_CPO_45Gd_Gd_NOEX_DIRA, impx=0)

    CPO_45Gd_Gd2_autop6_EXTR = lcm.new('LCM_INP', name_CPO_45Gd_Gd2_autop6_EXTR, impx=0)
    CPO_45Gd_Gd2_autop6_NOEX = lcm.new('LCM_INP', name_CPO_45Gd_Gd2_autop6_NOEX, impx=0)
    CPO_45Gd_Gd2_autop6_EXTR_DIRA = lcm.new('LCM_INP', name_CPO_45Gd_Gd2_autop6_EXTR_DIRA, impx=0)
    CPO_45Gd_Gd2_autop6_NOEX_DIRA = lcm.new('LCM_INP', name_CPO_45Gd_Gd2_autop6_NOEX_DIRA, impx=0)

    CPO_45Gd_Gd3_autop6_NOEX_NODI_GLOB = lcm.new('LCM_INP', name_CPO_Gd3_autop6_NOEX_NODI_GLOB, impx=0)
    CPO_45Gd_Gd3_autop6_EXTR_NODI_GLOB = lcm.new('LCM_INP', name_CPO_Gd3_autop6_EXTR_NODI_GLOB, impx=0)
    CPO_45Gd_Gd3_autop6_NOEX_NODI_NOGL = lcm.new('LCM_INP', name_CPO_Gd3_autop6_NOEX_NODI_NOGL, impx=0)
    CPO_45Gd_Gd3_autop6_EXTR_NODI_NOGL = lcm.new('LCM_INP', name_CPO_Gd3_autop6_EXTR_NODI_NOGL, impx=0)

    CPO_45Gd_Gd4_autop6_NOEX_NODI_GLOB = lcm.new('LCM_INP', name_CPO_Gd4_autop6_NOEX_NODI_GLOB, impx=0)
    CPO_45Gd_Gd4_autop6_EXTR_NODI_GLOB = lcm.new('LCM_INP', name_CPO_Gd4_autop6_EXTR_NODI_GLOB, impx=0)
    CPO_45Gd_Gd4_autop6_NOEX_NODI_NOGL = lcm.new('LCM_INP', name_CPO_Gd4_autop6_NOEX_NODI_NOGL, impx=0)
    CPO_45Gd_Gd4_autop6_EXTR_NODI_NOGL = lcm.new('LCM_INP', name_CPO_Gd4_autop6_EXTR_NODI_NOGL, impx=0)

    os.chdir(cwd_path)
    # Load S2 results with set fission Q-values, edepmode 0, pcc 0, 1 and 2

    # reminder of S2_case object contrustructor (case_name, lib_name, edep_id, areQfissSet, isEcaptSet, pcc_id, specific_power, tracked_nuclides, save_dir)
    S2_edep0_setQfiss_pcc0 = S2_case(case_name = "AT10_45Gd", 
                                    lib_name = "endfb8r1_pynjoy2012_kerma", 
                                    edep_id = 0, areQfissSet = True, isEcaptSet = False, 
                                    pcc_id = 0, specific_power = 26.5, tracked_nuclides = tracked_nuclides, save_dir = save_dir_AT10_45Gd)
    
    S2_edep0_setQfiss_pcc1 = S2_case(case_name = "AT10_45Gd",
                                    lib_name = "endfb8r1_pynjoy2012_kerma",
                                    edep_id = 0, areQfissSet = True, isEcaptSet = False,
                                    pcc_id = 1, specific_power = 26.5, tracked_nuclides = tracked_nuclides, save_dir = save_dir_AT10_45Gd)
    
    S2_edep0_setQfiss_pcc2 = S2_case(case_name = "AT10_45Gd",
                                    lib_name = "endfb8r1_pynjoy2012_kerma",
                                    edep_id = 0, areQfissSet = True, isEcaptSet = False,
                                    pcc_id = 2, specific_power = 26.5, tracked_nuclides = tracked_nuclides, save_dir = save_dir_AT10_45Gd)
    
    # Create DRAGON5 cases without (n,gamma) energy deposition

    # reminder of D5_case object contrustructor (pyCOMPO, dlib_name, bu_points, ssh_opt, correlation, sat, depl_sol, tracked_nuclides, BU_lists, save_dir)
    D5_45Gd_Gdautop3_EXTR = D5_case(pyCOMPO = CPO_45Gd_Gdautop3_EXTR,
                                dlib_name = "endfb8r1_295_NG0",
                                bu_points = "Gd_autop3",
                                ssh_opt = "RSE",
                                correlation = "CORR",
                                sat = "",
                                depl_sol = "RUNG",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("Gd_autop3"),
                                save_dir = save_dir_AT10_45Gd)
    D5_45Gd_Gdautop3_NOEX = D5_case(pyCOMPO = CPO_45Gd_Gdautop3_NOEX,
                                dlib_name = "endfb8r1_295_NG0",
                                bu_points = "Gd_autop3",
                                ssh_opt = "RSE",
                                correlation = "CORR",
                                sat = "",
                                depl_sol = "RUNG",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("Gd_autop3"),
                                save_dir = save_dir_AT10_45Gd)
    
    ## 16/04/2025 : added D5 cases are : 

    D5_45Gd_Gd_EXTR_NODI = D5_case(pyCOMPO = CPO_45Gd_Gd_EXTR,
                                dlib_name = "endfb8r1_295_NG0",
                                bu_points = "Gd",
                                ssh_opt = "RSE",
                                correlation = "CORR",
                                sat = "",
                                depl_sol = "RUNG",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("Gd"),
                                save_dir = save_dir_AT10_45Gd)
    D5_45Gd_Gd_NOEX_NODI = D5_case(pyCOMPO = CPO_45Gd_Gd_NOEX,
                                dlib_name = "endfb8r1_295_NG0",
                                bu_points = "Gd",
                                ssh_opt = "RSE",
                                correlation = "CORR",
                                sat = "",
                                depl_sol = "RUNG",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("Gd"),
                                save_dir = save_dir_AT10_45Gd)
    
    D5_45Gd_Gd_EXTR_DIRA = D5_case(pyCOMPO = CPO_45Gd_Gd_EXTR_DIRA,
                                dlib_name = "endfb8r1_295_NG0",
                                bu_points = "Gd",
                                ssh_opt = "RSE",
                                correlation = "CORR",
                                sat = "",
                                depl_sol = "RUNG",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("Gd"),
                                save_dir = save_dir_AT10_45Gd)
    
    D5_45Gd_Gd_NOEX_DIRA = D5_case(pyCOMPO = CPO_45Gd_Gd_NOEX_DIRA,
                                dlib_name = "endfb8r1_295_NG0",
                                bu_points = "Gd",
                                ssh_opt = "RSE",
                                correlation = "CORR",
                                sat = "",
                                depl_sol = "RUNG",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("Gd"),
                                save_dir = save_dir_AT10_45Gd)
    
    D5_45Gd_Gd2_autop6_EXTR_NODI = D5_case(pyCOMPO = CPO_45Gd_Gd2_autop6_EXTR,
                                dlib_name = "endfb8r1_295_NG0",
                                bu_points = "Gd2_autop6",
                                ssh_opt = "RSE",
                                correlation = "CORR",
                                sat = "",
                                depl_sol = "RUNG",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("Gd2_autop6"),
                                save_dir = save_dir_AT10_45Gd)
    
    D5_45Gd_Gd2_autop6_NOEX_NODI = D5_case(pyCOMPO = CPO_45Gd_Gd2_autop6_NOEX,
                                dlib_name = "endfb8r1_295_NG0",
                                bu_points = "Gd2_autop6",
                                ssh_opt = "RSE",
                                correlation = "CORR",
                                sat = "",
                                depl_sol = "RUNG",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("Gd2_autop6"),
                                save_dir = save_dir_AT10_45Gd)
    
    D5_45Gd_Gd2_autop6_EXTR_DIRA = D5_case(pyCOMPO = CPO_45Gd_Gd2_autop6_EXTR_DIRA,
                                dlib_name = "endfb8r1_295_NG0",
                                bu_points = "Gd2_autop6",
                                ssh_opt = "RSE",
                                correlation = "CORR",
                                sat = "",
                                depl_sol = "RUNG",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("Gd2_autop6"),
                                save_dir = save_dir_AT10_45Gd)
    
    D5_45Gd_Gd2_autop6_NOEX_DIRA = D5_case(pyCOMPO = CPO_45Gd_Gd2_autop6_NOEX_DIRA,
                                dlib_name = "endfb8r1_295_NG0",
                                bu_points = "Gd2_autop6",
                                ssh_opt = "RSE",
                                correlation = "CORR",
                                sat = "",
                                depl_sol = "RUNG",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("Gd2_autop6"),
                                save_dir = save_dir_AT10_45Gd)

    ## 23/04/2025 : added D5 cases are :
    D5_45Gd_Gd3_autop6_NOEX_NODI_GLOB = D5_case(pyCOMPO = CPO_45Gd_Gd3_autop6_NOEX_NODI_GLOB,
                                dlib_name = "endfb8r1_295_NG0",
                                bu_points = "Gd3_autop6",
                                ssh_opt = "RSE",
                                correlation = "NOCORR",
                                sat = "",
                                depl_sol = "KAPS",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("Gd3_autop6"),
                                save_dir = save_dir_AT10_45Gd)
   
    D5_45Gd_Gd3_autop6_EXTR_NODI_GLOB = D5_case(pyCOMPO = CPO_45Gd_Gd3_autop6_EXTR_NODI_GLOB,
                                dlib_name = "endfb8r1_295_NG0",
                                bu_points = "Gd3_autop6",
                                ssh_opt = "RSE",
                                correlation = "NOCORR",
                                sat = "",
                                depl_sol = "KAPS",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("Gd3_autop6"),
                                save_dir = save_dir_AT10_45Gd)
    
    D5_45Gd_Gd3_autop6_NOEX_NODI_NOGL = D5_case(pyCOMPO = CPO_45Gd_Gd3_autop6_NOEX_NODI_NOGL,
                                dlib_name = "endfb8r1_295_NG0",
                                bu_points = "Gd3_autop6",
                                ssh_opt = "RSE",
                                correlation = "NOCORR",
                                sat = "",
                                depl_sol = "KAPS",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("Gd3_autop6"),
                                save_dir = save_dir_AT10_45Gd)
    
    D5_45Gd_Gd3_autop6_EXTR_NODI_NOGL = D5_case(pyCOMPO = CPO_45Gd_Gd3_autop6_EXTR_NODI_NOGL,
                                dlib_name = "endfb8r1_295_NG0",
                                bu_points = "Gd3_autop6",
                                ssh_opt = "RSE",
                                correlation = "NOCORR",
                                sat = "",
                                depl_sol = "KAPS",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("Gd3_autop6"),
                                save_dir = save_dir_AT10_45Gd)
   
    D5_45Gd_Gd4_autop6_NOEX_NODI_GLOB = D5_case(pyCOMPO = CPO_45Gd_Gd4_autop6_NOEX_NODI_GLOB,
                                dlib_name = "endfb8r1_295_NG0",
                                bu_points = "Gd4_autop6",
                                ssh_opt = "RSE",
                                correlation = "NOCORR",
                                sat = "",
                                depl_sol = "KAPS",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("Gd4_autop6"),
                                save_dir = save_dir_AT10_45Gd)
    
    D5_45Gd_Gd4_autop6_EXTR_NODI_GLOB = D5_case(pyCOMPO = CPO_45Gd_Gd4_autop6_EXTR_NODI_GLOB,
                                dlib_name = "endfb8r1_295_NG0",
                                bu_points = "Gd4_autop6",
                                ssh_opt = "RSE",
                                correlation = "NOCORR",
                                sat = "",
                                depl_sol = "KAPS",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("Gd4_autop6"),
                                save_dir = save_dir_AT10_45Gd)
    
    D5_45Gd_Gd4_autop6_NOEX_NODI_NOGL = D5_case(pyCOMPO = CPO_45Gd_Gd4_autop6_NOEX_NODI_NOGL,
                                dlib_name = "endfb8r1_295_NG0",
                                bu_points = "Gd4_autop6",
                                ssh_opt = "RSE",
                                correlation = "NOCORR",
                                sat = "",
                                depl_sol = "KAPS",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("Gd4_autop6"),
                                save_dir = save_dir_AT10_45Gd)
    
    D5_45Gd_Gd4_autop6_EXTR_NODI_NOGL = D5_case(pyCOMPO = CPO_45Gd_Gd4_autop6_EXTR_NODI_NOGL,
                                dlib_name = "endfb8r1_295_NG0",
                                bu_points = "Gd4_autop6",
                                ssh_opt = "RSE",
                                correlation = "NOCORR",
                                sat = "",
                                depl_sol = "KAPS",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("Gd4_autop6"),
                                save_dir = save_dir_AT10_45Gd)

    

    
    # Compare the results of DRAGON5 and SERPENT2
    # NOEX and EXTR cases vs pcc 0
    delta_keff_45Gd_Gdautop3_NOEX_pcc0 = (D5_45Gd_Gdautop3_NOEX.keff - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm
    delta_keff_45Gd_Gdautop3_EXTR_pcc0 = (D5_45Gd_Gdautop3_EXTR.keff - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm

    # NOEX and EXTR cases vs pcc 1
    delta_keff_45Gd_Gdautop3_NOEX_pcc1 = (D5_45Gd_Gdautop3_NOEX.keff - S2_edep0_setQfiss_pcc1.keff) * 1e5 # pcm
    delta_keff_45Gd_Gdautop3_EXTR_pcc1 = (D5_45Gd_Gdautop3_EXTR.keff - S2_edep0_setQfiss_pcc1.keff) * 1e5 # pcm

    # NOEX and EXTR cases vs pcc 2
    delta_keff_45Gd_Gdautop3_NOEX_pcc2 = (D5_45Gd_Gdautop3_NOEX.keff - S2_edep0_setQfiss_pcc2.keff) * 1e5 # pcm
    delta_keff_45Gd_Gdautop3_EXTR_pcc2 = (D5_45Gd_Gdautop3_EXTR.keff - S2_edep0_setQfiss_pcc2.keff) * 1e5 # pcm


    # Plot the results
    # D5 vs pcc 0
    plt.figure(figsize=(10, 6))
    plt.plot(D5_45Gd_Gdautop3_NOEX.BU, delta_keff_45Gd_Gdautop3_NOEX_pcc0, label="D5 NOEX - S2 edepmode 0 - pcc 0", color='blue', linestyle='--', marker='x')
    plt.plot(D5_45Gd_Gdautop3_EXTR.BU, delta_keff_45Gd_Gdautop3_EXTR_pcc0, label="D5 EXTR - S2 edepmode 0 - pcc 0", color='green', linestyle='--', marker='D')
    plt.axhline(y=300, color='red', linestyle='--')
    plt.axhline(y=-300, color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("$\\Delta$ keff (pcm)")
    plt.title("$\\Delta$ keff between D5 and S2 edepmode 0 pcc 0")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_AT10_45Gd}/delta_keff_D5_Gd_autop3_NOEX_EXTR_vs_S2_edepmode0_pcc0.png")
    plt.close()
    
    # D5 vs pcc 1 
    plt.figure(figsize=(10, 6))
    plt.plot(D5_45Gd_Gdautop3_NOEX.BU, delta_keff_45Gd_Gdautop3_NOEX_pcc1, label="D5 NOEX - S2 edepmode 0 - pcc 1", color='red', linestyle='--', marker='x')
    plt.plot(D5_45Gd_Gdautop3_EXTR.BU, delta_keff_45Gd_Gdautop3_EXTR_pcc1, label="D5 EXTR - S2 edepmode 0 - pcc 1", color='purple', linestyle='--', marker='D')
    plt.axhline(y=300, color='red', linestyle='--')
    plt.axhline(y=-300, color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("$\\Delta$ keff (pcm)")
    plt.title("$\\Delta$ keff between D5 and S2 edepmode 0 pcc 1")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_AT10_45Gd}/delta_keff_D5_Gd_autop3_NOEX_EXTR_vs_S2_edepmode0_pcc1.png")
    plt.close()

    # D5 vs pcc 2
    plt.figure(figsize=(10, 6))
    plt.plot(D5_45Gd_Gdautop3_NOEX.BU, delta_keff_45Gd_Gdautop3_NOEX_pcc2, label="D5 NOEX - S2 edepmode 0 - pcc 2", color='blue', linestyle='--', marker='x')
    plt.plot(D5_45Gd_Gdautop3_EXTR.BU, delta_keff_45Gd_Gdautop3_EXTR_pcc2, label="D5 EXTR - S2 edepmode 0 - pcc 2", color='green', linestyle='--', marker='D')
    plt.axhline(y=300, color='red', linestyle='--')
    plt.axhline(y=-300, color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("$\\Delta$ keff (pcm)")
    plt.title("$\\Delta$ keff between D5 and S2 edepmode 0 pcc 2")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_AT10_45Gd}/delta_keff_D5_NOEX_EXTR_vs_S2_edepmode0_pcc2.png")
    plt.close()

    ## Convergence analysis of D5 solution : Compare Gd_autop3, Gd, Gd2_autop6
    # take Gd2_autop6 as reference --> hypothetically most converged

    # D5 Gd_autop3 vs Gd2_autop6

    delta_keff_D5_Gd_autop3_vs_Gd2_autop6 = (D5_45Gd_Gdautop3_EXTR.keff - D5_45Gd_Gd2_autop6_EXTR_NODI.keff) * 1e5 # pcm
    delta_keff_D5_Gd_autop3_NOEX_vs_Gd2_autop6 = (D5_45Gd_Gdautop3_NOEX.keff - D5_45Gd_Gd2_autop6_NOEX_NODI.keff) * 1e5 # pcm
    
    # D5 Gd vs Gd2_autop6
    delta_keff_D5_Gd_vs_Gd2_autop6 = (D5_45Gd_Gd_EXTR_DIRA.keff - D5_45Gd_Gd2_autop6_EXTR_NODI.keff) * 1e5 # pcm
    delta_keff_D5_Gd_NOEX_vs_Gd2_autop6 = (D5_45Gd_Gd_NOEX_DIRA.keff - D5_45Gd_Gd2_autop6_NOEX_NODI.keff) * 1e5 # pcm

    # plot the results :
    plt.figure(figsize=(10, 6))
    plt.plot(D5_45Gd_Gdautop3_EXTR.BU, delta_keff_D5_Gd_autop3_vs_Gd2_autop6, label="D5 Gd_autop3 - Gd2_autop6", color='blue', linestyle='--', marker='x')
    plt.plot(D5_45Gd_Gdautop3_NOEX.BU, delta_keff_D5_Gd_autop3_NOEX_vs_Gd2_autop6, label="D5 Gd_autop3 NOEX - Gd2_autop6", color='green', linestyle='--', marker='D')
    plt.plot(D5_45Gd_Gd_EXTR_DIRA.BU, delta_keff_D5_Gd_vs_Gd2_autop6, label="D5 Gd - Gd2_autop6", color='red', linestyle='--', marker='x')
    plt.plot(D5_45Gd_Gd_NOEX_DIRA.BU, delta_keff_D5_Gd_NOEX_vs_Gd2_autop6, label="D5 Gd NOEX - Gd2_autop6", color='purple', linestyle='--', marker='D')
    plt.axhline(y=300, color='red', linestyle='--')
    plt.axhline(y=-300, color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("$\\Delta$ keff (pcm)")
    plt.title("$\\Delta$ keff between D5 Gd_autop3, Gd and Gd2_autop6")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_AT10_45Gd}/delta_keff_D5_Gd_autop3_vs_Gd2_autop6.png")
    plt.close()


    ## Look at D5_45Gd_Gd2_autop6 cases vs S2 :

    # D5 Gd2_autop6 vs S2 edepmode 0 pcc 0
    delta_keff_D5_Gd2_autop6_NOEX_DIRA_vs_S2_pcc0 = (D5_45Gd_Gd2_autop6_NOEX_DIRA.keff - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm
    delta_keff_D5_Gd2_autop6_EXTR_DIRA_vs_S2_pcc0 = (D5_45Gd_Gd2_autop6_EXTR_DIRA.keff - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm
    delta_keff_D5_Gd2_autop6_NOEX_NODI_vs_S2_pcc0 = (D5_45Gd_Gd2_autop6_NOEX_NODI.keff - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm
    delta_keff_D5_Gd2_autop6_EXTR_NODI_vs_S2_pcc0 = (D5_45Gd_Gd2_autop6_EXTR_NODI.keff - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm

    # D5 Gd2_autop6 vs S2 edepmode 0 pcc 1
    delta_keff_D5_Gd2_autop6_NOEX_DIRA_vs_S2_pcc1 = (D5_45Gd_Gd2_autop6_NOEX_DIRA.keff - S2_edep0_setQfiss_pcc1.keff) * 1e5 # pcm
    delta_keff_D5_Gd2_autop6_EXTR_DIRA_vs_S2_pcc1 = (D5_45Gd_Gd2_autop6_EXTR_DIRA.keff - S2_edep0_setQfiss_pcc1.keff) * 1e5 # pcm
    delta_keff_D5_Gd2_autop6_NOEX_NODI_vs_S2_pcc1 = (D5_45Gd_Gd2_autop6_NOEX_NODI.keff - S2_edep0_setQfiss_pcc1.keff) * 1e5 # pcm
    delta_keff_D5_Gd2_autop6_EXTR_NODI_vs_S2_pcc1 = (D5_45Gd_Gd2_autop6_EXTR_NODI.keff - S2_edep0_setQfiss_pcc1.keff) * 1e5 # pcm
    # D5 Gd2_autop6 vs S2 edepmode 0 pcc 2
    delta_keff_D5_Gd2_autop6_NOEX_DIRA_vs_S2_pcc2 = (D5_45Gd_Gd2_autop6_NOEX_DIRA.keff - S2_edep0_setQfiss_pcc2.keff) * 1e5 # pcm
    delta_keff_D5_Gd2_autop6_EXTR_DIRA_vs_S2_pcc2 = (D5_45Gd_Gd2_autop6_EXTR_DIRA.keff - S2_edep0_setQfiss_pcc2.keff) * 1e5 # pcm
    delta_keff_D5_Gd2_autop6_NOEX_NODI_vs_S2_pcc2 = (D5_45Gd_Gd2_autop6_NOEX_NODI.keff - S2_edep0_setQfiss_pcc2.keff) * 1e5 # pcm
    delta_keff_D5_Gd2_autop6_EXTR_NODI_vs_S2_pcc2 = (D5_45Gd_Gd2_autop6_EXTR_NODI.keff - S2_edep0_setQfiss_pcc2.keff) * 1e5 # pcm


    # Plot the results
    # vs pcc 0
    plt.figure(figsize=(10, 6))
    plt.plot(D5_45Gd_Gd2_autop6_NOEX_DIRA.BU, delta_keff_D5_Gd2_autop6_NOEX_DIRA_vs_S2_pcc0, label="D5 Gd2_autop6 NOEX DIRA - S2 edep 0 - pcc 0", color='blue', linestyle='--', marker='x')
    plt.plot(D5_45Gd_Gd2_autop6_EXTR_DIRA.BU, delta_keff_D5_Gd2_autop6_EXTR_DIRA_vs_S2_pcc0, label="D5 Gd2_autop6 EXTR DIRA - S2 edep 0 - pcc 0", color='green', linestyle='--', marker='D')
    plt.plot(D5_45Gd_Gd2_autop6_NOEX_NODI.BU, delta_keff_D5_Gd2_autop6_NOEX_NODI_vs_S2_pcc0, label="D5 Gd2_autop6 NOEX NODI - S2 edep 0 - pcc 0", color='red', linestyle='--', marker='x')
    plt.plot(D5_45Gd_Gd2_autop6_EXTR_NODI.BU, delta_keff_D5_Gd2_autop6_EXTR_NODI_vs_S2_pcc0, label="D5 Gd2_autop6 EXTR NODI - S2 edep 0 - pcc 0", color='purple', linestyle='--', marker='D')
    plt.axhline(y=300, color='red', linestyle='--')
    plt.axhline(y=-300, color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("$\\Delta$ keff (pcm)")
    plt.title("$\\Delta$ keff between D5 Gd2_autop6 and S2 edepmode 0 pcc 0")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_AT10_45Gd}/delta_keff_D5_Gd2_autop6_vs_S2_edepmode0_pcc0.png")
    plt.close()

    # vs pcc 1
    plt.figure(figsize=(10, 6))
    plt.plot(D5_45Gd_Gd2_autop6_NOEX_DIRA.BU, delta_keff_D5_Gd2_autop6_NOEX_DIRA_vs_S2_pcc1, label="D5 Gd2_autop6 NOEX DIRA - S2 edep 0 - pcc 1", color='blue', linestyle='--', marker='x')
    plt.plot(D5_45Gd_Gd2_autop6_EXTR_DIRA.BU, delta_keff_D5_Gd2_autop6_EXTR_DIRA_vs_S2_pcc1, label="D5 Gd2_autop6 EXTR DIRA - S2 edep 0 - pcc 1", color='green', linestyle='--', marker='D')
    plt.plot(D5_45Gd_Gd2_autop6_NOEX_NODI.BU, delta_keff_D5_Gd2_autop6_NOEX_NODI_vs_S2_pcc1, label="D5 Gd2_autop6 NOEX NODI - S2 edep 0 - pcc 1", color='red', linestyle='--', marker='x')
    plt.plot(D5_45Gd_Gd2_autop6_EXTR_NODI.BU, delta_keff_D5_Gd2_autop6_EXTR_NODI_vs_S2_pcc1, label="D5 Gd2_autop6 EXTR NODI - S2 edep 0 - pcc 1", color='purple', linestyle='--', marker='D')
    plt.axhline(y=300, color='red', linestyle='--')
    plt.axhline(y=-300, color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("$\\Delta$ keff (pcm)")
    plt.title("$\\Delta$ keff between D5 Gd2_autop6 and S2 edepmode 0 pcc 1")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_AT10_45Gd}/delta_keff_D5_Gd2_autop6_vs_S2_edepmode0_pcc1.png")
    plt.close()
    
    # vs pcc 2
    plt.figure(figsize=(10, 6))
    plt.plot(D5_45Gd_Gd2_autop6_NOEX_DIRA.BU, delta_keff_D5_Gd2_autop6_NOEX_DIRA_vs_S2_pcc2, label="D5 Gd2_autop6 NOEX DIRA - S2 edep 0 - pcc 2", color='blue', linestyle='--', marker='x')
    plt.plot(D5_45Gd_Gd2_autop6_EXTR_DIRA.BU, delta_keff_D5_Gd2_autop6_EXTR_DIRA_vs_S2_pcc2, label="D5 Gd2_autop6 EXTR DIRA - S2 edep 0 - pcc 2", color='green', linestyle='--', marker='D')
    plt.plot(D5_45Gd_Gd2_autop6_NOEX_NODI.BU, delta_keff_D5_Gd2_autop6_NOEX_NODI_vs_S2_pcc2, label="D5 Gd2_autop6 NOEX NODI - S2 edep 0 - pcc 2", color='red', linestyle='--', marker='x')
    plt.plot(D5_45Gd_Gd2_autop6_EXTR_NODI.BU, delta_keff_D5_Gd2_autop6_EXTR_NODI_vs_S2_pcc2, label="D5 Gd2_autop6 EXTR NODI - S2 edep 0 - pcc 2", color='purple', linestyle='--', marker='D')
    plt.axhline(y=300, color='red', linestyle='--')
    plt.axhline(y=-300, color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("$\\Delta$ keff (pcm)")
    plt.title("$\\Delta$ keff between D5 Gd2_autop6 and S2 edepmode 0 pcc 2")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_AT10_45Gd}/delta_keff_D5_Gd2_autop6_vs_S2_edepmode0_pcc2.png")
    plt.close()

    ### D5 BU lists : convergence analysis
    # AT10_45Gd on BOC : [0, 20000] MWd/t
    # 17/04/2025 : Defined Gd0_BOC, Gdf_BOC, Gd1_BOC, Gd2_BOC, Gd3_BOC, Gd4_BOC, Gd5_BOC

    # corresponding COMPO names :
    name_45Gd_Gd0 = "CPO_endfb8r1_295_RSE_N_SALT_KAPS_NODI_GLOB_EXTR_Gd0_BOC"
    name_45Gd_Gdf = "CPO_endfb8r1_295_RSE_N_SALT_KAPS_NODI_GLOB_EXTR_Gdf_BOC"
    name_45Gd_Gd1 = "CPO_endfb8r1_295_RSE_N_SALT_KAPS_NODI_GLOB_EXTR_Gd1_BOC"
    name_45Gd_Gd2 = "CPO_endfb8r1_295_RSE_N_SALT_KAPS_NODI_GLOB_EXTR_Gd2_BOC"
    name_45Gd_Gd3 = "CPO_endfb8r1_295_RSE_N_SALT_KAPS_NODI_GLOB_EXTR_Gd3_BOC"
    name_45Gd_Gd4 = "CPO_endfb8r1_295_RSE_N_SALT_KAPS_NODI_GLOB_EXTR_Gd4_BOC"
    name_45Gd_Gd5 = "CPO_endfb8r1_295_RSE_N_SALT_KAPS_NODI_GLOB_EXTR_Gd5_BOC"
    # Load the data
    os.chdir(path_to_PYGAN_results)
    CPO_45Gd_Gd0 = lcm.new('LCM_INP', name_45Gd_Gd0, impx=0)
    CPO_45Gd_Gdf = lcm.new('LCM_INP', name_45Gd_Gdf, impx=0)
    CPO_45Gd_Gd1 = lcm.new('LCM_INP', name_45Gd_Gd1, impx=0)
    CPO_45Gd_Gd2 = lcm.new('LCM_INP', name_45Gd_Gd2, impx=0)
    CPO_45Gd_Gd3 = lcm.new('LCM_INP', name_45Gd_Gd3, impx=0)
    CPO_45Gd_Gd4 = lcm.new('LCM_INP', name_45Gd_Gd4, impx=0)
    CPO_45Gd_Gd5 = lcm.new('LCM_INP', name_45Gd_Gd5, impx=0)
    
    os.chdir(cwd_path)
    
    # Create D5 cases
    D5_45Gd_Gd0 = D5_case(pyCOMPO = CPO_45Gd_Gd0,
                                dlib_name = "endfb8r1_295",
                                bu_points = "Gd0_BOC",
                                ssh_opt = "RSE",
                                correlation = "NOCORR",
                                sat = "",
                                depl_sol = "KAPS",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("Gd0_BOC"),
                                save_dir = save_dir_AT10_45Gd)
    
    D5_45Gd_Gdf = D5_case(pyCOMPO = CPO_45Gd_Gdf,
                                dlib_name = "endfb8r1_295",
                                bu_points = "Gdf_BOC",
                                ssh_opt = "RSE",
                                correlation = "NOCORR",
                                sat = "",
                                depl_sol = "KAPS",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("Gdf_BOC"),
                                save_dir = save_dir_AT10_45Gd)
    D5_45Gd_Gd1 = D5_case(pyCOMPO = CPO_45Gd_Gd1,
                                dlib_name = "endfb8r1_295",
                                bu_points = "Gd1_BOC",
                                ssh_opt = "RSE",
                                correlation = "NOCORR",
                                sat = "",
                                depl_sol = "KAPS",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("Gd1_BOC"),
                                save_dir = save_dir_AT10_45Gd)
    D5_45Gd_Gd2 = D5_case(pyCOMPO = CPO_45Gd_Gd2,
                                dlib_name = "endfb8r1_295",
                                bu_points = "Gd2_BOC",
                                ssh_opt = "RSE",
                                correlation = "NOCORR",
                                sat = "",
                                depl_sol = "KAPS",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("Gd2_BOC"),
                                save_dir = save_dir_AT10_45Gd)
    
    D5_45Gd_Gd3 = D5_case(pyCOMPO = CPO_45Gd_Gd3,
                                dlib_name = "endfb8r1_295",
                                bu_points = "Gd3_BOC",
                                ssh_opt = "RSE",
                                correlation = "NOCORR",
                                sat = "",
                                depl_sol = "KAPS",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("Gd3_BOC"),
                                save_dir = save_dir_AT10_45Gd)
    
    D5_45Gd_Gd4 = D5_case(pyCOMPO = CPO_45Gd_Gd4,
                                dlib_name = "endfb8r1_295",
                                bu_points = "Gd4",
                                ssh_opt = "RSE",
                                correlation = "NOCORR",
                                sat = "",
                                depl_sol = "KAPS",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("Gd4_BOC"),
                                save_dir = save_dir_AT10_45Gd)

    D5_45Gd_Gd5 = D5_case(pyCOMPO = CPO_45Gd_Gd5,
                                dlib_name = "endfb8r1_295",
                                bu_points = "Gd5",
                                ssh_opt = "RSE",
                                correlation = "NOCORR",
                                sat = "",
                                depl_sol = "KAPS",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists("Gd5_BOC"),
                                save_dir = save_dir_AT10_45Gd)


    # Convergence analysis on D5 results :

    # D5 Gd1 vs Gd0
    delta_keff_D5_Gd1_vs_Gd0 = (D5_45Gd_Gd1.keff - D5_45Gd_Gd0.keff) * 1e5 # pcm
    # D5 Gd2 vs Gd0
    delta_keff_D5_Gd2_vs_Gd1 = (D5_45Gd_Gd2.keff - D5_45Gd_Gd0.keff) * 1e5 # pcm
    # D5 Gd3 vs Gd0
    delta_keff_D5_Gd3_vs_Gd0 = (D5_45Gd_Gd3.keff - D5_45Gd_Gd0.keff) * 1e5 # pcm
    # D5 Gd4 vs Gd0
    delta_keff_D5_Gd4_vs_Gd0 = (D5_45Gd_Gd4.keff - D5_45Gd_Gd0.keff) * 1e5 # pcm
    # D5 Gd5 vs Gd0
    delta_keff_D5_Gd5_vs_Gd0 = (D5_45Gd_Gd5.keff - D5_45Gd_Gd0.keff) * 1e5 # pcm
    # D5 Gdf vs Gd0
    delta_keff_D5_Gdf_vs_Gd0 = (D5_45Gd_Gdf.keff - D5_45Gd_Gd0.keff) * 1e5 # pcm
    

    # Plot the results
    plt.figure(figsize=(10, 6))
    plt.plot(D5_45Gd_Gd1.BU, delta_keff_D5_Gd1_vs_Gd0, label="D5 Gd1 - Gd0", color='blue', linestyle='--', marker='x')
    plt.plot(D5_45Gd_Gd2.BU, delta_keff_D5_Gd2_vs_Gd1, label="D5 Gd2 - Gd0", color='green', linestyle='--', marker='D')
    plt.plot(D5_45Gd_Gd3.BU, delta_keff_D5_Gd3_vs_Gd0, label="D5 Gd3 - Gd0", color='red', linestyle='--', marker='x')
    plt.plot(D5_45Gd_Gd4.BU, delta_keff_D5_Gd4_vs_Gd0, label="D5 Gd4 - Gd0", color='purple', linestyle='--', marker='D')
    plt.plot(D5_45Gd_Gd5.BU, delta_keff_D5_Gd5_vs_Gd0, label="D5 Gd5 - Gd0", color='orange', linestyle='--', marker='x')
    plt.plot(D5_45Gd_Gdf.BU, delta_keff_D5_Gdf_vs_Gd0, label="D5 Gdf - Gd0", color='brown', linestyle='--', marker='D')
    plt.axhline(y=300, color='red', linestyle='--')
    plt.axhline(y=-300, color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("$\\Delta$ keff (pcm)")
    plt.title("$\\Delta$ keff between D5 Gd1, Gd2, Gd3, Gd4, Gd5 and Gd0")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_AT10_45Gd}/delta_keff_D5_Gd1_Gd2_Gd3_Gd4_Gd5_vs_Gd0.png")
    plt.close()
    

    ### Compare Gd2_autop6 cases :
                                         

    S2_edep0_setQfiss_pcc1_Gd2 = S2_case(case_name = "AT10_45Gd_BUGd2",
                                    lib_name = "endfb8r1_pynjoy2012_kerma",
                                    edep_id = 0, areQfissSet = True, isEcaptSet = False,
                                    pcc_id = 1, specific_power = 26.5, tracked_nuclides = tracked_nuclides, save_dir = save_dir_AT10_45Gd)

    print("BU points S2 edep0 setQfiss pcc1 Gd2 : ", S2_edep0_setQfiss_pcc1_Gd2.BU)
    
    BU_points_Gd2_CPO = D5_45Gd_Gd2_autop6_NOEX_NODI.BU
    print("BU_points_Gd2_CPO : ", BU_points_Gd2_CPO)
    interp_keff_S2_edep0_setQfiss_pcc1_Gd2 = np.interp(BU_points_Gd2_CPO, S2_edep0_setQfiss_pcc1_Gd2.BU, S2_edep0_setQfiss_pcc1_Gd2.keff)
    delta_keff_D5_S2_Gd2_autop6_NOEX_NODI = (D5_45Gd_Gd2_autop6_NOEX_NODI.keff - interp_keff_S2_edep0_setQfiss_pcc1_Gd2) * 1e5 # pcm
    delta_keff_D5_S2_Gd2_autop6_EXTR_NODI = (D5_45Gd_Gd2_autop6_EXTR_NODI.keff - interp_keff_S2_edep0_setQfiss_pcc1_Gd2) * 1e5 # pcm
    delta_keff_D5_S2_Gd2_autop6_NOEX_DIRA = (D5_45Gd_Gd2_autop6_NOEX_DIRA.keff - interp_keff_S2_edep0_setQfiss_pcc1_Gd2) * 1e5 # pcm
    delta_keff_D5_S2_Gd2_autop6_EXTR_DIRA = (D5_45Gd_Gd2_autop6_EXTR_DIRA.keff - interp_keff_S2_edep0_setQfiss_pcc1_Gd2) * 1e5 # pcm

    # Plot the results
    plt.figure(figsize=(10, 6))
    plt.plot(D5_45Gd_Gd2_autop6_NOEX_NODI.BU, delta_keff_D5_S2_Gd2_autop6_NOEX_NODI, label="D5 Gd2_autop6 NOEX NODI - S2 edep 0 - pcc 1", color='blue', linestyle='--', marker='x')
    plt.plot(D5_45Gd_Gd2_autop6_EXTR_NODI.BU, delta_keff_D5_S2_Gd2_autop6_EXTR_NODI, label="D5 Gd2_autop6 EXTR NODI - S2 edep 0 - pcc 1", color='green', linestyle='--', marker='D')
    plt.plot(D5_45Gd_Gd2_autop6_NOEX_DIRA.BU, delta_keff_D5_S2_Gd2_autop6_NOEX_DIRA, label="D5 Gd2_autop6 NOEX DIRA - S2 edep 0 - pcc 1", color='red', linestyle='--', marker='x')
    plt.plot(D5_45Gd_Gd2_autop6_EXTR_DIRA.BU, delta_keff_D5_S2_Gd2_autop6_EXTR_DIRA, label="D5 Gd2_autop6 EXTR DIRA - S2 edep 0 - pcc 1", color='purple', linestyle='--', marker='D')
    plt.axhline(y=300, color='red', linestyle='--')
    plt.axhline(y=-300, color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("$\\Delta$ keff (pcm)")
    plt.title("$\\Delta$ keff between D5 Gd2_autop6 and S2 edepmode 0 pcc 1")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_AT10_45Gd}/delta_keff_D5_Gd2_autop6_vs_S2_edepmode0_pcc1_new.png")
    plt.close()

    # Plot keffs
    plt.figure(figsize=(10, 6))
    plt.plot(D5_45Gd_Gd2_autop6_NOEX_NODI.BU, D5_45Gd_Gd2_autop6_NOEX_NODI.keff, label="D5 Gd2_autop6 NOEX NODI", color='blue', linestyle='--', marker='x')
    plt.plot(S2_edep0_setQfiss_pcc1_Gd2.BU, S2_edep0_setQfiss_pcc1_Gd2.keff, label="S2_edep0 setQfiss pcc1 Gd2", color='green', linestyle='--', marker='D')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("keff")
    plt.title("keff between D5 Gd2_autop6 and S2 edepmode 0 pcc 1")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_AT10_45Gd}/keff_D5_Gd2_autop6_vs_S2_edepmode0_pcc1_new.png")
    plt.close()

    

    # Look at Gd3_autop6 and Gd4_autop6 cases vs S2 Gd
    # interpolate D5 keffs to S2 BU points
    interp_D5_keff_Gd3_autop6_NOEX_NODI_GLOB = np.interp(S2_edep0_setQfiss_pcc0.BU, D5_45Gd_Gd3_autop6_NOEX_NODI_GLOB.BU, D5_45Gd_Gd3_autop6_NOEX_NODI_GLOB.keff)
    interp_D5_keff_Gd3_autop6_EXTR_NODI_GLOB = np.interp(S2_edep0_setQfiss_pcc0.BU, D5_45Gd_Gd3_autop6_EXTR_NODI_GLOB.BU, D5_45Gd_Gd3_autop6_EXTR_NODI_GLOB.keff)
    interp_D5_keff_Gd3_autop6_NOEX_NODI_NOGL = np.interp(S2_edep0_setQfiss_pcc0.BU, D5_45Gd_Gd3_autop6_NOEX_NODI_NOGL.BU, D5_45Gd_Gd3_autop6_NOEX_NODI_NOGL.keff)
    interp_D5_keff_Gd3_autop6_EXTR_NODI_NOGL = np.interp(S2_edep0_setQfiss_pcc0.BU, D5_45Gd_Gd3_autop6_EXTR_NODI_NOGL.BU, D5_45Gd_Gd3_autop6_EXTR_NODI_NOGL.keff)
    # D5 Gd3_autop6 vs S2 edepmode 0 pcc 0
    delta_keff_Gd3_autop6_NOEX_NODI_GLOB_vs_S2_pcc0 = (interp_D5_keff_Gd3_autop6_NOEX_NODI_GLOB - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm
    delta_keff_Gd3_autop6_EXTR_NODI_GLOB_vs_S2_pcc0 = (interp_D5_keff_Gd3_autop6_EXTR_NODI_GLOB - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm
    delta_keff_Gd3_autop6_NOEX_NODI_NOGL_vs_S2_pcc0 = (interp_D5_keff_Gd3_autop6_NOEX_NODI_NOGL - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm
    delta_keff_Gd3_autop6_EXTR_NODI_NOGL_vs_S2_pcc0 = (interp_D5_keff_Gd3_autop6_EXTR_NODI_NOGL - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm
    # D5 Gd3_autop6 vs S2 edepmode 0 pcc 1
    delta_keff_Gd3_autop6_NOEX_NODI_GLOB_vs_S2_pcc1 = (interp_D5_keff_Gd3_autop6_NOEX_NODI_GLOB - S2_edep0_setQfiss_pcc1.keff) * 1e5 # pcm
    delta_keff_Gd3_autop6_EXTR_NODI_GLOB_vs_S2_pcc1 = (interp_D5_keff_Gd3_autop6_EXTR_NODI_GLOB - S2_edep0_setQfiss_pcc1.keff) * 1e5 # pcm
    delta_keff_Gd3_autop6_NOEX_NODI_NOGL_vs_S2_pcc1 = (interp_D5_keff_Gd3_autop6_NOEX_NODI_NOGL - S2_edep0_setQfiss_pcc1.keff) * 1e5 # pcm
    delta_keff_Gd3_autop6_EXTR_NODI_NOGL_vs_S2_pcc1 = (interp_D5_keff_Gd3_autop6_EXTR_NODI_NOGL - S2_edep0_setQfiss_pcc1.keff) * 1e5 # pcm
    # D5 Gd3_autop6 vs S2 edepmode 0 pcc 2
    delta_keff_Gd3_autop6_NOEX_NODI_GLOB_vs_S2_pcc2 = (interp_D5_keff_Gd3_autop6_NOEX_NODI_GLOB - S2_edep0_setQfiss_pcc2.keff) * 1e5 # pcm
    delta_keff_Gd3_autop6_EXTR_NODI_GLOB_vs_S2_pcc2 = (interp_D5_keff_Gd3_autop6_EXTR_NODI_GLOB - S2_edep0_setQfiss_pcc2.keff) * 1e5 # pcm
    delta_keff_Gd3_autop6_NOEX_NODI_NOGL_vs_S2_pcc2 = (interp_D5_keff_Gd3_autop6_NOEX_NODI_NOGL - S2_edep0_setQfiss_pcc2.keff) * 1e5 # pcm
    delta_keff_Gd3_autop6_EXTR_NODI_NOGL_vs_S2_pcc2 = (interp_D5_keff_Gd3_autop6_EXTR_NODI_NOGL - S2_edep0_setQfiss_pcc2.keff) * 1e5 # pcm
    
    
    # D5 Gd4_autop6 vs S2 edepmode 0 pcc 0
    interp_D5_keff_Gd4_autop6_NOEX_NODI_GLOB = np.interp(S2_edep0_setQfiss_pcc0.BU, D5_45Gd_Gd4_autop6_NOEX_NODI_GLOB.BU, D5_45Gd_Gd4_autop6_NOEX_NODI_GLOB.keff)
    interp_D5_keff_Gd4_autop6_EXTR_NODI_GLOB = np.interp(S2_edep0_setQfiss_pcc0.BU, D5_45Gd_Gd4_autop6_EXTR_NODI_GLOB.BU, D5_45Gd_Gd4_autop6_EXTR_NODI_GLOB.keff)
    interp_D5_keff_Gd4_autop6_NOEX_NODI_NOGL = np.interp(S2_edep0_setQfiss_pcc0.BU, D5_45Gd_Gd4_autop6_NOEX_NODI_NOGL.BU, D5_45Gd_Gd4_autop6_NOEX_NODI_NOGL.keff)
    interp_D5_keff_Gd4_autop6_EXTR_NODI_NOGL = np.interp(S2_edep0_setQfiss_pcc0.BU, D5_45Gd_Gd4_autop6_EXTR_NODI_NOGL.BU, D5_45Gd_Gd4_autop6_EXTR_NODI_NOGL.keff)
    # D5 Gd4_autop6 vs S2 edepmode 0 pcc 0
    delta_keff_Gd4_autop6_NOEX_NODI_GLOB_vs_S2_pcc0 = (interp_D5_keff_Gd4_autop6_NOEX_NODI_GLOB - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm
    delta_keff_Gd4_autop6_EXTR_NODI_GLOB_vs_S2_pcc0 = (interp_D5_keff_Gd4_autop6_EXTR_NODI_GLOB - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm
    delta_keff_Gd4_autop6_NOEX_NODI_NOGL_vs_S2_pcc0 = (interp_D5_keff_Gd4_autop6_NOEX_NODI_NOGL - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm
    delta_keff_Gd4_autop6_EXTR_NODI_NOGL_vs_S2_pcc0 = (interp_D5_keff_Gd4_autop6_EXTR_NODI_NOGL - S2_edep0_setQfiss_pcc0.keff) * 1e5 # pcm
    # D5 Gd4_autop6 vs S2 edepmode 0 pcc 1
    delta_keff_Gd4_autop6_NOEX_NODI_GLOB_vs_S2_pcc1 = (interp_D5_keff_Gd4_autop6_NOEX_NODI_GLOB - S2_edep0_setQfiss_pcc1.keff) * 1e5 # pcm
    delta_keff_Gd4_autop6_EXTR_NODI_GLOB_vs_S2_pcc1 = (interp_D5_keff_Gd4_autop6_EXTR_NODI_GLOB - S2_edep0_setQfiss_pcc1.keff) * 1e5 # pcm
    delta_keff_Gd4_autop6_NOEX_NODI_NOGL_vs_S2_pcc1 = (interp_D5_keff_Gd4_autop6_NOEX_NODI_NOGL - S2_edep0_setQfiss_pcc1.keff) * 1e5 # pcm
    delta_keff_Gd4_autop6_EXTR_NODI_NOGL_vs_S2_pcc1 = (interp_D5_keff_Gd4_autop6_EXTR_NODI_NOGL - S2_edep0_setQfiss_pcc1.keff) * 1e5 # pcm
    # D5 Gd4_autop6 vs S2 edepmode 0 pcc 2
    delta_keff_Gd4_autop6_NOEX_NODI_GLOB_vs_S2_pcc2 = (interp_D5_keff_Gd4_autop6_NOEX_NODI_GLOB - S2_edep0_setQfiss_pcc2.keff) * 1e5 # pcm
    delta_keff_Gd4_autop6_EXTR_NODI_GLOB_vs_S2_pcc2 = (interp_D5_keff_Gd4_autop6_EXTR_NODI_GLOB - S2_edep0_setQfiss_pcc2.keff) * 1e5 # pcm
    delta_keff_Gd4_autop6_NOEX_NODI_NOGL_vs_S2_pcc2 = (interp_D5_keff_Gd4_autop6_NOEX_NODI_NOGL - S2_edep0_setQfiss_pcc2.keff) * 1e5 # pcm
    delta_keff_Gd4_autop6_EXTR_NODI_NOGL_vs_S2_pcc2 = (interp_D5_keff_Gd4_autop6_EXTR_NODI_NOGL - S2_edep0_setQfiss_pcc2.keff) * 1e5 # pcm

    # plot results : Gd3_autop6/Gd4_autop3 vs S2 edepmode 0 pcc 0 NOEX vs EXTR / GLOB
    plt.figure(figsize=(10, 6))
    plt.plot(S2_edep0_setQfiss_pcc0.BU, delta_keff_Gd3_autop6_EXTR_NODI_GLOB_vs_S2_pcc0, label="D5 Gd3_autop6 EXTR GLOB - S2 edep 0 - pcc 0", color='blue', linestyle='--', marker='x')
    plt.plot(S2_edep0_setQfiss_pcc0.BU, delta_keff_Gd3_autop6_NOEX_NODI_GLOB_vs_S2_pcc0, label="D5 Gd3_autop6 NOEX GLOB - S2 edep 0 - pcc 0", color='green', linestyle='--', marker='D')
    plt.plot(S2_edep0_setQfiss_pcc0.BU, delta_keff_Gd4_autop6_EXTR_NODI_GLOB_vs_S2_pcc0, label="D5 Gd4_autop6 EXTR GLOB - S2 edep 0 - pcc 0", color='red', linestyle='--', marker='x')
    plt.plot(S2_edep0_setQfiss_pcc0.BU, delta_keff_Gd4_autop6_NOEX_NODI_GLOB_vs_S2_pcc0, label="D5 Gd4_autop6 NOEX GLOB - S2 edep 0 - pcc 0", color='purple', linestyle='--', marker='D')
    plt.axhline(y=300, color='red', linestyle='--')
    plt.axhline(y=-300, color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("$\\Delta$ keff (pcm)")
    plt.title("$\\Delta$ keff between D5 Gd3_autop6/Gd4_autop6 and S2 edepmode 0 pcc 0")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_AT10_45Gd}/delta_keff_D5_Gd3_autop6_Gd4_autop6_vs_S2_edepmode0_pcc0.png")
    plt.close()

    # plot results : Gd3_autop6/Gd4_autop3 vs S2 edepmode 0 pcc 1 NOEX vs EXTR / GLOB
    plt.figure(figsize=(10, 6))
    plt.plot(S2_edep0_setQfiss_pcc1.BU, delta_keff_Gd3_autop6_EXTR_NODI_GLOB_vs_S2_pcc1, label="D5 Gd3_autop6 EXTR GLOB - S2 edep 0 - pcc 1", color='blue', linestyle='--', marker='x')
    plt.plot(S2_edep0_setQfiss_pcc1.BU, delta_keff_Gd3_autop6_NOEX_NODI_GLOB_vs_S2_pcc1, label="D5 Gd3_autop6 NOEX GLOB - S2 edep 0 - pcc 1", color='green', linestyle='--', marker='D')
    plt.plot(S2_edep0_setQfiss_pcc1.BU, delta_keff_Gd4_autop6_EXTR_NODI_GLOB_vs_S2_pcc1, label="D5 Gd4_autop6 EXTR GLOB - S2 edep 0 - pcc 1", color='red', linestyle='--', marker='x')
    plt.plot(S2_edep0_setQfiss_pcc1.BU, delta_keff_Gd4_autop6_NOEX_NODI_GLOB_vs_S2_pcc1, label="D5 Gd4_autop6 NOEX GLOB - S2 edep 0 - pcc 1", color='purple', linestyle='--', marker='D')
    plt.axhline(y=300, color='red', linestyle='--')
    plt.axhline(y=-300, color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("$\\Delta$ keff (pcm)")
    plt.title("$\\Delta$ keff between D5 Gd3_autop6/Gd4_autop6 and S2 edepmode 0 pcc 1")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_AT10_45Gd}/delta_keff_D5_Gd3_autop6_Gd4_autop6_vs_S2_edepmode0_pcc1.png")
    plt.close()

    # plot results : Gd3_autop6/Gd4_autop3 vs S2 edepmode 0 pcc 2 NOEX vs EXTR / GLOB
    plt.figure(figsize=(10, 6))
    plt.plot(S2_edep0_setQfiss_pcc2.BU, delta_keff_Gd3_autop6_EXTR_NODI_GLOB_vs_S2_pcc2, label="D5 Gd3_autop6 EXTR GLOB - S2 edep 0 - pcc 2", color='blue', linestyle='--', marker='x')
    plt.plot(S2_edep0_setQfiss_pcc2.BU, delta_keff_Gd3_autop6_NOEX_NODI_GLOB_vs_S2_pcc2, label="D5 Gd3_autop6 NOEX GLOB - S2 edep 0 - pcc 2", color='green', linestyle='--', marker='D')
    plt.plot(S2_edep0_setQfiss_pcc2.BU, delta_keff_Gd4_autop6_EXTR_NODI_GLOB_vs_S2_pcc2, label="D5 Gd4_autop6 EXTR GLOB - S2 edep 0 - pcc 2", color='red', linestyle='--', marker='x')
    plt.plot(S2_edep0_setQfiss_pcc2.BU, delta_keff_Gd4_autop6_NOEX_NODI_GLOB_vs_S2_pcc2, label="D5 Gd4_autop6 NOEX GLOB - S2 edep 0 - pcc 2", color='purple', linestyle='--', marker='D')
    plt.axhline(y=300, color='red', linestyle='--')
    plt.axhline(y=-300, color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("$\\Delta$ keff (pcm)")
    plt.title("$\\Delta$ keff between D5 Gd3_autop6/Gd4_autop6 and S2 edepmode 0 pcc 2")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_AT10_45Gd}/delta_keff_D5_Gd3_autop6_Gd4_autop6_vs_S2_edepmode0_pcc2.png")
    plt.close()

    # Look at Gd3_autop6 and Gd4_autop6 cases vs S2 Gd2 : pcc1
    # interpolate D5 keffs to S2 BU points
    interp_D5_keff_Gd3_autop6_NOEX_NODI_GLOB = np.interp(S2_edep0_setQfiss_pcc1_Gd2.BU, D5_45Gd_Gd3_autop6_NOEX_NODI_GLOB.BU, D5_45Gd_Gd3_autop6_NOEX_NODI_GLOB.keff)
    interp_D5_keff_Gd3_autop6_EXTR_NODI_GLOB = np.interp(S2_edep0_setQfiss_pcc1_Gd2.BU, D5_45Gd_Gd3_autop6_EXTR_NODI_GLOB.BU, D5_45Gd_Gd3_autop6_EXTR_NODI_GLOB.keff)
    interp_D5_keff_Gd3_autop6_NOEX_NODI_NOGL = np.interp(S2_edep0_setQfiss_pcc1_Gd2.BU, D5_45Gd_Gd3_autop6_NOEX_NODI_NOGL.BU, D5_45Gd_Gd3_autop6_NOEX_NODI_NOGL.keff)
    interp_D5_keff_Gd3_autop6_EXTR_NODI_NOGL = np.interp(S2_edep0_setQfiss_pcc1_Gd2.BU, D5_45Gd_Gd3_autop6_EXTR_NODI_NOGL.BU, D5_45Gd_Gd3_autop6_EXTR_NODI_NOGL.keff)

    # D5 Gd3_autop6 vs S2 edepmode 0 pcc 1
    delta_keff_Gd3_autop6_NOEX_NODI_GLOB_vs_S2_pcc1 = (interp_D5_keff_Gd3_autop6_NOEX_NODI_GLOB - S2_edep0_setQfiss_pcc1_Gd2.keff) * 1e5 # pcm
    delta_keff_Gd3_autop6_EXTR_NODI_GLOB_vs_S2_pcc1 = (interp_D5_keff_Gd3_autop6_EXTR_NODI_GLOB - S2_edep0_setQfiss_pcc1_Gd2.keff) * 1e5 # pcm
    delta_keff_Gd3_autop6_NOEX_NODI_NOGL_vs_S2_pcc1 = (interp_D5_keff_Gd3_autop6_NOEX_NODI_NOGL - S2_edep0_setQfiss_pcc1_Gd2.keff) * 1e5 # pcm
    delta_keff_Gd3_autop6_EXTR_NODI_NOGL_vs_S2_pcc1 = (interp_D5_keff_Gd3_autop6_EXTR_NODI_NOGL - S2_edep0_setQfiss_pcc1_Gd2.keff) * 1e5 # pcm
    
    # D5 Gd4_autop6 vs S2 edepmode 0 pcc 1
    interp_D5_keff_Gd4_autop6_NOEX_NODI_GLOB = np.interp(S2_edep0_setQfiss_pcc1_Gd2.BU, D5_45Gd_Gd4_autop6_NOEX_NODI_GLOB.BU, D5_45Gd_Gd4_autop6_NOEX_NODI_GLOB.keff)
    interp_D5_keff_Gd4_autop6_EXTR_NODI_GLOB = np.interp(S2_edep0_setQfiss_pcc1_Gd2.BU, D5_45Gd_Gd4_autop6_EXTR_NODI_GLOB.BU, D5_45Gd_Gd4_autop6_EXTR_NODI_GLOB.keff)
    interp_D5_keff_Gd4_autop6_NOEX_NODI_NOGL = np.interp(S2_edep0_setQfiss_pcc1_Gd2.BU, D5_45Gd_Gd4_autop6_NOEX_NODI_NOGL.BU, D5_45Gd_Gd4_autop6_NOEX_NODI_NOGL.keff)
    interp_D5_keff_Gd4_autop6_EXTR_NODI_NOGL = np.interp(S2_edep0_setQfiss_pcc1_Gd2.BU, D5_45Gd_Gd4_autop6_EXTR_NODI_NOGL.BU, D5_45Gd_Gd4_autop6_EXTR_NODI_NOGL.keff)
    
    # D5 Gd4_autop6 vs S2 edepmode 0 pcc 1
    delta_keff_Gd4_autop6_NOEX_NODI_GLOB_vs_S2_pcc1 = (interp_D5_keff_Gd4_autop6_NOEX_NODI_GLOB - S2_edep0_setQfiss_pcc1_Gd2.keff) * 1e5 # pcm
    delta_keff_Gd4_autop6_EXTR_NODI_GLOB_vs_S2_pcc1 = (interp_D5_keff_Gd4_autop6_EXTR_NODI_GLOB - S2_edep0_setQfiss_pcc1_Gd2.keff) * 1e5 # pcm
    delta_keff_Gd4_autop6_NOEX_NODI_NOGL_vs_S2_pcc1 = (interp_D5_keff_Gd4_autop6_NOEX_NODI_NOGL - S2_edep0_setQfiss_pcc1_Gd2.keff) * 1e5 # pcm
    delta_keff_Gd4_autop6_EXTR_NODI_NOGL_vs_S2_pcc1 = (interp_D5_keff_Gd4_autop6_EXTR_NODI_NOGL - S2_edep0_setQfiss_pcc1_Gd2.keff) * 1e5 # pcm
    # D5 Gd4_autop6 vs S2 edepmode 0 pcc 2
    delta_keff_Gd4_autop6_NOEX_NODI_GLOB_vs_S2_pcc2 = (interp_D5_keff_Gd4_autop6_NOEX_NODI_GLOB - S2_edep0_setQfiss_pcc1_Gd2.keff) * 1e5 # pcm
    delta_keff_Gd4_autop6_EXTR_NODI_GLOB_vs_S2_pcc2 = (interp_D5_keff_Gd4_autop6_EXTR_NODI_GLOB - S2_edep0_setQfiss_pcc1_Gd2.keff) * 1e5 # pcm
    delta_keff_Gd4_autop6_NOEX_NODI_NOGL_vs_S2_pcc2 = (interp_D5_keff_Gd4_autop6_NOEX_NODI_NOGL - S2_edep0_setQfiss_pcc1_Gd2.keff) * 1e5 # pcm
    delta_keff_Gd4_autop6_EXTR_NODI_NOGL_vs_S2_pcc2 = (interp_D5_keff_Gd4_autop6_EXTR_NODI_NOGL - S2_edep0_setQfiss_pcc1_Gd2.keff) * 1e5 # pcm

    # plot results : Gd3_autop6/Gd4_autop3 vs S2 edepmode 0 pcc 1 NOEX vs EXTR / GLOB
    plt.figure(figsize=(10, 6))
    plt.plot(S2_edep0_setQfiss_pcc1_Gd2.BU, delta_keff_Gd3_autop6_EXTR_NODI_GLOB_vs_S2_pcc1, label="D5 Gd3_autop6 EXTR GLOB - S2 edep 0 - pcc 1", color='blue', linestyle='--', marker='x')
    plt.plot(S2_edep0_setQfiss_pcc1_Gd2.BU, delta_keff_Gd3_autop6_NOEX_NODI_GLOB_vs_S2_pcc1, label="D5 Gd3_autop6 NOEX GLOB - S2 edep 0 - pcc 1", color='green', linestyle='--', marker='D')
    plt.plot(S2_edep0_setQfiss_pcc1_Gd2.BU, delta_keff_Gd4_autop6_EXTR_NODI_GLOB_vs_S2_pcc1, label="D5 Gd4_autop6 EXTR GLOB - S2 edep 0 - pcc 1", color='red', linestyle='--', marker='x')
    plt.plot(S2_edep0_setQfiss_pcc1_Gd2.BU, delta_keff_Gd4_autop6_NOEX_NODI_GLOB_vs_S2_pcc1, label="D5 Gd4_autop6 NOEX GLOB - S2 edep 0 - pcc 1", color='purple', linestyle='--', marker='D')
    plt.axhline(y=300, color='red', linestyle='--')
    plt.axhline(y=-300, color='red', linestyle='--')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("$\\Delta$ keff (pcm)")
    plt.title("$\\Delta$ keff between D5 Gd3_autop6/Gd4_autop6 and S2 edepmode 0 pcc 1 / Gd2")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_AT10_45Gd}/delta_keff_D5_Gd3_autop6_Gd4_autop6_vs_S2_edepmode0_pcc1_Gd2.png")
    plt.close()


if post_treat_case3_1:
    """
    Post treat case 3.1 :
    Compare AT10_45Gd for pcc BU in Dragon5
    """

    path_to_PYGAN_results = f"{os.getcwd()}/PYGAN_COMPOS_path/AT10_45Gd_Cst_pow_evol_results/"
    cwd_path = os.getcwd()
    os.chdir(path_to_PYGAN_results)

    list_of_extr_options = [("NOEX","NOEX"), ("EXTR","EXTR"), ("NOEX","EXTR"), ("EXTR","NOEX")] # type of exrtapolation used at (predictor,corrector) step
    BU_lists_to_test = ["Gd0_BOC", "Gd1_BOC", "Gd2_BOC", "Gd3_BOC", "Gd4_BOC", "Gd5_BOC", "Gdf_BOC"]

    list_D5_PT_obj = []
    # Loop over the options
    for name_BUlist in BU_lists_to_test:
        for extr_opt in list_of_extr_options:
            name_compo = f"CPO_endfb8r1_295_NG0_PT_N_SALT_{name_BUlist}_KAPS_pcc{extr_opt[0]}_{extr_opt[1]}_GLOB"
            print(f"Processing {name_compo}")
            pyCOMPO = lcm.new('LCM_INP', name_compo, impx=0)
            case = D5_case(pyCOMPO,
                            dlib_name = "endfb8r1_295_NG0",
                            bu_points = name_BUlist,
                            ssh_opt = "PT",
                            correlation = "NOCORR",
                            sat = "",
                            depl_sol = "KAPS",
                            tracked_nuclides = tracked_nuclides,
                            BU_lists = getLists(name_BUlist),
                            save_dir = save_dir_AT10_45Gd
                        )
            print(f"for BU points : {name_BUlist}, compo points are {getLists(name_BUlist)['COMPO']} with length {len(getLists(name_BUlist)['COMPO'])}")
            case.set_BUscheme("predictor-corrector", f"{extr_opt[0]}-{extr_opt[1]}")
            list_D5_PT_obj.append(case)
    os.chdir(cwd_path)

    plt.figure(figsize=(10, 6))
    for D5_PCC_case in list_D5_PT_obj:
        plt.plot(D5_PCC_case.BU, D5_PCC_case.keff, label=f"{D5_PCC_case.BUScheme} {D5_PCC_case.extrapolation_type}", linestyle='--', marker='x')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("keff")
    plt.title("D5 keff Predictor-Corrector")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_AT10_45Gd}/D5_PCC_keff.png")
    plt.close()

    plt.figure(figsize=(10, 6))
    for D5_PCC_case in list_D5_PT_obj:
        plt.plot(D5_PCC_case.BU, D5_PCC_case.DRAGON_ISOTOPESDENS["Gd157"], label=f"{D5_PCC_case.BUScheme} {D5_PCC_case.extrapolation_type}", linestyle='--', marker='x')
    plt.xlabel("Burnup (MWd/tU)")
    plt.ylabel("Gd157 density (#/b-cm)")
    plt.title("D5 Gd157 density Predictor-Corrector")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_AT10_45Gd}/D5_PCC_Gd157_density.png")
    plt.close()

    # compute relative difference to Gd0_BOC
    # Gd0_BOC is the reference case, for each EXTRAPOLATION type
    # NOEX + NOEX
    for extr_type in ["NOEX-NOEX", "EXTR-EXTR", "NOEX-EXTR", "EXTR-NOEX"]:
        Gd0_BOC_case = [case for case in list_D5_PT_obj if case.bu_points == "Gd0_BOC" and case.extrapolation_type == extr_type][0]
        print(f"Gd0_BOC case : {Gd0_BOC_case}")
        
        plt.figure(figsize=(10, 6))
        # cases to compare : all other GdX_BOC with NOEX-NOEX
        for BU_pts in ["Gd1_BOC", "Gd2_BOC", "Gd3_BOC", "Gd4_BOC", "Gd5_BOC", "Gdf_BOC"]:
            
            GdX_BOC_case = [case for case in list_D5_PT_obj if case.bu_points == BU_pts and case.extrapolation_type == extr_type][0]

            delta_keff = (GdX_BOC_case.keff - Gd0_BOC_case.keff) * 1e5 # pcm

            plt.plot(GdX_BOC_case.BU, delta_keff, label=f"{GdX_BOC_case.bu_points} {GdX_BOC_case.extrapolation_type}", linestyle='--', marker='x')
        plt.axhline(y=300, color='red', linestyle='--')
        plt.axhline(y=-300, color='red', linestyle='--')
        plt.xlabel("Burnup (MWd/tU)")
        plt.ylabel("$\\Delta$ keff (pcm)")
        plt.title(f"$\\Delta$ keff between D5 GdX_BOC and D5 Gd0_BOC {extr_type}")
        plt.legend()
        plt.grid()
        plt.savefig(f"{save_dir_AT10_45Gd}/delta_keff_D5_GdX_BOC_vs_D5_Gd0_BOC_{extr_type}.png")
        plt.close()

        
    for extr_type in ["NOEX-NOEX", "EXTR-EXTR", "NOEX-EXTR", "EXTR-NOEX"]:
        Gd0_BOC_case = [case for case in list_D5_PT_obj if case.bu_points == "Gd0_BOC" and case.extrapolation_type == extr_type][0]
        print(f"Gd0_BOC case : {Gd0_BOC_case}")
        
        plt.figure(figsize=(10, 6))
        # cases to compare : all other GdX_BOC with NOEX-NOEX
        for BU_pts in ["Gd1_BOC", "Gd2_BOC", "Gd3_BOC", "Gd4_BOC", "Gd5_BOC", "Gdf_BOC"]:
            
            GdX_BOC_case = [case for case in list_D5_PT_obj if case.bu_points == BU_pts and case.extrapolation_type == extr_type][0]

            delta_NGd157 = (GdX_BOC_case.DRAGON_ISOTOPESDENS["Gd157"] - Gd0_BOC_case.DRAGON_ISOTOPESDENS["Gd157"]) * 100 / Gd0_BOC_case.DRAGON_ISOTOPESDENS["Gd157"] 

            plt.plot(GdX_BOC_case.BU, delta_NGd157, label=f"{GdX_BOC_case.bu_points} {GdX_BOC_case.extrapolation_type}", linestyle='--', marker='x')
        plt.axhline(y=2, color='red', linestyle='--')
        plt.axhline(y=-2, color='red', linestyle='--')
        plt.xlabel("Burnup (MWd/tU)")
        plt.ylabel("$\\Delta$ Gd157 density (%)")
        plt.title(f"$\\Delta$ Gd157 density between D5 GdX_BOC and D5 Gd0_BOC {extr_type}")
        plt.legend()
        plt.grid()
        plt.savefig(f"{save_dir_AT10_45Gd}/delta_Gd157_D5_GdX_BOC_vs_D5_Gd0_BOC_{extr_type}.png")
        plt.close()


    ## look at pcc results on Gd2
    list_D5_PT_obj_Gd2 = []
    os.chdir(path_to_PYGAN_results)
    for list_bu in ["Gd2_autop6", "Gd2_autop7"]:
        for extr_opt in list_of_extr_options:
                name_compo = f"CPO_endfb8r1_295_NG0_PT_N_SALT_{list_bu}_KAPS_pcc{extr_opt[0]}_{extr_opt[1]}_GLOB"
                print(f"Processing {name_compo}")
                pyCOMPO = lcm.new('LCM_INP', name_compo, impx=0)
                case = D5_case(pyCOMPO,
                                dlib_name = "endfb8r1_295_NG0",
                                bu_points = list_bu,
                                ssh_opt = "PT",
                                correlation = "NOCORR",
                                sat = "",
                                depl_sol = "KAPS",
                                tracked_nuclides = tracked_nuclides,
                                BU_lists = getLists(list_bu),
                                save_dir = save_dir_AT10_45Gd
                            )
                #print(f"for BU points : {name_BUlist}, compo points are {getLists(name_BUlist)['COMPO']} with length {len(getLists(name_BUlist)['COMPO'])}")
                case.set_BUscheme("predictor-corrector", f"{extr_opt[0]}-{extr_opt[1]}")
                list_D5_PT_obj_Gd2.append(case)
    os.chdir(cwd_path)

    S2_edep0_setQfiss_pcc0 = S2_case(case_name = "AT10_45Gd_BUGd2",
                                    lib_name = "endfb8r1_pynjoy2012_kerma",
                                    edep_id = 0, areQfissSet = True, isEcaptSet = False,
                                    pcc_id = 0, specific_power = 26.5, tracked_nuclides = tracked_nuclides, save_dir = save_dir_AT10_45Gd)
    
    S2_edep0_setQfiss_pcc1 = S2_case(case_name = "AT10_45Gd_BUGd2",
                                lib_name = "endfb8r1_pynjoy2012_kerma",
                                edep_id = 0, areQfissSet = True, isEcaptSet = False,
                                pcc_id = 1, specific_power = 26.5, tracked_nuclides = tracked_nuclides, save_dir = save_dir_AT10_45Gd)
    
    S2_edep0_setQfiss_pcc2 = S2_case(case_name = "AT10_45Gd_BUGd2",
                                lib_name = "endfb8r1_pynjoy2012_kerma",
                                edep_id = 0, areQfissSet = True, isEcaptSet = False,
                                pcc_id = 2, specific_power = 26.5, tracked_nuclides = tracked_nuclides, save_dir = save_dir_AT10_45Gd)
    
    for D5_test_case in list_D5_PT_obj_Gd2:
        # interpolate D5 keffs to S2 BU points
        interp_D5_keff_pcc0 = np.interp(S2_edep0_setQfiss_pcc0.BU, D5_test_case.BU, D5_test_case.keff)
        delta_keff_pcc0 = (interp_D5_keff_pcc0 - S2_edep0_setQfiss_pcc0.keff) * 1e5
        print("Sanity check on BU points")
        print((S2_edep0_setQfiss_pcc0.BU-D5_test_case.EVO_BU_steps)*100/D5_test_case.EVO_BU_steps)
        print((S2_edep0_setQfiss_pcc1.BU-D5_test_case.EVO_BU_steps)*100/D5_test_case.EVO_BU_steps)
        print((S2_edep0_setQfiss_pcc2.BU-D5_test_case.EVO_BU_steps)*100/D5_test_case.EVO_BU_steps)

        interp_D5_keff_pcc1 = np.interp(S2_edep0_setQfiss_pcc1.BU, D5_test_case.BU, D5_test_case.keff)
        delta_keff_pcc1 = (interp_D5_keff_pcc1 - S2_edep0_setQfiss_pcc1.keff) * 1e5
        
        interp_D5_keff_pcc2 = np.interp(S2_edep0_setQfiss_pcc2.BU, D5_test_case.BU, D5_test_case.keff)
        delta_keff_pcc2 = (interp_D5_keff_pcc2 - S2_edep0_setQfiss_pcc2.keff) * 1e5

        

        # plot results : delta keff
        plt.figure(figsize=(10, 6))
        plt.plot(S2_edep0_setQfiss_pcc0.BU, delta_keff_pcc0, label=f"{D5_test_case.extrapolation_type} - S2 edep 0 - pcc 0", linestyle='--', marker='x')
        plt.plot(S2_edep0_setQfiss_pcc1.BU, delta_keff_pcc1, label=f"{D5_test_case.extrapolation_type} - S2 edep 0 - pcc 1", linestyle='--', marker='D')
        plt.plot(S2_edep0_setQfiss_pcc2.BU, delta_keff_pcc2, label=f"{D5_test_case.extrapolation_type} - S2 edep 0 - pcc 2", linestyle='--', marker='o')
        plt.axhline(y=300, color='red', linestyle='--')
        plt.axhline(y=-300, color='red', linestyle='--')
        plt.xlabel("Burnup (MWd/tU)")
        plt.ylabel("$\\Delta$ keff (pcm)")
        plt.title(f"$\\Delta$ keff between D5 PCC {D5_test_case.extrapolation_type} Gd2_autop6 and S2 edepmode 0")
        plt.legend()
        plt.grid()
        plt.savefig(f"{save_dir_AT10_45Gd}/delta_keff_D5_pcc_{D5_test_case.bu_points}_vs_S2_edepmode0_{D5_test_case.extrapolation_type}.png")
        plt.close()


        for iso in tracked_nuclides:

            interp_D5_Niso_pcc0 = np.interp(S2_edep0_setQfiss_pcc0.BU, D5_test_case.BU, D5_test_case.DRAGON_ISOTOPESDENS[iso])
            delta_Niso_pcc0 = [(interp_D5_Niso_pcc0[i] - S2_edep0_setQfiss_pcc0.Ni[iso][i]) * 100 / S2_edep0_setQfiss_pcc0.Ni[iso][i]
                               if S2_edep0_setQfiss_pcc0.Ni[iso][i] != 0 else 0
                               for i in range(len(S2_edep0_setQfiss_pcc0.Ni[iso]))]

            interp_D5_Niso_pcc1 = np.interp(S2_edep0_setQfiss_pcc1.BU, D5_test_case.BU, D5_test_case.DRAGON_ISOTOPESDENS[iso])
            delta_Niso_pcc1 = [(interp_D5_Niso_pcc1[i] - S2_edep0_setQfiss_pcc1.Ni[iso][i]) * 100 / S2_edep0_setQfiss_pcc1.Ni[iso][i]
                                 if S2_edep0_setQfiss_pcc1.Ni[iso][i] != 0 else 0
                                 for i in range(len(S2_edep0_setQfiss_pcc1.Ni[iso]))]

            interp_D5_Niso_pcc2 = np.interp(S2_edep0_setQfiss_pcc2.BU, D5_test_case.BU, D5_test_case.DRAGON_ISOTOPESDENS[iso])
            delta_Niso_pcc2 = [(interp_D5_Niso_pcc2[i] - S2_edep0_setQfiss_pcc2.Ni[iso][i]) * 100 / S2_edep0_setQfiss_pcc2.Ni[iso][i]
                                    if S2_edep0_setQfiss_pcc2.Ni[iso][i] != 0 else 0
                                    for i in range(len(S2_edep0_setQfiss_pcc2.Ni[iso]))]
            # plot results : delta NGd157
            plt.figure(figsize=(10, 6))
            plt.plot(S2_edep0_setQfiss_pcc0.BU, delta_Niso_pcc0, label=f"{D5_test_case.extrapolation_type} - S2 edep 0 - pcc 0", linestyle='--', marker='x')
            plt.plot(S2_edep0_setQfiss_pcc1.BU, delta_Niso_pcc1, label=f"{D5_test_case.extrapolation_type} - S2 edep 0 - pcc 1", linestyle='--', marker='D')
            plt.plot(S2_edep0_setQfiss_pcc2.BU, delta_Niso_pcc2, label=f"{D5_test_case.extrapolation_type} - S2 edep 0 - pcc 2", linestyle='--', marker='o')
            plt.axhline(y=2, color='red', linestyle='--')
            plt.axhline(y=-2, color='red', linestyle='--')
            plt.xlabel("Burnup (MWd/tU)")
            plt.ylabel(f"$\\Delta$ {iso} density (%)")
            plt.title(f"$\\Delta$ {iso} density between D5 PCC {D5_test_case.extrapolation_type} Gd2_autop6 and S2 edepmode 0")
            plt.legend()
            plt.grid()
            plt.savefig(f"{save_dir_AT10_45Gd}/delta_{iso}_D5_pcc_{D5_test_case.bu_points}_vs_S2_edepmode0_{D5_test_case.extrapolation_type}.png")
            plt.close()



