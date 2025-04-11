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
post_treat_case0 = False
# Case 1 : HOM_Gd157_VBOC, focus on this
post_treat_case1 = False
# Case 2 : HOM_UOX_Gd157, focus on this
post_treat_case2 = False
# Case 3 : AT10_45Gd, focus on this
post_treat_case3 = True
# Case 4 : gduo2_295_kec1


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
    plt.plot(BU_points_plot, delta_keff_NG0_NOEX_DIRA_HOM_Gd157_VBOC_pcc0, label="D5 NG0 NOEX DIRA - S2 edepmode 0 - pcc 0", color='red', linestyle='--', marker='x')
    plt.plot(BU_points_plot, delta_keff_NG0_EXTR_NODI_HOM_Gd157_VBOC_pcc2, label="D5 NG0 EXTR NODI - S2 edepmode 0 - pcc 2", color='green', linestyle='--', marker='D')
    plt.plot(BU_points_plot, delta_keff_NG0_EXTR_DIRA_HOM_Gd157_VBOC_pcc2, label="D5 NG0 EXTR DIRA - S2 edepmode 0 - pcc 2", color='orange', linestyle='--', marker='D')
    plt.plot(BU_points_plot, delta_keff_NG0_NOEX_SAT_HOM_Gd157_VBOC_pcc0, label="D5 NG0 NOEX SAT - S2 edepmode 0 - pcc 0", color='purple', linestyle='--', marker='x')
    plt.plot(BU_points_plot, delta_keff_NG0_EXTR_SAT_HOM_Gd157_VBOC_pcc2, label="D5 NG0 EXTR SAT - S2 edepmode 0 - pcc 2", color='brown', linestyle='--', marker='D')
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
        delta_Niso_NOEX_NODI = [(D5_HOM_Gd157_VBOC_NG0_NOEX_NODI.DRAGON_ISOTOPESDENS[iso][idx] - S2_edep0_setQfiss_pcc0.Ni[iso][idx]) * 100 / S2_edep0_setQfiss_pcc0.Ni[iso][idx]
                        if S2_edep0_setQfiss_pcc0.Ni[iso][idx] != 0 else 0
                        for idx in range(len(S2_edep0_setQfiss_pcc0.Ni[iso]))]
        delta_Niso_EXTR_NODI = [(D5_HOM_Gd157_VBOC_NG0_EXTR_NODI.DRAGON_ISOTOPESDENS[iso][idx] - S2_edep0_setQfiss_pcc2.Ni[iso][idx]) * 100 / S2_edep0_setQfiss_pcc2.Ni[iso][idx]
                        if S2_edep0_setQfiss_pcc2.Ni[iso][idx] != 0 else 0
                        for idx in range(len(S2_edep0_setQfiss_pcc2.Ni[iso]))]
        plt.figure(figsize=(10, 6))
        plt.plot(BU_points_plot, delta_Niso_NOEX_NODI, label=f"D5 NG0 NOEX NODI - S2 edepmode 0 - pcc 0 - {iso}", color='blue', linestyle='--', marker='x')
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
    plt.plot(BU_points_plot, delta_keff_D5_NG0_EXTR_vs_NOEX_Case14, label="D5 NG0 EXTR vs NOEX : VBOC_finer2", color='orange', linestyle='--', marker='D')
    plt.plot(BU_points_plot, delta_keff_S2_pcc2_vs_pcc0, label="S2 edepmode 0 pcc2 vs pcc0 : VBOC", color='red', linestyle='--', marker='o')
    plt.plot(BU_points_plot_finerBU, delta_keff_S2_pcc2_vs_pcc0_finerBU, label="S2 edepmode 0 pcc2 vs pcc0 : VBOC_finerBU", color='purple', linestyle='--', marker='o')
    plt.plot(BU_points_plot_finer2, delta_keff_S2_pcc2_vs_pcc0_finer2, label="S2 edepmode 0 pcc2 vs pcc0 : VBOC_finer2", color='brown', linestyle='--', marker='o')
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
    plt.plot(BU_points_plot_finer2, delta_keff_VBOC_finer2_VBOC, label="S2 edepmode 0 pcc0 VBOC - VBOC_finer2", color='blue', linestyle='--', marker='x')
    plt.plot(BU_points_plot_finer2, delta_keff_VBOC_finer2_VBOC_finerBU, label="S2 edepmode 0 pcc0 VBOC_finerBU - VBOC_finer2", color='green', linestyle='--', marker='D')
    plt.plot(BU_points_plot_finer2, delta_keff_VBOC_finer2_VBOC_finerBU2, label="S2 edepmode 0 pcc0 VBOC_finerBU2 - VBOC_finer2", color='purple', linestyle='--', marker='o')
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
    plt.plot(BU_points_plot, delta_keff_VBOC_VBOC_finer2_pcc0, label="S2 VBOC - VBOC_finer 2 (pcc0)", color='blue', linestyle='--', marker='x')
    plt.plot(BU_points_plot, delta_keff_VBOC_finerBU_VBOC_finer2_pcc0, label="S2 VBOC_finerBU - VBOC_finer 2 (pcc0)", color='red', linestyle='--', marker='x')
    plt.plot(BU_points_plot, delta_keff_VBOC_finerBU2_VBOC_finer2_pcc0, label="S2 VBOC_finerBU2 - VBOC_finer 2 (pcc0)", color='green', linestyle='--', marker='D')
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
    plt.plot(BU_points_plot, delta_keff_VBOC_VBOC_finer2_pcc1, label="S2 VBOC - VBOC_finer 2 (pcc1)", color='green', linestyle='--', marker='D')
    plt.plot(BU_points_plot, delta_keff_VBOC_finerBU_VBOC_finer2_pcc1, label="S2 VBOC_finerBU - VBOC_finer 2 (pcc1)", color='purple', linestyle='--', marker='D')
    plt.plot(BU_points_plot, delta_keff_VBOC_finerBU2_VBOC_finer2_pcc1, label="S2 VBOC_finerBU2 - VBOC_finer 2 (pcc1)", color='orange', linestyle='--', marker='D')
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
    plt.plot(BU_points_plot, delta_keff_VBOC_VBOC_finer2_pcc2, label="S2 VBOC - VBOC_finer2 (pcc2)", color='orange', linestyle='--', marker='D')
    plt.plot(BU_points_plot, delta_keff_VBOC_finerBU_VBOC_finer2_pcc2, label="S2 VBOC_finerBU - VBOC_finer2 (pcc2)", color='brown', linestyle='--', marker='D')
    plt.plot(BU_points_plot, delta_keff_VBOC_finerBU2_VBOC_finer2_pcc2, label="S2 VBOC_finerBU2 - VBOC_finer2 (pcc2)", color='purple', linestyle='--', marker='D')
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


if post_treat_case3:
    print("Post-treatment of case 3 : AT10_45Gd_Cst_pow_evol")

    name_CPO_45Gd_Gdautop3_EXTR = "CPO_endfb8r1_295_NG0_RSE_SALT_Gd_autop3_RUNG_NODI_EXTR_GLOB"
    name_CPO_45Gd_Gdautop3_NOEX = "CPO_endfb8r1_295_NG0_RSE_SALT_Gd_autop3_RUNG_NODI_NOEX_GLOB"

    path_to_PYGAN_results = f"{os.getcwd()}/PYGAN_COMPOS_path/AT10_45Gd_Cst_pow_evol_results/"
    cwd_path = os.getcwd()
    os.chdir(path_to_PYGAN_results)
    # Load the data
    CPO_45Gd_Gdautop3_EXTR = lcm.new('LCM_INP', name_CPO_45Gd_Gdautop3_EXTR, impx=0)
    CPO_45Gd_Gdautop3_NOEX = lcm.new('LCM_INP', name_CPO_45Gd_Gdautop3_NOEX, impx=0)
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
    
    print(f"Gd_autop3 BU : {getLists('Gd_autop3')['BU']}, with len : {len(getLists('Gd_autop3')['BU'])}")
    print(f"Gd_autop3 COMPO : {getLists('Gd_autop3')['COMPO']}, with len : {len(getLists('Gd_autop3')['COMPO'])}")
    print(f"delta COMPO-BU = {np.array(getLists('Gd_autop3')['COMPO']) - np.array(getLists('Gd_autop3')['BU'])}")
    
    
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
    plt.savefig(f"{save_dir_AT10_45Gd}/delta_keff_D5_NOEX_EXTR_vs_S2_edepmode0_pcc0.png")
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
    plt.savefig(f"{save_dir_AT10_45Gd}/delta_keff_D5_NOEX_EXTR_vs_S2_edepmode0_pcc1.png")
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
