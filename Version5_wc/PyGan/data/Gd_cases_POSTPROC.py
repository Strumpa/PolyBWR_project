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



def post_treat_AT10_45Gd_vs_S2(evaluation, draglib_name, self_shielding_method, correlation, evo_solver, time_integrators, burnup_points, tracked_nuclides, save_dir_AT10_45Gd):
    """
    Post-treatment for the AT10_45Gd case
    compare DRAGON5 and SERPENT2 results for :
    - time integrator in time_integrators
    - calculation performed on burnup_points
    using the following options :
    - evaluation : endfb8r1
    - draglib_name : evaluation + draglib specific (ie mesh, or kermas or whatever)
    - self-shielding method
    - correlation
    - evolution solver

    - tracked_nuclides : list of isotopes to be tracked

    - save_dir_AT10_45Gd : directory where the results will be saved

    HARD CODING : origin of S2 ace file is pynjoy2012_kerma
    """
    pcc0 = True
    pcc1 = True
    pcc2 = True
    origin_of_S2_data = "pynjoy2012_kerma"
    if correlation == "NOCORR":
        correlation = "N"
    elif correlation == "CORR":
        correlation = "C"

    if draglib_name == "endfb8r1_295":
        D5_kerma_used = False
    else:
        D5_kerma_used = True
    
    if burnup_points == "Gd_autop3" or burnup_points == "Gd_autop4" or burnup_points == "Gd":
        S2_edep0_qfiss_pcc0 = S2_case(case_name = "AT10_45Gd",
                                lib_name = f"{evaluation}_{origin_of_S2_data}",
                                edep_id = 0,
                                areQfissSet = True,
                                isEcaptSet = False,
                                pcc_id = 0,
                                specific_power = 26.5,
                                tracked_nuclides = tracked_nuclides,
                                save_dir = save_dir_AT10_45Gd)
        
        S2_edep0_qfiss_pcc1 = S2_case(case_name = "AT10_45Gd",
                                lib_name = f"{evaluation}_{origin_of_S2_data}",
                                edep_id = 0,
                                areQfissSet = True,
                                isEcaptSet = False,
                                pcc_id = 1,
                                specific_power = 26.5,
                                tracked_nuclides = tracked_nuclides,
                                save_dir = save_dir_AT10_45Gd)
        
        S2_edep0_qfiss_pcc2 = S2_case(case_name = "AT10_45Gd",
                                lib_name =  f"{evaluation}_{origin_of_S2_data}",
                                edep_id = 0,
                                areQfissSet = True,
                                isEcaptSet = False,
                                pcc_id = 2,
                                specific_power = 26.5,
                                tracked_nuclides = tracked_nuclides,
                                save_dir = save_dir_AT10_45Gd)
        
        S2_edep2_pcc1 = S2_case(case_name = "AT10_45Gd",
                                lib_name = f"{evaluation}_{origin_of_S2_data}",
                                edep_id = 2,
                                areQfissSet = False,
                                isEcaptSet = False,
                                pcc_id = 1,
                                specific_power = 26.5,
                                tracked_nuclides = tracked_nuclides,
                                save_dir = save_dir_AT10_45Gd)
    else:
        s2_bu_id = burnup_points.split("_")[0]
        if pcc0:
            S2_edep0_qfiss_pcc0 = S2_case(case_name = f"AT10_45Gd_BU{s2_bu_id}",
                                    lib_name = f"{evaluation}_{origin_of_S2_data}",
                                    edep_id = 0,
                                    areQfissSet = True,
                                    isEcaptSet = False,
                                    pcc_id = 0,
                                    specific_power = 26.5,
                                    tracked_nuclides = tracked_nuclides,
                                    save_dir = save_dir_AT10_45Gd)
        if pcc1:
            S2_edep0_qfiss_pcc1 = S2_case(case_name = f"AT10_45Gd_BU{s2_bu_id}",
                                    lib_name = f"{evaluation}_{origin_of_S2_data}",
                                    edep_id = 0,
                                    areQfissSet = True,
                                    isEcaptSet = False,
                                    pcc_id = 1,
                                    specific_power = 26.5,
                                    tracked_nuclides = tracked_nuclides,
                                    save_dir = save_dir_AT10_45Gd)
        if pcc2:
            S2_edep0_qfiss_pcc2 = S2_case(case_name = f"AT10_45Gd_BU{s2_bu_id}",
                                    lib_name = f"{evaluation}_{origin_of_S2_data}",
                                    edep_id = 0,
                                    areQfissSet = True,
                                    isEcaptSet = False,
                                    pcc_id = 2,
                                    specific_power = 26.5,
                                    tracked_nuclides = tracked_nuclides,
                                    save_dir = save_dir_AT10_45Gd)
        """
        S2_edep2_pcc1 = S2_case(case_name = f"AT10_45Gd_BU{s2_bu_id}",
                                lib_name = f"{evaluation}_{origin_of_S2_data}",
                                edep_id = 2,
                                areQfissSet = False,
                                isEcaptSet = False,
                                pcc_id = 1,
                                specific_power = 26.5,
                                tracked_nuclides = tracked_nuclides,
                                save_dir = save_dir_AT10_45Gd)
        """
    # It was not noticed that saturation option played a huge role in the results
    path = os.getcwd()
    path_to_PyGan_results = f"PYGAN_COMPOS_path/AT10_45Gd_Cst_pow_evol_results/"
    os.chdir(path_to_PyGan_results)
    D5_result_cases = []
    D5_NG0_result_cases = []
    for time_integrator in time_integrators:
        if D5_kerma_used:
            CPO_NAME = f"CPO_{draglib_name}_{self_shielding_method}_{correlation}_SALT_{evo_solver}_GLOB_{time_integrator}_{burnup_points}"
            CPO = lcm.new('LCM_INP', CPO_NAME, impx=0)
            D5_NOM_case = D5_case(pyCOMPO = CPO,
                            dlib_name = draglib_name,
                            bu_points = burnup_points,
                            ssh_opt = self_shielding_method,
                            correlation = correlation,
                            sat = "NODI",
                            depl_sol = evo_solver,
                            tracked_nuclides = tracked_nuclides,
                            BU_lists = getLists(burnup_points),
                            save_dir = save_dir_AT10_45Gd)
            if time_integrator == "CECM" or time_integrator == "CECE":
                scheme = "predictor-corrector" 
            else:
                scheme = "predictor"
            D5_NOM_case.set_BUscheme(scheme, time_integrator)
            D5_result_cases.append(D5_case)
        else:
            CPO_NG0_NAME =  f"CPO_{draglib_name}_NG0_{self_shielding_method}_{correlation}_SALT_{burnup_points}_{evo_solver}_NODI_{time_integrator}_GLOB"
            #CPO_NAME =  f"CPO_{draglib_name}_{self_shielding_method}_{correlation}_SALT_{evo_solver}_NODI_GLOB_{time_integrator}_{burnup_points}"
            print(CPO_NG0_NAME)
            CPO_NG0 = lcm.new('LCM_INP', CPO_NG0_NAME, impx=0)
            #CPO = lcm.new('LCM_INP', CPO_NAME, impx=0)
            D5_NG0_case = D5_case(pyCOMPO = CPO_NG0,
                            dlib_name = f"{draglib_name}_NG0",
                            bu_points = burnup_points,
                            ssh_opt = self_shielding_method,
                            correlation = correlation,
                            sat = "NODI",
                            depl_sol = evo_solver,
                            tracked_nuclides = tracked_nuclides,
                            BU_lists = getLists(burnup_points),
                            save_dir = save_dir_AT10_45Gd)
            """
            D5_NOM_case = D5_case(pyCOMPO = CPO,
                            dlib_name = draglib_name,
                            bu_points = burnup_points,
                            ssh_opt = self_shielding_method,
                            correlation = correlation,
                            sat = "NODI",
                            depl_sol = evo_solver,
                            tracked_nuclides = tracked_nuclides,
                            BU_lists = getLists(burnup_points),
                            save_dir = save_dir_AT10_45Gd)
            """
            
            if time_integrator == "CECM" or time_integrator == "CECE":
                scheme = "predictor-corrector" 
            else:
                scheme = "predictor"
            #D5_NOM_case.set_BUscheme(scheme, time_integrator)
            D5_NG0_case.set_BUscheme(scheme, time_integrator)
            #D5_result_cases.append(D5_NOM_case)
            D5_NG0_result_cases.append(D5_NG0_case)
    os.chdir(path)
    if D5_NG0_result_cases:
        comparisonD5_NG0_S2_edep0_qfiss_pcc0 = multiD5S2(f"D5 {self_shielding_method}-{correlation} - S2 : modified DEPL vs edep0 with set qfiss, pcc0", D5_NG0_result_cases, S2_edep0_qfiss_pcc0, burnup_points, tracked_nuclides, f"{save_dir_AT10_45Gd}/{self_shielding_method}_{correlation}")
        comparisonD5_NG0_S2_edep0_qfiss_pcc1 = multiD5S2(f"D5 {self_shielding_method}-{correlation} - S2 : modified DEPL vs edep0 with set qfiss, pcc1", D5_NG0_result_cases, S2_edep0_qfiss_pcc1, burnup_points, tracked_nuclides, f"{save_dir_AT10_45Gd}/{self_shielding_method}_{correlation}")
        comparisonD5_NG0_S2_edep0_qfiss_pcc2 = multiD5S2(f"D5 {self_shielding_method}-{correlation} - S2 : modified DEPL vs edep0 with set qfiss, pcc2", D5_NG0_result_cases, S2_edep0_qfiss_pcc2, burnup_points, tracked_nuclides, f"{save_dir_AT10_45Gd}/{self_shielding_method}_{correlation}")


        comparisonD5_NG0_S2_edep0_qfiss_pcc0.compare_keffs()
        comparisonD5_NG0_S2_edep0_qfiss_pcc1.compare_keffs()
        comparisonD5_NG0_S2_edep0_qfiss_pcc2.compare_keffs()

        comparisonD5_NG0_S2_edep0_qfiss_pcc0.plot_delta_keff()
        comparisonD5_NG0_S2_edep0_qfiss_pcc1.plot_delta_keff()
        comparisonD5_NG0_S2_edep0_qfiss_pcc2.plot_delta_keff()

        comparisonD5_NG0_S2_edep0_qfiss_pcc0.compare_Ni()
        comparisonD5_NG0_S2_edep0_qfiss_pcc1.compare_Ni()
        comparisonD5_NG0_S2_edep0_qfiss_pcc2.compare_Ni()

        comparisonD5_NG0_S2_edep0_qfiss_pcc0.plot_delta_Ni()
        comparisonD5_NG0_S2_edep0_qfiss_pcc1.plot_delta_Ni()
        comparisonD5_NG0_S2_edep0_qfiss_pcc2.plot_delta_Ni()

    


if __name__ == "__main__":
    # evaluation and origin of S2 data

    evaluation = "endfb8r1"
    origin_of_S2_data = "pynjoy2012_kerma"

    # isotopes to be tracked 
    tracked_nuclides = ["U235","U238","Pu239","Pu240","Pu241","Pu242","Am241","Xe135","Sm149","Gd155","Gd157"]
    # time integrators to be used
    time_integrators = ["EXTR", "NOEX", "EXTR2"] # ,"CECM"]
    # burnup points to be used
    burnup_points = "Gd_autop3"

    # build the path to the results directory
    path=os.getcwd()
    save_dir_HOM_Gd157_VBOC = f"{path}/Gd_cases_POSTPROC/HOM_Gd157_VBOC"
    save_dir_HOM_UOX_Gd157 = f"{path}/Gd_cases_POSTPROC/HOM_UOX_Gd157"
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
    post_treat_case4 = False

    #post_treat_AT10_45Gd_vs_S2(evaluation, "endfb8r1_295", "RSE", "C", "KAPS", ["CECM"], "free", tracked_nuclides, save_dir_AT10_45Gd)
    #post_treat_AT10_45Gd_vs_S2(evaluation, "endfb8r1_295", "RSE", "C", "KAPS", ["CECM"], "free", tracked_nuclides, save_dir_AT10_45Gd)
    #post_treat_AT10_45Gd_vs_S2(evaluation, "endfb8r1_295", "PT", "N", "KAPS", ["CECM"], "free", tracked_nuclides, save_dir_AT10_45Gd)
    #post_treat_AT10_45Gd_vs_S2(evaluation, "endfb8r1_295", "PT", "N", "KAPS", ["CECM"], "free", tracked_nuclides, save_dir_AT10_45Gd)

    if post_treat_case3:
        ## --- DRAGON5-SERPENT2 comparison : post treatment
        # Gd_autop3 BU
        post_treat_AT10_45Gd_vs_S2(evaluation, "endfb8r1_295", "RSE", "C", "KAPS", time_integrators, "Gd_autop3", tracked_nuclides, save_dir_AT10_45Gd)
        post_treat_AT10_45Gd_vs_S2(evaluation, "endfb8r1_295", "RSE", "N", "KAPS", time_integrators, "Gd_autop3", tracked_nuclides, save_dir_AT10_45Gd)
        post_treat_AT10_45Gd_vs_S2(evaluation, "endfb8r1_295", "PT", "C", "KAPS", time_integrators, "Gd_autop3", tracked_nuclides, save_dir_AT10_45Gd)
        post_treat_AT10_45Gd_vs_S2(evaluation, "endfb8r1_295", "PT", "N", "KAPS", time_integrators, "Gd_autop3", tracked_nuclides, save_dir_AT10_45Gd)

        # Gd2_autop6 BU
        post_treat_AT10_45Gd_vs_S2(evaluation, "endfb8r1_295", "RSE", "C", "KAPS", time_integrators, "Gd2_autop6", tracked_nuclides, save_dir_AT10_45Gd)
        post_treat_AT10_45Gd_vs_S2(evaluation, "endfb8r1_295", "RSE", "N", "KAPS", time_integrators, "Gd2_autop6", tracked_nuclides, save_dir_AT10_45Gd)
        post_treat_AT10_45Gd_vs_S2(evaluation, "endfb8r1_295", "PT", "C", "KAPS", time_integrators, "Gd2_autop6", tracked_nuclides, save_dir_AT10_45Gd)
        post_treat_AT10_45Gd_vs_S2(evaluation, "endfb8r1_295", "PT", "N", "KAPS", time_integrators, "Gd2_autop6", tracked_nuclides, save_dir_AT10_45Gd)