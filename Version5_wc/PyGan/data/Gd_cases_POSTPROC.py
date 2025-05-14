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



def post_treat_Gd_case_vs_S2_edep_pcc(case, evaluation, draglib_name, self_shielding_method, correlation, 
                               evo_solver, time_integrators, burnup_points, tracked_nuclides, 
                               S2_edep, S2_pcc,
                               save_dir_case):
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

    - S2_edep : edepmode parameter for S2 comparison
    - S2_pcc : pcc parameter for S2 comparison

    - save_dir_AT10_45Gd : directory where the results will be saved

    HARD CODING : origin of S2 ace file is pynjoy2012_kerma
    """
    if case == "AT10_45Gd_Cst_pow_evol":
        S2_case_name = "AT10_45Gd"
    elif case == "gduo2_295_kec1_evol":
        S2_case_name = "gduo2_295_kec1"

    if S2_edep == 0:
        set_qfiss = True
    else:
        set_qfiss = False
    origin_of_S2_data = "pynjoy2012_kerma"
    if correlation == "NOCORR":
        correlation = "N"
    elif correlation == "CORR":
        correlation = "C"

    if draglib_name == "endfb8r1_295":
        D5_kerma_used = False
    else:
        D5_kerma_used = True
    
    if burnup_points == "Gd_autop3" or burnup_points == "Gd_autop4" or burnup_points == "Gd" or burnup_points == "Gd_autop5":
        S2_ref_case = S2_case(case_name = S2_case_name,
                                lib_name = f"{evaluation}_{origin_of_S2_data}",
                                edep_id = S2_edep,
                                areQfissSet = set_qfiss,
                                isEcaptSet = False,
                                pcc_id = S2_pcc,
                                specific_power = 26.5,
                                tracked_nuclides = tracked_nuclides,
                                save_dir = save_dir_case)
    
        
    else:
        s2_bu_id = burnup_points.split("_")[0]
        
        S2_ref_case = S2_case(case_name = f"{S2_case_name}_BU{s2_bu_id}",
                                lib_name = f"{evaluation}_{origin_of_S2_data}",
                                edep_id = S2_edep,
                                areQfissSet = set_qfiss,
                                isEcaptSet = False,
                                pcc_id = S2_pcc,
                                specific_power = 26.5,
                                tracked_nuclides = tracked_nuclides,
                                save_dir = save_dir_case)
    
    # It was noticed that DIRA / NODI saturation didnt heavily impact the results for AT10_45Gd
    path = os.getcwd()
    path_to_PyGan_results = f"PYGAN_COMPOS_path/{case}_results/"
    os.chdir(path_to_PyGan_results)
    D5_result_cases = []
    D5_NG0_result_cases = []
    for time_integrator in time_integrators:
        if D5_kerma_used:
            CPO_NAME = f"CPO_{draglib_name}_{self_shielding_method}_{correlation}_SALT_{evo_solver}_NODI_GLOB_{time_integrator}_{burnup_points}"
            print(CPO_NAME)
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
                            save_dir = save_dir_case)
            if time_integrator == "CECM" or time_integrator == "CECE":
                scheme = "predictor-corrector" 
            else:
                scheme = "predictor"
            D5_NOM_case.set_BUscheme(scheme, time_integrator)
            D5_result_cases.append(D5_NOM_case)
        else:
            CPO_NG0_NAME =  f"CPO_{draglib_name}_NG0_{self_shielding_method}_{correlation}_SALT_{burnup_points}_{evo_solver}_NODI_{time_integrator}_GLOB"
            CPO_NAME =  f"CPO_{draglib_name}_{self_shielding_method}_{correlation}_SALT_{burnup_points}_{evo_solver}_NODI_{time_integrator}_GLOB"
            print(CPO_NG0_NAME)
            CPO_NG0 = lcm.new('LCM_INP', CPO_NG0_NAME, impx=0)
            CPO = lcm.new('LCM_INP', CPO_NAME, impx=0)
            D5_NG0_case = D5_case(pyCOMPO = CPO_NG0,
                            dlib_name = f"{draglib_name}_NG0",
                            bu_points = burnup_points,
                            ssh_opt = self_shielding_method,
                            correlation = correlation,
                            sat = "NODI",
                            depl_sol = evo_solver,
                            tracked_nuclides = tracked_nuclides,
                            BU_lists = getLists(burnup_points),
                            save_dir = save_dir_case)
            
            D5_NOM_case = D5_case(pyCOMPO = CPO,
                            dlib_name = draglib_name,
                            bu_points = burnup_points,
                            ssh_opt = self_shielding_method,
                            correlation = correlation,
                            sat = "NODI",
                            depl_sol = evo_solver,
                            tracked_nuclides = tracked_nuclides,
                            BU_lists = getLists(burnup_points),
                            save_dir = save_dir_case)

            
            if time_integrator == "CECM" or time_integrator == "CECE":
                scheme = "predictor-corrector" 
            else:
                scheme = "predictor"
            D5_NOM_case.set_BUscheme(scheme, time_integrator)
            D5_NG0_case.set_BUscheme(scheme, time_integrator)
            D5_result_cases.append(D5_NOM_case)
            D5_NG0_result_cases.append(D5_NG0_case)
    os.chdir(path)
    if D5_NG0_result_cases:
        comparisonD5_NG0_S2 = multiD5S2(f"D5 {self_shielding_method}-{correlation} - S2 : modified DEPL vs edep0 with set qfiss, pcc {S2_pcc}", D5_NG0_result_cases, S2_ref_case, burnup_points, tracked_nuclides, f"{save_dir_case}/{draglib_name}/{self_shielding_method}_{correlation}")

        comparisonD5_NG0_S2.compare_keffs()
        comparisonD5_NG0_S2.plot_delta_keff()

        comparisonD5_NG0_S2.compare_Ni()
        comparisonD5_NG0_S2.plot_delta_Ni()
    else:
        comparisonD5_S2 = multiD5S2(f"D5 {self_shielding_method}-{correlation} - S2 : edep {S2_edep}, pcc {S2_pcc}", D5_result_cases, S2_ref_case, burnup_points, tracked_nuclides, f"{save_dir_case}/{draglib_name}/{self_shielding_method}_{correlation}")
        comparisonD5_S2.compare_keffs()
        comparisonD5_S2.plot_delta_keff()
        comparisonD5_S2.compare_Ni()
        comparisonD5_S2.plot_delta_Ni()

def custom_plots(save_dir_case):
    """
    Custom plots for the post-treatment
    """
    case = "gduo2_295_kec1_evol"
    burnup_points = "Gd_autop3"
    tracked_nuclides = ["U235","U238","Pu239","Pu240","Pu241","Pu242","Am241","Xe135","Sm149","Gd155","Gd157"]
    ssh_method = "PT"
    correlation = "C"
    time_integrator = "EXTR"
    draglib_name = "endfb8r1_295_NG0"
    evo_solver = "KAPS"
    path = os.getcwd()
    path_to_PyGan_results = f"PYGAN_COMPOS_path/{case}_results/"
    os.chdir(path_to_PyGan_results)
    cpo_name = f"CPO_{draglib_name}_{ssh_method}_{correlation}_SALT_{burnup_points}_{evo_solver}_NODI_{time_integrator}_GLOB"
    print(f"cpo_name : {cpo_name}")
    CPO = lcm.new('LCM_INP', cpo_name, impx=0)
    os.chdir(path)

    D5_case_custom = D5_case(pyCOMPO = CPO,
                            dlib_name = draglib_name,
                            bu_points = burnup_points,
                            ssh_opt = ssh_method,
                            correlation = correlation,
                            sat = "NODI",
                            depl_sol = evo_solver,
                            tracked_nuclides = tracked_nuclides,
                            BU_lists = getLists(burnup_points),
                            save_dir = save_dir_case)
    D5_case_custom.set_BUscheme("Predicor", time_integrator)
    
    S2_edep0_pcc0 = S2_case(case_name = "gduo2_295_kec1",
                                lib_name = "endfb8r1_pynjoy2012_kerma",
                                edep_id = 0,
                                areQfissSet = True,
                                isEcaptSet = False,
                                pcc_id = 0,
                                specific_power = 26.5,
                                tracked_nuclides = tracked_nuclides,
                                save_dir = save_dir_case)
    

    S2_edep0_pcc1 = S2_case(case_name = "gduo2_295_kec1",
                                lib_name = "endfb8r1_pynjoy2012_kerma",
                                edep_id = 0,
                                areQfissSet = True,
                                isEcaptSet = False,
                                pcc_id = 1,
                                specific_power = 26.5,
                                tracked_nuclides = tracked_nuclides,
                                save_dir = save_dir_case)
    
    S2_edep0_pcc2 = S2_case(case_name = "gduo2_295_kec1",
                                lib_name = "endfb8r1_pynjoy2012_kerma",
                                edep_id = 0,
                                areQfissSet = True,
                                isEcaptSet = False,
                                pcc_id = 2,
                                specific_power = 26.5,
                                tracked_nuclides = tracked_nuclides,
                                save_dir = save_dir_case)
    
    ### Plot D5_custom_case vs S2_edep0_pcc0

    comparisonD5_S2 = multiD5S2(f"gduo2_295_kec1_NG0_EXTR_S2_edep0_pcc0", [D5_case_custom], S2_edep0_pcc0, burnup_points, tracked_nuclides, f"{save_dir_case}/{draglib_name}/{ssh_method}_{correlation}")
    comparisonD5_S2.compare_keffs()
    comparisonD5_S2.plot_delta_keff()
    comparisonD5_S2.compare_Ni()
    comparisonD5_S2.plot_delta_Ni()

    ### Plot D5_custom_case vs S2_edep0_pcc1
    comparisonD5_S2 = multiD5S2(f"gduo2_295_kec1_NG0_EXTR_S2_edep0_pcc1", [D5_case_custom], S2_edep0_pcc1, burnup_points, tracked_nuclides, f"{save_dir_case}/{draglib_name}/{ssh_method}_{correlation}")
    comparisonD5_S2.compare_keffs()
    comparisonD5_S2.plot_delta_keff()
    comparisonD5_S2.compare_Ni()
    comparisonD5_S2.plot_delta_Ni()

    ### Plot D5_custom_case vs S2_edep0_pcc2
    comparisonD5_S2 = multiD5S2(f"gduo2_295_kec1_NG0_EXTR_S2_edep0_pcc2", [D5_case_custom], S2_edep0_pcc2, burnup_points, tracked_nuclides, f"{save_dir_case}/{draglib_name}/{ssh_method}_{correlation}")
    comparisonD5_S2.compare_keffs()
    comparisonD5_S2.plot_delta_keff()
    comparisonD5_S2.compare_Ni()
    comparisonD5_S2.plot_delta_Ni()


    


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

    ssh_methods = ["RSE","PT"]
    correlations = ["N", "C"]

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
    post_treat_case3 = False
    # Case 4 : gduo2_295_kec1
    post_treat_case4 = False
    custom_plots(save_dir_gduo2_295_kec1)

    if post_treat_case3:

        ## --- DRAGON5-SERPENT2 comparison : post treatment
        
        for pcc in [0, 1, 2]:
           # Gd BU
            post_treat_Gd_case_vs_S2_edep_pcc("AT10_45Gd_Cst_pow_evol", evaluation, "endfb8r1_295", "PT", "N", "KAPS", time_integrators, "Gd", tracked_nuclides,  0, pcc, save_dir_AT10_45Gd)
           
            # Gd_autop3 BU
           
            for ssh_method in ssh_methods:
                for correlation in correlations:
                    post_treat_Gd_case_vs_S2_edep_pcc("AT10_45Gd_Cst_pow_evol",evaluation, "endfb8r1_295", ssh_method, correlation, "KAPS", time_integrators, "Gd_autop3", tracked_nuclides, 0, pcc, save_dir_AT10_45Gd)

            # Gd2_autop6 BU
            for ssh_method in ssh_methods:
                for correlation in correlations:
                    post_treat_Gd_case_vs_S2_edep_pcc("AT10_45Gd_Cst_pow_evol", evaluation, "endfb8r1_295", ssh_method, correlation, "KAPS", time_integrators, "Gd2_autop6", tracked_nuclides, 0, pcc, save_dir_AT10_45Gd)
    
        # Gd BU
        post_treat_Gd_case_vs_S2_edep_pcc("AT10_45Gd_Cst_pow_evol", evaluation, "endfb81295K", "PT", "N", "KAPS", time_integrators, "Gd", tracked_nuclides,  2, 1, save_dir_AT10_45Gd)
        
        # Gd_autop3 BU
        #post_treat_AT10_45Gd_vs_S2_edep_pcc(evaluation, "endfb81295K", "RSE", "C", "KAPS", time_integrators, "Gd_autop3", tracked_nuclides, 2, 1, save_dir_AT10_45Gd)
        #post_treat_AT10_45Gd_vs_S2_edep_pcc(evaluation, "endfb81295K", "RSE", "N", "KAPS", time_integrators, "Gd_autop3", tracked_nuclides,  2, 1, save_dir_AT10_45Gd)
        #post_treat_AT10_45Gd_vs_S2_edep_pcc(evaluation, "endfb81295K", "PT", "C", "KAPS", time_integrators, "Gd_autop3", tracked_nuclides,  2, 1, save_dir_AT10_45Gd)
        post_treat_Gd_case_vs_S2_edep_pcc("AT10_45Gd_Cst_pow_evol", evaluation, "endfb81295K", "PT", "N", "KAPS", time_integrators, "Gd_autop3", tracked_nuclides,  2, 1, save_dir_AT10_45Gd)

        # Gd2_autop6 BU
        #post_treat_AT10_45Gd_vs_S2_edep_pcc(evaluation, "endfb81295K", "RSE", "C", "KAPS", time_integrators, "Gd2_autop6", tracked_nuclides, 2, 1, save_dir_AT10_45Gd)
        #post_treat_AT10_45Gd_vs_S2_edep_pcc(evaluation, "endfb81295K", "RSE", "N", "KAPS", time_integrators, "Gd2_autop6", tracked_nuclides, 2, 1, save_dir_AT10_45Gd)
        #post_treat_AT10_45Gd_vs_S2_edep_pcc(evaluation, "endfb81295K", "PT", "C", "KAPS", time_integrators, "Gd2_autop6", tracked_nuclides, 2, 1, save_dir_AT10_45Gd)
        post_treat_Gd_case_vs_S2_edep_pcc("AT10_45Gd_Cst_pow_evol", evaluation, "endfb81295K", "PT", "N", "KAPS", time_integrators, "Gd2_autop6", tracked_nuclides, 2, 1, save_dir_AT10_45Gd)

        
        # Gd BU
        post_treat_Gd_case_vs_S2_edep_pcc("AT10_45Gd_Cst_pow_evol", evaluation, "endfb81295K2", "PT", "N", "KAPS", time_integrators, "Gd_autop3", tracked_nuclides,  2, 1, save_dir_AT10_45Gd)

        # Gd_autop3 BU
        #post_treat_AT10_45Gd_vs_S2_edep_pcc(evaluation, "endfb81295K2", "RSE", "C", "KAPS", time_integrators, "Gd_autop3", tracked_nuclides, 2, 1, save_dir_AT10_45Gd)
        #post_treat_AT10_45Gd_vs_S2_edep_pcc(evaluation, "endfb81295K2", "RSE", "N", "KAPS", time_integrators, "Gd_autop3", tracked_nuclides,  2, 1, save_dir_AT10_45Gd)
        #post_treat_AT10_45Gd_vs_S2_edep_pcc(evaluation, "endfb81295K2", "PT", "C", "KAPS", time_integrators, "Gd_autop3", tracked_nuclides,  2, 1, save_dir_AT10_45Gd)
        post_treat_Gd_case_vs_S2_edep_pcc("AT10_45Gd_Cst_pow_evol", evaluation, "endfb81295K2", "PT", "N", "KAPS", time_integrators, "Gd_autop3", tracked_nuclides,  2, 1, save_dir_AT10_45Gd)

        # Gd2_autop6 BU
        #post_treat_AT10_45Gd_vs_S2_edep_pcc(evaluation, "endfb81295K2", "RSE", "C", "KAPS", time_integrators, "Gd2_autop6", tracked_nuclides, 2, 1, save_dir_AT10_45Gd)
        #post_treat_AT10_45Gd_vs_S2_edep_pcc(evaluation, "endfb81295K2", "RSE", "N", "KAPS", time_integrators, "Gd2_autop6", tracked_nuclides, 2, 1, save_dir_AT10_45Gd)
        #post_treat_AT10_45Gd_vs_S2_edep_pcc(evaluation, "endfb81295K2", "PT", "C", "KAPS", time_integrators, "Gd2_autop6", tracked_nuclides, 2, 1, save_dir_AT10_45Gd)
        post_treat_Gd_case_vs_S2_edep_pcc("AT10_45Gd_Cst_pow_evol", evaluation, "endfb81295K2", "PT", "N", "KAPS", time_integrators, "Gd2_autop6", tracked_nuclides, 2, 1, save_dir_AT10_45Gd)

    if post_treat_case4:
        
        # gduo2_295_kec1
        # edepmode = 0
        time_integrators = ["EXTR"]
        for ssh_method in ["RSE", "PT"]:
            for correlation in ["N", "C"]:
                for bu_pts in ["Gd_autop3", "Gd_autop4", "Gd_autop5"]:
                    post_treat_Gd_case_vs_S2_edep_pcc("gduo2_295_kec1_evol", evaluation, "endfb8r1_295", ssh_method, correlation, "KAPS", time_integrators, bu_pts, tracked_nuclides, 0, 1, save_dir_gduo2_295_kec1)
                    post_treat_Gd_case_vs_S2_edep_pcc("gduo2_295_kec1_evol", evaluation, "endfb8r1_295", ssh_method, correlation, "KAPS", time_integrators, bu_pts, tracked_nuclides, 0, 2, save_dir_gduo2_295_kec1)

        # edepmode = 2 
        for ssh_method in ["RSE", "PT"]:
            for correlation in ["N", "C"]:
                for bu_pts in ["Gd_autop3", "Gd_autop4", "Gd_autop5"]:

                    post_treat_Gd_case_vs_S2_edep_pcc("gduo2_295_kec1_evol", evaluation, "endfb81295K", ssh_method, correlation, "KAPS", time_integrators, bu_pts, tracked_nuclides, 2, 1, save_dir_gduo2_295_kec1)
                    post_treat_Gd_case_vs_S2_edep_pcc("gduo2_295_kec1_evol", evaluation, "endfb81295K", ssh_method, correlation, "KAPS", time_integrators, bu_pts, tracked_nuclides, 2, 2, save_dir_gduo2_295_kec1) 
        
        

                    post_treat_Gd_case_vs_S2_edep_pcc("gduo2_295_kec1_evol", evaluation, "endfb81295K2", ssh_method, correlation, "KAPS", time_integrators, bu_pts, tracked_nuclides, 2, 1, save_dir_gduo2_295_kec1)
                    post_treat_Gd_case_vs_S2_edep_pcc("gduo2_295_kec1_evol", evaluation, "endfb81295K2", ssh_method, correlation, "KAPS", time_integrators, bu_pts, tracked_nuclides, 2, 2, save_dir_gduo2_295_kec1) 
        
        
