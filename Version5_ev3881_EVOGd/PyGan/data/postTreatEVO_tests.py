## PyGan script to post treat cases in burnup evolution
# The DRAGON5 EVO: module is responsible for computing burnup evolution.
# System considered : Homogeenized pincell with Gd157
# Comparing DRAGON5 cases : before and after specific normalization modifications

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import lifo 
import lcm
import cle2000 
from getBULists import getBULists
from postproc_cst_pow_evol import DRAGON_case as D5_case
from postproc_cst_pow_evol import Serpent2_case as S2_case



def create_D5_comparison(modif_case, ref_case):
    """
    Create a comparison between a modified case and a reference case in DRAGON5.
    """
    keff_modif = modif_case.keff
    keff_ref = ref_case.keff

    # Create a DataFrame to hold the comparison results
    comparison_df = pd.DataFrame({
        'BU Points': modif_case.BU,
        'delta_keff': (keff_modif - keff_ref)*1e5, # Convert to pcm
    })
    for nuclide in modif_case.tracked_nuclides:
        comparison_df[f'delta {nuclide}'] = modif_case.DRAGON_ISOTOPESDENS[nuclide] - ref_case.DRAGON_ISOTOPESDENS[nuclide]
        comparison_df[f'delta {nuclide} (%)'] = [(modif_case.DRAGON_ISOTOPESDENS[nuclide][i] - ref_case.DRAGON_ISOTOPESDENS[nuclide][i])*100/ref_case.DRAGON_ISOTOPESDENS[nuclide][i] 
                                                if ref_case.DRAGON_ISOTOPESDENS[nuclide][i] != 0 else 0
                                                for i in range(len(modif_case.DRAGON_ISOTOPESDENS[nuclide]))]
        
    return comparison_df

def create_D5_S2_comparison(D5_case, S2_case):
    """
    Create a comparison between DRAGON5 cases and a Serpent2 case.
    """
    # Initialize a DataFrame to hold the comparison results
    if len(S2_case.BU) == len(D5_case.BU):
        comparison_df = pd.DataFrame({
            'BU Points': D5_case.BU,
            'delta_keff': (D5_case.keff - S2_case.keff) * 1e5 ,
        })

        # Add isotopic comparisons
        for nuclide in D5_case.tracked_nuclides:
            comparison_df[f'delta {nuclide}'] = D5_case.DRAGON_ISOTOPESDENS[nuclide] - S2_case.Ni[nuclide]
            comparison_df[f'delta {nuclide} (%)'] = [(D5_case.DRAGON_ISOTOPESDENS[nuclide][i] - S2_case.Ni[nuclide][i])*100/S2_case.Ni[nuclide][i] 
                                                    if S2_case.Ni[nuclide][i] != 0 else 0
                                                    for i in range(len(S2_case.Ni[nuclide]))]
            
    else:
        # interpolate D5 values on available S2 burnup points
        interp_keff = np.interp(S2_case.BU, D5_case.BU, D5_case.keff)
        comparison_df = pd.DataFrame({
            'BU Points': S2_case.BU,
            'delta_keff': (interp_keff - S2_case.keff) * 1e5 ,
        })
        # Add isotopic comparisons
        for nuclide in D5_case.tracked_nuclides:
            if nuclide in D5_case.DRAGON_ISOTOPESDENS.keys():
                interp_isotope = np.interp(S2_case.BU, D5_case.BU, D5_case.DRAGON_ISOTOPESDENS[nuclide])
                comparison_df[f'delta {nuclide}'] = interp_isotope - S2_case.Ni[nuclide]
                comparison_df[f'delta {nuclide} (%)'] = [(interp_isotope[i] - S2_case.Ni[nuclide][i])*100/S2_case.Ni[nuclide][i] 
                                                        if S2_case.Ni[nuclide][i] != 0 else 0
                                                        for i in range(len(S2_case.Ni[nuclide]))]
            else:
                print(f"Nuclide {nuclide} not found in DRAGON5 case.")

    return comparison_df

def create_D5_OpenMC_comparison(D5_case, OpenMC_case, OpenMC_tracked_nuclides):
    """
    Create a comparison between DRAGON5 cases and an OpenMC case.
    OpenMC case is already a DataFrame.

    """

    comparison_df = pd.DataFrame({
        'BU Points': D5_case.BU,
        'delta_keff': (D5_case.keff - OpenMC_case['Keff']) * 1e5 ,
    })

    for nuclide in OpenMC_tracked_nuclides:
        if nuclide in D5_case.DRAGON_ISOTOPESDENS:
            comparison_df[f'delta {nuclide}'] = D5_case.DRAGON_ISOTOPESDENS[nuclide] - OpenMC_case[nuclide]
            comparison_df[f'delta {nuclide} (%)'] = [(D5_case.DRAGON_ISOTOPESDENS[nuclide][i] - OpenMC_case[nuclide][i])*100/OpenMC_case[nuclide][i] 
                                                    if OpenMC_case[nuclide][i] != 0 else 0
                                                    for i in range(len(OpenMC_case[nuclide]))]
        else:
            print(f"Nuclide {nuclide} not found in DRAGON5 case.")

    return comparison_df


def post_treat_HOM_Gd157_VBOC():
    name_D5_case = "HOM_Gd157_VBOC_evo_tests"
    name_S2_case = "HOM_Gd157_VBOC_EVO"
    # define path to NOGL_HM data
    path_to_EVOmodif_data =  "EVOGd_modif_path"
    path_to_ref_data = "DRAGON_COMPOS_path"
    
    # Define calculation options for the CPOs name reconstruction
    time_integrator = "EXTR"
    BU_points = "VBOC_evo_tests" # 25 burnup points to be used
    tracking_option = "SALT"
    draglib = "endfb8r1_295" # "endfb8r1_295_NG0" #"endfb8r1_295" # "endfb8r1_295", "endfb81295K", "endfb81295K2"
    ssh_option = "PT"
    correlation = "N"
    tracked_nuclides = ["U235", "U238", "Pu239", "Pu240", "Pu241", "Pu242", "Am241", "Gd155", "Gd157", "Gd158", "Xe135", "Sm149"]


    ## 
    list_CPOS_EVO_modif = ["COMPO_HOM_Gd157_VBOC_endfb8r1_295_EXTR_NODI_default",
                           "COMPO_HOM_Gd157_VBOC_endfb8r1_295_NG0_EXTR_NODI_default",
                           "COMPO_HOM_Gd157_VBOC_endfb8r1_295_EXTR_NODI_fissionE_decay",
                           "COMPO_HOM_Gd157_VBOC_endfb8r1_295_EXTR_NODI_fissonE_nodecay",
                           "COMPO_HOM_Gd157_VBOC_endfb8r1_295_EXTR_NODI_noNG_decay",
                           "COMPO_HOM_Gd157_VBOC_endfb8r1_295_EXTR_NODI_noNG_nodecay"]
    
    list_CPOS_wc = ["COMPO_HOM_Gd157_VBOC_endfb8r1_EXTR_NODI",
                    "COMPO_HOM_Gd157_VBOC_endfb8r1_NG0_EXTR_NODI"]

    path = os.getcwd()
    save_dir_case = f"{path}/postTreatEVO_tests_results/{name_D5_case}/{draglib}_D5"
    if not os.path.exists(save_dir_case):
        os.makedirs(save_dir_case)
    save_dir_D5_vs_D5 = f"{save_dir_case}/D5_vs_D5"
    if not os.path.exists(save_dir_D5_vs_D5):
        os.makedirs(save_dir_D5_vs_D5)
    save_dir_D5_vs_S2 = f"{save_dir_case}/D5_vs_S2"
    if not os.path.exists(save_dir_D5_vs_S2):
        os.makedirs(save_dir_D5_vs_S2)
    
    # Load the "reference" wc data, before modifications to source code
    os.chdir(path_to_ref_data)
    CPO_wc_ref = lcm.new('LCM_INP', "COMPO_HOM_Gd157_VBOC_endfb8r1_EXTR_NODI", impx=0)
    CPO_wc_NG0 = lcm.new('LCM_INP', "COMPO_HOM_Gd157_VBOC_endfb8r1_NG0_EXTR_NODI", impx=0)
    os.chdir(path)

    # Create D5 wc cases : not modified in EVO
    D5_case_wc_ref = D5_case(pyCOMPO=CPO_wc_ref, dlib_name=draglib, bu_points=BU_points, 
                                ssh_opt=ssh_option, correlation=correlation, sat="NODI", depl_sol="EXTR", 
                                tracked_nuclides=tracked_nuclides, BU_lists=getBULists(BU_points), save_dir=save_dir_case)
    
    D5_case_wc_NG0 = D5_case(pyCOMPO=CPO_wc_NG0, dlib_name=draglib, bu_points=BU_points, 
                                ssh_opt=ssh_option, correlation=correlation, sat="NODI", depl_sol="EXTR", 
                                tracked_nuclides=tracked_nuclides, BU_lists=getBULists(BU_points), save_dir=save_dir_case)
    
    # Load EVOGd data : testing modifications to the normalisation procedure in EVOSIG
    os.chdir(path_to_EVOmodif_data)
    print(os.getcwd())
    print(os.listdir())
    # _ref : reference case, no modifications to the normalisation procedure
    # _NG0 : no modifications to the normalisation procedure, but using the NG0 library
    # _fissionE_decay : modifications to the normalisation procedure, using fission energy decay energy
    # _fissionE_nodecay : modifications to the normalisation procedure, using fission energy without decay energy
    # _noNG_decay : modifications to the normalisation procedure, no NG reactions contributions, using decay energy
    # _noNG_nodecay : modifications to the normalisation procedure, no NG reactions contributions, no decay energy
    CPO_EVOmodif_ref = lcm.new('LCM_INP', "COMPO_HOM_Gd157_VBOC_endfb8r1_295_EXTR_NODI_default", impx=0)
    CPO_EVOmodif_NG0 = lcm.new('LCM_INP', "COMPO_HOM_Gd157_VBOC_endfb8r1_295_NG0_EXTR_NODI_default", impx=0)

    CPO_EVOmodif_fissionE_decay = lcm.new('LCM_INP', "COMPO_HOM_Gd157_VBOC_endfb8r1_295_EXTR_NODI_fissionE_decay", impx=0)
    CPO_EVOmodif_fissionE_nodecay = lcm.new('LCM_INP', "COMPO_HOM_Gd157_VBOC_endfb8r1_295_EXTR_NODI_fissionE_nodecay", impx=0)
    CPO_EVOmodif_noNG_decay = lcm.new('LCM_INP', "COMPO_HOM_Gd157_VBOC_endfb8r1_295_EXTR_NODI_noNG_decay", impx=0)
    CPO_EVOmodif_noNG_nodecay = lcm.new('LCM_INP', "COMPO_HOM_Gd157_VBOC_endfb8r1_295_EXTR_NODI_noNG_nodecay", impx=0)

    os.chdir(path)

    # Create D5 EVOGd cases : test different normalisations
    # reference case : no modifications, should be the same as the wc case
    D5_case_modif_ref = D5_case(pyCOMPO=CPO_EVOmodif_ref, dlib_name=draglib, bu_points=BU_points, 
                                ssh_opt=ssh_option, correlation=correlation, sat="NODI", depl_sol="EXTR", 
                                tracked_nuclides=tracked_nuclides, BU_lists=getBULists(BU_points), save_dir=save_dir_case)
    
    # NG0 case : no modifications, should be the same as the wc NG0 case
    D5_case_modif_NG0 = D5_case(pyCOMPO=CPO_EVOmodif_NG0, dlib_name=draglib, bu_points=BU_points, 
                                ssh_opt=ssh_option, correlation=correlation, sat="NODI", depl_sol="EXTR", 
                                tracked_nuclides=tracked_nuclides, BU_lists=getBULists(BU_points), save_dir=save_dir_case)
    
    # fissionE_decay case : modifications to the normalisation procedure, using fission energy decay energy
    D5_case_modif_fissionE_decay = D5_case(pyCOMPO=CPO_EVOmodif_fissionE_decay, dlib_name=draglib, bu_points=BU_points, 
                                ssh_opt=ssh_option, correlation=correlation, sat="NODI", depl_sol="EXTR", 
                                tracked_nuclides=tracked_nuclides, BU_lists=getBULists(BU_points), save_dir=save_dir_case)
    
    # fissionE_nodecay case : modifications to the normalisation procedure, using fission energy without decay energy
    D5_case_modif_fissionE_nodecay = D5_case(pyCOMPO=CPO_EVOmodif_fissionE_nodecay, dlib_name=draglib, bu_points=BU_points,
                                ssh_opt=ssh_option, correlation=correlation, sat="NODI", depl_sol="EXTR",
                                tracked_nuclides=tracked_nuclides, BU_lists=getBULists(BU_points), save_dir=save_dir_case)
    
    # noNG_decay case : modifications to the normalisation procedure, no NG reactions contributions, using decay energy
    D5_case_modif_noNG_decay = D5_case(pyCOMPO=CPO_EVOmodif_noNG_decay, dlib_name=draglib, bu_points=BU_points,
                                ssh_opt=ssh_option, correlation=correlation, sat="NODI", depl_sol="EXTR",
                                tracked_nuclides=tracked_nuclides, BU_lists=getBULists(BU_points), save_dir=save_dir_case)
    
    # noNG_nodecay case : modifications to the normalisation procedure, no NG reactions contributions, no decay energy
    D5_case_modif_noNG_nodecay = D5_case(pyCOMPO=CPO_EVOmodif_noNG_nodecay, dlib_name=draglib, bu_points=BU_points,
                                ssh_opt=ssh_option, correlation=correlation, sat="NODI", depl_sol="EXTR",
                                tracked_nuclides=tracked_nuclides, BU_lists=getBULists(BU_points), save_dir=save_dir_case)
    
    
    # Sanity check : if modif_ref and wc_ref are the same, if modif_NG0 and wc_NG0 are the same
    # Create data frames for the comparisons
    comparison_modif_ref_wc = create_D5_comparison(D5_case_modif_ref, D5_case_wc_ref)
    comparison_modif_NG0_wc = create_D5_comparison(D5_case_modif_NG0, D5_case_wc_NG0)

    ## Plot sanity check.
    # Plot delta keff
    plt.figure(figsize=(10, 6))
    plt.plot(comparison_modif_ref_wc['BU Points'], comparison_modif_ref_wc['delta_keff'], label='Modif Ref - WC Ref', marker='o')
    plt.plot(comparison_modif_NG0_wc['BU Points'], comparison_modif_NG0_wc['delta_keff'], label='Modif NG0 - WC NG0', marker='o')
    plt.xlabel('Burnup Points (MWd/kgU)')
    plt.ylabel('Delta keff (pcm)')
    plt.title('Delta keff Comparison')
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_D5_vs_D5}/sanity_check_delta_keff.png")
    plt.close()

    # Compare NG0 with ref for both wc and modif cases
    comparison_wc_NG0_ref = create_D5_comparison(D5_case_wc_NG0, D5_case_wc_ref)
    comparison_modif_NG0_ref = create_D5_comparison(D5_case_modif_NG0, D5_case_modif_ref)

    # Plot delta keff for NG0 vs ref
    plt.figure(figsize=(10, 6))
    plt.plot(comparison_wc_NG0_ref['BU Points'], comparison_wc_NG0_ref['delta_keff'], label='WC NG0 - WC Ref', marker='o')
    plt.plot(comparison_modif_NG0_ref['BU Points'], comparison_modif_NG0_ref['delta_keff'], label='Modif NG0 - Modif Ref', marker='o')
    plt.xlabel('Burnup Points (MWd/kgU)')
    plt.ylabel('Delta keff (pcm)')
    plt.title('NG0 vs Ref Comparison')
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_D5_vs_D5}/NG0_vs_ref_delta_keff.png")
    plt.close()


    # Create Serpent2 cases :
    #S2_case_edep0_no_qfiss = S2_case(case_name=name_S2_case, lib_name="endfb8r1_pynjoy2012_kerma", 
    #                                edep_id=0, areQfissSet=False, isEcaptSet=False, pcc_id=2, 
    #                                specific_power=25.6, tracked_nuclides=tracked_nuclides, save_dir=save_dir_case)
    S2_case_edep0_qfiss = S2_case(case_name=name_S2_case, lib_name="endfb8r1_pynjoy2012_kerma",
                                edep_id=0, areQfissSet=True, isEcaptSet=False, pcc_id=2, 
                                specific_power=25.6, tracked_nuclides=tracked_nuclides, save_dir=save_dir_case)
    

    ## Compare DRAGON5 cases with Serpent2 case
    comp_D5_ref_wc_S2 = create_D5_S2_comparison(D5_case_wc_ref, S2_case_edep0_qfiss)
    comp_D5_NG0_wc_S2 = create_D5_S2_comparison(D5_case_wc_NG0, S2_case_edep0_qfiss)
    comp_D5_modif_ref_S2 = create_D5_S2_comparison(D5_case_modif_ref, S2_case_edep0_qfiss)
    comp_D5_modif_NG0_S2 = create_D5_S2_comparison(D5_case_modif_NG0, S2_case_edep0_qfiss)

    # Plot delta keff for D5 wc ref vs S2
    plt.figure(figsize=(10, 6))
    plt.plot(comp_D5_ref_wc_S2['BU Points'], comp_D5_ref_wc_S2['delta_keff'], label='D5 WC Ref - S2', marker='o')
    plt.plot(comp_D5_NG0_wc_S2['BU Points'], comp_D5_NG0_wc_S2['delta_keff'], label='D5 WC NG0 - S2', marker='o')
    plt.plot(comp_D5_modif_ref_S2['BU Points'], comp_D5_modif_ref_S2['delta_keff'], label='D5 Modif Ref - S2', marker='o')
    plt.plot(comp_D5_modif_NG0_S2['BU Points'], comp_D5_modif_NG0_S2['delta_keff'], label='D5 Modif NG0 - S2', marker='o')
    plt.xlabel('Burnup Points (MWd/kgU)')
    plt.ylabel('Delta keff (pcm)')
    plt.title('DRAGON5 vs Serpent2 Delta keff Comparison')
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_D5_vs_S2}/delta_keff_default_ref_NG0.png")
    plt.close()


    ## Compare DRAGON5 NG0 with DRAGON5_modif NO NG
    comp_D5_modif_noNG_decay_NG0 = create_D5_comparison(D5_case_modif_noNG_decay, D5_case_modif_NG0)
    comp_D5_modif_noNG_nodecay_NG0 = create_D5_comparison(D5_case_modif_noNG_nodecay, D5_case_modif_NG0)

    # Plot delta keff for D5 modif noNG decay vs D5 modif NG0
    plt.figure(figsize=(10, 6))
    plt.plot(comp_D5_modif_noNG_decay_NG0['BU Points'], comp_D5_modif_noNG_decay_NG0['delta_keff'], label='D5 Modif noNG decay - D5 Modif NG0', marker='o')
    plt.plot(comp_D5_modif_noNG_nodecay_NG0['BU Points'], comp_D5_modif_noNG_nodecay_NG0['delta_keff'], label='D5 Modif noNG nodecay - D5 Modif NG0', marker='o') 
    plt.xlabel('Burnup Points (MWd/kgU)')
    plt.ylabel('Delta keff (pcm)')
    plt.title('DRAGON5 Modif noNG vs Modif NG0 Delta keff Comparison')
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_D5_vs_D5}/modif_noNG_vs_modif_NG0_delta_keff.png")
    plt.close()

    ## Compare DRAGON5 fissionE with and without decay energy to NG0
    comp_D5_modif_fissionE_decay_NG0 = create_D5_comparison(D5_case_modif_fissionE_decay, D5_case_modif_NG0)
    comp_D5_modif_fissionE_nodecay_NG0 = create_D5_comparison(D5_case_modif_fissionE_nodecay, D5_case_modif_NG0)
    # Plot delta keff for D5 modif fissionE decay vs D5 modif NG0
    plt.figure(figsize=(10, 6))
    plt.plot(comp_D5_modif_fissionE_decay_NG0['BU Points'], comp_D5_modif_fissionE_decay_NG0['delta_keff'], label='D5 Modif fissionE decay - D5 Modif NG0', marker='o')
    plt.plot(comp_D5_modif_fissionE_nodecay_NG0['BU Points'], comp_D5_modif_fissionE_nodecay_NG0['delta_keff'], label='D5 Modif fissionE nodecay - D5 Modif NG0', marker='o')
    plt.xlabel('Burnup Points (MWd/kgU)')
    plt.ylabel('Delta keff (pcm)')
    plt.title('DRAGON5 Modif fissionE vs Modif NG0 Delta keff Comparison')
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_D5_vs_D5}/modif_fissionE_vs_modif_NG0_delta_keff.png")
    plt.close()

    ## Compare DRAGON5 fissionE and noNG with and without decay energy included to Serpent2 case
    comp_D5_modif_fissionE_decay_S2 = create_D5_S2_comparison(D5_case_modif_fissionE_decay, S2_case_edep0_qfiss)
    comp_D5_modif_fissionE_nodecay_S2 = create_D5_S2_comparison(D5_case_modif_fissionE_nodecay, S2_case_edep0_qfiss)
    comp_D5_modif_noNG_decay_S2 = create_D5_S2_comparison(D5_case_modif_noNG_decay, S2_case_edep0_qfiss)
    comp_D5_modif_noNG_nodecay_S2 = create_D5_S2_comparison(D5_case_modif_noNG_nodecay, S2_case_edep0_qfiss)

    # Plot delta keff for D5 modif fissionE decay vs S2
    plt.figure(figsize=(10, 6))
    plt.plot(comp_D5_modif_fissionE_decay_S2['BU Points'], comp_D5_modif_fissionE_decay_S2['delta_keff'], label='D5 Modif fission energy, decay energy - S2', marker='o')
    plt.plot(comp_D5_modif_fissionE_nodecay_S2['BU Points'], comp_D5_modif_fissionE_nodecay_S2['delta_keff'], label='D5 Modif fission energy, no decay energy - S2', marker='o')
    plt.plot(comp_D5_modif_noNG_decay_S2['BU Points'], comp_D5_modif_noNG_decay_S2['delta_keff'], label='D5 Modif no NG energy, decay energy - S2', marker='o')
    plt.plot(comp_D5_modif_noNG_nodecay_S2['BU Points'], comp_D5_modif_noNG_nodecay_S2['delta_keff'], label='D5 Modif no NG energy, no decay energy - S2', marker='o')
    plt.xlabel('Burnup Points (MWd/kgU)')
    plt.ylabel('Delta keff (pcm)')
    plt.title('DRAGON5 Modif fissionE/noNG vs Serpent2 Delta keff Comparison')
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_D5_vs_S2}/modif_fissionE_noNG_vs_S2_delta_keff.png")
    plt.close()


def post_treat_AT10_45Gd():
    """
    Post-treat the AT10_45Gd case.
    """
    name_D5_case = "AT10_45Gd_Cst_pow_evol"
    # define path to NOGL_HM data
    path_to_EVOmodif_data =  "EVOGd_modif_path_PyGan"
    # Define calculation options for the CPOs name reconstruction
    time_integrator = "EXTR"
    BU_points = "Gd_autop3" # 25 burnup points to be used
    name_S2_case = "AT10_45Gd"
    if BU_points == "Gd2_autop6":
        name_S2_case = "AT10_45Gd_BUGd2"

    tracking_option = "SALT"
    draglib = "endfb8r1_295" # "endfb8r1_295_NG0" #"endfb8r1_295" # "endfb8r1_295", "endfb81295K", "endfb81295K2"
    ssh_option = "PT"
    correlation = "N"
    tracked_nuclides = ["U235", "U238", "Pu239", "Pu240", "Pu241", "Pu242", "Am241", "Gd155", "Gd157", "Xe135", "Sm149"]

    name_CPO_EDEP0_decay = f"CPO_endfb8r1_295_{ssh_option}_{correlation}_{tracking_option}_{BU_points}_KAPS_NODI_{time_integrator}_NOGL"
    #name_CPO_EDEP0_no_decay = f"CPO_endfb8r1_295_{ssh_option}_{correlation}_{tracking_option}_{BU_points}_KAPS_NODI_{time_integrator}_NOGL_no_decay"
    path = os.getcwd()
    save_dir_case = f"{path}/postTreatEVO_tests_results/{name_D5_case}/{draglib}_D5"
    if not os.path.exists(save_dir_case):
        os.makedirs(save_dir_case)
    save_dir_D5_vs_D5 = f"{save_dir_case}/D5_vs_D5_{BU_points}"
    if not os.path.exists(save_dir_D5_vs_D5):
        os.makedirs(save_dir_D5_vs_D5)
    save_dir_D5_vs_S2 = f"{save_dir_case}/D5_vs_S2_{BU_points}"
    if not os.path.exists(save_dir_D5_vs_S2):
        os.makedirs(save_dir_D5_vs_S2)
    
    # Load the "reference" wc data, before modifications to source code
    os.chdir(f"{path_to_EVOmodif_data}/{name_D5_case}_results")
    CPO_EDEP0_decay = lcm.new('LCM_INP', name_CPO_EDEP0_decay, impx=0)
    #CPO_EDEP0_no_decay = lcm.new('LCM_INP', name_CPO_EDEP0_no_decay, impx=0)
    os.chdir(path)


    # Create D5 case 
    D5_case_EDEP0_decay = D5_case(pyCOMPO=CPO_EDEP0_decay, dlib_name=draglib, bu_points=BU_points, 
                                ssh_opt=ssh_option, correlation=correlation, sat="NODI", depl_sol="EXTR", 
                                tracked_nuclides=tracked_nuclides, BU_lists=getBULists(BU_points), save_dir=save_dir_case)
    #D5_case_EDEP0_no_decay = D5_case(pyCOMPO=CPO_EDEP0_no_decay, dlib_name=draglib, bu_points=BU_points, 
    #                            ssh_opt=ssh_option, correlation=correlation, sat="NODI", depl_sol="EXTR", 
    #                            tracked_nuclides=tracked_nuclides, BU_lists=getBULists(BU_points), save_dir=save_dir_case)

    print(f"D5 BU points are : {D5_case_EDEP0_decay.BU}")
    # Create S2 case
    S2_case_edep0_pcc0 = S2_case(case_name=name_S2_case, lib_name="endfb8r1_pynjoy2012_kerma",
                                edep_id=0, areQfissSet=False, isEcaptSet=False, pcc_id=0, 
                                specific_power=26.5, tracked_nuclides=tracked_nuclides, save_dir=save_dir_case)
    S2_case_edep0_pcc1 = S2_case(case_name=name_S2_case, lib_name="endfb8r1_pynjoy2012_kerma",
                                edep_id=0, areQfissSet=False, isEcaptSet=False, pcc_id=1, 
                                specific_power=26.5, tracked_nuclides=tracked_nuclides, save_dir=save_dir_case)
    print(f"S2 edep0 pcc1 BU points are : {S2_case_edep0_pcc1.BU}")
    S2_case_edep0_pcc2 = S2_case(case_name=name_S2_case, lib_name="endfb8r1_pynjoy2012_kerma",
                                edep_id=0, areQfissSet=False, isEcaptSet=False, pcc_id=2, 
                                specific_power=26.5, tracked_nuclides=tracked_nuclides, save_dir=save_dir_case)
    print(f"S2 edep0 pcc2 BU points are : {S2_case_edep0_pcc2.BU}")
    
    ## Check BU points 
    if len(D5_case_EDEP0_decay.BU) != len(S2_case_edep0_pcc1.BU):
        print("Warning: Different number of burnup points in D5 and S2 cases.")
    else:
        print("Same number of burnup points in D5 and S2 cases.")
        difference_BU = np.array(D5_case_EDEP0_decay.BU) - np.array(S2_case_edep0_pcc1.BU)
        if np.all(difference_BU == 0):
            print("Burnup points are the same in D5 and S2 pcc 1 case.")
        else:
            print(difference_BU)
            print("Burnup points are different in D5 and S2 pcc 1 case.")

    if len(D5_case_EDEP0_decay.BU) != len(S2_case_edep0_pcc2.BU):
        print("Warning: Different number of burnup points in D5 and S2 cases.")
    else:
        print("Same number of burnup points in D5 and S2 cases.")
        difference_BU = np.array(D5_case_EDEP0_decay.BU) - np.array(S2_case_edep0_pcc2.BU)
        if np.all(difference_BU == 0):
            print("Burnup points are the same in D5 and S2 pcc 2 case.")
        else:
            print(difference_BU)
            print("Burnup points are different in D5 and S2 pcc 2 case.")

    ## Compare DRAGON5 cases with Serpent2 case
    comp_D5_EDEP0_decay_S2_pcc0 = create_D5_S2_comparison(D5_case_EDEP0_decay, S2_case_edep0_pcc0)
    comp_D5_EDEP0_decay_S2_pcc1 = create_D5_S2_comparison(D5_case_EDEP0_decay, S2_case_edep0_pcc1)
    comp_D5_EDEP0_decay_S2_pcc2 = create_D5_S2_comparison(D5_case_EDEP0_decay, S2_case_edep0_pcc2)
    #comp_D5_EDEP0_no_decay_S2_pcc1 = create_D5_S2_comparison(D5_case_EDEP0_no_decay, S2_case_edep0_pcc1)
    #comp_D5_EDEP0_no_decay_S2_pcc2 = create_D5_S2_comparison(D5_case_EDEP0_no_decay, S2_case_edep0_pcc2)
    # Plot delta keff for D5 EDEP0 decay vs S2
    plt.figure(figsize=(10, 6))
    plt.plot(comp_D5_EDEP0_decay_S2_pcc0['BU Points'], comp_D5_EDEP0_decay_S2_pcc0['delta_keff'], label='D5 EDEP0 + decay energy - S2 edep 0 pcc 0', marker='o')
    plt.plot(comp_D5_EDEP0_decay_S2_pcc1['BU Points'], comp_D5_EDEP0_decay_S2_pcc1['delta_keff'], label='D5 EDEP0 + decay energy - S2 edep 0 pcc 1', marker='x')
    plt.plot(comp_D5_EDEP0_decay_S2_pcc2['BU Points'], comp_D5_EDEP0_decay_S2_pcc2['delta_keff'], label='D5 EDEP0 + decay energy - S2 edep 0 pcc 2', marker='x')
    # plot +/- 300 pcm lines
    plt.axhline(y=300, color='r', linestyle='--', label='+300 pcm')
    plt.axhline(y=-300, color='r', linestyle='--', label='-300 pcm')
    plt.xlabel('Burnup Points (MWd/kgU)')
    plt.ylabel('Delta keff (pcm)')
    plt.title('DRAGON5 EDEP0 + decay energy vs Serpent2 Delta keff Comparison')
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_D5_vs_S2}/D5_EDEP0_decay_vs_S2_delta_keff_{ssh_option}_{correlation}_{time_integrator}.png")
    plt.close()

    for iso in tracked_nuclides:
        plt.figure(figsize=(10, 6))
        plt.plot(comp_D5_EDEP0_decay_S2_pcc0['BU Points'], comp_D5_EDEP0_decay_S2_pcc0[f'delta {iso} (%)'], label=f'D5 EDEP0 + decay energy - S2 edep 0 pcc 0 {iso}', marker='o')
        plt.plot(comp_D5_EDEP0_decay_S2_pcc1['BU Points'], comp_D5_EDEP0_decay_S2_pcc1[f'delta {iso} (%)'], label=f'D5 EDEP0 + decay energy - S2 edep 0 pcc 1 {iso}', marker='x')
        plt.plot(comp_D5_EDEP0_decay_S2_pcc2['BU Points'], comp_D5_EDEP0_decay_S2_pcc2[f'delta {iso} (%)'], label=f'D5 EDEP0 + decay energy - S2 edep 0 pcc 2 {iso}', marker='x')
        # plot +/- 2 % lines
        plt.axhline(y=2, color='r', linestyle='--', label='+2 %')
        plt.axhline(y=-2, color='r', linestyle='--', label='-2 %')
        plt.xlabel('Burnup Points (MWd/kgU)')
        plt.ylabel(f'Delta {iso} (%)')
        plt.title(f'DRAGON5 EDEP0 + decay energy vs Serpent2 Delta {iso} Comparison')
        plt.legend()
        plt.grid()
        plt.savefig(f"{save_dir_D5_vs_S2}/D5_EDEP0_decay_vs_S2_delta_{iso}_{ssh_option}_{correlation}_{time_integrator}.png")
        plt.close()


    


if __name__ == "__main__":

    post_treat_HOM_Gd157_VBOC()
    post_treat_AT10_45Gd()
    
    

    