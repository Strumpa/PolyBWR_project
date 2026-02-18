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


def plot_comparison(comparison_df_dict, tracked_nuclides, title, ssh_option, correlation, time_integrator, save_name, save_path):
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    # Plot delta keff
    plt.figure(figsize=(10, 6))
    for key in comparison_df_dict:
        plt.plot(comparison_df_dict[key]['BU Points'], comparison_df_dict[key]['delta_keff'], label=key, marker='x')
    # plot +/- 300 pcm lines
    plt.axhline(y=300, color='r', linestyle='--', label='+300 pcm')
    plt.axhline(y=-300, color='r', linestyle='--', label='-300 pcm')
    plt.xlabel('Burnup Points (MWd/kgU)')
    plt.ylabel('Delta keff (pcm)')
    plt.title(f"{title} : delta Keff")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_path}/delta_keff_{save_name}.png")
    plt.close()
    
    # Plot delta isotopes (%)
    for iso in tracked_nuclides:
        plt.figure(figsize=(10, 6))
        for key in comparison_df_dict:
            plt.plot(comparison_df_dict[key]['BU Points'], comparison_df_dict[key][f'delta {iso} (%)'], label=key, marker='x')
        # plot +/- 2 % lines
        plt.axhline(y=2, color='r', linestyle='--', label='+2 %')
        plt.axhline(y=-2, color='r', linestyle='--', label='-2 %')
        plt.xlabel('Burnup Points (MWd/kgU)')
        plt.ylabel(f'Delta {iso} (%)')
        plt.title(f'{title} : delta {iso} (%)')
        plt.legend()
        plt.grid()
        plt.savefig(f"{save_path}/delta_{iso}_{save_name}.png")
        plt.close()

    
def post_treat_AT10_cell(cell_name, bu_points):
    """
    Post-treat a AT10 pincell case.
    """
    name_D5_case = f"{cell_name}_Cst_pow_evol"
    # define path to NOGL_HM data
    path_to_Version5_wc_data =  "Version5_EVOGd_data_path"
    # Define calculation options for the CPOs name reconstruction
    time_integrator = "EXTR"
    BU_points = bu_points
    name_S2_case = f"{cell_name}"
    if BU_points == "Gd2_autop6":
        name_S2_case = f"{cell_name}_BUGd2"

    tracking_option = "SALT"
    evaluation = "ENDFB8R1"
#    draglib = "endfb8r1dec" # "endfb8r1_295_NG0" #"endfb8r1_295" # "endfb8r1_295", "endfb81295K", "endfb81295K2"
    ssh_option = "PT"
    correlation = "N"
    tracked_nuclides = ["U235", "U238", "Pu239", "Pu240", "Pu241", "Pu242", "Am241", "Gd155", "Gd157", "Xe135", "Sm149"]

    name_CPO_NOGL = f"CPO_endfb8r1_295_{ssh_option}_{correlation}_{tracking_option}_{BU_points}_KAPS_NODI_{time_integrator}_NOGL"
    name_CPO_EDP0 = f"CPO_endfb8r1_295_{ssh_option}_{correlation}_{tracking_option}_{BU_points}_KAPS_NODI_{time_integrator}_EDP0"
    #name_CPO_NOGL_dec = f"CPO_endfb8r1dec_{ssh_option}_{correlation}_{tracking_option}_{BU_points}_KAPS_NODI_{time_integrator}_NOGL"
    #name_CPO_EDP0_dec = f"CPO_endfb8r1dec_{ssh_option}_{correlation}_{tracking_option}_{BU_points}_KAPS_NODI_{time_integrator}_EDP0"
    path = os.getcwd()
    save_dir_case = f"{path}/postTreatEVO_tests_results/{name_D5_case}/{evaluation}_D5"
    if not os.path.exists(save_dir_case):
        os.makedirs(save_dir_case)
    save_dir_D5_vs_D5 = f"{save_dir_case}/D5_vs_D5_{BU_points}"
    if not os.path.exists(save_dir_D5_vs_D5):
        os.makedirs(save_dir_D5_vs_D5)
    save_dir_D5_vs_S2 = f"{save_dir_case}/D5_vs_S2_{BU_points}"
    if not os.path.exists(save_dir_D5_vs_S2):
        os.makedirs(save_dir_D5_vs_S2)
    
    # Load wc data, NOGL case
    os.chdir(f"{path_to_Version5_wc_data}/{name_D5_case}_results")
    print(os.listdir())
    print(name_CPO_NOGL)
    CPO_NOGL = lcm.new('LCM_INP', name_CPO_NOGL, impx=0)
    CPO_EDP0 = lcm.new('LCM_INP', name_CPO_EDP0, impx=0)

    os.chdir(path)



    # Create D5 cases
    D5_case_EDP0 = D5_case(pyCOMPO=CPO_EDP0, dlib_name="endfb8r1v5p1", bu_points=BU_points, 
                                ssh_opt=ssh_option, correlation=correlation, sat="NODI", depl_sol="EXTR", 
                                tracked_nuclides=tracked_nuclides, BU_lists=getBULists(BU_points), save_dir=save_dir_case)

    D5_case_NOGL = D5_case(pyCOMPO=CPO_NOGL, dlib_name="endfb8r1v5p1", bu_points=BU_points, 
                                ssh_opt=ssh_option, correlation=correlation, sat="NODI", depl_sol="EXTR", 
                                tracked_nuclides=tracked_nuclides, BU_lists=getBULists(BU_points), save_dir=save_dir_case)
    


    # Create S2 case
    #S2_case_edep0_pcc0 = S2_case(case_name=name_S2_case, lib_name="J311_pynjoy2016",
    #                            edep_id=0, areQfissSet=False, isEcaptSet=False, pcc_id=0, 
    #                            specific_power=26.5, tracked_nuclides=tracked_nuclides, save_dir=save_dir_case)
    S2_case_edep0_pcc1 = S2_case(case_name=name_S2_case, lib_name="endfb8r1_pynjoy2012_kerma",
                                edep_id=0, areQfissSet=False, isEcaptSet=False, pcc_id=1, 
                                specific_power=26.5, tracked_nuclides=tracked_nuclides, save_dir=save_dir_case)
    print(f"S2 edep0 pcc1 BU points are : {S2_case_edep0_pcc1.BU}")
    #S2_case_edep0_pcc2 = S2_case(case_name=name_S2_case, lib_name="endfb8r1_pynjoy2012_kerma",
    #                            edep_id=0, areQfissSet=False, isEcaptSet=False, pcc_id=2, 
    #                            specific_power=26.5, tracked_nuclides=tracked_nuclides, save_dir=save_dir_case)
    #print(f"S2 edep0 pcc2 BU points are : {S2_case_edep0_pcc2.BU}")
    S2_case_edep2_pcc1 = S2_case(case_name=name_S2_case, lib_name="endfb8r1_pynjoy2012_kerma",
                                edep_id=2, areQfissSet=False, isEcaptSet=False, pcc_id=1, 
                                specific_power=26.5, tracked_nuclides=tracked_nuclides, save_dir=save_dir_case)
    #S2_case_edep2_pcc2 = S2_case(case_name=name_S2_case, lib_name="J311_pynjoy2016",
    #                            edep_id=2, areQfissSet=False, isEcaptSet=False, pcc_id=2, 
    #                            specific_power=26.5, tracked_nuclides=tracked_nuclides, save_dir=save_dir_case)
    
   

    ## Compare DRAGON5 cases with Serpent2 edep0 pcc1 and pcc2 cases
    comp_D5_NOGL_S2_pcc1 = create_D5_S2_comparison(D5_case_NOGL, S2_case_edep0_pcc1)
    #comp_D5_NOGL_S2_pcc2 = create_D5_S2_comparison(D5_case_NOGL, S2_case_edep0_pcc2)
    comp_D5_EDP0_S2_pcc1 = create_D5_S2_comparison(D5_case_EDP0, S2_case_edep0_pcc1)
    #comp_D5_EDP0_S2_pcc2 = create_D5_S2_comparison(D5_case_EDP0, S2_case_edep0_pcc2)

    # Plot comparisons
    
    comparison_dict_NOGL = {
        'D5 NOGL - S2 edep 0 pcc 1': comp_D5_NOGL_S2_pcc1,
        #'D5 NOGL - S2 edep 0 pcc 2': comp_D5_NOGL_S2_pcc2,
    }
    plot_comparison(comparison_dict_NOGL, tracked_nuclides, title="D5 NOGL vs S2 edep 0 pcc 1", 
                    ssh_option=ssh_option, correlation=correlation, time_integrator=time_integrator, 
                    save_name="D5_NOGL_vs_S2_edep0_pcc1", save_path=save_dir_D5_vs_S2)
    
    comparison_dict_EDP0 = {
        'D5 EDP0 - S2 edep 0 pcc 1': comp_D5_EDP0_S2_pcc1,
        #'D5 EDP0 - S2 edep 0 pcc 2': comp_D5_EDP0_S2_pcc2,
    }
    plot_comparison(comparison_dict_EDP0, tracked_nuclides, title="D5 EDEP0 vs S2 edep 0 pcc 1", 
                    ssh_option=ssh_option, correlation=correlation, time_integrator=time_integrator, 
                    save_name="D5_EDEP0_vs_S2_edep0_pcc1", save_path=save_dir_D5_vs_S2)


def post_treat_AT10_24UOX_schemes():
    """
    Post treat the AT10_24UOX case for different depletion schemes.
    checking if 1L and 2L schemes give comparable results :
    it was noticed that in the case of 2L scheme fission Q values were not present in the RER variable, presenting a potential issue in the normalisation procedure.
    1L scheme : 1 level MOC calculation with 295 groups
    2L scheme : 2 levels MOC : PIJ on 295 groups, then 1L MOC on 26 groups.
    
    
    """
    
    name_D5_case = "ATRIUM10_evol"
    # define path to NOGL_HM data
    path_to_Version5_1_data =  "Version5_1_results"
    # Define calculation options for the CPOs name reconstruction
    path_to_data_2L = f"2L_PIJ_MOC_norm_NONE/AT10_24UOX/finest_geom_AT10_void_0_J311_295_PT_NOCORR_region_num"
    path_to_data_1L = f"1L_MOC_norm_NONE/AT10_24UOX/finest_geom_AT10_void_0_J311_295_PT_NOCORR_region_num"
    CPO_name = "CPO_n24_ld150_n8_ld25_TSPC_4_PIJ_MOC_GAUS_4"
    time_integrator = "EXTR"
    BU_points = "ATRIUM_10_S2_BU"
    name_S2_case = "AT10_24UOX"
    ssh_option = "PT"
    correlation = "N"
    if BU_points == "Gd2_autop6":
        name_S2_case = "AT10_45Gd_BUGd2"

    tracking_option = "SALT"
    draglib = "J311_295" # "endfb8r1_295_NG0" #"endfb8r1_295" # "endfb8r1_295", "endfb81295K", "endfb81295K2"
    tracked_nuclides = ["U235", "U238", "Pu239", "Pu240", "Pu241", "Pu242", "Am241", "Gd155", "Gd157", "Xe135", "Sm149"]

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
    
    # Load the 1L scheme 5.1 data
    os.chdir(f"{path_to_Version5_1_data}/{name_D5_case}_results/{path_to_data_1L}")
    print(os.listdir())
    print(CPO_name)
    CPO_1L = lcm.new('LCM_INP', CPO_name, impx=0)
    os.chdir(path)

    # Load the 2L scheme 5.1 data
    os.chdir(f"{path_to_Version5_1_data}/{name_D5_case}_results/{path_to_data_2L}")
    print(os.listdir())
    CPO_2L = lcm.new('LCM_INP', CPO_name, impx=0)
    os.chdir(path)


    # Create D5 case 
    D5_case_1L = D5_case(pyCOMPO=CPO_1L, dlib_name=draglib, bu_points=BU_points, 
                                ssh_opt=ssh_option, correlation=correlation, sat="NODI", depl_sol="EXTR", 
                                tracked_nuclides=tracked_nuclides, BU_lists=getBULists(BU_points), save_dir=save_dir_case)
    D5_case_2L = D5_case(pyCOMPO=CPO_2L, dlib_name=draglib, bu_points=BU_points, 
                                ssh_opt=ssh_option, correlation=correlation, sat="NODI", depl_sol="EXTR", 
                                tracked_nuclides=tracked_nuclides, BU_lists=getBULists(BU_points), save_dir=save_dir_case)
    
    
    # Create D5 comparison
    comparison_1L_2L = create_D5_comparison(D5_case_1L, D5_case_2L)
    comparison_dict_1L_2L = {
        'D5 1L scheme - D5 2L scheme': comparison_1L_2L,
    }
    plot_comparison(comparison_dict_1L_2L, tracked_nuclides, title="D5 1L vs D5 2L scheme", 
                    ssh_option=ssh_option, correlation=correlation, time_integrator=time_integrator, 
                    save_name="D5_1L_vs_D5_2L_scheme", save_path=save_dir_D5_vs_D5)
    
    # Create S2 case problem : no J311_pynjoy2016_kerma runs left ? only endfb8r1
    """
    S2_case_edep0_pcc1 = S2_case(case_name=name_S2_case, lib_name="J311_pynjoy2016",
                                edep_id=0, areQfissSet=False, isEcaptSet=False, pcc_id=1, 
                                specific_power=26.5, tracked_nuclides=tracked_nuclides, save_dir=save_dir_case)
    print(f"S2 edep0 pcc1 BU points are : {S2_case_edep0_pcc1.BU}")
    S2_case_edep0_pcc2 = S2_case(case_name=name_S2_case, lib_name="J311_pynjoy2016",
                                edep_id=0, areQfissSet=False, isEcaptSet=False, pcc_id=2, 
                                specific_power=26.5, tracked_nuclides=tracked_nuclides, save_dir=save_dir_case)
    """

    
def post_treat_AT10_45Gd_schemes():
    name_D5_case = "ATRIUM10_evol"
    # define path to NOGL_HM data
    path_to_Version5_1_data =  "Version5_1_results"
    path_to_Version5_1_edp0_data =  "Version5_1_edp0_mixpwr_results"
    # Define calculation options for the CPOs name reconstruction
    #CPO_name = "CPO_n24_ld35_n8_ld25_TSPC_4_PIJ_MOC_GAUS_4"
    
    time_integrator = "EXTR"
    BU_points = "Gd_autop3" #
    name_S2_case = "AT10_45Gd"
    ssh_option = "PT"
    correlation = "N"
    if BU_points == "Gd2_autop6":
        name_S2_case = "AT10_45Gd_BUGd2"
    CPO_name = f"CPO_AT10_45Gd_{BU_points}"
    tracking_option = "SALT"
    draglib = "J311_295E0" # "endfb8r1_295_NG0" #"endfb8r1_295" # "endfb8r1_295", "endfb81295K", "endfb81295K2"
    tracked_nuclides = ["U235", "U238", "Pu239", "Pu240", "Pu241", "Pu242", "Am241", "Gd155", "Gd157", "Xe135", "Sm149"]

    path_to_data_2L = f"2L_PIJ_MOC_norm_NONE/{name_S2_case}/finest_geom_AT10_void_0_{draglib}_{ssh_option}_NOCORR_region_num"
    path_to_data_1L = f"1L_MOC_norm_NONE/{name_S2_case}/finest_geom_AT10_void_0_{draglib}_{ssh_option}_NOCORR_region_num"

    path_to_data_2L_QFIS = f"2L_PIJ_MOC_norm_QFIS/{name_S2_case}/finest_geom_AT10_void_0_{draglib}_{ssh_option}_NOCORR_region_num"
    path_to_data_1L_QFIS = f"1L_MOC_norm_QFIS/{name_S2_case}/finest_geom_AT10_void_0_{draglib}_{ssh_option}_NOCORR_region_num"

    path_to_data_2L_EDP0 = f"2L_PIJ_MOC_norm_EDP0/{name_S2_case}/finest_geom_AT10_void_0_{draglib}_{ssh_option}_NOCORR_region_num"
    path_to_data_1L_EDP0 = f"1L_MOC_norm_EDP0/{name_S2_case}/finest_geom_AT10_void_0_{draglib}_{ssh_option}_NOCORR_region_num"

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
    
    # Load the 1L scheme 5.1 data, normalisation without Qfiss correction
    os.chdir(f"{path_to_Version5_1_data}/{name_D5_case}_results/{path_to_data_1L}")
    print(os.listdir())
    CPO_1L = lcm.new('LCM_INP', CPO_name, impx=0)
    os.chdir(path)
    # Load the 1L scheme 5.1 data, normalisation with Qfiss correction
    os.chdir(f"{path_to_Version5_1_data}/{name_D5_case}_results/{path_to_data_1L_QFIS}")
    print(os.listdir())
    CPO_1L_QFIS = lcm.new('LCM_INP', CPO_name, impx=0)
    os.chdir(path)

    # Load the 2L scheme 5.1 data, normalisation without Qfiss correction
    os.chdir(f"{path_to_Version5_1_data}/{name_D5_case}_results/{path_to_data_2L}")
    print(os.listdir())
    CPO_2L = lcm.new('LCM_INP', CPO_name, impx=0)
    os.chdir(path)
    
    # Load the 2L scheme 5.1 data, normalisation with Qfiss correction
    os.chdir(f"{path_to_Version5_1_data}/{name_D5_case}_results/{path_to_data_2L_QFIS}")
    print(os.listdir())
    CPO_2L_QFIS = lcm.new('LCM_INP', CPO_name, impx=0)
    os.chdir(path)


    # Load the 1L scheme 5.1 edp0 + MIXP data, normalisation Qfiss
    os.chdir(f"{path_to_Version5_1_edp0_data}/{name_D5_case}_results/{path_to_data_1L_QFIS}")
    print(os.listdir())
    CPO_1L_QFIS_MIXP = lcm.new('LCM_INP', CPO_name, impx=0)
    os.chdir(path)
    
    # Load the 1L scheme 5.1 data, normalisation with Qfiss correction
    if "E0" in draglib:
        os.chdir(f"{path_to_Version5_1_edp0_data}/{name_D5_case}_results/{path_to_data_1L_EDP0}")
        print(os.listdir())
        CPO_1L_EDP0_MIXP = lcm.new('LCM_INP', CPO_name, impx=0)
        os.chdir(path)

    # Load the 2L scheme 5.1 edp0 + MIXP data, normalisation Qfiss
    os.chdir(f"{path_to_Version5_1_edp0_data}/{name_D5_case}_results/{path_to_data_2L_QFIS}")
    print(os.listdir())
    CPO_2L_QFIS_MIXP = lcm.new('LCM_INP', CPO_name, impx=0)
    os.chdir(path)
    
    # Load the 2L scheme 5.1 data, normalisation with Qfiss correction
    os.chdir(f"{path_to_Version5_1_edp0_data}/{name_D5_case}_results/{path_to_data_2L_EDP0}")
    print(os.listdir())
    CPO_2L_EDP0_MIXP = lcm.new('LCM_INP', CPO_name, impx=0)
    os.chdir(path)



    # Create D5 case 
    D5_case_1L = D5_case(pyCOMPO=CPO_1L, dlib_name=draglib, bu_points=BU_points, 
                                ssh_opt=ssh_option, correlation=correlation, sat="NODI", depl_sol="EXTR", 
                                tracked_nuclides=tracked_nuclides, BU_lists=getBULists(BU_points), save_dir=save_dir_case)
    D5_case_2L = D5_case(pyCOMPO=CPO_2L, dlib_name=draglib, bu_points=BU_points, 
                                ssh_opt=ssh_option, correlation=correlation, sat="NODI", depl_sol="EXTR", 
                                tracked_nuclides=tracked_nuclides, BU_lists=getBULists(BU_points), save_dir=save_dir_case)
    
    D5_case_1L_QFIS = D5_case(pyCOMPO=CPO_1L_QFIS, dlib_name=draglib, bu_points=BU_points, 
                                ssh_opt=ssh_option, correlation=correlation, sat="NODI", depl_sol="EXTR", 
                                tracked_nuclides=tracked_nuclides, BU_lists=getBULists(BU_points), save_dir=save_dir_case)
    D5_case_2L_QFIS = D5_case(pyCOMPO=CPO_2L_QFIS, dlib_name=draglib, bu_points=BU_points, 
                                ssh_opt=ssh_option, correlation=correlation, sat="NODI", depl_sol="EXTR", 
                                tracked_nuclides=tracked_nuclides, BU_lists=getBULists(BU_points), save_dir=save_dir_case)
    

    ### Cases from edp0_debug : using MIXPWR keyword
    # Qfiss normalisation :
    D5_case_1L_QFIS_MIXP = D5_case(pyCOMPO=CPO_1L_QFIS_MIXP, dlib_name=draglib, bu_points=BU_points, 
                                ssh_opt=ssh_option, correlation=correlation, sat="NODI", depl_sol="EXTR", 
                                tracked_nuclides=tracked_nuclides, BU_lists=getBULists(BU_points), save_dir=save_dir_case)
    D5_case_2L_QFIS_MIXP = D5_case(pyCOMPO=CPO_2L_QFIS_MIXP, dlib_name=draglib, bu_points=BU_points, 
                                ssh_opt=ssh_option, correlation=correlation, sat="NODI", depl_sol="EXTR", 
                                tracked_nuclides=tracked_nuclides, BU_lists=getBULists(BU_points), save_dir=save_dir_case)

    # EDP0 normalisation :
    if "E0" in draglib:
        D5_case_1L_EDP0_MIXP = D5_case(pyCOMPO=CPO_1L_EDP0_MIXP, dlib_name=draglib, bu_points=BU_points, 
                                    ssh_opt=ssh_option, correlation=correlation, sat="NODI", depl_sol="EXTR", 
                                    tracked_nuclides=tracked_nuclides, BU_lists=getBULists(BU_points), save_dir=save_dir_case)
        
    D5_case_2L_EDP0_MIXP = D5_case(pyCOMPO=CPO_2L_EDP0_MIXP, dlib_name=draglib, bu_points=BU_points, 
                                ssh_opt=ssh_option, correlation=correlation, sat="NODI", depl_sol="EXTR", 
                                tracked_nuclides=tracked_nuclides, BU_lists=getBULists(BU_points), save_dir=save_dir_case)

    # Compare 2L and 1L schemes for same normalisation option
    # Default NOGL normalisation
    comparison_1L_2L = create_D5_comparison(D5_case_1L, D5_case_2L)
    comparison_dict_1L_2L = {
        'D5 1L scheme - D5 2L scheme': comparison_1L_2L,
    }
    plot_comparison(comparison_dict_1L_2L, tracked_nuclides, title="D5 1L vs D5 2L scheme", 
                    ssh_option=ssh_option, correlation=correlation, time_integrator=time_integrator, 
                    save_name="D5_1L_vs_D5_2L_scheme", save_path=save_dir_D5_vs_D5)
    
    # QFIS normalisation
    comparison_1L_QFIS_2L_QFIS = create_D5_comparison(D5_case_1L_QFIS, D5_case_2L_QFIS)
    comparison_dict_1L_QFIS_2L_QFIS = {
        'D5 1L scheme QFIS - D5 2L scheme QFIS': comparison_1L_QFIS_2L_QFIS,
    }
    plot_comparison(comparison_dict_1L_QFIS_2L_QFIS, tracked_nuclides, title="D5 1L vs 2L, QFIS normalisation", 
                    ssh_option=ssh_option, correlation=correlation, time_integrator=time_integrator, 
                    save_name="D5_1L_QFIS_vs_D5_2L_QFIS", save_path=save_dir_D5_vs_D5)
    
    # QFIS normalisation with MIXP keyword
    comparison_1L_QFIS_MIXP_2L_QFIS_MIXP = create_D5_comparison(D5_case_1L_QFIS_MIXP, D5_case_2L_QFIS_MIXP)
    comparison_dict_1L_QFIS_MIXP_2L_QFIS_MIXP = {
        'D5 1L scheme QFIS MIXP - D5 2L scheme QFIS MIXP': comparison_1L_QFIS_MIXP_2L_QFIS_MIXP,
    }
    plot_comparison(comparison_dict_1L_QFIS_MIXP_2L_QFIS_MIXP, tracked_nuclides, title="D5 1L vs 2L, QFIS normalisation, MIXP", 
                    ssh_option=ssh_option, correlation=correlation, time_integrator=time_integrator, 
                    save_name="D5_1L_QFIS_MIXP_vs_D5_2L_QFIS_MIXP", save_path=save_dir_D5_vs_D5)
    
    # EDP0 normalisation with MIXP keyword
    if "E0" in draglib:
        comparison_1L_EDP0_MIXP_2L_EDP0_MIXP = create_D5_comparison(D5_case_1L_EDP0_MIXP, D5_case_2L_EDP0_MIXP)
        comparison_dict_1L_EDP0_MIXP_2L_EDP0_MIXP = {
            'D5 1L scheme EDP0 MIXP - D5 2L scheme EDP0 MIXP': comparison_1L_EDP0_MIXP_2L_EDP0_MIXP,
        }
        plot_comparison(comparison_dict_1L_EDP0_MIXP_2L_EDP0_MIXP, tracked_nuclides, title="D5 1L vs 2L, EDP0 normalisation, MIXP keyword", 
                        ssh_option=ssh_option, correlation=correlation, time_integrator=time_integrator, 
                        save_name="D5_1L_EDP0_MIXP_vs_D5_2L_EDP0_MIXP", save_path=save_dir_D5_vs_D5)
        
    # Create S2 case
    S2_case_edep0_pcc2 = S2_case(case_name=name_S2_case, lib_name="J311_pynjoy2016",
                                edep_id=0, areQfissSet=False, isEcaptSet=False, pcc_id=2, 
                                specific_power=26.5, tracked_nuclides=tracked_nuclides, save_dir=save_dir_case)
    print(f"S2 edep0 pcc2 BU points are : {S2_case_edep0_pcc2.BU}")
    
    
    # Compare D5 cases with S2 case
    # Default NOGL normalisation, 1L and 2L schemes
    comp_D5_1L_S2_pcc2 = create_D5_S2_comparison(D5_case_1L, S2_case_edep0_pcc2)
    comp_D5_2L_S2_pcc2 = create_D5_S2_comparison(D5_case_2L, S2_case_edep0_pcc2)
    comparison_dict_1L_S2 = {
        'D5 1L scheme - S2 edep 0 pcc 2': comp_D5_1L_S2_pcc2,
        'D5 2L scheme - S2 edep 0 pcc 2': comp_D5_2L_S2_pcc2,
    }
    plot_comparison(comparison_dict_1L_S2, tracked_nuclides, title="D5 1L and 2L vs S2 edep 0 pcc 2", 
                    ssh_option=ssh_option, correlation=correlation, time_integrator=time_integrator, 
                    save_name="D5_1L_2L_vs_S2_edep0_pcc2", save_path=save_dir_D5_vs_S2)
    
    # QFIS normalisation, 1L and 2L schemes
    comp_D5_1L_QFIS_S2_pcc2 = create_D5_S2_comparison(D5_case_1L_QFIS, S2_case_edep0_pcc2)
    comp_D5_2L_QFIS_S2_pcc2 = create_D5_S2_comparison(D5_case_2L_QFIS, S2_case_edep0_pcc2)
    comparison_dict_1L_QFIS_S2 = {
        'D5 1L scheme QFIS - S2 edep 0 pcc 2': comp_D5_1L_QFIS_S2_pcc2,
        'D5 2L scheme QFIS - S2 edep 0 pcc 2': comp_D5_2L_QFIS_S2_pcc2,
    }
    plot_comparison(comparison_dict_1L_QFIS_S2, tracked_nuclides, title="D5 1L and 2L QFIS vs S2 edep 0 pcc 2", 
                    ssh_option=ssh_option, correlation=correlation, time_integrator=time_integrator, 
                    save_name="D5_1L_2L_QFIS_vs_S2_edep0_pcc2", save_path=save_dir_D5_vs_S2)
    
    # QFIS normalisation with MIXP keyword, 1L and 2L schemes
    comp_D5_1L_QFIS_MIXP_S2_pcc2 = create_D5_S2_comparison(D5_case_1L_QFIS_MIXP, S2_case_edep0_pcc2)
    comp_D5_2L_QFIS_MIXP_S2_pcc2 = create_D5_S2_comparison(D5_case_2L_QFIS_MIXP, S2_case_edep0_pcc2)
    comparison_dict_1L_QFIS_MIXP_S2 = {
        'D5 1L scheme QFIS MIXP - S2 edep 0 pcc 2': comp_D5_1L_QFIS_MIXP_S2_pcc2,
        'D5 2L scheme QFIS MIXP - S2 edep 0 pcc 2': comp_D5_2L_QFIS_MIXP_S2_pcc2,
    }
    plot_comparison(comparison_dict_1L_QFIS_MIXP_S2, tracked_nuclides, title="D5 1L and 2L QFIS MIXP vs S2 edep 0 pcc 2", 
                    ssh_option=ssh_option, correlation=correlation, time_integrator=time_integrator,
                    save_name="D5_1L_2L_QFIS_MIXP_vs_S2_edep0_pcc2", save_path=save_dir_D5_vs_S2)
    
    # EDP0 normalisation with MIXP keyword, 1L and 2L schemes
    if "E0" in draglib:
        comp_D5_1L_EDP0_MIXP_S2_pcc2 = create_D5_S2_comparison(D5_case_1L_EDP0_MIXP, S2_case_edep0_pcc2)
        comp_D5_2L_EDP0_MIXP_S2_pcc2 = create_D5_S2_comparison(D5_case_2L_EDP0_MIXP, S2_case_edep0_pcc2)
        comparison_dict_1L_EDP0_MIXP_S2 = {
            'D5 1L scheme EDP0 MIXP - S2 edep 0 pcc 2': comp_D5_1L_EDP0_MIXP_S2_pcc2,
            'D5 2L scheme EDP0 MIXP - S2 edep 0 pcc 2': comp_D5_2L_EDP0_MIXP_S2_pcc2,
        }
        plot_comparison(comparison_dict_1L_EDP0_MIXP_S2, tracked_nuclides, title="D5 1L and 2L EDP0 MIXP vs S2 edep 0 pcc 2", 
                        ssh_option=ssh_option, correlation=correlation, time_integrator=time_integrator,
                        save_name="D5_1L_2L_EDP0_MIXP_vs_S2_edep0_pcc2", save_path=save_dir_D5_vs_S2)
        

def post_treat_EDP0_QFIS_tests(case_name, BU_points):
    """
    Post-treatment for EDP0 tests.
    """
    path = os.getcwd()
    name_D5_case = "ATRIUM10_evol"
    path_to_Version5_1_data =  "Version5_1_results"
    path_to_Version5_1_edp0_data =  "Version5_1_edp0_mixpwr_results"
    # Define calculation options for the CPOs name / D5_cases parsing and reconstruction
    time_integrator = "EXTR"
    name_S2_case = case_name #"AT10_45Gd", "AT10_ASSBLY"
    if case_name == "AT10_ASSBLY":
        name_S2_case_path = f"{os.environ['SERPENT_RESULTS']}/ATRIUM10_pin_jeff311_pynjoy2016_edep0_pcc2_BU_inp"
    else: 
        name_S2_case_path = None
    ssh_option = "RSE"
    ssh_solution_door = "IC"
    correlation = "N"
    if BU_points == "Gd2_autop6":
        name_S2_case = "AT10_45Gd_BUGd2"
    CPO_name = f"CPO_{case_name}_{BU_points}"
    tracking_option = "SALT"
    draglib = "J311_295" # "endfb8r1_295_NG0" #"endfb8r1_295" # "endfb8r1_295", "endfb81295K", "endfb81295K2"
    tracked_nuclides = ["U235", "U238", "Pu239", "Pu240", "Pu241", "Pu242", "Am241", "Gd155", "Gd157", "Xe135", "Sm149"]
    calculation_schemes = ["2L_IC_MOC", "2L_PIJ_MOC"]#, "2L_PIJ_MOC"] #["1L_MOC", "2L_PIJ_MOC"]
    normalisation_options = ["EDP0","QFIS"]

    modif_Q_values_options = [#("std_Q_fiss", "std_NG_Q"), 
                            #("modif_Q_fiss", "std_NG_Q"),
                            #("modif_Q_fiss", "no_NG"),
                            #("modif_Q_fiss", "no_NG_ALL"),
                            ("modif_Q_fiss", "no_NG_ALL_ADDED_Pu_Am_Np"),]
                            #("XS_modif_Q_fiss", "no_NG_ALL")]

    save_dir_case = f"{path}/postTreatEVO_tests_results/{draglib}/{case_name}"
    if not os.path.exists(save_dir_case):
        os.makedirs(save_dir_case)
    save_dir_D5_vs_D5 = f"{save_dir_case}/D5_vs_D5_{BU_points}"
    if not os.path.exists(save_dir_D5_vs_D5):
        os.makedirs(save_dir_D5_vs_D5)
    save_dir_D5_vs_S2 = f"{save_dir_case}/D5_vs_S2_{BU_points}"
    if not os.path.exists(save_dir_D5_vs_S2):
        os.makedirs(save_dir_D5_vs_S2)   

    DRAGON_cases = {} # Dictionary to store results from cases.

    for calc_scheme in calculation_schemes:
        for norm_opt in normalisation_options:
            for combination in modif_Q_values_options:
                path_to_data = f"{ssh_solution_door}_{calc_scheme}_norm_{norm_opt}_{combination[0]}_{combination[1]}/{name_S2_case}/finest_geom_AT10_void_0_{draglib}_{ssh_option}_NOCORR_region_num"
                os.chdir(f"{path_to_Version5_1_edp0_data}/{name_D5_case}_results/{path_to_data}")
                py_CPO = lcm.new('LCM_INP', CPO_name, impx=0)
                os.chdir(path)
                D5_case_opts = D5_case(pyCOMPO=py_CPO, dlib_name=draglib, bu_points=BU_points,
                    ssh_opt=ssh_option, correlation=correlation, sat="NODI", depl_sol="EXTR",
                    tracked_nuclides=tracked_nuclides, BU_lists=getBULists(BU_points), save_dir=save_dir_case)
                D5_case_opts.calculation_scheme = calc_scheme
                D5_case_opts.rates_normalisation = norm_opt
                D5_case_opts.fission_Q_opt = combination[0]
                D5_case_opts.NG_Q_opt = combination[1]
                DRAGON_cases[f"{calc_scheme} : {norm_opt}, {combination[0]} {combination[1]}"] = D5_case_opts

    #


    # Serpent2 case : edep 0 ppc 2
    S2_case_edep0_pcc2 = S2_case(case_name=name_S2_case, lib_name="J311_pynjoy2016",
                    edep_id=0, areQfissSet=False, isEcaptSet=False, pcc_id=2, 
                    specific_power=26.5, tracked_nuclides=tracked_nuclides, save_dir=save_dir_case, path=name_S2_case_path)

    ## Create comparisons based on normaliation opt. 
    # Filter out specific options to be included in comparisons
    # 
    comparisons_dict_D5_vs_S2_EDP0 = {}
    comparisons_dict_D5_1L_vs_S2_EDP0 = {}
    comparisons_dict_D5_2L_vs_S2_EDP0 = {}

    comparisons_dict_D5_vs_S2_QFIS = {}
    comparisons_dict_D5_1L_vs_S2_QFIS = {}
    comparisons_dict_D5_2L_vs_S2_QFIS = {}

    D5_cases_desc = DRAGON_cases.keys()
    for desc in D5_cases_desc:
        DRAGON_case = DRAGON_cases[desc]
        
        comp_D5_S2 = create_D5_S2_comparison(DRAGON_case, S2_case_edep0_pcc2)
        scheme = DRAGON_case.calculation_scheme
        norm = DRAGON_case.rates_normalisation
        fission_Q_opt = DRAGON_case.fission_Q_opt
        NG_Q_opt = DRAGON_case.NG_Q_opt
        if norm == "EDP0":
            comparisons_dict_D5_vs_S2_EDP0[f"D5 {scheme} EDP0 {fission_Q_opt} {NG_Q_opt}"] = comp_D5_S2
            if scheme == "1L_MOC":
                comparisons_dict_D5_1L_vs_S2_EDP0[f"D5 {scheme} EDP0 {fission_Q_opt} {NG_Q_opt}"] = comp_D5_S2
            elif scheme == "2L_PIJ_MOC":
                comparisons_dict_D5_2L_vs_S2_EDP0[f"D5 {scheme} EDP0 {fission_Q_opt} {NG_Q_opt}"] = comp_D5_S2
        if norm == "QFIS":
            comparisons_dict_D5_vs_S2_QFIS[f"D5 {scheme} QFIS {fission_Q_opt} {NG_Q_opt}"] = comp_D5_S2
            if scheme == "1L_MOC":
                comparisons_dict_D5_1L_vs_S2_QFIS[f"D5 {scheme} QFIS {fission_Q_opt} {NG_Q_opt}"] = comp_D5_S2
            elif scheme == "2L_PIJ_MOC" or "2L_IC_MOC":
                comparisons_dict_D5_2L_vs_S2_QFIS[f"D5 {scheme} QFIS {fission_Q_opt} {NG_Q_opt}"] = comp_D5_S2
    
    # All EDP0 cases                
    plot_comparison(comparisons_dict_D5_vs_S2_EDP0, tracked_nuclides, title="D5, EDP0 normalisation vs S2 edep 0 pcc 2",
                    ssh_option=ssh_option, correlation=correlation, time_integrator=time_integrator,
                    save_name=f"{name_S2_case}_D5_vs_S2_EDP0", save_path=f"{save_dir_D5_vs_S2}/EDP0")
    
    plot_comparison(comparisons_dict_D5_1L_vs_S2_EDP0, tracked_nuclides, title="D5 1L, EDP0 normalisation vs S2 edep 0 pcc 2",
                    ssh_option=ssh_option, correlation=correlation, time_integrator=time_integrator,
                    save_name=f"1L_{name_S2_case}_D5_vs_S2_EDP0", save_path=f"{save_dir_D5_vs_S2}/EDP0/1L")
    
    plot_comparison(comparisons_dict_D5_2L_vs_S2_EDP0, tracked_nuclides, title="D5 2L, EDP0 normalisation vs S2 edep 0 pcc 2",
                    ssh_option=ssh_option, correlation=correlation, time_integrator=time_integrator,
                    save_name=f"2L_{name_S2_case}_D5_vs_S2_EDP0", save_path=f"{save_dir_D5_vs_S2}/EDP0/2L")

    plot_comparison(comparisons_dict_D5_vs_S2_QFIS, tracked_nuclides, title="D5, QFIS normalisation vs S2 edep 0 pcc 2",
                    ssh_option=ssh_option, correlation=correlation, time_integrator=time_integrator,
                    save_name=f"{name_S2_case}_D5_vs_S2_QFIS", save_path=f"{save_dir_D5_vs_S2}/QFIS")

    plot_comparison(comparisons_dict_D5_1L_vs_S2_QFIS, tracked_nuclides, title="D5 1L, QFIS normalisation vs S2 edep 0 pcc 2",
                    ssh_option=ssh_option, correlation=correlation, time_integrator=time_integrator,
                    save_name=f"1L_{name_S2_case}_D5_vs_S2_QFIS", save_path=f"{save_dir_D5_vs_S2}/QFIS/1L")

    plot_comparison(comparisons_dict_D5_2L_vs_S2_QFIS, tracked_nuclides, title="D5 2L, QFIS normalisation vs S2 edep 0 pcc 2",
                    ssh_option=ssh_option, correlation=correlation, time_integrator=time_integrator,
                    save_name=f"2L_{name_S2_case}_D5_vs_S2_QFIS", save_path=f"{save_dir_D5_vs_S2}/QFIS/2L")






        
    


if __name__ == "__main__":

    #post_treat_HOM_Gd157_VBOC()
    post_treat_AT10_cell(cell_name="AT10_24UOX", bu_points="UOx2_autop5")
    post_treat_AT10_cell(cell_name="AT10_24UOX", bu_points="UOx4_autop5")
    post_treat_AT10_cell(cell_name="AT10_45Gd", bu_points="Gd_autop3")
    post_treat_AT10_cell(cell_name="AT10_45Gd", bu_points="Gd_autop5")
    #post_treat_AT10_24UOX_schemes()
    #post_treat_AT10_45Gd_schemes()
    #post_treat_EDP0_QFIS_tests(case_name="AT10_ASSBLY", BU_points="ATRIUM_10_S2_BU")
    
    

    