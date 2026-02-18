import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import lifo 
import lcm
import cle2000 
from getLists import getLists
from postproc_cst_pow_evol import DRAGON_case as D5_case
from postproc_cst_pow_evol import Serpent2_case as S2_case

## Python3 script to post treat results from modified EVO: vs standard EVO: solvers, 
# In Dragon5 EVO: the "NOGL" keyword was previously used to specify that reaction rates normalization 
# should be performed using reactions in the fuel only.
# with "GLOB" keyword, the normalization is done using reactino from all isotopes present in the global geometry.

# This induces differences in the case when exothermic reactions happen outside of regions filled with fuel.

# In a new approach A. Hébert proposed a modification of the "NOGL" keyword such that it only "Heavy Metals" are used for normalization.

# This script considers 3 types of normalization:
# 1. "NOGL" - normalization using reactions in the fuel only
# 2. "GLOB" - normalization using reactions in all isotopes present in the global geometry
# 3. "NOGL_HM" - normalization using reactions in heavy metals only

# Heating factors can be computed through 2 main approaches:
# 1. Using the "Q-value" historical approach where the heating factors are computed through the mass difference equations for each reaction,
# 2. Using the "KERMA" approach where the heating factors are computed through the KERMA values of each reaction.

# Author : R. Guasch
# Date : 2025-06-02



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


if __name__ == "__main__":
    name_D5_case = "AT10_45Gd_Cst_pow_evol"
    name_S2_case = "AT10_45Gd"
    # define path to NOGL_HM data
    path_to_EVOmodif_data =  f"EVO_modif_path/{name_D5_case}_results"
    path_to_ref_data = f"PYGAN_COMPOS_path/{name_D5_case}_results"
    

    # path to OpenMC data
    path_to_OpenMC_data = f"{os.environ['HOME']}/working_dir/OpenMC_cases/post_treat"
    # OpenMC tracked nuclides :
    OpenMC_tracked_nuclides = ["Gd155", "Gd157", "Gd158", "U235", "U238", "Pu239", "Xe135", "Sm149"]
    # Define calculation options for the CPOs name reconstruction
    time_integrator = "EXTR"
    BU_points = "Gd_autop3" #"UOx2_autop5"  # "Gd_autop3"
    tracking_option = "SALT"
    draglib = "endfb81295K2" # "endfb8r1_295_NG0" #"endfb8r1_295" # "endfb8r1_295", "endfb81295K", "endfb81295K2"
    ssh_option = "PT"
    correlation = "N"
    tracked_nuclides = ["U235","U238","Pu239","Pu240","Pu241","Pu242","Am241","Xe135","Sm149","Gd155","Gd157"]

    name_cpo_NOGL = f"CPO_{draglib}_{ssh_option}_{correlation}_{tracking_option}_{BU_points}_KAPS_NODI_EXTR_NOGL"
    name_cpo_GLOB = f"CPO_{draglib}_{ssh_option}_{correlation}_{tracking_option}_{BU_points}_KAPS_NODI_EXTR_GLOB"

    path = os.getcwd()
    save_dir_case = f"{path}/postTreatGLOB_opt_results/{name_D5_case}/{draglib}_D5"
    if not os.path.exists(save_dir_case):
        os.makedirs(save_dir_case)
    
    # Load the NOGL_HM data
    os.chdir(path_to_EVOmodif_data)
    
    CPO_NOGL_modif = lcm.new('LCM_INP', name_cpo_NOGL, impx=0)
    CPO_GLOB_modif = lcm.new('LCM_INP', name_cpo_GLOB, impx=0)
    os.chdir(path)
    # Load the reference data
    os.chdir(path_to_ref_data)
    print(os.listdir())
    CPO_NOGL_ref = lcm.new('LCM_INP', name_cpo_NOGL, impx=0)
    CPO_GLOB_ref = lcm.new('LCM_INP', name_cpo_GLOB, impx=0)
    os.chdir(path)

    # Create D5 cases :
    D5_case_NOGL_ref = D5_case(pyCOMPO=CPO_NOGL_ref, dlib_name=draglib, bu_points=BU_points, 
                                ssh_opt="PT", correlation=correlation, sat="NODI", depl_sol="EXTR", 
                                tracked_nuclides=tracked_nuclides, BU_lists=getLists(BU_points), save_dir=save_dir_case)
    D5_case_GLOB_ref = D5_case(pyCOMPO=CPO_GLOB_ref, dlib_name=draglib, bu_points=BU_points,
                                ssh_opt="PT", correlation=correlation, sat="NODI", depl_sol="EXTR", 
                                tracked_nuclides=tracked_nuclides, BU_lists=getLists(BU_points), save_dir=save_dir_case)
    D5_case_NOGL_modif = D5_case(pyCOMPO=CPO_NOGL_modif, dlib_name=draglib, bu_points=BU_points,
                                ssh_opt="PT", correlation=correlation, sat="NODI", depl_sol="EXTR", 
                                tracked_nuclides=tracked_nuclides, BU_lists=getLists(BU_points), save_dir=save_dir_case)
    D5_case_GLOB_modif = D5_case(pyCOMPO=CPO_GLOB_modif, dlib_name=draglib, bu_points=BU_points,
                                ssh_opt="PT", correlation=correlation, sat="NODI", depl_sol="EXTR", 
                                tracked_nuclides=tracked_nuclides, BU_lists=getLists(BU_points), save_dir=save_dir_case)
    

    # Create D5 comparisons :

    comparison_NOGL = create_D5_comparison(D5_case_NOGL_modif, D5_case_NOGL_ref)
    comparison_GLOB = create_D5_comparison(D5_case_GLOB_modif, D5_case_GLOB_ref)

    # plot the results
    plt.figure(figsize=(12, 6))
    plt.plot(comparison_NOGL['BU Points'], comparison_NOGL['delta_keff'], marker='o', label='NOGL', linestyle='--', linewidth=.5)
    plt.plot(comparison_GLOB['BU Points'], comparison_GLOB['delta_keff'], marker='x', label='GLOB', linestyle='--', linewidth=.5)
    plt.title('Comparison NOGL - Modif vs Ref')
    plt.xlabel('BU Points')
    plt.ylabel('Delta keff (pcm)')
    # plot +/- 300 pcm lines
    plt.axhline(y=300, color='r', linestyle='--', label='+300 pcm')
    plt.axhline(y=-300, color='r', linestyle='--', label='-300 pcm')
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_case}/comparison_keff_NOGL_GLOB_{draglib}.png")
    plt.close()
    #plot_NOGL = comparison_NOGL.plot(x='BU Points', y='delta_keff', kind='bar', title='Comparison NOGL - Modif vs Ref')

    for isotope in tracked_nuclides:
        plt.figure(figsize=(12, 6))
        plt.plot(comparison_NOGL['BU Points'], comparison_NOGL[f'delta {isotope} (%)'], marker='o', label=f'NOGL {isotope}', linestyle='--', linewidth=.5)
        plt.plot(comparison_GLOB['BU Points'], comparison_GLOB[f'delta {isotope} (%)'], marker='x', label=f'GLOB {isotope}', linestyle='--', linewidth=.5)
        plt.title(f'Comparison NOGL - Modif vs Ref for {isotope}')
        plt.xlabel('BU Points')
        plt.ylabel(f'Delta {isotope} (%)')
        # plot +2/-2% lines
        plt.axhline(y=2, color='r', linestyle='--', label='+2%')
        plt.axhline(y=-2, color='r', linestyle='--', label='-2%')
        plt.legend()
        plt.grid()
        plt.savefig(f"{save_dir_case}/comparison_{isotope}_NOGL_GLOB_{draglib}.png")
        plt.close()

    # Create Serpent2 cases :
    S2_case_edep0_no_qfiss = S2_case(case_name=name_S2_case, lib_name="endfb8r1_pynjoy2012_kerma", 
                                    edep_id=0, areQfissSet=False, isEcaptSet=False, pcc_id=1, 
                                    specific_power=26.5, tracked_nuclides=tracked_nuclides, save_dir=save_dir_case)
    S2_case_edep0_qfiss = S2_case(case_name=name_S2_case, lib_name="endfb8r1_pynjoy2012_kerma",
                                edep_id=0, areQfissSet=True, isEcaptSet=False, pcc_id=1, 
                                specific_power=26.5, tracked_nuclides=tracked_nuclides, save_dir=save_dir_case)
    
    S2_case_edep2 = S2_case(case_name=name_S2_case, lib_name="endfb8r1_pynjoy2012_kerma",
                            edep_id=2, areQfissSet=False, isEcaptSet=False, pcc_id=1, 
                            specific_power=26.5, tracked_nuclides=tracked_nuclides, save_dir=save_dir_case)
    
    # Create OpenMC cases :
    OpenMC_case_edep0 = pd.read_csv(f"{path_to_OpenMC_data}/df_openmc_deplete_results_AT10_45Gd_default_fissq_CELI.csv")
    OpenMC_case_edep0_qfiss = pd.read_csv(f"{path_to_OpenMC_data}/df_openmc_deplete_results_AT10_45Gd_set_fissq_CELI.csv")
    OpenMC_case_edep2 = pd.read_csv(f"{path_to_OpenMC_data}/df_openmc_deplete_results_AT10_45Gd_energy_deposition_CELI.csv")

    

    # Create S2 comparisons :
    # Compare to S2 edep0 with default Qfiss
    comparison_NOGL_modif_S2 = create_D5_S2_comparison(D5_case_NOGL_modif, S2_case_edep0_no_qfiss)
    comparison_NOGL_S2 = create_D5_S2_comparison(D5_case_NOGL_ref, S2_case_edep0_no_qfiss)
    comparison_GLOB_modif_S2 = create_D5_S2_comparison(D5_case_GLOB_modif, S2_case_edep0_no_qfiss)
    comparison_GLOB_S2 = create_D5_S2_comparison(D5_case_GLOB_ref, S2_case_edep0_no_qfiss)
    print(f"delta_keff vs S2_case_edep0: {comparison_NOGL_modif_S2['delta_keff'].values}")
    print(f"Initial delta keff for modif NOGL vs S2 edep0 : {comparison_NOGL_modif_S2['delta_keff'][0]} pcm")
    print(f"Initial delta keff for original NOGL vs S2 edep0 : {comparison_NOGL_S2['delta_keff'][0]} pcm")
    print(f"Initial delta keff for modif GLOB vs S2 edep0 : {comparison_GLOB_modif_S2['delta_keff'][0]} pcm")
    print(f"Initial delta keff for original GLOB vs S2 edep0 : {comparison_GLOB_S2['delta_keff'][0]} pcm")

    # create D5 vs OpenMC comparisons
    comparison_NOGL_modif_OpenMC = create_D5_OpenMC_comparison(D5_case_NOGL_modif, OpenMC_case_edep0, OpenMC_tracked_nuclides)
    comparison_NOGL_OpenMC = create_D5_OpenMC_comparison(D5_case_NOGL_ref, OpenMC_case_edep0, OpenMC_tracked_nuclides)

    # plot the results
    plt.figure(figsize=(12, 6))
    plt.plot(comparison_NOGL_modif_S2['BU Points'], comparison_NOGL_modif_S2['delta_keff'], marker='o', label='NOGL, modified vs S2', linestyle='--', linewidth=.5)
    plt.plot(comparison_NOGL_S2['BU Points'], comparison_NOGL_S2['delta_keff'], marker='o', label='NOGL vs S2', linestyle='--', linewidth=.5)
    plt.plot(comparison_NOGL_modif_OpenMC['BU Points'], comparison_NOGL_modif_OpenMC['delta_keff'], marker='x', label='NOGL, modified vs OpenMC', linestyle='--', linewidth=.5)
    plt.plot(comparison_NOGL_OpenMC['BU Points'], comparison_NOGL_OpenMC['delta_keff'], marker='x', label='NOGL vs OpenMC', linestyle='--', linewidth=.5)
    #plt.plot(comparison_GLOB_S2['BU Points'], comparison_GLOB_S2['delta_keff'], marker='x', label='GLOB', linestyle='--', linewidth=.5)
    plt.title('Keff comparison D5 vs S2/OpenMC default fission Q-values')
    plt.xlabel('BU Points')
    plt.ylabel('Delta keff (pcm)')
    # plot +/- 300 pcm lines
    plt.axhline(y=300, color='r', linestyle='--', label='+300 pcm')
    plt.axhline(y=-300, color='r', linestyle='--', label='-300 pcm')
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_case}/comparison_keff_NOGL_GLOB_{draglib}_S2_edep0_OpenMC_default_qfiss.png")
    plt.close()
    
    for isotope in tracked_nuclides:
        plt.figure(figsize=(12, 6))
        plt.plot(comparison_NOGL_modif_S2['BU Points'], comparison_NOGL_modif_S2[f'delta {isotope} (%)'], marker='o', label=f'NOGL modif {isotope}', linestyle='--', linewidth=.5)
        plt.plot(comparison_NOGL_S2['BU Points'], comparison_NOGL_S2[f'delta {isotope} (%)'], marker='o', label=f'NOGL {isotope}', linestyle='--', linewidth=.5)
        if isotope in OpenMC_tracked_nuclides:
            plt.plot(comparison_NOGL_modif_OpenMC['BU Points'], comparison_NOGL_modif_OpenMC[f'delta {isotope} (%)'], marker='x', label=f'NOGL modif {isotope} vs OpenMC', linestyle='--', linewidth=.5)
            plt.plot(comparison_NOGL_OpenMC['BU Points'], comparison_NOGL_OpenMC[f'delta {isotope} (%)'], marker='x', label=f'NOGL {isotope} vs OpenMC', linestyle='--', linewidth=.5)
        #plt.plot(comparison_GLOB_modif_S2['BU Points'], comparison_GLOB_modif_S2[f'delta {isotope} (%)'], marker='x', label=f'GLOB modif {isotope}', linestyle='--', linewidth=.5)
        #plt.plot(comparison_GLOB_S2['BU Points'], comparison_GLOB_S2[f'delta {isotope} (%)'], marker='x', label=f'GLOB {isotope}', linestyle='--', linewidth=.5)
        plt.title(f'$\\Delta$ N{isotope} D5 vs S2 / OpenMC, edep=0 default fission Q-values')
        plt.xlabel('BU Points')
        plt.ylabel(f'Delta {isotope} (%)')
        # plot +2/-2% lines
        plt.axhline(y=2, color='r', linestyle='--', label='+2%')
        plt.axhline(y=-2, color='r', linestyle='--', label='-2%')
        plt.legend()
        plt.grid()
        plt.savefig(f"{save_dir_case}/comparison_{isotope}_NOGL_GLOB_{draglib}_S2_edep0_OpenMC_default_qfiss.png")
        plt.close()

    # Compare to S2 edep0 with Qfiss set
    comparison_NOGL_modif_S2 = create_D5_S2_comparison(D5_case_NOGL_modif, S2_case_edep0_qfiss)
    comparison_NOGL_S2 = create_D5_S2_comparison(D5_case_NOGL_ref, S2_case_edep0_qfiss)
    comparison_GLOB_modif_S2 = create_D5_S2_comparison(D5_case_GLOB_modif, S2_case_edep0_qfiss)
    comparison_GLOB_S2 = create_D5_S2_comparison(D5_case_GLOB_ref, S2_case_edep0_qfiss)
    
    print(f"Initial delta keff for modif NOGL vs S2 edep0 with setQfiss : {comparison_NOGL_modif_S2['delta_keff'][0]} pcm")
    print(f"Initial delta keff for original NOGL vs S2 edep0 setQfiss : {comparison_NOGL_S2['delta_keff'][0]} pcm")
    print(f"Initial delta keff for modif GLOB vs S2 edep0 setQfiss : {comparison_GLOB_modif_S2['delta_keff'][0]} pcm")
    print(f"Initial delta keff for original GLOB vs S2 edep0 setQfiss : {comparison_GLOB_S2['delta_keff'][0]} pcm")

    # Create OpenMC cases with Qfiss set:
    comparison_NOGL_modif_OpenMC = create_D5_OpenMC_comparison(D5_case_NOGL_modif, OpenMC_case_edep0_qfiss, OpenMC_tracked_nuclides)
    comparison_NOGL_OpenMC = create_D5_OpenMC_comparison(D5_case_NOGL_ref, OpenMC_case_edep0_qfiss, OpenMC_tracked_nuclides)

    # plot the results
    print(f"delta_keff vs S2_case_edep0_qfiss : {comparison_NOGL_modif_S2['delta_keff'].values}")
    plt.figure(figsize=(12, 6))
    plt.plot(comparison_NOGL_modif_S2['BU Points'], comparison_NOGL_modif_S2['delta_keff'], marker='o', label='NOGL modified vs S2', linestyle='--', linewidth=.5)
    plt.plot(comparison_NOGL_S2['BU Points'], comparison_NOGL_S2['delta_keff'], marker='o', label='NOGL vs S2', linestyle='--', linewidth=.5)
    plt.plot(comparison_NOGL_modif_OpenMC['BU Points'], comparison_NOGL_modif_OpenMC['delta_keff'], marker='x', label='NOGL modified vs OpenMC', linestyle='--', linewidth=.5)
    plt.plot(comparison_NOGL_OpenMC['BU Points'], comparison_NOGL_OpenMC['delta_keff'], marker='x', label='NOGL vs OpenMC', linestyle='--', linewidth=.5)
    plt.title('Keff comparisons D5 vs S2/OpenMC, edep=0 with set fission Q-values')
    plt.xlabel('BU Points')
    plt.ylabel('Delta keff (pcm)')
    # plot +/- 300 pcm lines
    plt.axhline(y=300, color='r', linestyle='--', label='+300 pcm')
    plt.axhline(y=-300, color='r', linestyle='--', label='-300 pcm')
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_case}/comparison_keff_NOGL_GLOB_{draglib}_S2_edep0_setQfiss_OpenMC_set_fissq.png")
    plt.close()
    
    for isotope in tracked_nuclides:
        plt.figure(figsize=(12, 6))
        plt.plot(comparison_NOGL_modif_S2['BU Points'], comparison_NOGL_modif_S2[f'delta {isotope} (%)'], marker='o', label=f'NOGL modif {isotope}', linestyle='--', linewidth=.5)
        plt.plot(comparison_NOGL_S2['BU Points'], comparison_NOGL_S2[f'delta {isotope} (%)'], marker='o', label=f'NOGL {isotope}', linestyle='--', linewidth=.5)
        if isotope in OpenMC_tracked_nuclides:
            plt.plot(comparison_NOGL_modif_OpenMC['BU Points'], comparison_NOGL_modif_OpenMC[f'delta {isotope} (%)'], marker='x', label=f'NOGL modif {isotope} vs OpenMC', linestyle='--', linewidth=.5)
            plt.plot(comparison_NOGL_OpenMC['BU Points'], comparison_NOGL_OpenMC[f'delta {isotope} (%)'], marker='x', label=f'NOGL {isotope} vs OpenMC', linestyle='--', linewidth=.5)
        plt.title(f'$\\Delta$ N{isotope} D5 vs S2/OpenMC, edep=0 with set fission Q-values')
        plt.xlabel('BU Points')
        plt.ylabel(f'Delta {isotope} (%)')
        # plot +2/-2% lines
        plt.axhline(y=2, color='r', linestyle='--', label='+2%')
        plt.axhline(y=-2, color='r', linestyle='--', label='-2%')
        plt.legend()
        plt.grid()
        plt.savefig(f"{save_dir_case}/comparison_{isotope}_NOGL_GLOB_{draglib}_S2_edep0_setQfiss_OpenMC_set_fissq.png")
        plt.close()

    # Create S2 comparisons with edep=2:
    comparison_NOGL_modif_S2 = create_D5_S2_comparison(D5_case_NOGL_modif, S2_case_edep2)
    comparison_NOGL_S2 = create_D5_S2_comparison(D5_case_NOGL_ref, S2_case_edep2)
    comparison_GLOB_modif_S2 = create_D5_S2_comparison(D5_case_GLOB_modif, S2_case_edep2)
    comparison_GLOB_S2 = create_D5_S2_comparison(D5_case_GLOB_ref, S2_case_edep2)
    # Create OpenMC cases with edep=2:
    #OpenMC_case_edep2 = pd.read_csv(f"{path_to_OpenMC_data}/df_openmc_deplete_results_AT10_45Gd_energy_deposition_CELI.csv")
    comparison_NOGL_modif_OpenMC = create_D5_OpenMC_comparison(D5_case_NOGL_modif, OpenMC_case_edep2, OpenMC_tracked_nuclides)
    comparison_NOGL_OpenMC = create_D5_OpenMC_comparison(D5_case_NOGL_ref, OpenMC_case_edep2, OpenMC_tracked_nuclides)
    
    print(f"delta_keff vs S2_case_edep2 : {comparison_NOGL_modif_S2['delta_keff'].values}")
    print(f"Initial delta keff for modif NOGL vs S2 edep2 : {comparison_NOGL_modif_S2['delta_keff'][0]} pcm")
    print(f"Initial delta keff for original NOGL vs S2 edep2 : {comparison_NOGL_S2['delta_keff'][0]} pcm")
    print(f"Initial delta keff for modif GLOB vs S2 edep2 : {comparison_GLOB_modif_S2['delta_keff'][0]} pcm")
    print(f"Initial delta keff for original GLOB vs S2 edep2 : {comparison_GLOB_S2['delta_keff'][0]} pcm")

    # plot the results
    plt.figure(figsize=(12, 6))
    plt.plot(comparison_NOGL_modif_S2['BU Points'], comparison_NOGL_modif_S2['delta_keff'], marker='o', label='NOGL, modified', linestyle='--', linewidth=.5)
    plt.plot(comparison_NOGL_S2['BU Points'], comparison_NOGL_S2['delta_keff'], marker='o', label='NOGL', linestyle='--', linewidth=.5)
    plt.plot(comparison_NOGL_modif_OpenMC['BU Points'], comparison_NOGL_modif_OpenMC['delta_keff'], marker='x', label='NOGL, modified vs OpenMC', linestyle='--', linewidth=.5)
    plt.plot(comparison_NOGL_OpenMC['BU Points'], comparison_NOGL_OpenMC['delta_keff'], marker='x', label='NOGL vs OpenMC', linestyle='--', linewidth=.5)
    # plot +/- 300 pcm lines
    plt.axhline(y=300, color='r', linestyle='--', label='+300 pcm')
    plt.axhline(y=-300, color='r', linestyle='--', label='-300 pcm')
    plt.title('Comparison NOGL - GLOB vs S2, edep=2')
    plt.xlabel('BU Points')
    plt.ylabel('Delta keff (pcm)')
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_case}/comparison_keff_NOGL_GLOB_{draglib}_S2_edep2_OpenMC_energy_deposition.png")
    plt.close()
    
    for isotope in tracked_nuclides:
        plt.figure(figsize=(12, 6))
        plt.plot(comparison_NOGL_modif_S2['BU Points'], comparison_NOGL_modif_S2[f'delta {isotope} (%)'], marker='o', label=f'NOGL modif {isotope}', linestyle='--', linewidth=.5)
        plt.plot(comparison_NOGL_S2['BU Points'], comparison_NOGL_S2[f'delta {isotope} (%)'], marker='o', label=f'NOGL {isotope}', linestyle='--', linewidth=.5)
        if isotope in OpenMC_tracked_nuclides: 
            plt.plot(comparison_NOGL_modif_OpenMC['BU Points'], comparison_NOGL_modif_OpenMC[f'delta {isotope} (%)'], marker='x', label=f'NOGL modif {isotope} vs OpenMC', linestyle='--', linewidth=.5)
            plt.plot(comparison_NOGL_OpenMC['BU Points'], comparison_NOGL_OpenMC[f'delta {isotope} (%)'], marker='x', label=f'NOGL {isotope} vs OpenMC', linestyle='--', linewidth=.5)
        plt.title(f'Comparison NOGL - GLOB vs S2 for {isotope}, edep=2')
        plt.xlabel('BU Points')
        plt.ylabel(f'Delta {isotope} (%)')
        # plot +2/-2% lines
        plt.axhline(y=2, color='r', linestyle='--', label='+2%')
        plt.axhline(y=-2, color='r', linestyle='--', label='-2%')
        plt.legend()
        plt.grid()
        plt.savefig(f"{save_dir_case}/comparison_{isotope}_NOGL_GLOB_{draglib}_S2_edep2.png")
        plt.close()

    