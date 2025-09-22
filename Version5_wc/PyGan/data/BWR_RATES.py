##############################################################################
#                                                                            #
# Description : Python3 scritp for BWR test cases                            #
# Author      : R. Guasch, adapted from Lucas Fede and Mathias Francois      #
# Date        : August 2024                                                  #
#                                                                            #
##############################################################################
#
# ----------------------------------------------------------------------------
#                                  IMPORTS                                   |
# ----------------------------------------------------------------------------

import os, shutil, sys
import pickle
import lifo
import lcm
import cle2000
from collections import defaultdict
import random as rd
import numpy as np
import matplotlib
import matplotlib.cm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from typing import List, Union
import serpentTools as st
from serpentTools.settings import rc
from Serpent2_parsing import parse_S2_ASSBLY_rates, parse_Serpent2_material_det, parse_Serpent2_lattice_det, parse_S2_ASSBLY_rates_lat_det, parse_S2_pin_mat_det
from Dragon5_cpo_parsing import parse_DRAGON_rates_enrich_num, parse_DRAGON_DIAG_rates_regi_num, parse_DRAGON_SCHEME
from plotters import plot_errors_BWR_assembly,  plot_BWR_assembly, plot_pinwise_errors_BWR_assembly

matplotlib.use('Agg')


def renormalize_rates(rates):
    """
    rates is the list of rates values to renormalize.
    They correspond to reaction rates vs cells in the geometry in a single energy group.
    returns rates normalized to nCells == length of rates.
    """
    print(f"nCells = {len(rates)}")
    print(f"rates = {rates}")
    rates = np.array(rates)
    return rates/np.sum(rates)


def renormalize_rates_dict(rates_dict):
    """
    rates_dict is a dictionary with the following structure:
        "mix1": {
            [gr1_value, gr2_value, ...],
            ...
        },
        ...
    },

    This function renormalizes the values for each isotope and mix/group.
    Returns a dictionary of the same structure with renormalized values.
    """
    normalized_dict = {}

    for mix, group_values in rates_dict.items():
            total = sum(group_values)
            if total != 0:
                normalized_values = [val / total for val in group_values]
            else:
                normalized_values = group_values[:]  # Leave unchanged if total is zero
            normalized_dict[mix] = normalized_values

    return normalized_dict


def BWR_CLUSTER(name_case, name_compo, reaction_type, n_groups, bu):

    ################## FILE MANAGEMENT AND COMPOS ##################

    # --- Current directory path
    path = os.getcwd()
    # List of all the files in the current directory
    files = os.listdir(path)
    #print(files)
    #print('$----------------------------'+path)
    
    # --- Creating the results storage directory
    name_DIR = f'DRAGON_RATES_{name_case}'
    a = os.path.exists(name_DIR)
    if not a:
        os.mkdir(name_DIR)
    
    # Access paths to multicompo and serpent results
    SERPENT_path = f'{os.environ["SERPENT_RESULTS"]}'
    # Recovering multi-compos
    #name_compo = 'COMPO_' + name_case + name_case 
    print(f"name_compo = {name_compo}")
    os.chdir("Linux_aarch64/")
    #print(name_compo)
    pyCOMPO = lcm.new('LCM_INP', name_compo, impx=0)
    os.chdir(path)

    if name_case == "AT10_2x2_UOX":
        """
        C1 C2
        C2 C4
        """
        N_col=2
        nCell = 4
        MIXES_idx = [0,1,2]
        cell_S2_to_mix_indx = {0:0, 1:1, 2:1, 3:2}
        MIXES = ["C1", "C2", "C4"]
        number_of_each_mix = {"C1":1, "C2":2, "C4":1}
    elif name_case == "bench_3x3_UOX":
        """
            C1 C2 C3    cell1 cell2 cell3
            C2 C4 C1 =  cell4 cell5 cell6
            C3 C1 C6    cell7 cell8 cell9
        """
        N_col=3
        nCell = 9   # 1,2,3,4,5 = C1, C2, C3, C4, C6
        MIXES_idx = [0,1,2,3,4]
        cell_S2_to_mix_indx = {0:0, 1:1, 2:2, 3:1, 4:3, 5:0, 6:2, 7:0, 8:4}
        reg_idx_to_cell = {0 : "C1", 1: "C2", 2 : "C3", 3: "C2", 4: "C4", 5: "C1", 6: "C3", 7: "C1", 8: "C6"}
        MIXES = ["C1", "C2", "C3", "C4", "C6"]    
        number_of_each_mix = {"C1":3, "C2":2, "C3":2, "C4":1, "C6":1}
    elif name_case == "AT10_3x3_UOX_Gd":
        """
            C1 C2 C3    cell1 cell2 cell3
            C2 C4 C7 =  cell4 cell5 cell6
            C3 C7 C6    cell7 cell8 cell9
        """
        N_col=3
        nCell = 9   # C1, C2, C3, C4, C6, C7
        MIXES_idx = [0,1,2,3,4,5]
        cell_S2_to_mix_indx = {0:0, 1:1, 2:2, 3:1, 4:3, 5:5, 6:2, 7:5, 8:4}
        reg_idx_to_cell = {0 : "C1", 1: "C2", 2 : "C3", 3: "C2", 4: "C4", 5: "C7", 6: "C3", 7: "C7", 8: "C6"}
        MIXES = ["C1", "C2", "C3", "C4", "C6", "C7"]
        number_of_each_mix = {"C1":1, "C2":2, "C3":2, "C4":1, "C6":1, "C7":2}

    # Isotopes souhaites et nb total d'isotopes
    iso_study = ['U235 ', 'U238 ', 'Pu239', 'Pu241']
    reaction_name = 'NFTOT'
    n_reac = 6
    rates_names = 'fission rates'
    #pyCOMPO['EDIBU'].lib()
    #print(pyCOMPO.keys())
    len_isotot = np.shape(pyCOMPO['EDIBU_HOM']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESDENS'])[0] - 1
    ########## CALCULATIONS ##########
    # Retrieve keff from pyCOMPO
    keff_D5 = pyCOMPO['EDIBU_HOM']['MIXTURES'][0]['CALCULATIONS'][0]['K-EFFECTIVE']
    # Calcul du facteur de normalisation
    prodD5 = 0.0
    prodS2 = 0.0
    ListS2 = []    
    

    # Create index to ALIAS dictionary
    Iso_index_to_ALIAS = {}
    U235_fiss_rate = {}
    for iso in range(len_isotot):
        isotope = pyCOMPO['EDIBU_HOM']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][iso]['ALIAS'][0:5]
        #print(f"isotope = {isotope}, iso number = {iso}")
        Iso_index_to_ALIAS[iso] = isotope
    for iso in range(len_isotot):
        #print(f"iso index = {iso}, isotope = {Iso_index_to_ALIAS[iso]}")
        isotope = pyCOMPO['EDIBU_HOM']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][iso]['ALIAS'][0:5]
        #print(f"isotope = {isotope}")
        if isotope in ['U235 ']: #, 'U238 ']: #, 'Pu239', 'Pu241']:
            for mix in MIXES_idx:
                U235_fiss_rate[f"mix{mix+1}"] = {}
                print(f"mix index is = {mix} with corresponding C name = {MIXES[mix]}")
                print(f"isotope is = {isotope}")
                print("burnup is = ", bu)
                print(pyCOMPO['HOM2g']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso].keys())
                NWT0 = pyCOMPO['HOM2g']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NWT0']
                N = pyCOMPO['HOM1g']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESDENS'][iso]
                vol = pyCOMPO['HOM1g']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESVOL'][iso]
                print(f"for mix {mix+1}, N = {N}, vol = {vol}")
                NFTOT = pyCOMPO['HOM2g']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NFTOT']
                #NGAMMA = pyCOMPO['EDIBU_2gr']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NG']
                print(f"mix index is = {mix} with corresponding C name = {MIXES[mix]}")
                print(f"Number of each mix = {number_of_each_mix[MIXES[mix]]}")
                for gr in range(n_groups):
                    U235_fiss_rate[f"mix{mix+1}"][f"gr{gr+1}"] = NFTOT[gr]*NWT0[gr]*N*vol/number_of_each_mix[MIXES[mix]]

    print(f"U235_fiss_rate = {U235_fiss_rate}")
    U235_fiss_rates_regi = {}
    # duplicate rates from individual mixes to form the cells lattice
    for i in range(nCell):
        U235_fiss_rates_regi[f"reg{i+1}"] = {}
        for gr in range(n_groups):
            U235_fiss_rates_regi[f"reg{i+1}"][f"G{gr+1}"] = U235_fiss_rate[f"mix{cell_S2_to_mix_indx[i]+1}"][f"gr{gr+1}"]


    
    
    
    # Parse Serpent2 results
    keff_S2, U235_material_det_fission_Rates_S2 = parse_Serpent2_material_det(SERPENT_path, name_case, XS_lib_S2="endfb8r1_pynjoy2012", bcond=2, edepmode=0, pcc=1, bu=0)
    
    print(f"U235_material_det_fission_Rates_S2 = {U235_material_det_fission_Rates_S2}")
    # normalize the rates : sum of rates = nCells
    # group 1 : U235 fast fissions : note that group numbers are in increasing lethargy in Dragon and in increasing energy in Serpent
    U235_fiss_rate_norm_1 = renormalize_rates([U235_fiss_rates_regi[reg]["G1"] for reg in U235_fiss_rates_regi.keys()])

    # group 2 : U235 thermal fissions
    U235_fiss_rate_norm_2 = renormalize_rates([U235_fiss_rates_regi[reg]["G2"] for reg in U235_fiss_rates_regi.keys()])

    print(f"for bcond = 2, name_compo = {name_compo}, delta keff = {(keff_D5 - keff_S2)*1e5}")


    ### Material detector
    # group 1 : U235 fast fissions : note that group numbers are in increasing lethargy in Dragon and in increasing energy in Serpent
    U235_fission_rates_S2_mat_det_norm_1 = renormalize_rates([U235_material_det_fission_Rates_S2[cell]["G2"] for cell in U235_material_det_fission_Rates_S2.keys()])
    # group 2 : U235 thermal fissions
    U235_fission_rates_S2_mat_det_norm_2 = renormalize_rates([U235_material_det_fission_Rates_S2[cell]["G1"] for cell in U235_material_det_fission_Rates_S2.keys()])

    # relative difference
    relative_difference1_mat_det = [(U235_fiss_rate_norm_1[i] - U235_fission_rates_S2_mat_det_norm_1[i])*100/U235_fission_rates_S2_mat_det_norm_1[i] for i in range(nCell)]
    relative_difference2_mat_det = [(U235_fiss_rate_norm_2[i] - U235_fission_rates_S2_mat_det_norm_2[i])*100/U235_fission_rates_S2_mat_det_norm_2[i] for i in range(nCell)]

    plot_error_grid_from_list(f"{name_compo}_bc2_error_grid_U235_fission_rates_g1_vs_mat_det", name_case, name_compo, 1, relative_difference1_mat_det, N_col)
    plot_error_grid_from_list(f"{name_compo}_bc2_error_grid_U235_fission_rates_g2_vs_mat_det", name_case, name_compo, 2, relative_difference2_mat_det, N_col)

    return


def compute_diff_rates(D5_rates, S2_rates):
    """
    D5_rates (dict) : Dragon rates, mix as key, array of ngroups rates as values
    S2_rates (dict) : Serpent rates, mix as key, array of ngroups rates as values

    !! Warning : S2 : increasing in energy, D5 : increasing in lethargy ==> decreasing in energy
    Invert the order of the groups in S2 rates to compare with D5 rates 

    Return:
        List of dictionaries. Index in the list = group number.
        Each dictionary maps mix name to the relative difference (%) in that group.
    """

    delta_rates_dict = {}
    key1 = next(iter(D5_rates))  # Get any mix to determine number of groups
    n_groups = len(D5_rates[key1])

    # Compute per-mix relative difference arrays
    for mix in D5_rates:
        d5 = np.array(D5_rates[mix])
        s2 = np.array(S2_rates[mix][::-1])  # Reverse group order in Serpent data
        delta_rates_dict[mix] = (d5 - s2) * 100 / s2

    # Convert to list of dictionaries by group
    delta_rates = []
    for i in range(n_groups):
        group_dict = {}
        for mix, delta_array in delta_rates_dict.items():
            group_dict[mix] = delta_array[i]
        delta_rates.append(group_dict)

    return delta_rates
# Execute post treatment procedure

def analyze_main_flux_geometry_refinement():
    """
    Analyze the impact of geometry refinement on flux and rates accuracy for a BWR assembly.
    Compares Dragon5 results with Serpent2 reference.
    """
#### ATRIUM-10 assembly
    name_case = "ATRIUM10"
    # geometry_type : choice of discretisation for flux calculation geometry.
    #geometry_type = "default" # "cool_ring_SECT40" # "default", "fine1", "corners1"
    geometry_refinement_options = ["finest_geom"]#["default", "fine1", "fine2", "coolant_ring", "finest_on_Gd_coolant_ring", "finest_on_Gd", "cool_ring_SECT40"] # "finest_geom", 
    
    composition_option = "AT10_void_0"
    evaluation = "J311_295"
    ssh_methods = ["RSE"]#, "PT"]
    correlation_options = ["NOCORR"]#, "CORR"]
    anisotropy_level = 2

    for geometry_type in geometry_refinement_options:
        if geometry_type in ["default", "fine1", "fine2", "coolant_ring"]:
            refinement_opt_name = geometry_type
            num_angles = 24
            line_density = 140.0
            batch = 750
            
        elif geometry_type in ["finest_on_Gd", "finest_on_Gd_coolant_ring"]:
            refinement_opt_name = geometry_type
            num_angles = 24
            line_density = 140.0
            if geometry_type == "finest_on_Gd_coolant_ring":
                batch = 2000
            elif geometry_type == "finest_on_Gd":
                batch = 1000
        elif geometry_type in ["finest_geom", "cool_ring_SECT40"]:
            refinement_opt_name = geometry_type
            num_angles = 24
            line_density = 150.0
            batch = 6000
            if geometry_type == "cool_ring_SECT40":
                ssh_methods = ["PT"]
                batch = 3000
        
        CPO_name = f"CPO_n{num_angles}_ld{int(line_density)}_n8_ld25_TSPC_{anisotropy_level}_MOC_GAUS_4_{batch}_200"
        CPO_name = "CPO_n24_ld150_n8_ld25_TSPC_2_MOC_GAUS_4_6000_200"
        print(f"geometry_type : {geometry_type}, refinement_opt_name : {refinement_opt_name}, CPO_name = {CPO_name}")
        deltas = {}
        for ssh_method in ssh_methods:
            if ssh_method not in deltas.keys():
                deltas[ssh_method] = {}
            for correlation_option in correlation_options:
                if correlation_option not in deltas[ssh_method].keys():
                    deltas[ssh_method][correlation_option] = {}
                keff_D5, fiss_rates_D5, n_gamma_rates_D5, FLUX_295groups_D5 = parse_DRAGON_DIAG_rates_regi_num(name_case, CPO_name, composition_option, evaluation, ssh_method, correlation_option, geometry_type, fission_isotopes = ["U235", "U238"] , n_gamma_isotopes = ["U238", "Gd155", "Gd157"], bu=0)
                #S2_results_test = parse_S2_ASSBLY_rates_lat_det(name_case="AT10_ASSBLY_t0_Rates", XS_lib_S2="J311_pynjoy2016", fission_isotopes=["U235", "U238", "U234"], n_gamma_isotopes=["U238", "Gd155", "Gd157"], bu=0)
                keff_S2, fission_rates_S2, ngamma_rates_S2 = parse_S2_pin_mat_det(name_case="ATRIUM10_pin", XS_lib_S2="jeff311_pynjoy2016", fission_isotopes=["U235", "U238"], ngamma_isotopes=["U238", "Gd155", "Gd157"], bu=0)


                print(fission_rates_S2)
                # Normalise to nCells = 49 cells with fissile material.
                fiss_rates_D5[0] =  fiss_rates_D5[0] * 49 / np.sum(fiss_rates_D5[0])
                fission_rates_S2[0] = fission_rates_S2[0] * 49 / np.sum(fission_rates_S2[0])
                fiss_rates_D5[1] =  fiss_rates_D5[1] * 49 / np.sum(fiss_rates_D5[1])
                fission_rates_S2[1] = fission_rates_S2[1] * 49 / np.sum(fission_rates_S2[1])

                delta_therm = fiss_rates_D5[0] - fission_rates_S2[0]
                delta_fast = fiss_rates_D5[1] - fission_rates_S2[1]
                
                delta_rel_therm = [(fiss_rates_D5[0][i] - fission_rates_S2[0][i])*100 / fission_rates_S2[0][i] if fission_rates_S2[0][i]!=0 else 0.0 for i in range(len(fission_rates_S2[0]))]
                delta_rel_fast = [(fiss_rates_D5[1][i] - fission_rates_S2[1][i])*100 / fission_rates_S2[1][i] if fission_rates_S2[1][i]!=0 else 0.0 for i in range(len(fission_rates_S2[1]))]

                CPOS_to_treat = [] # ["CPO_AT10_ASSBLY_RSE_fine_J311", "CPO_AT10_ASSBLY_RSE_CORR_fine_J311"]
                print(f"For geometry name : {geometry_type}")
                print(f"keff_D5 = {keff_D5}, keff_S2 = {keff_S2}")
                delta_keff = (keff_D5 - keff_S2)*1e5
                print(f"error (D5-S2) on keff = {delta_keff:.1f}")
                deltas[ssh_method][correlation_option]["keff"] = delta_keff
                
                # Compute RMS, MAX and AVG errors on fission rates
                rms_therm = np.sqrt(np.mean(np.array(delta_rel_therm)**2))
                rms_fast = np.sqrt(np.mean(np.array(delta_rel_fast)**2))
                avg_therm = np.mean(np.abs(np.array(delta_rel_therm)))
                avg_fast = np.mean(np.abs(np.array(delta_rel_fast)))
                max_therm = np.max(np.abs(np.array(delta_rel_therm)))
                max_fast = np.max(np.abs(np.array(delta_rel_fast)))
                print(f"RMS errors on fission rates (thermal, fast) = {rms_therm:.2f} %, {rms_fast:.2f} %")
                print(f"AVG errors on fission rates (thermal, fast) = {avg_therm:.2f} %, {avg_fast:.2f} %")
                print(f"MAX errors on fission rates (thermal, fast) = {max_therm:.2f} %, {max_fast:.2f} %")
                deltas[ssh_method][correlation_option]["rms_therm"] = rms_therm
                deltas[ssh_method][correlation_option]["avg_therm"] = avg_therm
                deltas[ssh_method][correlation_option]["max_therm"] = max_therm
                deltas[ssh_method][correlation_option]["rms_fast"] = rms_fast
                deltas[ssh_method][correlation_option]["avg_fast"] = avg_fast
                deltas[ssh_method][correlation_option]["max_fast"] = max_fast
                

                plot_pinwise_errors_BWR_assembly(np.array([delta_rel_therm,delta_rel_fast]), name_case, CPO_name, calculation_opt = f"{geometry_type}_{ssh_method}_{correlation_option}", fig_name="Total_fission_rates_diff")


        print(deltas)
        
def analyze_computational_schemes():
    """
    Analyze the impact of different computational schemes on flux and rates accuracy for a BWR assembly.
    Compares Dragon5 results with Serpent2 reference.
    """
    fission_isotopes = ["U235", "U238"]
    #### ATRIUM-10 assembly
    computational_schemes = ["1L_MOC", "2L_PIJ_MOC", "2L_IC_MOC"]
    name_case = "ATRIUM10"
    geometry_type = "finest_geom" # "cool_ring_SECT40" # "default", "fine1", "corners1"
    composition_option = "AT10_void_0"
    evaluation = "J311_295"
    ssh_methods = ["PT"]
    ssh_sol_doors = ["IC", "PIJ"]
    correlation_option = "NOCORR"
    anisotropy_level = 4
    num_angles = 24
    line_density = 150.0
    batch = 3000
    refinement_opt_name = geometry_type
    deltas = {}
    
    for computational_scheme in computational_schemes:
        if computational_scheme not in deltas.keys():
            deltas[computational_scheme] = {}
        for ssh_method in ssh_methods:
            for ssh_sol in ssh_sol_doors:
                if f"{ssh_method}_{ssh_sol}" not in deltas[computational_scheme].keys():
                    deltas[computational_scheme][f"{ssh_method}_{ssh_sol}"] = {}
                CPO_name = f"CPO_n{num_angles}_ld{int(line_density)}_n8_ld25_TSPC_{anisotropy_level}_{ssh_sol}_MOC_GAUS_4_{batch}_200"
                keff_D5, fiss_rates_D5, n_gamma_rates_D5, FLUX_295groups_D5 = parse_DRAGON_SCHEME(name_case, CPO_name, composition_option, evaluation, ssh_method, correlation_option, geometry_type, fission_isotopes , n_gamma_isotopes = ["U238", "Gd155", "Gd157"], bu=0, computational_scheme=computational_scheme)
                keff_S2, fission_rates_S2, ngamma_rates_S2 = parse_S2_pin_mat_det(name_case="ATRIUM10_pin", XS_lib_S2="jeff311_pynjoy2016", fission_isotopes=fission_isotopes, ngamma_isotopes=["U238", "Gd155", "Gd157"], bu=0)


                print(fission_rates_S2)
                # Normalise to nCells = 49 cells with fissile material.
                fiss_rates_D5[0] =  fiss_rates_D5[0] * 49 / np.sum(fiss_rates_D5[0])
                fission_rates_S2[0] = fission_rates_S2[0] * 49 / np.sum(fission_rates_S2[0])
                fiss_rates_D5[1] =  fiss_rates_D5[1] * 49 / np.sum(fiss_rates_D5[1])
                fission_rates_S2[1] = fission_rates_S2[1] * 49 / np.sum(fission_rates_S2[1])

                delta_therm = fiss_rates_D5[0] - fission_rates_S2[0]
                delta_fast = fiss_rates_D5[1] - fission_rates_S2[1]
                
                delta_rel_therm = [(fiss_rates_D5[0][i] - fission_rates_S2[0][i])*100 / fission_rates_S2[0][i] if fission_rates_S2[0][i]!=0 else 0.0 for i in range(len(fission_rates_S2[0]))]
                delta_rel_fast = [(fiss_rates_D5[1][i] - fission_rates_S2[1][i])*100 / fission_rates_S2[1][i] if fission_rates_S2[1][i]!=0 else 0.0 for i in range(len(fission_rates_S2[1]))]

                print(f"For geometry name : {geometry_type}")
                print(f"For computational scheme : {computational_scheme}, ssh solution door : {ssh_sol}")
                print(f"keff_D5 = {keff_D5}, keff_S2 = {keff_S2}")
                delta_keff = (keff_D5 - keff_S2)*1e5
                print(f"error (D5-S2) on keff = {delta_keff:.1f}")
                deltas[computational_scheme][f"{ssh_method}_{ssh_sol}"]["keff"] = delta_keff

                # Compute RMS, MAX and AVG errors on fission rates
                rms_therm = np.sqrt(np.mean(np.array(delta_rel_therm)**2))
                rms_fast = np.sqrt(np.mean(np.array(delta_rel_fast)**2))
                avg_therm = np.mean(np.abs(np.array(delta_rel_therm)))
                avg_fast = np.mean(np.abs(np.array(delta_rel_fast)))
                max_therm = np.max(np.abs(np.array(delta_rel_therm)))
                max_fast = np.max(np.abs(np.array(delta_rel_fast)))
                print(f"RMS errors on fission rates (thermal, fast) = {rms_therm:.2f} %, {rms_fast:.2f} %")
                print(f"AVG errors on fission rates (thermal, fast) = {avg_therm:.2f} %, {avg_fast:.2f} %")
                print(f"MAX errors on fission rates (thermal, fast) = {max_therm:.2f} %, {max_fast:.2f} %")
                deltas[computational_scheme][f"{ssh_method}_{ssh_sol}"]["rms_therm"] = rms_therm
                deltas[computational_scheme][f"{ssh_method}_{ssh_sol}"]["avg_therm"] = avg_therm
                deltas[computational_scheme][f"{ssh_method}_{ssh_sol}"]["max_therm"] = max_therm
                deltas[computational_scheme][f"{ssh_method}_{ssh_sol}"]["rms_fast"] = rms_fast
                deltas[computational_scheme][f"{ssh_method}_{ssh_sol}"]["avg_fast"] = avg_fast
                deltas[computational_scheme][f"{ssh_method}_{ssh_sol}"]["max_fast"] = max_fast
                

                plot_pinwise_errors_BWR_assembly(np.array([delta_rel_therm,delta_rel_fast]), name_case, CPO_name, calculation_opt = f"{computational_scheme}_{geometry_type}_{ssh_method}_{ssh_sol}_{correlation_option}", fig_name="Total_fission_rates_diff")
    print(deltas)

if __name__ == "__main__":


    #analyze_main_flux_geometry_refinement()
    analyze_computational_schemes()