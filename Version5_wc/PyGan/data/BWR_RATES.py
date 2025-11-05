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
    # Test 2L_PIJ_noMM_MOC_SC : should be similar to 2L_PIJ_MOC_SC but longer as MERGMIX keyword not used.
    computational_schemes = ["2L_IC_MOC_SC"]#["2L_IC_MOC_SC", "2L_PIJ_MOC_SC"]#["2L_IC_MOC", "2L_IC_MOC_SC_GRMAX15", "2L_IC_MOC_SC_GRMAX16", "2L_IC_MOC_SC_GRMAX17", "2L_IC_MOC_SC_GRMAX18","2L_IC_MOC_SC_GRMAX19", "2L_IC_MOC_SC_GRMAX20", "2L_IC_MOC_SC_GRMAX21", "2L_IC_MOC_SC_GRMAX22", "2L_IC_MOC_SC_GRMAX23"]
    analyze_first_level = False # False, flag to indicate if first levelof calculation should be analysed then another COMPO object is opened and 
    name_case = "ATRIUM10"
    geometry_type = "finest_geom" # "cool_ring_SECT40" # "default", "fine1", "corners1"
    composition_option = "AT10_void_0"
    evaluation = "J311_295"
    ssh_methods = ["RSE"]#, "RSE"]#, "PT"]
    ssh_sol_doors = ["IC"]#["PIJ", "IC"]
    correlation_option = "NOCORR"
    anisotropy_level = 4
    transport_correction = "NONE" # "APOL"
    num_angles = 24
    line_density = 150.0
    batch = 2000
    refinement_opt_name = geometry_type
    deltas = {}
    
    for computational_scheme in computational_schemes:
        if computational_scheme not in deltas.keys():
            deltas[computational_scheme] = {}
        for ssh_method in ssh_methods:
            for ssh_sol in ssh_sol_doors:
                if f"{ssh_method}_{ssh_sol}" not in deltas[computational_scheme].keys():
                    deltas[computational_scheme][f"{ssh_method}_{ssh_sol}"] = {}
                
                if transport_correction == "NONE":
                    CPO_name = f"CPO_n{num_angles}_ld{int(line_density)}_n8_ld25_TSPC_{anisotropy_level}_{ssh_sol}_MOC_GAUS_4_{batch}_200"
                else:
                    CPO_name = f"CPO_n{num_angles}_ld{int(line_density)}_n8_ld25_TSPC_{anisotropy_level}_{ssh_sol}_MOC_GAUS_4_{batch}_200_{transport_correction}"
                print(f"computational_scheme : {computational_scheme}, CPO_name = {CPO_name}")
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
                print(f"{delta_keff:.0f} & {rms_fast:.2f}\\% & {rms_therm:.2f}\\% & {max_fast:.2f}\\% &  {max_therm:.2f}\\%")
                

                plot_pinwise_errors_BWR_assembly(np.array([delta_rel_therm,delta_rel_fast]), name_case, CPO_name, calculation_opt = f"{computational_scheme}_{geometry_type}_{ssh_method}_{ssh_sol}_{correlation_option}", fig_name="Total_fission_rates_diff")
                
                ## Analysis of 1st level in the calculation scheme : PIJ / IC computed keff and rates.
                if analyze_first_level and "2L" in computational_scheme:
                    print(f"Begin analysis of level 1 results")
                    CPO_name_level1 = f"LVL1__CPO_n{num_angles}_ld{int(line_density)}_n8_ld25_TSPC_{anisotropy_level}_{ssh_sol}_MOC_GAUS_4_{batch}_200"
                    keff_D5, fiss_rates_D5, n_gamma_rates_D5, FLUX_295groups_D5 = parse_DRAGON_SCHEME(name_case, CPO_name_level1, composition_option, evaluation, ssh_method, correlation_option, geometry_type, fission_isotopes , n_gamma_isotopes = ["U238", "Gd155", "Gd157"], bu=0, computational_scheme=computational_scheme)
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
                    

                    plot_pinwise_errors_BWR_assembly(np.array([delta_rel_therm,delta_rel_fast]), name_case, CPO_name_level1, calculation_opt = f"LVL1_{computational_scheme}_{geometry_type}_{ssh_method}_{ssh_sol}_{correlation_option}", fig_name=f"Total_fission_rates_diff")
    print(deltas)
    
    
    

if __name__ == "__main__":


    #analyze_main_flux_geometry_refinement()
    analyze_computational_schemes()
    