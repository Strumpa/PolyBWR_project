## Python3 PyGan script to post treat OECD NEA PHASE IIIB benchmark calculations performed with DRAGON5 and Serpent2.

# This script is used to post-process the results of the OECD NEA Phase IIIB benchmark calculations.
# It compares the results from DRAGON5 and Serpent2, specifically focusing : 
# At t=0, the following quantities are compared:
#  1) keff
#  2) flux spectrum
#  3) reaction rates
## Consider generalizing to BU evolution next.

import os
import glob
import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import lcm
import sys
import serpentTools as st
from Serpent2_parsing import parse_Serpent2_lattice_det, parse_Spectrum_detector, parse_eighth_lattice_S2, parse_material_det
from Dragon5_cpo_parsing import parse_Spectrum_D5, parse_D5_rates
from plotters import plot_spectrum_comparison, plot_eighth_lattice

def post_treat_CPO(COMPO_name, path_to_CPOs, draglib):
    """
    Args:
        COMPO_name (str): name of the multicompo object to open
        path_to_CPOs (str): path to read multicompo from
        drablib (str): choice of draglib used : specifies nuclear data evaluation and energy mesh used. options are "J311_295" and "endfb8r1_295"
    
    Purpose : post treat contents of multi compo generated from D5

    """
    path = os.getcwd()
    # read compo file
    os.chdir(path_to_CPOs)
    os.listdir()
    print(COMPO_name)
    pyCOMPO=lcm.new('LCM_INP',COMPO_name,impx=0)
    os.chdir(path)

    # Parse COMPO
    keff_D5, energy_mesh, FLUX_295groups = parse_Spectrum_D5(pyCOMPO)
    print(f"Energy bounds")
    print(energy_mesh[0])
    print(energy_mesh[-1])

    # Read Serpent results
    path_to_S2 = os.environ["SERPENT_RESULTS"]
    if draglib == "J311_295":
        name_case = "OECD_NEA_PHASE_IIIB_jeff311_inp"
    

    keff_S2, S2_spectrum = parse_Spectrum_detector(path_to_S2, name_case)

    delta_keff = (keff_D5 - keff_S2)*1e5
    print("keff_D5:", keff_D5)
    print("keff_S2:", keff_S2)
    print(f"delta_keff: {(keff_D5 - keff_S2)*1e5:.1f}")

    plot_spectrum_comparison(energy_mesh, FLUX_295groups, S2_spectrum, COMPO_name, results_dir)

    fiss_rates_D5, fiss_iso_dens = parse_D5_rates(pyCOMPO)
    keff_S2, fiss_rates_S2, n_gamma_rates_S2, mix_desc_S2 = parse_material_det(path_to_S2, name_case)

    POS_TO_ROD_IDX = {9:1, 8:1, 7:1, 6:5, 5:2, 4:2, 3:2, 2:3, 1:4}
    fissile_isotopes = ["U235", "U238"]
    for pos_idx in POS_TO_ROD_IDX.keys():
        mix_idx = POS_TO_ROD_IDX[pos_idx]
        print(f"For mix {mix_idx} at position {pos_idx}")
        for isotope in fissile_isotopes:
            print(f"{isotope} density S2= {mix_desc_S2[mix_idx][isotope]}")
            print(f"{isotope} density in CPO = {fiss_iso_dens[pos_idx][isotope]}")
            print(f"delta = {(fiss_iso_dens[pos_idx][isotope] - mix_desc_S2[mix_idx][isotope])*100/mix_desc_S2[mix_idx][isotope]}")

    print(f"fiss_rates_S2: {fiss_rates_S2}, with shape {fiss_rates_S2.shape}")
    print(f"fiss_rates_D5: {fiss_rates_D5}, with shape {fiss_rates_D5.shape}")

    print(f"keff_S2: {keff_S2}")
    print(f"Delta keff = {(keff_D5 - keff_S2)*1e5:.6f}")
    print(f"fiss_rates_S2: {fiss_rates_S2}")

    print(f"n_gamma_rates_S2: {n_gamma_rates_S2}")

    # normalize fission rates to a total of n_Cells
    thermal_fission_rates_S2 = fiss_rates_S2[0] / np.sum(fiss_rates_S2[0]) * fiss_rates_S2.shape[1]
    fast_fission_rates_S2 = fiss_rates_S2[1] / np.sum(fiss_rates_S2[1]) * fiss_rates_S2.shape[1]
    print(f"fiss_rates_D5: {fiss_rates_D5}, with shape {fiss_rates_D5.shape}")

    thermal_fiss_rates_D5 = fiss_rates_D5[1] / np.sum(fiss_rates_D5[1]) * fiss_rates_D5.shape[1]
    fast_fiss_rates_D5 = fiss_rates_D5[0] / np.sum(fiss_rates_D5[0]) * fiss_rates_D5.shape[1]

    print(f"thermal fission rates S2 : {thermal_fission_rates_S2}")
    print(f"thermal fission rates D5 : {thermal_fiss_rates_D5}")

    print(f"fast fission rates S2 : {fast_fission_rates_S2}")
    print(f"fast fission rates D5 : {fast_fiss_rates_D5}")


    # Relative difference in fission rates
    delta_rel_therm_fission_rates = (thermal_fiss_rates_D5 - thermal_fission_rates_S2) * 100 / thermal_fission_rates_S2
    delta_rel_fast_fission_rates = (fast_fiss_rates_D5 - fast_fission_rates_S2) * 100 / fast_fission_rates_S2

    print(f"delta_rel_therm_fission_rates: {delta_rel_therm_fission_rates}")
    print(f"delta_rel_fast_fission_rates: {delta_rel_fast_fission_rates}")

    print(f"sum of delta_rel_therm_fission_rates: {np.sum(delta_rel_therm_fission_rates)}")
    print(f"sum of delta_rel_fast_fission_rates: {np.sum(delta_rel_fast_fission_rates)}")

    delta_therm_fiss_rates = (thermal_fiss_rates_D5 - thermal_fission_rates_S2)
    delta_fast_fiss_rates = (fast_fiss_rates_D5 - fast_fission_rates_S2)
    
    print(f"delta_therm_fission_rates: {delta_therm_fiss_rates}")
    print(f"delta_fast_fission_rates: {delta_fast_fiss_rates}")

    print(f"sum of delta_therm_fission_rates: {np.sum(delta_therm_fiss_rates)}")
    print(f"sum of delta_fast_fission_rates: {np.sum(delta_fast_fiss_rates)}")

    plot_eighth_lattice(COMPO_name, [delta_rel_therm_fission_rates, delta_rel_fast_fission_rates], results_dir)

    return delta_keff


if __name__ == "__main__":
    # Create results directory if it doesn't exist
    results_dir = "OECD_NEA_IIIB_results"
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)

    path_to_CPOs = "DRAGON_RESULTS"
    case_to_treat = "IIIB" # "OECD_NEA_PHASE_IIIB"
    draglib = "J311_295" # "endfb8r1_295"
    # find all compos with format *case_to_treat*draglib*
    cwd_path = path = os.getcwd()
    os.chdir(path_to_CPOs)
    print(glob.glob(f"*{case_to_treat}*{draglib}*"))
    list_of_cpos_matching = glob.glob(f"*{case_to_treat}*{draglib}*")
    os.chdir(cwd_path)
    name_to_delta_keff = {}
    for cpo_name in list_of_cpos_matching:
        delta_keff = post_treat_CPO(cpo_name[1::], path_to_CPOs, draglib)
        name_to_delta_keff[cpo_name] = abs(delta_keff)

    min_delta_keff_case = min(name_to_delta_keff, key=name_to_delta_keff.get)
    print(f"Min delta keff (D5-S2) is {name_to_delta_keff[min_delta_keff_case]} for case {min_delta_keff_case}")

    print(name_to_delta_keff)




