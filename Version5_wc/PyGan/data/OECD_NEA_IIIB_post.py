## Python3 PyGan script to post treat OECD NEA PHASE IIIB benchmark calculations performed with DRAGON5 and Serpent2.

# This script is used to post-process the results of the OECD NEA Phase IIIB benchmark calculations.
# It compares the results from DRAGON5 and Serpent2, specifically focusing : 
# At t=0, the following quantities are compared:
#  1) keff
#  2) flux spectrum
#  3) reaction rates
## Consider generalizing to BU evolution next.

import os
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

# Create results directory if it doesn't exist
results_dir = "OECD_NEA_IIIB_results"
if not os.path.exists(results_dir):
    os.makedirs(results_dir)

# Add the path to the PyGan module

COMPO_name = "CPO_OECD_NEA_PHASE_IIIB_J311_295_WR_SPLITR_NOSECT_inrs1_eighth" # "CPO_OECD_NEA_PHASE_IIIB_SECT_WR" #"CPO_OECD_NEA_PHASE_IIIB_SECT_SPLIT_WR_indiv_inrs"#
path_to_CPO = "DRAGON_RESULTS"
path = os.getcwd()

# read compo file
os.chdir(path_to_CPO)
pyCOMPO=lcm.new('LCM_INP',COMPO_name,impx=0)
os.chdir(path)

# Parse COMPO
keff_D5, energy_mesh, FLUX_295groups = parse_Spectrum_D5(pyCOMPO)

# Read Serpent results
path_to_S2 = os.environ["SERPENT_RESULTS"]
name_case = "OECD_NEA_PHASE_IIIB_jeff311_inp"

keff_S2, S2_spectrum = parse_Spectrum_detector(path_to_S2, name_case)

print("keff_D5:", keff_D5)
print("keff_S2:", keff_S2)
print(f"delta_keff: {(keff_D5 - keff_S2)*1e5:.6f}")

plot_spectrum_comparison(energy_mesh, FLUX_295groups, S2_spectrum, results_dir)

# Parse Serpent2 lattice detector
fiss_rates_S2, mix_desc_S2 = parse_eighth_lattice_S2(path_to_S2, name_case)
# Parse DRAGON5 CPO for fission rates
fiss_rates_D5, fiss_iso_dens = parse_D5_rates(pyCOMPO)

POS_TO_ROD_IDX = {9:1, 8:1, 7:1, 6:5, 5:2, 4:2, 3:2, 2:3, 1:4}
fissile_isotopes = ["U234", "U235", "U236", "U238"]
for pos_idx in POS_TO_ROD_IDX.keys():
    mix_idx = POS_TO_ROD_IDX[pos_idx]
    print(f"For mix {mix_idx} at position {pos_idx}")
    for isotope in fissile_isotopes:
        print(f"{isotope} density S2= {mix_desc_S2[mix_idx][isotope]}")
        print(f"{isotope} density in CPO = {fiss_iso_dens[pos_idx][isotope]}")
        print(f"delta = {(fiss_iso_dens[pos_idx][isotope] - mix_desc_S2[mix_idx][isotope])*100/mix_desc_S2[mix_idx][isotope]}")

print(f"fiss_rates_S2: {fiss_rates_S2}, with shape {fiss_rates_S2.shape}")
print(f"fiss_rates_D5: {fiss_rates_D5}, with shape {fiss_rates_D5.shape}")

## Normalize the fission rates to a total of n_Cells
thermal_fiss_rates_S2 = fiss_rates_S2[0] / np.sum(fiss_rates_S2[0]) * fiss_rates_S2.shape[1]
thermal_fiss_rates_D5 = fiss_rates_D5[1] / np.sum(fiss_rates_D5[1]) * fiss_rates_D5.shape[1]

fast_fiss_rates_S2 = fiss_rates_S2[1] / np.sum(fiss_rates_S2[1]) * fiss_rates_S2.shape[1]
fast_fiss_rates_D5 = fiss_rates_D5[0] / np.sum(fiss_rates_D5[0]) * fiss_rates_D5.shape[1]

# Relative difference in fission rates
delta_rel_therm_fiss_rates = (thermal_fiss_rates_D5 - thermal_fiss_rates_S2) * 100 / thermal_fiss_rates_S2
delta_rel_fast_fiss_rates = (fast_fiss_rates_D5 - fast_fiss_rates_S2) * 100 / fast_fiss_rates_S2

print(f"delta_rel_therm_fiss_rates: {delta_rel_therm_fiss_rates}")
print(f"delta_rel_fast_fiss_rates: {delta_rel_fast_fiss_rates}")

plot_eighth_lattice(COMPO_name, [delta_rel_therm_fiss_rates, delta_rel_fast_fiss_rates], results_dir)


keff_S2, fission_rates_S2, n_gamma_rates_S2 = parse_material_det(path_to_S2, name_case)
print(f"keff_S2: {keff_S2}")
print(f"fission_rates: {fission_rates_S2}")
print(f"n_gamma_rates: {n_gamma_rates_S2}")

# normalize fission rates to a total of n_Cells
thermal_fission_rates_S2 = fission_rates_S2[0] / np.sum(fission_rates_S2[0]) * fission_rates_S2.shape[1]
fast_fission_rates_S2 = fission_rates_S2[1] / np.sum(fission_rates_S2[1]) * fission_rates_S2.shape[1]


# Relative difference in fission rates
delta_rel_therm_fission_rates = (thermal_fiss_rates_D5 - thermal_fission_rates_S2) * 100 / thermal_fission_rates_S2
delta_rel_fast_fission_rates = (fast_fiss_rates_D5 - fast_fission_rates_S2) * 100 / fast_fission_rates_S2

print(f"delta_rel_therm_fission_rates: {delta_rel_therm_fission_rates}")
print(f"delta_rel_fast_fission_rates: {delta_rel_fast_fission_rates}")

plot_eighth_lattice(COMPO_name, [delta_rel_therm_fission_rates, delta_rel_fast_fission_rates], results_dir)

