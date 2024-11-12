# Python3 script for post treatment of Gd157 XS.
# Author : R. Guasch
# Date : 2024-11-12
# Purpose : Post treatment of Gd157 XS and abs rates, Serpent2 PyNjoy2016 vs Serpent2 oldlib comparison

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import serpentTools as st
from energyMeshHandler import energyMeshHandler


def plot_histograms_XS(values_lists, iso, colors, grmin, grmax, legends=None, save_dir=None, save_name=None):
    # Check that the number of values is 295
    for values in values_lists:
        if len(values) != grmax - grmin + 1:
            raise ValueError(f"All input rates lists must contain exactly {grmax - grmin + 1} values.")
    
    # Create bins from 1 to 295
    energy_groups = list(range(grmin, grmax + 1))
    #colors = ['skyblue', 'red']
    # Plot the histogram
    
    plt.figure(figsize=(10, 6))
    for i in range(len(values_lists)):
        plt.bar(energy_groups, values_lists[i], edgecolor='black', alpha=0.5, color=colors[i], label=legends[i])
    plt.ylabel('Relative difference on XS (%)')
    plt.xticks(range(grmin, grmax+1, 20))  # Set x-ticks with an interval of 20 for readability
    plt.xlabel('Energy group')
    plt.title(f'{iso} absorption and n,gamma XS difference between PyNjoy2016 and oldlib')
    plt.grid()
    plt.legend()
    if save_dir:
        plt.savefig(f'{save_dir}/{save_name}.png')
    plt.show()
    plt.close()

def plot_histograms_rates(values_lists, iso, colors, grmin, grmax, legends=None, save_dir=None, save_name=None):
    # Check that the number of values is 295
    for values in values_lists:
        if len(values) != grmax - grmin + 1:
            raise ValueError("All input rates lists must contain exactly 295 values.")
    
    # Create bins from 1 to 295
    energy_groups = list(range(grmin, grmax + 1))
    
    plt.figure(figsize=(10, 6))
    for i in range(len(values_lists)):
        plt.bar(energy_groups, values_lists[i], edgecolor='black', alpha=0.5, color=colors[i], label=legends[i])
    if "error" in legends[0]:
        plt.ylabel('Relative error (%)')
    elif "ratio" in save_name:
        plt.ylabel('Ratio of absorption rates')
    else:
        plt.ylabel('Absorption rate')
    plt.xticks(range(grmin, grmax+1, 20))  # Set x-ticks with an interval of 20 for readability
    plt.xlabel('Energy group')
    if "error" in legends[0]:
        plt.title(f'Relative errors on absorption rates for {iso}')
    else:
        plt.title(f'Absorption rates for {iso}')
    plt.grid()
    plt.legend()
    #plt.tight_layout()
    if save_dir:
        plt.savefig(f'{save_dir}/{save_name}.png')
    plt.show()
    plt.close()

def parse_Serpent_microdepletion(path_to_serpent_results, library, bu_step):
    # Read mdep file
    mdep = st.read(f"{path_to_serpent_results}/HOM_UOX_Gd157_XS_{library}_mc_mdx{bu_step}.m")

    mdep_scores = mdep.xsVal

    for key in mdep_scores['1'].keys():
        Gd157_XS_102 = mdep_scores['1'][key]
    for key in mdep_scores['2'].keys():
        Gd157_XS_101 = mdep_scores['2'][key]
    for key in mdep_scores['3'].keys():
        U8_XS_102 = mdep_scores['3'][key]
    for key in mdep_scores['4'].keys():
        U8_XS_101 = mdep_scores['4'][key]
    return  Gd157_XS_102, Gd157_XS_101, U8_XS_102, U8_XS_101

def parse_Serpent_detector(path_to_serpent_results, library, bu_step):
    # Read detector file
    det = st.read(f"{path_to_serpent_results}/HOM_UOX_Gd157_XS_{library}_mc_det{bu_step}.m")
    # Get detector names
    Gd_det = det.detectors["Gd_det"].tallies
    n_groups = Gd_det.shape[0]
    n_reactions = Gd_det.shape[1]
    Gd_nGamma_rates = []
    U8_nGamma_rates = []
    U8_abs_rates = []
    for group in range(n_groups):
        Gd_nGamma_rates.append(Gd_det[group, 0])
        U8_nGamma_rates.append(Gd_det[group, 1])
        U8_abs_rates.append(Gd_det[group, 2])
    Gd_nGamma_rates_reversed = np.array(Gd_nGamma_rates[::-1])
    U8_nGamma_rates_reversed = np.array(U8_nGamma_rates[::-1])
    U8_abs_rates_reversed = np.array(U8_abs_rates[::-1])
    return Gd_nGamma_rates_reversed, U8_nGamma_rates_reversed, U8_abs_rates_reversed


if __name__ == "__main__":
    # Parse Serpent microdepletion results
    input_case = "HOM_UOX_Gd157"
    path_S2 = f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/HOM_CELL_study/{input_case}/XS_study"
    #library = "oldlib"

    local_path = os.getcwd()
    libraries = ["PyNjoy2016", "oldlib"]
    bu_step = 0

    
    path_exists = os.path.exists(f"{local_path}/{input_case}_S2/XS_study")
    if not path_exists:
        os.makedirs(f"{local_path}/{input_case}_S2/XS_study")
    path_save = f"{local_path}/{input_case}_S2/XS_study"

    LETHARGY_MESH = energyMeshHandler('SHEM295.txt', E0=1.0E+07, energyUnits='eV').lethargyMesh

    # Parse Serpent microdepletion results for both libraries at bu_step = 0
    for library in libraries:
        Gd157_nGamma_rates, U8_nGamma_rates, U8_abs_rates = parse_Serpent_detector(path_S2, library, bu_step)
        Gd157_XS_ngamma, Gd157_XS_abstot, U8_XS_ngamma, U8_XS_abstot = parse_Serpent_microdepletion(path_S2, library, bu_step)
        if library == "PyNjoy2016":
            Gd157_rates_ngamma_PyNjoy2016 = Gd157_nGamma_rates
            U8_rates_ngamma_PyNjoy2016 = U8_nGamma_rates
            # Exclude first entry as it is N_iso
            N_Gd157_PyNjoy2016 = Gd157_XS_abstot[0]
            Gd157_XS_ngamma_PyNjoy2016 = Gd157_XS_ngamma[1:]
            Gd157_XS_abstot_PyNjoy2016 = Gd157_XS_abstot[1:]
            N_U8_PyNjoy2016 = U8_XS_abstot[0]
            U8_XS_ngamma_PyNjoy2016 = U8_XS_ngamma[1:]
            U8_XS_abstot_PyNjoy2016 = U8_XS_abstot[1:]
        elif library == "oldlib":
            Gd157_rates_ngamma_oldlib = Gd157_nGamma_rates
            U8_rates_ngamma_oldlib = U8_nGamma_rates
            N_Gd157_oldlib = Gd157_XS_abstot[0]
            Gd157_XS_ngamma_oldlib = Gd157_XS_ngamma[1:]
            Gd157_XS_abstot_oldlib = Gd157_XS_abstot[1:]
            N_U8_oldlib = U8_XS_abstot[0]
            U8_XS_ngamma_oldlib = U8_XS_ngamma[1:]
            U8_XS_abstot_oldlib = U8_XS_abstot[1:]
    
    # Check that N_Gd157 is the same for both libraries
    if N_Gd157_PyNjoy2016 != N_Gd157_oldlib:
        raise ValueError("N_Gd157 is not the same for both libraries.")
    if N_U8_PyNjoy2016 != N_U8_oldlib:
        raise ValueError("N_U8 is not the same for both libraries.")

    

    # plot XS comparing Serpent2 PyNjoy2016 vs oldlib
    plt.figure()
    plt.step(LETHARGY_MESH, Gd157_XS_ngamma_PyNjoy2016, label='PyNjoy2016')
    plt.step(LETHARGY_MESH, Gd157_XS_ngamma_oldlib, label='oldlib')
    plt.title('Gd157 nGamma XS')
    plt.legend()
    plt.yscale('log')
    plt.show()
    plt.savefig(f'{path_save}/SERPENT2_oldlib_PyNjoy2016_comparison_Gd157nGamma.png')
    plt.close()

    plt.figure()
    plt.step(LETHARGY_MESH, U8_XS_ngamma_PyNjoy2016, label='PyNjoy2016')
    plt.step(LETHARGY_MESH, U8_XS_ngamma_oldlib, label='oldlib')
    plt.title('U8 nGamma XS')
    plt.legend()
    plt.yscale('log')
    plt.show()
    plt.savefig(f'{path_save}/SERPENT2_oldlib_PyNjoy2016_comparison_U8nGamma.png')
    plt.close()


    plt.figure()
    plt.step(LETHARGY_MESH, Gd157_XS_abstot_PyNjoy2016, label='PyNjoy2016')
    plt.step(LETHARGY_MESH, Gd157_XS_abstot_oldlib, label='oldlib')
    plt.title('Gd157 total absorption XS')
    plt.legend()
    plt.yscale('log')
    plt.show()
    plt.savefig(f'{path_save}/SERPENT2_oldlib_PyNjoy2016_comparison_Gd157abstot.png')
    plt.close()

    plt.figure()
    plt.step(LETHARGY_MESH, U8_XS_abstot_PyNjoy2016, label='PyNjoy2016')
    plt.step(LETHARGY_MESH, U8_XS_abstot_oldlib, label='oldlib')
    plt.title('U238 total absorption XS')
    plt.legend()
    plt.yscale('log')
    plt.show()
    plt.savefig(f'{path_save}/SERPENT2_oldlib_PyNjoy2016_comparison_U8abstot.png')
    plt.close()


    # Compute relative errors on rates

    Gd157_error_ngamma_rates = (Gd157_rates_ngamma_PyNjoy2016 - Gd157_rates_ngamma_oldlib) * 100 / Gd157_rates_ngamma_oldlib
    U8_error_ngamma_rates = (U8_rates_ngamma_PyNjoy2016 - U8_rates_ngamma_oldlib) * 100 / U8_rates_ngamma_oldlib

    # Plot histograms for rates
    plot_histograms_rates([Gd157_error_ngamma_rates, ], "Gd157", ['skyblue'], 1, 295, legends=['Gd157 nGamma rates error'], save_dir=path_save, save_name='SERPENT2_oldlib_PyNjoy2016_comparison_diff_Gd157_nGamma_rates_error')
    plot_histograms_rates([U8_error_ngamma_rates], "U238", ['red'], 1, 295, legends=['U238 nGamma rates error'], save_dir=path_save, save_name='SERPENT2_oldlib_PyNjoy2016_comparison_diff_U8_nGamma_rates_error')


    # Compute relative errors on XS
    Gd157_error_abstot = (Gd157_XS_abstot_PyNjoy2016 - Gd157_XS_abstot_oldlib) * 100 / Gd157_XS_abstot_oldlib
    U8_error_abstot = (U8_XS_abstot_PyNjoy2016 - U8_XS_abstot_oldlib) * 100 / U8_XS_abstot_oldlib

    Gd157_error_ngamma = (Gd157_XS_ngamma_PyNjoy2016 - Gd157_XS_ngamma_oldlib) * 100 / Gd157_XS_ngamma_oldlib
    U8_error_ngamma = (U8_XS_ngamma_PyNjoy2016 - U8_XS_ngamma_oldlib) * 100 / U8_XS_ngamma_oldlib

    # Plot histograms for XS
    plot_histograms_XS([Gd157_error_abstot], "Gd157", ['skyblue'], 1, 295, legends=['Gd157 total absorption'], save_dir=path_save, save_name='SERPENT2_oldlib_PyNjoy2016_comparison_diff_Gd157_abstot_XS')
    plot_histograms_XS([Gd157_error_ngamma], "Gd157", ['red'], 1, 295, legends=['Gd157 n,gamma'], save_dir=path_save, save_name='SERPENT2_oldlib_PyNjoy2016_comparison_diff_Gd157_ngamma_XS')
    plot_histograms_XS([U8_error_abstot], "U238", ['skyblue'], 1, 295, legends=['U238 total absorption'], save_dir=path_save, save_name='SERPENT2_oldlib_PyNjoy2016_comparison_diff_U238_abstot_XS')
    plot_histograms_XS([U8_error_ngamma], "U238", ['red'], 1, 295, legends=['U238 n,gamma'], save_dir=path_save, save_name='SERPENT2_oldlib_PyNjoy2016_comparison_diff_U238_ngamma_XS')