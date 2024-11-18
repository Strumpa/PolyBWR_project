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
    if "abstot" in save_name:
        plt.title(f'{iso} : relative difference in $\sigma_a$ between PyNjoy2016 and oldlib')
    elif "ngamma" in save_name:
        plt.title(f'{iso} : relative difference in $\sigma (n,\gamma)$ between PyNjoy2016 and oldlib')
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
        if "ngamma" in save_name:
            plt.title(f'Relative difference on $\\tau (n,\gamma)$ rates for {iso}')
        elif "abs" in save_name:
            plt.title(f'Relative difference on $\\tau_a$ rates for {iso}')
        #plt.title(f'Relative errors on absorption rates for {iso}')
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
    print(Gd_det.shape)
    n_groups = Gd_det.shape[0]
    n_reactions = Gd_det.shape[1]
    print(f"Number of energy groups: {n_groups}")
    print(f"Number of reactions: {n_reactions}")
    Gd_nGamma_rates = []
    Gd_abs_rates = []
    U8_nGamma_rates = []
    U8_abs_rates = []
    for group in range(n_groups):
        Gd_nGamma_rates.append(Gd_det[group, 0])
        Gd_abs_rates.append(Gd_det[group, 1])
        U8_nGamma_rates.append(Gd_det[group, 2])
        U8_abs_rates.append(Gd_det[group, 3])
    Gd_nGamma_rates_reversed = np.array(Gd_nGamma_rates[::-1])
    Gd_abs_rates_reversed = np.array(Gd_abs_rates[::-1])
    U8_nGamma_rates_reversed = np.array(U8_nGamma_rates[::-1])
    U8_abs_rates_reversed = np.array(U8_abs_rates[::-1])
    return Gd_nGamma_rates_reversed, Gd_abs_rates_reversed, U8_nGamma_rates_reversed, U8_abs_rates_reversed

def plot_piecewise_constant(boundaries, values):
    """
    Plot a piecewise constant distribution.
    
    Parameters:
    boundaries (list or array): List of n+1 boundary values defining the intervals.
    values (list or array): List of n values for each interval.
    """
    if len(boundaries) != len(values) + 1:
        raise ValueError("The number of boundaries must be one more than the number of values.")
    
    # Prepare data for step plot
    x = []
    y = []
    
    for i in range(len(values)):
        x.extend([boundaries[i], boundaries[i+1]])
        y.extend([values[i], values[i]])
    
    print(f"x = {x}")
    print(f"y = {y}")
    print(len(x))
    print(len(y))
    # Plot the piecewise constant function
    plt.step(x, y, where='post', color='b', linewidth=2)
    plt.xlabel('x')
    plt.ylabel('Value')
    plt.yscale('log')
    plt.title('Piecewise Constant Distribution')
    plt.grid(True)
    plt.savefig(f'piecewise_constant.png')
    plt.show()

if __name__ == "__main__":
    # Parse Serpent microdepletion results
    input_case = "HOM_UOX_Gd157"
    path_S2 = f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/HOM_CELL_study/{input_case}/XS_study"
    #library = "oldlib"
    delta_NGd157 = []
    delta_NU238 = []
    delta_Gd157_ngamma_rates = {"avg": [], "rms": [], "max": []}
    delta_U238_ngamma_rates = {"avg": [], "rms": [], "max": []}
    delta_Gd157_ngamma_XS = {"avg": [], "rms": [], "max": []}
    delta_U238_ngamma_XS = {"avg": [], "rms": [], "max": []}

    local_path = os.getcwd()
    libraries = ["PyNjoy2016", "oldlib"]
    bu_steps_to_treat = range(0, 70)
    bu_steps_to_print = [0, 69]

    BU_list = [0.0, 0.015, 0.03, 0.05, 0.075, 0.1125, 0.15, 0.2, 0.25, 0.375, 0.5, 
            0.625, 0.75, 0.875, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5,
            2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5, 4.75, 5.0,
            5.25, 5.5, 5.75, 6.0, 6.25, 6.5, 6.75, 7.0, 7.25, 7.5,
            7.75, 8.0, 8.25, 8.5, 8.75, 9.0, 9.25, 9.5, 9.75, 10.0,
            10.25, 10.5, 10.75, 11.0, 11.25, 11.5, 11.75, 12.0, 12.5, 13.0, 
            13.5, 14.0, 14.5, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0]
    BU_list = np.array(BU_list)*1000 # Convert to MWd/t
    BU_axis = BU_list[bu_steps_to_treat]
    print(f"BU steps to treat: {BU_axis}")
    print(len(BU_axis))
    print(len(bu_steps_to_treat))

    
    path_exists = os.path.exists(f"{local_path}/{input_case}_S2/XS_study")
    if not path_exists:
        os.makedirs(f"{local_path}/{input_case}_S2/XS_study")
    path_save = f"{local_path}/{input_case}_S2/XS_study"

    LETHARGY_MESH = energyMeshHandler('SHEM295.txt', E0=1.0E+07, energyUnits='eV').lethargyMesh
    print(len(LETHARGY_MESH))
    ENERGY_mesh = energyMeshHandler('SHEM295.txt', E0=1.0E+07, energyUnits='eV').energyMesh
    for bu_step in bu_steps_to_treat:
        # Parse Serpent microdepletion results for both libraries at bu_step = 0
        for library in libraries:
            Gd157_nGamma_rates, Gd157_abs_rates, U8_nGamma_rates, U8_abs_rates = parse_Serpent_detector(path_S2, library, bu_step)
            Gd157_XS_ngamma, Gd157_XS_abstot, U8_XS_ngamma, U8_XS_abstot = parse_Serpent_microdepletion(path_S2, library, bu_step)
            if library == "PyNjoy2016":
                Gd157_rates_ngamma_PyNjoy2016 = Gd157_nGamma_rates
                Gd157_abs_rates_PyNjoy2016 = Gd157_abs_rates
                U8_rates_ngamma_PyNjoy2016 = U8_nGamma_rates
                U8_rates_abs_PyNjoy2016 = U8_abs_rates
                # Exclude first entry as it is N_iso
                N_Gd157_PyNjoy2016 = Gd157_XS_abstot[0]
                Gd157_XS_ngamma_PyNjoy2016 = Gd157_XS_ngamma[1:]
                Gd157_XS_abstot_PyNjoy2016 = Gd157_XS_abstot[1:]
                N_U8_PyNjoy2016 = U8_XS_abstot[0]
                U8_XS_ngamma_PyNjoy2016 = U8_XS_ngamma[1:]
                U8_XS_abstot_PyNjoy2016 = U8_XS_abstot[1:]
            elif library == "oldlib":
                Gd157_rates_ngamma_oldlib = Gd157_nGamma_rates
                Gd157_rates_abs_oldlib = Gd157_abs_rates
                U8_rates_ngamma_oldlib = U8_nGamma_rates
                U8_rates_abs_oldlib = U8_abs_rates
                N_Gd157_oldlib = Gd157_XS_abstot[0]
                Gd157_XS_ngamma_oldlib = Gd157_XS_ngamma[1:]
                Gd157_XS_abstot_oldlib = Gd157_XS_abstot[1:]
                N_U8_oldlib = U8_XS_abstot[0]
                U8_XS_ngamma_oldlib = U8_XS_ngamma[1:]
                U8_XS_abstot_oldlib = U8_XS_abstot[1:]
        
        if bu_step == 0:
            # Check that N_Gd157 is the same for both libraries
            if N_Gd157_PyNjoy2016 != N_Gd157_oldlib:
                raise ValueError("N_Gd157 is not the same for both libraries.")
            if N_U8_PyNjoy2016 != N_U8_oldlib:
                raise ValueError("N_U8 is not the same for both libraries.")
            N_Gd157_error = 0.0
            N_U8_error = 0.0
        else: 
            # Check relative differences on N_Gd157 and N_U8
            N_Gd157_error = (N_Gd157_PyNjoy2016 - N_Gd157_oldlib) * 100 / N_Gd157_oldlib
            N_U8_error = (N_U8_PyNjoy2016 - N_U8_oldlib) * 100 / N_U8_oldlib
            print(f"Errors on evolution of Gd157 and U238 between PyNjoy2016 and oldlib at bu_step = {bu_step}")
            print(f"Relative difference on N_Gd157: {N_Gd157_error}%")
            print(f"Relative difference on N_U8: {N_U8_error}%")
        delta_NGd157.append(N_Gd157_error)
        delta_NU238.append(N_U8_error)

        #plot_piecewise_constant(LETHARGY_MESH, Gd157_XS_ngamma_PyNjoy2016)
        U_mesh = []
        Gd157_abs_rates_PyNjoy2016_to_plot = []
        Gd157_abs_rates_oldlib_to_plot = []
        Gd157_ngamma_rates_PyNjoy2016_to_plot = []
        Gd157_ngamma_rates_oldlib_to_plot = []

        U8_abs_rates_PyNjoy2016_to_plot = []
        U8_abs_rates_oldlib_to_plot = []
        U8_ngamma_rates_PyNjoy2016_to_plot = []
        U8_ngamma_rates_oldlib_to_plot = []
        
        Gd157_XS_ngamma_PyNjoy2016_to_plot = []
        Gd157_XS_ngamma_oldlib_to_plot = []

        Gd157_XS_abstot_PyNjoy2016_to_plot = []
        Gd157_XS_abstot_oldlib_to_plot = []

        U8_XS_ngamma_PyNjoy2016_to_plot = []
        U8_XS_ngamma_oldlib_to_plot = []

        U8_XS_abstot_PyNjoy2016_to_plot = []
        U8_XS_abstot_oldlib_to_plot = []

        for i in range(len(Gd157_abs_rates_PyNjoy2016)):
            U_mesh.extend([LETHARGY_MESH[i], LETHARGY_MESH[i+1]])
            Gd157_abs_rates_PyNjoy2016_to_plot.extend([Gd157_abs_rates_PyNjoy2016[i], Gd157_abs_rates_PyNjoy2016[i]])
            Gd157_abs_rates_oldlib_to_plot.extend([Gd157_rates_abs_oldlib[i], Gd157_rates_abs_oldlib[i]])
            U8_abs_rates_PyNjoy2016_to_plot.extend([U8_rates_abs_PyNjoy2016[i], U8_rates_abs_PyNjoy2016[i]])
            U8_abs_rates_oldlib_to_plot.extend([U8_rates_abs_oldlib[i], U8_rates_abs_oldlib[i]])
            Gd157_ngamma_rates_PyNjoy2016_to_plot.extend([Gd157_rates_ngamma_PyNjoy2016[i], Gd157_rates_ngamma_PyNjoy2016[i]])
            Gd157_ngamma_rates_oldlib_to_plot.extend([Gd157_rates_ngamma_oldlib[i], Gd157_rates_ngamma_oldlib[i]])
            U8_ngamma_rates_PyNjoy2016_to_plot.extend([U8_rates_ngamma_PyNjoy2016[i], U8_rates_ngamma_PyNjoy2016[i]])
            U8_ngamma_rates_oldlib_to_plot.extend([U8_rates_ngamma_oldlib[i], U8_rates_ngamma_oldlib[i]])
            
            Gd157_XS_ngamma_PyNjoy2016_to_plot.extend([Gd157_XS_ngamma_PyNjoy2016[i], Gd157_XS_ngamma_PyNjoy2016[i]])
            Gd157_XS_ngamma_oldlib_to_plot.extend([Gd157_XS_ngamma_oldlib[i], Gd157_XS_ngamma_oldlib[i]])
            U8_XS_ngamma_PyNjoy2016_to_plot.extend([U8_XS_ngamma_PyNjoy2016[i], U8_XS_ngamma_PyNjoy2016[i]])
            U8_XS_ngamma_oldlib_to_plot.extend([U8_XS_ngamma_oldlib[i], U8_XS_ngamma_oldlib[i]])
            Gd157_XS_abstot_PyNjoy2016_to_plot.extend([Gd157_XS_abstot_PyNjoy2016[i], Gd157_XS_abstot_PyNjoy2016[i]])
            Gd157_XS_abstot_oldlib_to_plot.extend([Gd157_XS_abstot_oldlib[i], Gd157_XS_abstot_oldlib[i]])
            U8_XS_abstot_PyNjoy2016_to_plot.extend([U8_XS_abstot_PyNjoy2016[i], U8_XS_abstot_PyNjoy2016[i]])
            U8_XS_abstot_oldlib_to_plot.extend([U8_XS_abstot_oldlib[i], U8_XS_abstot_oldlib[i]])

            

        if bu_step in bu_steps_to_print:
            # plot XS comparing Serpent2 PyNjoy2016 vs oldlib
            plt.figure()
            plt.step(U_mesh, Gd157_XS_ngamma_PyNjoy2016_to_plot, where='post', label='PyNjoy2016')
            plt.step(U_mesh, Gd157_XS_ngamma_oldlib_to_plot, where='post', label='oldlib')
            plt.title('Gd157 $(n,\gamma)$ XS')
            plt.legend()
            plt.yscale('log')
            plt.show()
            plt.savefig(f'{path_save}/SERPENT2_oldlib_PyNjoy2016_comparison_Gd157nGamma_bu{bu_step}.png')
            plt.close()

            plt.figure()
            plt.step(U_mesh, U8_XS_ngamma_PyNjoy2016_to_plot, where='post', label='PyNjoy2016')
            plt.step(U_mesh, U8_XS_ngamma_oldlib_to_plot, where='post', label='oldlib')
            plt.title('U238 $(n,\gamma)$ XS')
            plt.legend()
            plt.yscale('log')
            plt.show()
            plt.savefig(f'{path_save}/SERPENT2_oldlib_PyNjoy2016_comparison_U8nGamma_bu{bu_step}.png')
            plt.close()


            plt.figure()
            plt.step(U_mesh, Gd157_XS_abstot_PyNjoy2016_to_plot, where='post', label='PyNjoy2016')
            plt.step(U_mesh, Gd157_XS_abstot_oldlib_to_plot, where='post', label='oldlib')
            plt.title('Gd157 total absorption XS')
            plt.legend()
            plt.yscale('log')
            plt.show()
            plt.savefig(f'{path_save}/SERPENT2_oldlib_PyNjoy2016_comparison_Gd157abstot_bu{bu_step}.png')
            plt.close()

            plt.figure()
            plt.step(U_mesh, U8_XS_abstot_PyNjoy2016_to_plot, where='post', label='PyNjoy2016')
            plt.step(U_mesh, U8_XS_abstot_oldlib_to_plot, where='post', label='oldlib')
            plt.title('U238 total absorption XS')
            plt.legend()
            plt.yscale('log')
            plt.show()
            plt.savefig(f'{path_save}/SERPENT2_oldlib_PyNjoy2016_comparison_U8abstot_bu{bu_step}.png')
            plt.close()

            # plot rates comparing Serpent2 PyNjoy2016 vs oldlib
            plt.figure()
            plt.step(U_mesh, Gd157_ngamma_rates_PyNjoy2016_to_plot, where='post', label='PyNjoy2016')
            plt.step(U_mesh, Gd157_ngamma_rates_oldlib_to_plot, where='post', label='oldlib')
            plt.title('Gd157 $(n,\gamma)$ rates')
            plt.legend()
            plt.yscale('log')
            plt.show()
            plt.savefig(f'{path_save}/SERPENT2_oldlib_PyNjoy2016_comparison_Gd157nGamma_rates_bu{bu_step}.png')
            plt.close()

            plt.figure()
            plt.step(U_mesh, U8_ngamma_rates_PyNjoy2016_to_plot, where='post', label='PyNjoy2016')
            plt.step(U_mesh, U8_ngamma_rates_oldlib_to_plot, where='post', label='oldlib')
            plt.title('U238 $(n,\gamma)$ rates')
            plt.legend()
            plt.yscale('log')
            plt.show()
            plt.savefig(f'{path_save}/SERPENT2_oldlib_PyNjoy2016_comparison_U8nGamma_rates_bu{bu_step}.png')
            plt.close()

            plt.figure()
            plt.step(U_mesh, U8_abs_rates_PyNjoy2016_to_plot, where='post', label='PyNjoy2016')
            plt.step(U_mesh, U8_abs_rates_oldlib_to_plot, where='post', label='oldlib')
            plt.title('U238 total absorption rates')
            plt.legend()
            plt.yscale('log')
            plt.show()
            plt.savefig(f'{path_save}/SERPENT2_oldlib_PyNjoy2016_comparison_U8abstot_rates_bu{bu_step}.png')
            plt.close()



        # Compute relative errors on rates

        Gd157_error_ngamma_rates = (Gd157_rates_ngamma_PyNjoy2016 - Gd157_rates_ngamma_oldlib) * 100 / Gd157_rates_ngamma_oldlib
        U8_error_ngamma_rates = (U8_rates_ngamma_PyNjoy2016 - U8_rates_ngamma_oldlib) * 100 / U8_rates_ngamma_oldlib
        U8_error_abs_rates = (U8_rates_abs_PyNjoy2016 - U8_rates_abs_oldlib) * 100 / U8_rates_abs_oldlib

        delta_Gd157_ngamma_rates["avg"].append(np.mean(Gd157_error_ngamma_rates))
        delta_Gd157_ngamma_rates["rms"].append(np.sqrt(np.mean(Gd157_error_ngamma_rates**2)))
        delta_Gd157_ngamma_rates["max"].append(np.max(Gd157_error_ngamma_rates))
        delta_U238_ngamma_rates["avg"].append(np.mean(U8_error_ngamma_rates))
        delta_U238_ngamma_rates["rms"].append(np.sqrt(np.mean(U8_error_ngamma_rates**2)))
        delta_U238_ngamma_rates["max"].append(np.max(U8_error_ngamma_rates))

        if bu_step in bu_steps_to_print:
            # Plot histograms for rates
            plot_histograms_rates([Gd157_error_ngamma_rates, ], "Gd157", ['skyblue'], 1, 295, legends=['Gd157 $\\tau (n,\gamma)$ error'], save_dir=path_save, save_name=f'SERPENT2_oldlib_PyNjoy2016_comparison_diff_Gd157_ngamma_rates_error_bu{bu_step}')
            plot_histograms_rates([U8_error_ngamma_rates], "U238", ['red'], 1, 295, legends=['U238  $\\tau (n,\gamma)$ error'], save_dir=path_save, save_name=f'SERPENT2_oldlib_PyNjoy2016_comparison_diff_U8_ngamma_rates_error_bu{bu_step}')
            plot_histograms_rates([U8_error_abs_rates], "U238", ['green'], 1, 295, legends=['U238  $\\tau_a$ error'], save_dir=path_save, save_name=f'SERPENT2_oldlib_PyNjoy2016_comparison_diff_U8_abstot_rates_error_bu{bu_step}')

        # Compute relative errors on XS
        Gd157_error_abstot = (Gd157_XS_abstot_PyNjoy2016 - Gd157_XS_abstot_oldlib) * 100 / Gd157_XS_abstot_oldlib
        U8_error_abstot = (U8_XS_abstot_PyNjoy2016 - U8_XS_abstot_oldlib) * 100 / U8_XS_abstot_oldlib

        Gd157_error_ngamma = (Gd157_XS_ngamma_PyNjoy2016 - Gd157_XS_ngamma_oldlib) * 100 / Gd157_XS_ngamma_oldlib
        U8_error_ngamma = (U8_XS_ngamma_PyNjoy2016 - U8_XS_ngamma_oldlib) * 100 / U8_XS_ngamma_oldlib

        delta_Gd157_ngamma_XS["avg"].append(np.mean(Gd157_error_ngamma))
        delta_Gd157_ngamma_XS["rms"].append(np.sqrt(np.mean(Gd157_error_ngamma**2)))
        delta_Gd157_ngamma_XS["max"].append(np.max(Gd157_error_ngamma))
        delta_U238_ngamma_XS["avg"].append(np.mean(U8_error_ngamma))
        delta_U238_ngamma_XS["rms"].append(np.sqrt(np.mean(U8_error_ngamma**2)))
        delta_U238_ngamma_XS["max"].append(np.max(U8_error_ngamma))


        if bu_step in bu_steps_to_print:
            # Plot histograms for XS
            plot_histograms_XS([Gd157_error_abstot], "Gd157", ['skyblue'], 1, 295, legends=['Gd157 total absorption'], save_dir=path_save, save_name=f'SERPENT2_oldlib_PyNjoy2016_comparison_diff_Gd157_abstot_XS_bu{bu_step}')
            plot_histograms_XS([Gd157_error_ngamma], "Gd157", ['red'], 1, 295, legends=['Gd157 $(n,\gamma)$'], save_dir=path_save, save_name=f'SERPENT2_oldlib_PyNjoy2016_comparison_diff_Gd157_ngamma_XS_bu{bu_step}')
            plot_histograms_XS([U8_error_abstot], "U238", ['skyblue'], 1, 295, legends=['U238 total absorption'], save_dir=path_save, save_name=f'SERPENT2_oldlib_PyNjoy2016_comparison_diff_U238_abstot_XS_bu{bu_step}')
            plot_histograms_XS([U8_error_ngamma], "U238", ['red'], 1, 295, legends=['U238 $(n,\gamma)$'], save_dir=path_save, save_name=f'SERPENT2_oldlib_PyNjoy2016_comparison_diff_U238_ngamma_XS_bu{bu_step}')

    # Plot evolution of relative differences on N_Gd157 and N_U8
    plt.figure()
    plt.plot(BU_axis, delta_NGd157, label='Gd157')
    plt.plot(BU_axis, delta_NU238, label='U238')
    plt.xlabel('Burnup (MWd/t)')
    plt.ylabel('Relative difference (%)')
    plt.title('Evolution of relative differences on N_Gd157 and N_U8 between PyNjoy2016 and oldlib')
    plt.legend()
    plt.grid()
    plt.show()
    plt.savefig(f'{path_save}/SERPENT2_oldlib_PyNjoy2016_comparison_diff_NGd157_NU8.png')
    plt.close()

    # Plot evolution of relative differences on rates
    plt.figure()
    plt.plot(BU_axis, delta_Gd157_ngamma_rates["avg"], label='AVG $\\Delta$ Gd157')
    plt.plot(BU_axis, delta_Gd157_ngamma_rates["rms"], label='RMS $\\Delta$ Gd157')
    plt.plot(BU_axis, delta_Gd157_ngamma_rates["max"], label='MAX $\\Delta$ Gd157')

    plt.plot(BU_axis, delta_U238_ngamma_rates["avg"], label='AVG $\\Delta$ U238')
    plt.plot(BU_axis, delta_U238_ngamma_rates["rms"], label='RMS $\\Delta$ U238')
    plt.plot(BU_axis, delta_U238_ngamma_rates["max"], label='MAX $\\Delta$ U238')
    plt.xlabel('Burnup (MWd/t)')
    plt.ylabel('Relative difference (%)')
    plt.title('Evolution of relative differences on $\\tau (n,\gamma)$ (PyNjoy2016-oldlib)')
    plt.legend()
    plt.grid()
    plt.show()
    plt.savefig(f'{path_save}/SERPENT2_oldlib_PyNjoy2016_comparison_diff_ngamma_rates_vsBU.png')
    plt.close()

    # Plot evolution of relative differences on XS
    plt.figure()
    plt.plot(BU_axis, delta_Gd157_ngamma_XS["avg"], label='AVG $\\Delta$ Gd157')
    plt.plot(BU_axis, delta_Gd157_ngamma_XS["rms"], label='RMS $\\Delta$ Gd157')
    plt.plot(BU_axis, delta_Gd157_ngamma_XS["max"], label='MAX $\\Delta$ Gd157')

    plt.plot(BU_axis, delta_U238_ngamma_XS["avg"], label='AVG $\\Delta$ U238')
    plt.plot(BU_axis, delta_U238_ngamma_XS["rms"], label='RMS $\\Delta$ U238')
    plt.plot(BU_axis, delta_U238_ngamma_XS["max"], label='MAX $\\Delta$ U238')
    plt.xlabel('Burnup (MWd/t)')
    plt.ylabel('Relative difference (%)')
    plt.title('Evolution of relative differences on $\sigma (n,\gamma)$ (PyNjoy2016-oldlib)')
    plt.legend()
    plt.grid()
    plt.show()
    plt.savefig(f'{path_save}/SERPENT2_oldlib_PyNjoy2016_comparison_diff_ngamma_XS_vsBU.png')
    plt.close()
    print("End of post treatment of Gd157 XS and rates")

