# Python3 script for post treatment of Gd157 XS.
# Author : R. Guasch
# Date : 2024-11-04
# Purpose : Post treatment of Gd157 XS, Serpent2 vs DRAGON5 comparison, aimed at understanding the discrepancies between the two codes.
# Usage : python3 PostTreat_XS.py

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import serpentTools as st
from energyMeshHandler import energyMeshHandler


def plot_histograms(values_lists, iso, colors, grmin, grmax, legends=None, save_dir=None, save_name=None):
    # Check that the number of values is 295
    for values in values_lists:
        print(len(values))
        if len(values) != grmax - grmin + 1:
            raise ValueError(f"All input rates lists must contain exactly {grmax - grmin + 1} values.")
    
    # Create bins from 1 to 295
    energy_groups = list(range(grmin, grmax + 1))
    #colors = ['skyblue', 'red']
    # Plot the histogram
    
    plt.figure(figsize=(10, 6))
    for i in range(len(values_lists)):
        print(f"values_lists[{i}] = {values_lists[i]}")
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
        plt.title(f'Relative errors on XS for {iso}')
    else:
        plt.title(f'XS for {iso}')
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

def parse_DRAGON_data(case_name):
    """
    case_name associates to a specific DRAGON case and self-shielding method : HOM_UOX_Gd157_{noCORR/CORR}_{selfshielding_method} are supported.
    selfshielding_method : AUTO, RSE, PT, SUBG

    returns : reconstructed XS from DRAGON/VDG output
    """
    DRAGON_fluxes = []
    DRAGON_Gd_abs_rates = []
    AUTO_rates = []
    AUTO_fluxes = []

    filename = f"/home/p117902/working_dir/PolyBWR_project/Version5_wc/Dragon/Linux_aarch64/{case_name}.result"
    Gd157_counter = 0
    with open(filename) as f:
        lines = f.readlines()
    
    for i in range(len(lines)):
        if "VDG: PROCESSING ISOTOPE 'Gd157   0001'" in lines[i]:
            case_fluxes = []
            case_rates = []
            for line in lines[i+37:i+68]:
                case_fluxes.extend(line.strip().split())
            case_fluxes = [float(x) for x in case_fluxes]
            #print(f"Case fluxes are = {case_fluxes}")
            #print(f"Length of case fluxes = {len(case_fluxes)}")
            for line in lines[i+69:i+100]:
                case_rates.extend(line.strip().split())
            if Gd157_counter == 0:
                for line in lines[i+125:i+141]:
                    AUTO_fluxes.extend(line.strip().split())
                for line in lines[i+143:i+159]:
                    AUTO_rates.extend(line.strip().split())

            print(f"AUTO_fluxes = {AUTO_fluxes}")
            print(f"Length of AUTO_fluxes = {len(AUTO_fluxes)}")

            print(f"AUTO_rates = {AUTO_rates}")
            print(f"Length of AUTO_rates = {len(AUTO_rates)}")
            Gd157_counter += 1
            case_rates = [float(x) for x in case_rates]
            #print(f"Case rates are = {case_rates}")
            #print(f"Length of case rates = {len(case_rates)}")
            DRAGON_fluxes.append(np.array(case_fluxes))
            DRAGON_Gd_abs_rates.append(np.array(case_rates))
    AUTO_fluxes = [float(x) for x in AUTO_fluxes]
    AUTO_rates = [float(x) for x in AUTO_rates]
    

    return np.array(DRAGON_fluxes), np.array(DRAGON_Gd_abs_rates), np.array(AUTO_fluxes), np.array(AUTO_rates)


if __name__ == "__main__":
    # Parse Serpent microdepletion results
    input_case = "HOM_UOX_Gd157"
    correlation = "_noCORR"
    path_S2 = f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/HOM_CELL_study/{input_case}/XS_study"
    #library = "oldlib"

    local_path = os.getcwd()
    library = "PyNjoy2016"
    bu_step = 0
    GRMIN = 52
    GRMAX = 206

    
    path_exists = os.path.exists(f"{local_path}/{input_case}{correlation}/XS_study")
    if not path_exists:
        os.makedirs(f"{local_path}/{input_case}{correlation}/XS_study")
    path_save = f"{local_path}/{input_case}{correlation}/XS_study"

    LETHARGY_MESH = energyMeshHandler('SHEM295.txt', E0=1.0E+07, energyUnits='eV').lethargyMesh
    D5_fluxes, D5_Gd_abs_rates, AUTO_fluxes, AUTO_Gd_abs_rates = parse_DRAGON_data("HOM_UOX_Gd157_autop_noCORR")
    print(f"Length of D5_fluxes = {len(D5_fluxes)}")
    print(f"Length of D5_Gd_abs_rates = {len(D5_Gd_abs_rates)}")


    # Parse Serpent microdepletion results
    #U8_XS_102, Gd157_XS_102 = 
    Gd157_XS_ngamma, Gd157_XS_abstot, U8_XS_ngamma, U8_XS_abstot = parse_Serpent_microdepletion(path_S2, library, 0)
    N_Gd157 = Gd157_XS_abstot[0]
    N_U8 = U8_XS_abstot[0]

    # Exclude first entry as it is N_iso
    Gd157_XS_ngamma = Gd157_XS_ngamma[1:]
    Gd157_XS_abstot = Gd157_XS_abstot[1:]
    U8_XS_ngamma = U8_XS_ngamma[1:]
    U8_XS_abstot = U8_XS_abstot[1:]
    print(f"Length of Gd157_XS_ngamma = {len(Gd157_XS_ngamma)}")

    # Now recover XS from DRAGON D5_fluxes and Gd157 absorption rates, RR = Phi * Sigma_a = Phi * N * sigma_a

    print(f"D5_fluxes 0 = {D5_fluxes[0]}")
    micro_XS_a_Gd157_RSE = D5_Gd_abs_rates[0]/(D5_fluxes[0]*N_Gd157)
    micro_XS_a_Gd157_PT = D5_Gd_abs_rates[1]/(D5_fluxes[1]*N_Gd157)
    micro_XS_a_Gd157_SUBG = D5_Gd_abs_rates[2]/(D5_fluxes[2]*N_Gd157)

    AUTO_micro_XS_a_Gd157 = AUTO_Gd_abs_rates/(AUTO_fluxes*N_Gd157)

    # plot XS comparing Serpent2 and DRAGON
    print(f"Length of Gd157_XS_ngamma = {len(Gd157_XS_ngamma)}")
    print(f"len lethargy mesh = {len(LETHARGY_MESH)}")

    plt.figure()
    plt.step(LETHARGY_MESH[52:207], Gd157_XS_ngamma[52:207], label='Serpent2')
    plt.step(LETHARGY_MESH[52:207], AUTO_micro_XS_a_Gd157, label='AUTO')
    #plt.step(LETHARGY_MESH, micro_XS_a_Gd157_RSE, label='DRAGON RSE')
    #plt.step(LETHARGY_MESH, micro_XS_a_Gd157_PT, label='DRAGON PT')
    #plt.step(LETHARGY_MESH, micro_XS_a_Gd157_SUBG, label='DRAGON SUBG')
    plt.legend()
    plt.yscale('log')
    plt.show()
    plt.savefig(f'{path_save}/AUTO_S2_Gd157_XS_comparison_{input_case}{correlation}.png')
    plt.close()

    plt.figure()
    plt.step(LETHARGY_MESH[80:206], Gd157_XS_ngamma[80:206], label='Serpent2')
    plt.step(LETHARGY_MESH[80:206], micro_XS_a_Gd157_RSE[80:206], label='DRAGON RSE')
    plt.legend()
    plt.yscale('log')
    plt.show()
    plt.savefig(f'{path_save}/RSE_S2_Gd157_XS_comparison_{input_case}{correlation}.png')
    plt.close()


    plt.figure()
    plt.step(LETHARGY_MESH, Gd157_XS_ngamma, label='Serpent2')
    plt.step(LETHARGY_MESH, micro_XS_a_Gd157_PT, label='DRAGON PT')
    plt.legend()
    plt.yscale('log')
    plt.show()
    plt.savefig(f'{path_save}/PT_S2_Gd157_XS_comparison_{input_case}{correlation}.png')
    plt.close()

    plt.figure()
    plt.step(LETHARGY_MESH, Gd157_XS_ngamma, label='Serpent2')
    plt.step(LETHARGY_MESH, micro_XS_a_Gd157_SUBG, label='DRAGON SUBG')
    plt.legend()
    plt.yscale('log')
    plt.show()
    plt.savefig(f'{path_save}/SUBG_S2_Gd157_XS_comparison_{input_case}{correlation}.png')
    plt.close()

    # plot AUTO vs DRAGON fluxes
    plt.figure()
    plt.step(LETHARGY_MESH[52:207], AUTO_fluxes, label='AUTO')
    plt.step(LETHARGY_MESH, D5_fluxes[0], label='DRAGON RSE')
    plt.step(LETHARGY_MESH, D5_fluxes[1], label='DRAGON PT')
    plt.step(LETHARGY_MESH, D5_fluxes[2], label='DRAGON SUBG')
    plt.legend()
    plt.yscale('log')
    plt.show()
    plt.savefig(f'Flux_comparison_{input_case}{correlation}.png')
    plt.close()


    # Study groups 52 to 206
    # Compute errors on XS DRAGON vs Serpent2

    error_RSE_XS = (micro_XS_a_Gd157_RSE - Gd157_XS_ngamma)*100/Gd157_XS_ngamma
    error_PT_XS = (micro_XS_a_Gd157_PT - Gd157_XS_ngamma)*100/Gd157_XS_ngamma
    error_SUBG_XS = (micro_XS_a_Gd157_SUBG - Gd157_XS_ngamma)*100/Gd157_XS_ngamma
    error_AUTO_XS = (AUTO_micro_XS_a_Gd157 - Gd157_XS_ngamma[52:207])*100/Gd157_XS_ngamma[52:207]

    # Compute errors on fluxes DRAGON vs Serpent2 ? --> extract fluxes from Serpent2 output

    # Compute errors on XS DRAGON vs AUTO --> should be 0

    error_RSE_AUTO_XS = (micro_XS_a_Gd157_RSE[52:207] - AUTO_micro_XS_a_Gd157)
    error_PT_AUTO_XS = (micro_XS_a_Gd157_PT[52:207] - AUTO_micro_XS_a_Gd157)
    error_SUBG_AUTO_XS = (micro_XS_a_Gd157_SUBG[52:207] - AUTO_micro_XS_a_Gd157)
     
    # Compute errors on fluxes AUTO vs DRAGON --> should not be 0 !!! 
    error_RSE_AUTO_flx = (D5_fluxes[0][52:207] - AUTO_fluxes)
    error_PT_AUTO_flx = (D5_fluxes[1][52:207] - AUTO_fluxes)
    error_SUBG_AUTO_flx = (D5_fluxes[2][52:207] - AUTO_fluxes)
     


    #plot_histograms([ratio_D5, ratio_S2], f"{iso}/U8", ["skyblue", "red"], grmin, grmax, [f"D5 : {iso} abs rates / U8 abs rates", f"S2 : {iso} abs rates / U8 abs rates"], save_dir = path_save, save_name=f"hist_{iso}_U8_ratio_absorption_rates_{input_case}{correlation}_{SSH_method}")

    plot_histograms([error_AUTO_XS], f"Gd157", ["orange"], GRMIN, GRMAX, ["error on XS"], save_dir = path_save, save_name=f"AUTO_S2_Gd157_error_absorption_XS_{input_case}{correlation}")
    plot_histograms([error_RSE_XS[GRMIN:GRMAX+1], error_PT_XS[GRMIN:GRMAX+1], error_SUBG_XS[GRMIN:GRMAX+1]], f"Gd157", ["skyblue", "red", "green"], GRMIN, GRMAX, ["error on XS RSE", "error on XS PT", "error on XS SUBG"], save_dir = path_save, save_name=f"DRAGON_S2_Gd157_error_absorption_XS_{input_case}{correlation}")
    plot_histograms([error_RSE_XS[GRMIN:GRMAX+1]], f"Gd157", ["skyblue"], GRMIN, GRMAX, ["error on XS RSE"], save_dir = path_save, save_name=f"RSE_S2_Gd157_error_absorption_XS_{input_case}{correlation}")
    plot_histograms([error_PT_XS[GRMIN:GRMAX+1]], f"Gd157", ["red"], GRMIN, GRMAX, ["error on XS PT"], save_dir = path_save, save_name=f"PT_S2_Gd157_error_absorption_XS_{input_case}{correlation}")
    plot_histograms([error_SUBG_XS[GRMIN:GRMAX+1]], f"Gd157", ["green"], GRMIN, GRMAX, ["error on XS SUBG"], save_dir = path_save, save_name=f"SUBG_S2_Gd157_error_absorption_XS_{input_case}{correlation}")