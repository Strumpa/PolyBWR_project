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

matplotlib.use('Agg')

def parse_Serpent2_lattice_det(path_to_S2, name_case, XS_lib_S2, edepmode, pcc, bu):
    """
    Lattice detector post-treatment for Serpent2 simulations
    path_to_S2 : (str) path to the Serpent2 output files
    name_case : (str) name of the Serpent2 output file
    XS_lib_S2 : (str) name of the Serpent2 XS library
    edepmode : (int) energy deposition mode
    pcc : (int) predictor corrector option
    bu : (int) burnup step
    """

    # energy deposition mode 0 (default) : Constant energy deposition per fission.
    # at bu = 0
    if edepmode == 0:
        detector = st.read(f"{path_to_S2}/{name_case}_{XS_lib_S2}_mc_det{bu}.m")
    else:
        detector = st.read(f"{path_to_S2}/{name_case}_edep{edepmode}_mc_det{bu}.m")

    ### _pins_2G is a lattice detector
    # For 2x2 lattice :
    #print(detector.detectors["_pins_2G"].tallies.shape) # shape = (2, 4, 10) = (energy group, cell, repsonse to mt ?)
    # Requested 2G energy groups, 4 cells/pins in the 2x2 lattice and 10 tallies / reactions (MT) in the detector
    # Reactions : MT=102 U235, U238, Pu239, Pu241, Xe135, Sm149 (microscopic rates), MT = -6 : fission "macroscopic" rates for U235, U238, Pu239, Pu241
    # For 3x3 lattice :
    # Requested 2G energy groups, 9 cells/pins in the 3x3 lattice and 10 tallies / reactions (MT) in the detector
    # Reactions : MT=102 U235, U238, Pu239, Pu241, Xe135, Sm149 (microscopic rates), MT = -6 : fission "macroscopic" rates for U235, U238, Pu239, Pu241
    # Expected detector shape : (2, 9, 10) 
    ngroups = detector.detectors["_pins_2G"].tallies.shape[0]
    ncells = detector.detectors["_pins_2G"].tallies.shape[1]
    ntallies = detector.detectors["_pins_2G"].tallies.shape[2]
    tally_index_to_react = {0: "U235_ngamma", 1: "U238_ngamma", 2: "Pu239_ngamma", 3: "Pu241_ngamma", 4: "Xe135_ngamma", 5: "Sm149_ngamma", 6: "U235_fiss", 7: "U238_fiss", 8: "Pu239_fiss", 9: "Pu241_fiss"}
    if name_case == "AT10_2x2_UOX":
                #      C1           C2           C2             C4
        N_U235 = [5.67035E-04, 7.560370E-04, 7.560370E-04, 1.051340E-03] # U235 densities in each cell, in 1. positive x, 2. positive y order
    elif name_case == "bench_3x3_UOX":
                       
        N_U235 = [5.67035E-04, 7.560370E-04, 9.686590E-04, #C1, C2, C3  
                  7.560370E-04, 1.051340E-03, 5.67035E-04, #C2, C4, C1 
                  9.686590E-04, 5.67035E-04, 1.169460E-03] #C3, C1, C6
    reaction_rates = {}
    U235_fiss_rates_S2_lat_det = {}
    # Extracting the detector response
    # Strange that responses are not 0 for Pu isotopes fission rates at bu = 0 ?
    for i in range(ngroups):
        for j in range(ncells):
            if f"cell{j+1}" not in U235_fiss_rates_S2_lat_det.keys():
                U235_fiss_rates_S2_lat_det[f"cell{j+1}"] = {}
            for k in range(ntallies):
                print(f"Energy group {i+1}, Cell {j+1}, Reaction {tally_index_to_react[k]} : {detector.detectors['_pins_2G'].tallies[i,j,k]}")
                #print(detector.detectors['_pins_2G'].tallies[i,j,k])
                if tally_index_to_react[k] not in reaction_rates.keys():
                    reaction_rates[tally_index_to_react[k]] = {}
                reaction_rates[tally_index_to_react[k]][f"cell{j+1}_G{i+1}"] = detector.detectors['_pins_2G'].tallies[i,j,k]    
                # For now extracting only the fission rates of U235
                if k == 6:
                    #reaction_rates[f"U235_fiss_cell{j+1}_G{i+1}"] = detector.detectors['_pins_2G'].tallies[i,j,k]*N_U235[f"C{j+1}"]
                    U235_fiss_rates_S2_lat_det[f"cell{j+1}"][f"G{i+1}"] = detector.detectors['_pins_2G'].tallies[i,j,k]*N_U235[j]

    # Symmetrizing the fission rates for the 2x2 lattice
    if name_case == "AT10_2x2_UOX":
        U235_fiss_rates_S2_lat_det["cell2"]["G1"] = (U235_fiss_rates_S2_lat_det["cell2"]["G1"] + U235_fiss_rates_S2_lat_det["cell3"]["G1"]) / 2
        U235_fiss_rates_S2_lat_det["cell3"]["G1"] = U235_fiss_rates_S2_lat_det["cell2"]["G1"]
        U235_fiss_rates_S2_lat_det["cell2"]["G2"] = (U235_fiss_rates_S2_lat_det["cell2"]["G2"] + U235_fiss_rates_S2_lat_det["cell3"]["G2"]) / 2
        U235_fiss_rates_S2_lat_det["cell3"]["G2"] = U235_fiss_rates_S2_lat_det["cell2"]["G2"]
    
    return U235_fiss_rates_S2_lat_det

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

def plot_error_grid_from_list(name_case, name_compo, group, error_list, N, cmap="bwr", text_color="black"):
    """
    Plots an N x N grid where each cell is colored according to the error value.
    
    Parameters:
    - name_case (str): Name of the geometry.
    - group (int): Group number.
    - error_list (list or numpy.ndarray): 1D list of error values, ordered by increasing x first, then increasing y.
    - N (int): Size of the grid (assumed to be N x N).
    - cmap (str): Colormap name (default is 'coolwarm').
    - text_color (str): Color of the text inside the squares (default is black).
    
    Returns:
    - None (displays the plot)
    """
    # Convert list into an NxN array
    error_matrix = np.array(error_list).reshape(N, N)  # Reshape row-wise

    fig, ax = plt.subplots(figsize=(8, 8))

    # Normalize colors based on error values
    min_err, max_err = np.min(error_matrix), np.max(error_matrix)
    norm = plt.Normalize(min_err, max_err)
    color_map = plt.cm.get_cmap(cmap)

    # Plot each square
    for i in range(N):
        for j in range(N):
            value = error_matrix[i, j]
            color = color_map(norm(value))  # Get color from colormap
            
            # The input ordering assumes row-major format (increasing x first), 
            # but Matplotlib's (0,0) is at the top-left, so we flip vertically
            ax.add_patch(plt.Rectangle((j, i), 1, 1, color=color, edgecolor='black'))

            # Add text at the center of each square
            ax.text(j + 0.5, i + 0.5, f"{value:.2f} %", 
                    ha='center', va='center', fontsize=15, color=text_color)

    # Set axis limits
    ax.set_xlim(0, N)
    ax.set_ylim(0, N)

    # Remove axes ticks
    ax.set_xticks([])
    ax.set_yticks([])

    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=color_map, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Relative errors on \\tau_f for U235 (%)", fontsize=12)

    # Show the plot
    plt.title("Spatial Error (D5-S2) Distribution on U235 fission rates", fontsize=14)
    plt.savefig(f"DRAGON_RATES_{name_case}/{name_compo}_error_grid_U235_fission_rates_g{group}.png", dpi=300)


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

    # Parse Serpent2 results
    U235_fission_rates_S2 = parse_Serpent2_lattice_det(SERPENT_path, name_case, XS_lib_S2="j311_pynjoy2016" , edepmode=0, pcc=1, bu=0)

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
        MIXES = ["C1", "C2", "C3", "C4", "C6"]    
        number_of_each_mix = {"C1":3, "C2":2, "C3":2, "C4":1, "C6":1}
    # Isotopes souhaites et nb total d'isotopes
    iso_study = ['U235 ', 'U238 ', 'Pu239', 'Pu241']
    reaction_name = 'NFTOT'
    n_reac = 6
    rates_names = 'fission rates'
    #pyCOMPO['EDIBU'].lib()
    #print(pyCOMPO.keys())
    len_isotot = np.shape(pyCOMPO['EDIBU_HOM']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESDENS'])[0] - 1
    ########## CALCULATIONS ##########
    
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
                NFTOT = pyCOMPO['HOM2g']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NFTOT']
#                NGAMMA = pyCOMPO['EDIBU_2gr']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NG']
                print(f"mix index is = {mix} with corresponding C name = {MIXES[mix]}")
                print(f"Number of each mix = {number_of_each_mix[MIXES[mix]]}")
                for gr in range(n_groups):
                    U235_fiss_rate[f"mix{mix+1}"][f"gr{gr+1}"] = NFTOT[gr]*NWT0[gr]*N*vol/number_of_each_mix[MIXES[mix]]

    print(f"U235_fiss_rate = {U235_fiss_rate}")
    # duplicate rates from individual mixes to form the cells lattice
    U235_fiss_rates_cells = {}
    
    for cell in cell_S2_to_mix_indx.keys():
        U235_fiss_rates_cells[f"cell{cell+1}"] = {}
        for gr in range(n_groups):
            print(f"adding score of mix{cell_S2_to_mix_indx[cell]+1} to cell{cell+1}, in group {gr+1}")
            U235_fiss_rates_cells[f"cell{cell+1}"][f"G{gr+1}"] = U235_fiss_rate[f"mix{cell_S2_to_mix_indx[cell]+1}"][f"gr{gr+1}"]
    # normalize the rates : sum of rates = nCells
    # group 1 : U235 fast fissions : note that group numbers are in increasing lethargy in Dragon and in increasing energy in Serpent
    U235_fiss_rate_norm_1 = renormalize_rates([U235_fiss_rates_cells[cell]["G1"] for cell in U235_fiss_rates_cells.keys()])
    U235_fission_rates_S2_norm_1 = renormalize_rates([U235_fission_rates_S2[f"cell{mix_idx+1}"]["G2"] for mix_idx in cell_S2_to_mix_indx.keys()])

    # group 2 : U235 thermal fissions
    U235_fiss_rate_norm_2 = renormalize_rates([U235_fiss_rates_cells[cell]["G2"] for cell in U235_fiss_rates_cells.keys()])
    U235_fission_rates_S2_norm_2 = renormalize_rates([U235_fission_rates_S2[f"cell{mix_idx+1}"]["G1"] for mix_idx in cell_S2_to_mix_indx.keys()])

    relative_difference1 = [(U235_fiss_rate_norm_1[i] - U235_fission_rates_S2_norm_1[i])*100/U235_fission_rates_S2_norm_1[i] for i in range(nCell)]
    relative_difference2 = [(U235_fiss_rate_norm_2[i] - U235_fission_rates_S2_norm_2[i])*100/U235_fission_rates_S2_norm_2[i] for i in range(nCell)]

    #print(f"relative_difference1 = {relative_difference1}")
    plot_error_grid_from_list(name_case, name_compo, 1, relative_difference1, N_col) 
    plot_error_grid_from_list(name_case, name_compo, 2, relative_difference2, N_col) 
    return relative_difference1, relative_difference2

# Execute post treatment procedure


if __name__ == "__main__":
    name_case = "AT10_2x2_UOX" # "bench_3x3_UOX" #"AT10_2x2_UOX"
    name_compo = "COMPO_2x2_UOX_TRAN" #"COMPO_2x2_UOX_12032025_RSE" #"COMPO_bench_3x3_UOX_CALC" # "COMPO_bench_3x3_UOX_C1_inrs1" # "COMPO_2x2_UOX_12032025_RSE" 
    print(f"name_case = {name_case}, name_case = {name_case}, name_compo = {name_compo}")

    #diffg1,diffg2 = BWR_CLUSTER(name_case, name_compo, 'fission', 2, 0)
    #print(f"diffg1 = {diffg1}")
    #print(f"diffg2 = {diffg2}")


    #name_case = "AT10_2x2_UOX"
    #name_compo = "COMPO_2x2_UOX_12032025_RSE" 
    #print(f"name_case = {name_case}, name_case = {name_case}, name_compo = {name_compo}")

    #diffg1,diffg2 = BWR_CLUSTER(name_case, name_compo, 'fission', 2, 0)
    #print(f"diffg1 = {diffg1}")
    #print(f"diffg2 = {diffg2}")

    name_case = "bench_3x3_UOX"
    name_compo = "COMPO_bench_3x3_UOX_TRAN" # "COMPO_bench_3x3_UOX_TRAN" #"COMPO_bench_3x3_UOX_C1_inrs1" # "COMPO_2x2_UOX_12032025_RSE"
    diffg1,diffg2 = BWR_CLUSTER(name_case, name_compo, 'fission', 2, 0)
    print(f"name_case = {name_case}, name_case = {name_case}, name_compo = {name_compo}")
    print(f"diffg1 = {diffg1}")
    print(f"diffg2 = {diffg2}")


