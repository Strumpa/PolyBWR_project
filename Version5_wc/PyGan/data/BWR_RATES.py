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
from Serpent2_parsing import parse_S2_ASSBLY_rates, parse_Serpent2_material_det, parse_Serpent2_lattice_det, parse_S2_ASSBLY_rates_lat_det

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

def plot_error_grid_from_list(plot_title, name_case, name_compo, group, error_list, N, cmap="bwr", text_color="black"):
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
    cbar.set_label("Relative errors on $\\tau_f$ for U235 (%)", fontsize=12)

    # Show the plot
    plt.title("Spatial Error (D5-S2) Distribution on U235 fission rates", fontsize=14)
    plt.savefig(f"DRAGON_RATES_{name_case}/{plot_title}.png", dpi=300)
    plt.close(fig)
    return

def sum_rates_over_iso(rates_dict):
    """
    rates_dict is a nested dictionary with the structure:
    {
        "iso1": {
            "mix1": [val1, val2, ...],
            "mix2": [...],
            ...
        },
        ...
    }

    This function sums the values over isotopes while preserving the mix structure,
    returning:
    {
        "mix1": [sum_val1, sum_val2, ...],
        "mix2": [...],
        ...
    }
    """
    result = defaultdict(list)

    for isotope_data in rates_dict.values():
        for mix, values in isotope_data.items():
            if mix not in result:
                result[mix] = [0.0] * len(values)
            for i, val in enumerate(values):
                result[mix][i] += val

    return dict(result)

def plot_errors_BWR_assembly(errors_rates, name_case, name_compo, fig_name, unfold_symmetry, cmap="Spectral"):
    """
        Assembly map for ATRIUM-10

        C1 C2 C3 C5 C6 C5 C4 C3 C2 C1
        C2 C4 C7 C6 C7 C6 C6 C7 C4 C2
        C3 C7 C6 C6 C6 C7 C6 C6 C7 C3
        C5 C6 C6 C6 C6 C6 C7 C6 C5 C4
        C6 C7 C6 C6 WH WH WH C4 C6 C4
        C5 C6 C7 C6 WH WH WH C3 C7 C4
        C4 C6 C6 C7 WH WH WH C4 C4 C4
        C3 C7 C6 C6 C4 C3 C4 C4 C8 C3
        C2 C4 C7 C5 C6 C7 C4 C8 C4 C2
        C1 C2 C3 C4 C4 C4 C4 C3 C2 C1

        WH = Water Hole
        Cn = Fuel Cell number n

        Parameters:
        - errors_rates (list) : list containing dictionnaries with relative errors. Element in list corresponds to the group number.
            for each group number is a dictionary with key mix name and value the relative error.
        - name_case (str): Name of the geometry.
        - name_compo (str): Name of the MULTICOMPO obect with calculation parameter identifiers.
        - fig_name (str): Name of the figure to be saved.
        - unfold_symmetry (bool): If True, apply symmetry to the assembly map.
    
    """

    if unfold_symmetry:
            # Define the assembly map
            assembly_map = np.array([
                ["C1", "C2", "C3", "C5", "C6", "C5", "C4", "C3", "C2", "C1"],
                ["C2", "C4", "C7", "C6", "C7", "C6", "C6", "C7", "C4", "C2"],
                ["C3", "C7", "C6", "C6", "C6", "C7", "C6", "C6", "C7", "C3"],
                ["C5", "C6", "C6", "C6", "C6", "C6", "C7", "C6", "C5", "C4"],
                ["C6", "C7", "C6", "C6", "WH", "WH", "WH", "C4", "C6", "C4"],
                ["C5", "C6", "C7", "C6", "WH", "WH", "WH", "C3", "C7", "C4"],
                ["C4", "C6", "C6", "C7", "WH", "WH", "WH", "C4", "C4", "C4"],
                ["C3", "C7", "C6", "C6", "C4", "C3", "C4", "C4", "C8", "C3"],
                ["C2", "C4", "C7", "C5", "C6", "C7", "C4", "C8", "C4", "C2"],
                ["C1", "C2", "C3", "C4", "C4", "C4", "C4", "C3", "C2", "C1"]
            ])
    else:
            assembly_map = np.array([
                ["C1", "C2", "C3", "C5", "C6", "C5", "C4", "C3", "C2", "C1"],
                ["sy", "C4", "C7", "C6", "C7", "C6", "C6", "C7", "C4", "C2"],
                ["sy", "sy", "C6", "C6", "C6", "C7", "C6", "C6", "C7", "C3"],
                ["sy", "sy", "sy", "C6", "C6", "C6", "C7", "C6", "C5", "C4"],
                ["sy", "sy", "sy", "sy", "WH", "WH", "WH", "C4", "C6", "C4"],
                ["sy", "sy", "sy", "sy", "sy", "WH", "WH", "C3", "C7", "C4"],
                ["sy", "sy", "sy", "sy", "sy", "sy", "WH", "C4", "C4", "C4"],
                ["sy", "sy", "sy", "sy", "sy", "sy", "sy", "C4", "C8", "C3"],
                ["sy", "sy", "sy", "sy", "sy", "sy", "sy", "sy", "C4", "C2"],
                ["sy", "sy", "sy", "sy", "sy", "sy", "sy", "sy", "sy", "C1"]
                ])

    a = os.path.exists(f"DRAGON_RATES_{name_case}/{name_compo}")
    if not a:
        os.makedirs(f"DRAGON_RATES_{name_case}/{name_compo}")
    for group_num in range(len(errors_rates)):
        # Create a figure and axis
        fig, ax = plt.subplots(figsize=(10, 10))

        # Define the grid size
        grid_size = 10

        min_error = min(list(errors_rates[group_num].values()))
        max_error = max(list(errors_rates[group_num].values()))
        # Normalize colors based on error values
        norm = plt.Normalize(min_error, max_error)
        color_map = plt.cm.get_cmap(cmap)

        # Loop through the assembly map and plot each cell
        for row in range(assembly_map.shape[0]):
            for col in range(assembly_map.shape[1]):
                cell_name = assembly_map[row][col]
                print(f"cell_name = {cell_name}")
                
                # Get the error on reaction rate for the current cell
                if cell_name in errors_rates[group_num].keys():
                    error_value = errors_rates[group_num][cell_name] if cell_name in errors_rates[group_num].keys() else 0
                    color = color_map(norm(error_value))  # Get color from colormap
                else:
                    rate_g1 = cell_name
                # Get the reaction rate for the current cell
                # Plot the cell with the corresponding color
                if cell_name == "sy" or cell_name == "WH":
                    print("skip") 
                else:
                    rect = plt.Rectangle((col, grid_size - row - 1), 1, 1, facecolor=color, edgecolor='black', linewidth=2)
                    ax.add_patch(rect)
                    # Add the reaction rate text inside the cell
                    ax.text(col + 0.5, grid_size - row - 0.5, f"{error_value:.2f}", ha='center', va='center', fontsize=12, color='black')
        # Set the limits and aspect ratio
        ax.set_xlim(0, grid_size)
        ax.set_ylim(0, grid_size)
        ax.set_aspect('equal')
        # Set the title and labels
        ax.set_title(f'{fig_name.replace("_"," ")}, group {group_num+1}', fontsize=16)
        # Hide the axes
        ax.set_xticks([])
        ax.set_yticks([])
        # Show the plot
        plt.tight_layout()
        fig.savefig(f"DRAGON_RATES_{name_case}/{name_compo}/assembly_map_{fig_name}_g{group_num+1}.png", dpi=300)
        plt.close(fig)

    return

def plot_BWR_assembly(reaction_rates, name_case, name_compo, fig_name, unfold_symmetry):
        """
        Assembly map for ATRIUM-10

        C1 C2 C3 C5 C6 C5 C4 C3 C2 C1
        C2 C4 C7 C6 C7 C6 C6 C7 C4 C2
        C3 C7 C6 C6 C6 C7 C6 C6 C7 C3
        C5 C6 C6 C6 C6 C6 C7 C6 C5 C4
        C6 C7 C6 C6 WH WH WH C4 C6 C4
        C5 C6 C7 C6 WH WH WH C3 C7 C4
        C4 C6 C6 C7 WH WH WH C4 C4 C4
        C3 C7 C6 C6 C4 C3 C4 C4 C8 C3
        C2 C4 C7 C5 C6 C7 C4 C8 C4 C2
        C1 C2 C3 C4 C4 C4 C4 C3 C2 C1

        WH = Water Hole
        Cn = Fuel Cell number n

        Parameters:
        - reaction_rates (dict): Dictionary containing reaction rates for each cell.
        - name_case (str): Name of the geometry.
        - name_compo (str): Name of the MULTICOMPO obect with calculation parameter identifiers.
        - fig_name (str): Name of the figure to be saved.
        - unfold_symmetry (bool): If True, apply symmetry to the assembly map.
        """
        a = os.path.exists(f"DRAGON_RATES_{name_case}/{name_compo}")
        if not a:
            os.makedirs(f"DRAGON_RATES_{name_case}/{name_compo}")
        # Create a figure and axis
        fig1, ax1 = plt.subplots(figsize=(10, 10))
        fig2, ax2 = plt.subplots(figsize=(10, 10))

        # Define the grid size
        grid_size = 10


        # Define the colors for each cell
        colors = {
            "C1": "blue",
            "C2": "green",
            "C3": "yellow",
            "C4": "orange",
            "C5": "red",
            "C6": "purple",
            "C7": "pink",
            "C8": "cyan",
            "WH": "white",
            "sy": "lightgrey"
        }
        if unfold_symmetry:
            # Define the assembly map
            assembly_map = np.array([
                ["C1", "C2", "C3", "C5", "C6", "C5", "C4", "C3", "C2", "C1"],
                ["C2", "C4", "C7", "C6", "C7", "C6", "C6", "C7", "C4", "C2"],
                ["C3", "C7", "C6", "C6", "C6", "C7", "C6", "C6", "C7", "C3"],
                ["C5", "C6", "C6", "C6", "C6", "C6", "C7", "C6", "C5", "C4"],
                ["C6", "C7", "C6", "C6", "WH", "WH", "WH", "C4", "C6", "C4"],
                ["C5", "C6", "C7", "C6", "WH", "WH", "WH", "C3", "C7", "C4"],
                ["C4", "C6", "C6", "C7", "WH", "WH", "WH", "C4", "C4", "C4"],
                ["C3", "C7", "C6", "C6", "C4", "C3", "C4", "C4", "C8", "C3"],
                ["C2", "C4", "C7", "C5", "C6", "C7", "C4", "C8", "C4", "C2"],
                ["C1", "C2", "C3", "C4", "C4", "C4", "C4", "C3", "C2", "C1"]
            ])
        else:
            assembly_map = np.array([
                ["C1", "C2", "C3", "C5", "C6", "C5", "C4", "C3", "C2", "C1"],
                ["sy", "C4", "C7", "C6", "C7", "C6", "C6", "C7", "C4", "C2"],
                ["sy", "sy", "C6", "C6", "C6", "C7", "C6", "C6", "C7", "C3"],
                ["sy", "sy", "sy", "C6", "C6", "C6", "C7", "C6", "C5", "C4"],
                ["sy", "sy", "sy", "sy", "WH", "WH", "WH", "C4", "C6", "C4"],
                ["sy", "sy", "sy", "sy", "sy", "WH", "WH", "C3", "C7", "C4"],
                ["sy", "sy", "sy", "sy", "sy", "sy", "WH", "C4", "C4", "C4"],
                ["sy", "sy", "sy", "sy", "sy", "sy", "sy", "C4", "C8", "C3"],
                ["sy", "sy", "sy", "sy", "sy", "sy", "sy", "sy", "C4", "C2"],
                ["sy", "sy", "sy", "sy", "sy", "sy", "sy", "sy", "sy", "C1"]
                ])

        # Loop through the assembly map and plot each cell
        for row in range(assembly_map.shape[0]):
            for col in range(assembly_map.shape[1]):
                cell_name = assembly_map[row][col]
                if cell_name in colors.keys():
                    color = colors[cell_name]
                    # Get the reaction rate for the current cell
                    if cell_name in reaction_rates.keys():
                        rate_g1 = reaction_rates[cell_name][0]
                        rate_g2 = reaction_rates[cell_name][1]
                    else:
                        rate_g1 = cell_name
                        rate_g2 = cell_name
                    # Get the reaction rate for the current cell
                    # Plot the cell with the corresponding color
                    if cell_name == "sy" or cell_name == "WH":
                        print("skip") 
                    else:
                        rect1 = plt.Rectangle((col, grid_size - row - 1), 1, 1, color=color, edgecolor='black')
                        rect2 = plt.Rectangle((col, grid_size - row - 1), 1, 1, color=color, edgecolor='black')
                        ax1.add_patch(rect1)
                        ax2.add_patch(rect2)
                        # Add the reaction rate text inside the cell
                        ax1.text(col + 0.5, grid_size - row - 0.5, f"{rate_g1:.2f}", ha='center', va='center', fontsize=8, color='black')
                        ax2.text(col + 0.5, grid_size - row - 0.5, f"{rate_g2:.2f}", ha='center', va='center', fontsize=8, color='black')
        # Set the limits and aspect ratio
        ax1.set_xlim(0, grid_size)
        ax1.set_ylim(0, grid_size)
        ax1.set_aspect('equal')
        # Set the title and labels
        ax1.set_title(f"Assembly Map for {name_case} - {name_compo}, fast group", fontsize=16)
        ax1.set_xlabel("X-axis", fontsize=12)
        ax1.set_ylabel("Y-axis", fontsize=12)
        # Hide the axes
        ax1.set_xticks([])
        ax1.set_yticks([])
        # Show the plot
        plt.tight_layout()
        fig1.savefig(f"DRAGON_RATES_{name_case}/{name_compo}/assembly_map_{fig_name}_g1.png", dpi=300)
        plt.close(fig1)

        ax2.set_xlim(0, grid_size)
        ax2.set_ylim(0, grid_size)
        ax2.set_aspect('equal')
        # Set the title and labels
        ax2.set_title(f"Assembly Map for {name_case} - {name_compo},thermal group", fontsize=16)
        ax2.set_xlabel("X-axis", fontsize=12)
        ax2.set_ylabel("Y-axis", fontsize=12)
        # Hide the axes
        ax2.set_xticks([])
        ax2.set_yticks([])
        # Show the plot
        plt.tight_layout()
        fig2.savefig(f"DRAGON_RATES_{name_case}/{name_compo}/assembly_map_{fig_name}_g2.png", dpi=300)
        plt.close(fig2)
        return

def parse_DRAGON_rates_enrich_num(name_case, name_compo, fission_isotopes, n_gamma_isotopes, n_groups, bu, unfold_symmetry):
    """
    Parse DRAGON5 rates from the specified COMPO file.
    
    Parameters:
    - name_case (str): Name of the case.
    - name_compo (str): Name of the composition.
    - reaction_type (str): Type of reaction to parse.
    - n_groups (int): Number of energy groups.
    - bu (int): Burnup step.
    
    Returns:
    - None
    """
    if unfold_symmetry:
        sym_factor = 2
    else:
        sym_factor = 1
    # Load the DRAGON rates
    path = os.getcwd()
    os.chdir("Linux_aarch64")
    print(f"Loading {name_case} rates from {name_compo}")
    pyCOMPO = lcm.new('LCM_INP', name_compo, impx=0)
    os.chdir(path)
    # Retrieve the fission rates

    len_isotot = np.shape(pyCOMPO['EDIBU_HOM']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESDENS'])[0] - 1
    print(f"len_isotot = {len_isotot}")
    ########## CALCULATIONS ##########
    # Retrieve keff from pyCOMPO
    keff_D5 = pyCOMPO['EDIBU_HOM']['MIXTURES'][0]['CALCULATIONS'][0]['K-EFFECTIVE']
    print(f"keff_D5 = {keff_D5}")
    MIXES_idx = [0,1,2,3,4,5,6,7]
    MIXES = ["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"]
    """
        C1 C2 C3 C5 C6 C5 C4 C3 C2 C1
        C2 C4 C7 C6 C7 C6 C6 C7 C4 C2
        C3 C7 C6 C6 C6 C7 C6 C6 C7 C3
        C5 C6 C6 C6 C6 C6 C7 C6 C5 C4
        C6 C7 C6 C6 W1 WB W2 C4 C6 C4
        C5 C6 C7 C6 WL W0 WR C3 C7 C4
        C4 C6 C6 C7 W3 WT W4 C4 C4 C4
        C3 C7 C6 C6 C4 C3 C4 C4 C8 C3 
        C2 C4 C7 C5 C6 C7 C4 C8 C4 C2
        C1 C2 C3 C4 C4 C4 C4 C3 C2 C1
    
    """
    number_of_each_mix = {"C1":4, "C2":8, "C3":10, "C4":20, "C5":6, "C6":27, "C7":14, "C8":2}
    Iso_index_to_ALIAS = {}
    fiss_rates = {}
    n_gamma_rates = {}
    for iso in range(len_isotot):
        isotope = pyCOMPO['EDIBU_HOM']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][iso]['ALIAS'][0:5].strip()
        #print(f"isotope = {isotope}, iso number = {iso}")
        Iso_index_to_ALIAS[iso] = isotope
    for iso in range(len_isotot):
        #print(f"iso index = {iso}, isotope = {Iso_index_to_ALIAS[iso]}")
        isotope = pyCOMPO['EDIBU_HOM']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][iso]['ALIAS'][0:5].strip()
        #print(f"isotope = {isotope}")
        if isotope in fission_isotopes:
            isotope_fission_rate = {}
            for mix in MIXES_idx:
                isotope_fission_rate[f"C{mix+1}"] = {}
                NWT0 = pyCOMPO['HOM2g']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NWT0']
                N = pyCOMPO['HOM1g']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESDENS'][iso]
                vol = pyCOMPO['HOM1g']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESVOL'][iso]
                NFTOT = pyCOMPO['HOM2g']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NFTOT']
                isotope_fission_rate[f"C{mix+1}"] = np.array(NFTOT)*np.array(NWT0)*N*vol*sym_factor/number_of_each_mix[MIXES[mix]] # multiply volume by 2 to account for diagonal symmetry of the assembly
            fiss_rates[isotope] = isotope_fission_rate
        if isotope in n_gamma_isotopes:
            isotope_n_gamma_rate = {}
            for mix in MIXES_idx:
                isotope_n_gamma_rate[f"C{mix+1}"] = {}
                NWT0 = pyCOMPO['HOM2g']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NWT0']
                N = pyCOMPO['HOM1g']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESDENS'][iso]
                vol = pyCOMPO['HOM1g']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESVOL'][iso]
                NGAMMA = pyCOMPO['EDIBU_2gr']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NG']
                isotope_n_gamma_rate[f"C{mix+1}"] = np.array(NGAMMA)*np.array(NWT0)*N*vol*sym_factor/number_of_each_mix[MIXES[mix]] # multiply volume by 2 to account for diagonal symmetry of the assembly
            n_gamma_rates[isotope] = isotope_n_gamma_rate
    
    print(f"fiss_rates = {fiss_rates}")
    sum_summed_fiss_rates = sum_rates_over_iso(fiss_rates)

    sum_summed_n_gamma_rates = sum_rates_over_iso(n_gamma_rates)

    fiss_rates["TOT"] = sum_summed_fiss_rates
    #n_gamma_rates["TOT"] = sum_summed_n_gamma_rates

    return keff_D5, fiss_rates, n_gamma_rates

def parse_DRAGON_rates_regi_num(name_case, name_compo, fission_isotopes, n_gamma_isotopes, n_groups, bu, unfold_symmetry):
    """
    Parse DRAGON5 rates from the specified COMPO file.
    
    Parameters:
    - name_case (str): Name of the case.
    - name_compo (str): Name of the composition.
    - reaction_type (str): Type of reaction to parse.
    - n_groups (int): Number of energy groups.
    - bu (int): Burnup step.
    
    Returns:
    - None
    """
    if unfold_symmetry:
        sym_factor = 2
    else:
        sym_factor = 1
    # Load the DRAGON rates
    path = os.getcwd()
    os.chdir("Linux_aarch64")
    print(f"Loading {name_case} rates from {name_compo}")
    pyCOMPO = lcm.new('LCM_INP', name_compo, impx=0)
    os.chdir(path)
    # Retrieve the fission rates

    len_isotot = np.shape(pyCOMPO['EDIBU_HOM']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESDENS'])[0] - 1
    print(f"len_isotot = {len_isotot}")
    ########## CALCULATIONS ##########
    # Retrieve keff from pyCOMPO
    keff_D5 = pyCOMPO['EDIBU_HOM']['MIXTURES'][0]['CALCULATIONS'][0]['K-EFFECTIVE']
    print(f"keff_D5 = {keff_D5}")
    MIXES_idx = [0,1,2,3,4,5,6,7,8,9,10,
                11,12,13,14,15,16,17,18,19,
                20,21,22,23,24,25,26,27,
                28,29,30,31,32,33,34, 
                 # 35, 36, 37 are W1, WB, W2
                38, 39, 40,
                 # 41, 42 are  W0, WR
                43, 44, 45,
                # 46 is W4
                47, 48, 49,
                50, 51, 52, 
                53, 54, 
                55]
    MIXES = ["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"]
    """
    C1_1 C2_2 C3_3 C5_4 C6_5 C5_6 C4_7 C3_8 C2_9 C1_10
    C2_2 C4_11 C7_12 C6_13 C7_14 C6_15 C6_16 C7_17 C4_18 C2_19
    C3_3 C7_12 C6_20 C6_21 C6_22 C7_23 C6_24 C6_25 C7_26 C3_27
    C5_4 C6_13 C6_21 C6_28 C6_29 C6_30 C7_31 C6_32 C5_33 C4_34
    C6_5 C7_14 C6_22 C6_29 W1 WB W2 C4_38 C6_39 C4_40
    C5_6 C6_15 C7_23 C6_30 WL W0 WR C3_43 C7_44 C4_45
    C4_7 C6_16 C6_24 C7_31 W3 WT W4 C4_47 C4_48 C4_49
    C3_8 C7_17 C6_25 C6_32 C4_38 C3_43 C4_47 C4_50 C8_51 C3_52
    C2_9 C4_18 C7_26 C5_33 C6_39 C7_44 C4_48 C8_51 C4_53 C2_54
    C1_10 C2_19 C3_27 C4_34 C4_40 C4_45 C4_49 C3_52 C2_54 C1_55
    
    """
    number_of_each_mix = {"C1":4, "C2":8, "C3":10, "C4":20, "C5":6, "C6":27, "C7":14, "C8":2}
    Iso_index_to_ALIAS = {}
    fiss_rates = {}
    n_gamma_rates = {}
    for iso in range(len_isotot):
        #print(f"iso index = {iso}, isotope = {Iso_index_to_ALIAS[iso]}")
        isotope = pyCOMPO['EDIBU_HOM']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][iso]['ALIAS'][0:5].strip()
        #print(f"isotope = {isotope}")
        if isotope in fission_isotopes:
            isotope_fission_rate = {}
            for mix in MIXES_idx:
                isotope_fission_rate[f"C{mix+1}"] = {}
                NWT0 = pyCOMPO['HOM2g']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NWT0']
                N = pyCOMPO['HOM1g']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESDENS'][iso]
                vol = pyCOMPO['HOM1g']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESVOL'][iso]
                NFTOT = pyCOMPO['HOM2g']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NFTOT']
                isotope_fission_rate[f"C{mix+1}"] = np.array(NFTOT)*np.array(NWT0)*N*vol*sym_factor # multiply volume by 2 to account for diagonal symmetry of the assembly
            fiss_rates[isotope] = isotope_fission_rate
        """
        if isotope in n_gamma_isotopes:
            isotope_n_gamma_rate = {}
            for mix in MIXES_idx:
                isotope_n_gamma_rate[f"C{mix+1}"] = {}
                NWT0 = pyCOMPO['HOM2g']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NWT0']
                N = pyCOMPO['HOM1g']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESDENS'][iso]
                vol = pyCOMPO['HOM1g']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESVOL'][iso]
                NGAMMA = pyCOMPO['EDIBU_2gr']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NG']
                isotope_n_gamma_rate[f"C{mix+1}"] = np.array(NGAMMA)*np.array(NWT0)*N*vol*sym_factor/number_of_each_mix[MIXES[mix]] # multiply volume by 2 to account for diagonal symmetry of the assembly
            n_gamma_rates[isotope] = isotope_n_gamma_rate
        """
    print(f"fiss_rates = {fiss_rates}")
    sum_summed_fiss_rates = sum_rates_over_iso(fiss_rates)

    sum_summed_n_gamma_rates = sum_rates_over_iso(n_gamma_rates)

    fiss_rates["TOT"] = sum_summed_fiss_rates
    #n_gamma_rates["TOT"] = sum_summed_n_gamma_rates

    return keff_D5, fiss_rates, n_gamma_rates

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
#                NGAMMA = pyCOMPO['EDIBU_2gr']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NG']
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


if __name__ == "__main__":
    # AT10_2x2_UOX
    name_case = "AT10_2x2_UOX"
    name_compo = "COMPO_2x2_UOX_REFL_REFL" #"COMPO_2x2_UOX_12032025_RSE" #"COMPO_bench_3x3_UOX_CALC" # "COMPO_bench_3x3_UOX_C1_inrs1" # "COMPO_2x2_UOX_12032025_RSE" 
    print(f"name_case = {name_case}, name_case = {name_case}, name_compo = {name_compo}")
    #BWR_CLUSTER("AT10_2x2_UOX", "COMPO_2x2_UOX_REFL_REFL", 'fission', 2, 0)

    
    # bench_3x3_UOX
    name_case = "bench_3x3_UOX"
    name_compo = "COMPO_bench_3x3_UOX_REFL_REFL" #"COMPO_bench_3x3_UOX_CALC" # "COMPO_bench_3x3_UOX_C1_inrs1" # "COMPO_2x2_UOX_12032025_RSE"

    #BWR_CLUSTER("bench_3x3_UOX", "COMPO_bench_3x3_UOX_REFL_REFL", 'fission', 2, 0)

    # AT10_3x3_UOX_Gd
    name_case = "AT10_3x3_UOX_Gd"
    name_compo = "COMPO_AT10_3x3_UOX_Gd_REFL_REFL" #"COMPO_bench_3x3_UOX_CALC" # "COMPO_bench_3x3_UOX_C1_inrs1" # "COMPO_2x2_UOX_12032025_RSE"

    #BWR_CLUSTER("AT10_3x3_UOX_Gd", "COMPO_AT10_3x3_UOX_Gd_REFL_REFL_CORR_7inrsGd", 'fission', 2, 0)

    #BWR_CLUSTER("AT10_3x3_UOX_Gd", "COMPO_AT10_3x3_UOX_Gd_REFL_REFL_CORR_allinrs", 'fission', 2, 0)

    #BWR_CLUSTER("AT10_3x3_UOX_Gd", "COMPO_AT10_3x3_UOX_Gd_REFL_REFL_allinrs", 'fission', 2, 0)

    #BWR_CLUSTER("AT10_3x3_UOX_Gd", "COMPO_AT10_3x3_UOX_Gd_REFL_REFL_CTRA_allinrs_REGI_U238", 'fission', 2, 0)

    #BWR_CLUSTER("AT10_3x3_UOX_Gd", "COMPO_AT10_3x3_UOX_Gd_REFL_REFL_allinrs_SECT", 'fission', 2, 0)

    #BWR_CLUSTER("AT10_3x3_UOX_Gd", "COMPO_AT10_3x3_UOX_Gd_REFL_REFL_allinrs_SECT_REGI", 'fission', 2, 0)
    #COMPO_AT10_3x3_UOX_Gd_REFL_REFL_allinrs_SECT_fine
    #BWR_CLUSTER("AT10_3x3_UOX_Gd", "COMPO_AT10_3x3_UOX_Gd_REFL_REFL_allinrs_SECT_fine", 'fission', 2, 0)

    #BWR_CLUSTER("AT10_3x3_UOX_Gd", "COMPO_AT10_3x3_UOX_Gd_REFL_REFL_allinrs_SECT_REGI_fine", 'fission', 2, 0)


    #### ATRIUM-10 assembly
    name_case = "ATRIUM_10"
    CPOS_to_treat = ["CPO_AT10_ASSBLY_RSE_fine_J311", "CPO_AT10_ASSBLY_RSE_CORR_fine_J311"]
    #n_gamma_isotopes=["Gd155", "Gd157", "U238"]
    for name_compo in CPOS_to_treat:
        keff_D5, fiss_rates, ngamma_rates = parse_DRAGON_rates(name_case, name_compo, fission_isotopes=["U235", "U238"], n_gamma_isotopes=[], n_groups=2, bu=0, unfold_symmetry=True)

        renormed_tot_fiss_rate = renormalize_rates_dict(fiss_rates["TOT"])
        plot_BWR_assembly(renormed_tot_fiss_rate, name_case, name_compo, fig_name="Total_fission_rates_renorm", unfold_symmetry=True)
        
        S2keff, S2fission_rates, S2ngamma_rates = parse_S2_ASSBLY_rates(name_case="AT10_ASSBLY_t0_Rates", XS_lib_S2="J311_pynjoy2016", fission_isotopes=["U235", "U238", "U234"], n_gamma_isotopes=[], bu=0, unfold_symmetry=True)
        
        renormed_tot_fiss_rate_S2 = renormalize_rates_dict(S2fission_rates["TOT"])

        delta_rates = compute_diff_rates(renormed_tot_fiss_rate, renormed_tot_fiss_rate_S2)
        print(f"For case : {name_case}, and compo : {name_compo}")
        print(f"delta_keff = {(keff_D5 - S2keff)*1e5}")

        plot_errors_BWR_assembly(delta_rates, name_case, name_compo, fig_name="Total_fission_rates_diff", unfold_symmetry=True)
