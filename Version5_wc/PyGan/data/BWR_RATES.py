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

def parse_Serpent2_lattice_det(path_to_S2, name_case, edepmode, bu):
    """
    Lattice detector post-treatment for Serpent2 simulations
    path_to_S2 : (str) path to the Serpent2 output files
    name_case : (str) name of the Serpent2 output file
    edepmode : (int) energy deposition mode
    bu : (int) burnup step
    """

    # energy deposition mode 0 (default) : Constant energy deposition per fission.
    # at bu = 0
    if edepmode == 0:
        detector = st.read(f"{path_to_S2}/{name_case}_mc_det{bu}.m")
    else:
        detector = st.read(f"{path_to_S2}/{name_case}_edep{edepmode}_mc_det{bu}.m")

    print(detector.detectors)
    print(detector.detectors["_pins_2G"].tallies)
    print(detector.detectors["_pins_2G"].tallies.shape) # shape = (2, 4, 10) = (energy group, cell, repsonse to mt ?)
    # Requested 2G energy groups, 4 cells/pins in the 2x2 lattice and 10 tallies / reactions (MT) in the detector
    # Reactions : MT=102 U235, U238, Pu239, Pu241, Xe135, Sm149 (microscopic rates), MT = -6 : fission macroscopic rates for U235, U238, Pu239, Pu241
    ngroups = detector.detectors["_pins_2G"].tallies.shape[0]
    ncells = detector.detectors["_pins_2G"].tallies.shape[1]
    ntallies = detector.detectors["_pins_2G"].tallies.shape[2]
    tally_index_to_react = {0: "U235_ngamma", 1: "U238_ngamma", 2: "Pu239_ngamma", 3: "Pu241_ngamma", 4: "Xe135_ngamma", 5: "Sm149_ngamma", 6: "U235_fiss", 7: "U238_fiss", 8: "Pu239_fiss", 9: "Pu241_fiss"}

    reaction_rates = {}
    U235_fiss_rates_S2 = {}
    # Extracting the detector response
    # Strange that responses are not 0 for Pu isotopes fission rates at bu = 0 ?
    for i in range(ngroups):
        for j in range(ncells):
            if f"cell{j+1}" not in U235_fiss_rates_S2.keys():
                U235_fiss_rates_S2[f"cell{j+1}"] = {}
            for k in range(ntallies):
                #print(f"Energy group {i+1}, Cell {j+1}, Reaction {tally_index_to_react[k]} : {detector.detectors['_pins_2G'].tallies[i,j,k]}")
                print(detector.detectors['_pins_2G'].tallies[i,j,k])
                # For now extracting only the fission rates of U235
                if k == 6:
                    reaction_rates[f"U235_fiss_cell{j+1}_G{i+1}"] = detector.detectors['_pins_2G'].tallies[i,j,k]
                    
                    U235_fiss_rates_S2[f"cell{j+1}"][f"G{i+1}"] = detector.detectors['_pins_2G'].tallies[i,j,k]
    return U235_fiss_rates_S2

def plot_field_on_geometry(geom_name, field):
    if geom_name == "AT10_2x2":
        n_cells = 4
        pitch = 1.295
        r_fuel_rod = 0.5140

        # place origin at (0,0), bottom left corner of the lattice
        x = [0, pitch, 2*pitch]
        y = [0, pitch, 2*pitch]
        X, Y = np.meshgrid(x, y)
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        ax.set_xlim(-0.5, 2*pitch+0.5)
        ax.set_ylim(-0.5, 2*pitch+0.5)
        
        #check number of field entries == number of cells
        if len(field) != n_cells:
            raise ValueError("Number of field entries does not match number of cells in the lattice")
        else:
            for i in range(n_cells):
                # plot fuel rods
                fuel_rod = mpatches.Circle((X[i], Y[i]), r_fuel_rod, color='black', fill=True)
                ax.add_patch(fuel_rod)
                # plot field value
                ax.text(X[i], Y[i], f"{field[i]:.2e}", ha='center', va='center', color='white')
        plt.savefig(f"{geom_name}_field.png")
    else:
        raise ValueError("Geometry not supported")




def BWR_CLUSTER(name_geom, name_case, name_compo, reaction_type, n_groups, bu):

    # Graphical settings (create a scale from white to red - possible to go from bleu to red to differenciate negative differences from positive ones)
    cmap = matplotlib.cm.get_cmap('bwr')

    sym = 'diagonal'
    if name_geom == "AT10_2x2":
        sym = "none"

    ################## FILE MANAGEMENT AND COMPOS ##################

    # --- Current directory path
    path = os.getcwd()
    # List of all the files in the current directory
    files = os.listdir(path)
    print(files)
    #print('$----------------------------'+path)
    
    # --- Creating the results storage directory
    name_DIR = f'DRAGON_RATES_{name_geom}'
    a = os.path.exists(name_DIR)
    if not a:
        os.mkdir(name_DIR)
    
    # Access paths to multicompo and serpent results in Resultats_serpent4PyGan
    SERPENT_path = f'/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/PyNjoy2016_results/{name_geom}/'
    
    # Recovering multi-compos
    #name_compo = 'COMPO_' + name_geom + name_case 
    print(f"name_compo = {name_compo}")
    lenBU_D5 = 1
    os.chdir("BWRresults_PyGan_2x2_UOX/")
    pyCOMPO = lcm.new('LCM_INP', name_compo, impx=0)
    os.chdir(path)
    #print(pyCOMPO['EDIBU_2gr']['MIXTURES'][0]['CALCULATIONS'][0]['MIXTURESVOL'])

    # Parse Serpent2 results

    U235_fission_rates_S2 = parse_Serpent2_lattice_det(SERPENT_path, name_case, 2, 0)

    if name_geom == "AT10_2x2_UOX":
        MIXES = [0,1,2,3]
    
    # Isotopes souhaites et nb total d'isotopes
    iso_study = ['U235 ', 'U238 ', 'Pu239', 'Pu241']
    reaction_name = 'NFTOT'
    n_reac = 6
    rates_names = 'fission rates'
    #pyCOMPO['EDIBU'].lib()
    len_isotot = np.shape(pyCOMPO['EDIBU']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESDENS'])[0] - 1
    ########## CALCULATIONS ##########
    
    # Calcul du facteur de normalisation
    prodD5 = 0.0
    prodS2 = 0.0
    ListS2 = []    
    nCell = 4
    MIXES = [0,1,2,3]
    # Create index to ALIAS dictionary
    Iso_index_to_ALIAS = {}
    U235_fiss_rate = {}
    for iso in range(len_isotot):
        isotope = pyCOMPO['EDIBU']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][iso]['ALIAS'][0:5]
        print(f"isotope = {isotope}, iso number = {iso}")
        Iso_index_to_ALIAS[iso] = isotope
    for iso in range(len_isotot):
        print(f"iso index = {iso}, isotope = {Iso_index_to_ALIAS[iso]}")
        isotope = pyCOMPO['EDIBU']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][iso]['ALIAS'][0:5]
        print(f"isotope = {isotope}")
        if isotope in ['U235 ']: #, 'U238 ']: #, 'Pu239', 'Pu241']:
            for mix in MIXES:
                U235_fiss_rate[f"mix{mix+1}"] = {}
                NWT0 = pyCOMPO['EDIBU_2gr']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NWT0']
                N = pyCOMPO['EDIBU']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESDENS'][iso]
                vol = pyCOMPO['EDIBU']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESVOL'][iso]
                NFTOT = pyCOMPO['EDIBU_2gr']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NFTOT']
#                NGAMMA = pyCOMPO['EDIBU_2gr']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NG']
                
                for gr in range(n_groups):
#                    prodD5 += 6*NWT0[gr]*(NFTOT[gr]+NGAMMA[gr])*N*vol
                    prodD5 += NWT0[gr]*NFTOT[gr]*N*vol
                    U235_fiss_rate[f"mix{mix+1}"][f"gr{gr+1}"] = NFTOT[gr]*NWT0[gr]*N*vol
                    print(f"mix{mix+1}, gr{gr+1}, U235_fiss_rate = {U235_fiss_rate[f'mix{mix+1}'][f'gr{gr+1}']}")
    
    print('$$$$$$$$$$$$$$$$$$$$$$$$$$', prodS2,prodD5)
    print([U235_fiss_rate[f"mix{mix+1}"][f"gr1"] for mix in MIXES])
    print([U235_fission_rates_S2[f"cell{mix+1}"][f"G1"] for mix in MIXES])
    prodS2 = sum([U235_fission_rates_S2[f"cell{mix+1}"][f"G{gr+1}"] for mix in MIXES for gr in range(n_groups)])
    renorm_factor = prodD5/prodS2

    #U235_fiss_rate = {f"mix{mix}": {f"gr{gr}": U235_fiss_rate[f"mix{mix+1}"][f"gr{gr}"]*renorm_factor for gr in range(1,3)} for mix in MIXES}

    U235_fission_rates_S2 = {f"cell{mix+1}": {f"G{gr+1}": U235_fission_rates_S2[f'cell{mix+1}'][f'G{gr+1}']*renorm_factor for gr in range(n_groups)} for mix in MIXES}
    #plot_field_on_geometry(name_geom, [U235_fiss_rate[f"mix{mix}"][f"gr1"] for mix in MIXES])
    print([U235_fiss_rate[f"mix{mix+1}"][f"gr2"] for mix in MIXES])
    print([U235_fission_rates_S2[f"cell{mix+1}"][f"G1"] for mix in MIXES])

    print([U235_fiss_rate[f"mix{mix+1}"][f"gr1"] for mix in MIXES])
    print([U235_fission_rates_S2[f"cell{mix+1}"][f"G2"] for mix in MIXES])



# Execute post treatment procedure


AT10_2x2_CLUSTER=['AT10_2x2_UOX', '_UOx6_autop5_PT_MOC_noCalc']
AT10_2x2_CLUSTER_GD=['AT10_2x2_UOX_Gd', '_UOx6_autop5_PT_MOC_noCalc']

name_geom = "AT10_2x2"
name_case = "AT10_2x2_UOX"
name_compo = "COMPO_AT10_2x2_UOX_UOx6_autop5_PT_MOC_noCalc" 
print(f"name_geom = {name_geom}, name_case = {name_case}, name_compo = {name_compo}")
BWR_CLUSTER(name_geom, name_case, name_compo, 'fission', 2, 0)







