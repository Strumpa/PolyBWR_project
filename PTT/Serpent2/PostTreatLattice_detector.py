# Python3 script making use of Serpent tools to post-treat the detector output of a Serpent2 simulation
# The script is used to extract the detector response from a lattice detector
# ATRIUM10 : UOX_2x2 and UOX_Gd_2x2 lattices to start with
# Author : R. Guasch
# Date : 27/11/2024

# Importing libraries
import os
import numpy as np
import serpentTools as st

# Path to the Serpent2 output file
path_to_serpent_results = "/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/PyNjoy2016_results/AT10_2x2"
# Name of the Serpent2 output file
case_name = "AT10_2x2_UOX"

# energy deposition mode 0 (default) : Constant energy deposition per fission.
# at bu = 0
detector = st.read(f"{path_to_serpent_results}/{case_name}_mc_det0.m")

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
# Extracting the detector response
# Strange that responses are not 0 for Pu isotopes fission rates at bu = 0 ?
for i in range(ngroups):
    for j in range(ncells):
        for k in range(ntallies):
            #print(f"Energy group {i+1}, Cell {j+1}, Reaction {tally_index_to_react[k]} : {detector.detectors['_pins_2G'].tallies[i,j,k]}")
            print(detector.detectors['_pins_2G'].tallies[i,j,k])
            # For now extracting only the fission rates of U235
            if k == 6:
                reaction_rates[f"U235_fiss_cell{j+1}_G{i+1}"] = detector.detectors['_pins_2G'].tallies[i,j,k]
print(reaction_rates)

