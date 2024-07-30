# Python3 script to generate the pins and detectors for the ATRIUM-10 assembly.
# Purpose : automate the generation of Serpent2 reference case
# Author : R. Guasch
# Date : 2024-07-30

import DMLGInterface as DMLG
import GeometricTools as GeoT
import numpy as np
from Composition_handler import Compo_Processor



def generateS2_detectors(individual_pins, isotopes_for_fission_scores, isotopes_for_radiative_capture_scores):
    """
    Generate the pins and detectors for the ATRIUM-10 assembly

    individual_pins : list of pins, allowing to request detector scores for each pin.
    Detailed materials in a single pin are merged to collect pin-by-pin detector scores.
    requested reaction rates for :  - macroscopic total fission cross section : -6 reaction number for isotopes in isotopes_for_fission_scores
                                    - radiative capture : 102 MT reaction for isotopes in isotopes_for_radiative_capture_scores

    Assuming a 2 energy group structure for the detectors has already been defined, labeled 1.
    """
    print("% $$$ ---- Begin detectors definition ---- $$$")
    for pin_name in individual_pins:
        print(f"det _{pin_name}_2G")
        print(f"    de 1")
        print(f"    dt -4")
        mat, number = pin_name.split("_")
        if "GADO" in mat or "Gd" in mat:
            print(f"    dm {mat}_A_{number} dm {mat}_B_{number} dm {mat}_C_{number} dm {mat}_D_{number} dm {mat}_E_{number} dm {mat}_F_{number}")
        else:
            print(f"    dm {mat}_A_{number} dm {mat}_B_{number} dm {mat}_C_{number} dm {mat}_D_{number}")
        for iso in isotopes_for_fission_scores:
            print(f"    dr -6 {iso}")
        for iso in isotopes_for_radiative_capture_scores:
            print(f"    dr 102 {iso}")
        
        print(f"\n")


    return

def generateS2_pins(ATRIUM_10_NL_24_bundle):
    """
    Generate the pins for the ATRIUM-10 assembly
    """
    print(f"% $$$ ---- Begin pins definition ---- $$$")
    detailed_material_names = []
    individual_pins = []
    for pin in ATRIUM_10_NL_24_bundle.pins:
        pin_name = pin.name
        number_of_pins = ATRIUM_10_NL_24_bundle.numberOfPinsperType[pin_name]
        print(f"% ---> {number_of_pins} pins of type {pin_name}")
        for i in range(number_of_pins):
            print(f"pin {pin_name}_{i+1}")
            for vol in pin.volumes.keys():
                if "Gap" in vol or "gap" in vol:
                    print(f"    gap {pin.volume_names_to_radii[vol]:.06f}")
                elif "Clad" in vol or "clad" in vol:
                    print(f"    clad {pin.volume_names_to_radii[vol]:.06f}")
                else:
                    print(f"    {vol}_{i+1} {pin.volume_names_to_radii[vol]:.06f}")
                    detailed_material_names.append(f"{vol}_{i+1}")
                    [mat,crown]=vol.split("_")
                    if f"{mat}_{i+1}" not in individual_pins:
                        individual_pins.append(f"{mat}_{i+1}")
            print(f"    cool")
            print(f"\n")
    return detailed_material_names, individual_pins

def generateS2_materials(detailed_material_names, temperatures_dict, S2_lib, path_to_composition_files="../compo_data/"):
    """
    Generate the S2 materials for the ATRIUM-10 assembly, using the detailed material names, compositions and temperatures.
    detailed_material_names : list of detailed material names, to be used as S2 material names and associated to indiviual sub-pin volumes,
    temperatures_dict : dictionary associating the material names to the temperatures, the temperatures are floats in Kelvin.
    S2_lib : string, the Serpent2 xs library used, either "PyNjoy2016" or "sss_jeff311". Required to set the correct temperature suffix.
    """
    print("% $$$ ---- Begin material definition ---- $$$")
    print(f"% ---> Fuel Materials")
    for detailed_mat in detailed_material_names:
        material_label = detailed_mat.split("_")[0]
        print(f"mat {detailed_mat} sum tmp {temperatures_dict[material_label]} burn 1")
        if S2_lib == "PyNjoy2016":
            if temperatures_dict[material_label] == 750.0:
                temp_suffix = ".05c"
            elif temperatures_dict[material_label] == 559.0:
                temp_suffix = ".05c"
        elif S2_lib == "sss_jeff311":
            if temperatures_dict[material_label] == 750.0:
                temp_suffix = ".06c"
            elif temperatures_dict[material_label] == 559.0:
                temp_suffix = ".03c"
        # Use Composition_handler to parse compositions for the ATRIUM-10 assembly from the composition files
        material_compo = Compo_Processor(material_label, path_to_composition_files+f"{material_label}_compo.txt")
        composition_dict = material_compo.compo
        for isotope_code in composition_dict.keys():
            print(f"    {isotope_code}{temp_suffix} {composition_dict[isotope_code]}")
        print(f"\n")
       

    return


# Use GeometricTools to generate the pins and detectors for the ATRIUM-10 assembly from the analytical definition,

# Geometrical data for the ATRIUM-10 assembly :

path_to_composition_files = "../compo_data/"

Channel_box_out = 2.3975+1.1025
Channel_box_in = 2.3175+1.0225
pitch_A=7.62*2
Box_o = 6.87*2
Box_i = 6.7*2
Channel_box_xL_out = -1.1025
Channel_box_XR_out = 2.3975

pins_names=["24UOx", "32UOx", "42UOx", "45UOx", "48UOx", "50UOx", "45GADO", "42GADO"]
pin_radii =[0.4435,0.4520,0.5140] # Fuel, gap, clad radii

pins_number_dict={"24UOx": 4, "32UOx" : 8, "42UOx": 10, "45UOx": 21, "48UOx": 6, "50UOx":26, "45GADO": 14, "42GADO": 2}

# Temperature data for fuel compositions, in Kelvin,
material_temperatures_dict = {"24UOx": 750.0, "32UOx" : 750.0, "42UOx": 750.0, "45UOx": 750.0, "48UOx": 750.0, "50UOx": 750.0, "45GADO": 750.0, "42GADO": 750.0}


# Create the ATRIUM-10 assembly bundle, with the pins and the number of pins per type,
ATRIUM_10_NL_24_bundle = GeoT.geom_ASSBLY(pitch_A,pins_names, Box_o, Channel_box_out, Box_i, Channel_box_in)
ATRIUM_10_NL_24_bundle.setPins(pin_radii[0], pin_radii[1], pin_radii[2])
ATRIUM_10_NL_24_bundle.setNumberofPins(pins_number_dict)

# Choose the isotopes for the fission and radiative capture reaction rate scores,

isotopes_for_fission_scores = ["U234", "U235", "U238", "Pu239", "Pu241"]
isotopes_for_radiative_capture_scores = ["U234", "U235", "U238", "Pu239", "Pu241", "Gd155", "Gd157", "Xe135", "Sm149"]




# retrieve detailed material names and generate the pins for the ATRIUM-10 assembly, allowing for the generation of the S2 materials and S2 detectors,
detailed_material_names, individual_pins = generateS2_pins(ATRIUM_10_NL_24_bundle)

# Generate the S2 materials for the ATRIUM-10 assembly, using the detailed material names, compositions and temperatures,
generateS2_materials(detailed_material_names, material_temperatures_dict, "PyNjoy2016", path_to_composition_files)

# Generate the detectors for the ATRIUM-10 assembly, using the individual pins and the isotopes for the fission and gamma scores
generateS2_detectors(individual_pins, isotopes_for_fission_scores, isotopes_for_radiative_capture_scores)

