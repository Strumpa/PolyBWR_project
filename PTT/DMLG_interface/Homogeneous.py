# Python3 scritp used to homogenize BWR pin cell
# Author : R. Guasch
# Date : 17/07/2024
# Purpose : Homogenize BWR pin cell to check Serpent2 vs Dragon5 bateman evolution
# Part of the BWR project, assessing the differences between PyNjoy2016 and sss_jeff311 Serpent2 libraries.

import numpy as np
from GeometricTools import geom_PIN
from Composition_handler import Compo_Processor

def create_exact_BWR_pin(name, fuel_radius, gap_radius, clad_radius, isGd, height, pitch, path_to_composition_files):
    """
    Create BWR pin cell geometry
    Call the GeometricTools class to create a geom_PIN object.
    geom_PIN requires arguments: label, fuel_radius, gap_radius, clad_radius, isGd, height
    """
    exact_pin = geom_PIN(name, fuel_radius, gap_radius, clad_radius, isGd, height, pitch, sectorize=False)
    fuel_compo = Compo_Processor("45Gd_Fuel", path_to_composition_files+"FUEL_compo.txt")
    fuel_compo.check_consistency()
    fuel_compo.iso_code_to_nuclide()
    
    gap_compo = Compo_Processor("45Gd_gap", path_to_composition_files+"GAP_compo.txt")
    gap_compo.check_consistency()
    gap_compo.iso_code_to_nuclide()

    clad_compo = Compo_Processor("45Gd_clad", path_to_composition_files+"CLAD_compo.txt")
    clad_compo.check_consistency()
    clad_compo.iso_code_to_nuclide()

    coolant_compo = Compo_Processor("45Gd_coolant", path_to_composition_files+"COOLANT_compo.txt")
    coolant_compo.check_consistency()
    coolant_compo.iso_code_to_nuclide()

    mix_numbering_dict = {"1": fuel_compo.Ndens_isos, "2": gap_compo.Ndens_isos, "3": clad_compo.Ndens_isos, "4": coolant_compo.Ndens_isos}
    exact_pin.setMixes(mix_numbering_dict)

    return exact_pin

def homogenize_pin(exact_pin, pitch_for_homogenization):
    """
    Homogenize BWR pin cell to a single mix with the same pitch
    """
    homogeneous_mix = {}
    homogeneous_vol = pitch_for_homogenization**2*exact_pin.height
    print(exact_pin.volumes)
    print(exact_pin.mixes)
    for mix in exact_pin.mixes.keys():
        for iso in exact_pin.mixes[mix].keys():
            if iso not in homogeneous_mix.keys():
                homogeneous_mix[iso] = exact_pin.mixes[mix][iso]*exact_pin.volumes[exact_pin.name+"_"+exact_pin.mix_to_name_connectivity_dict[mix]]/homogeneous_vol
            else:
                homogeneous_mix[iso] += exact_pin.mixes[mix][iso]*exact_pin.volumes[exact_pin.name+"_"+exact_pin.mix_to_name_connectivity_dict[mix]]/homogeneous_vol
    return homogeneous_mix

def print_to_Dragon5_format(mix):
    """
    rewrite the info in Dragon iso = iso [atomic density] format
    isotopes as keys and atomic density as values
    """
    print("$$$-- Output in Dragon5 format --$$$")
    for iso in mix.keys():
        #print(str(iso)+"  = "+str(iso)+"   "+str((mix[iso])))
        print(f"{str(iso)}  = {str(iso)}   {str((mix[iso]))}")

def get_code_from_iso(dictionary, target_value):
    for key, value in dictionary.items():
        if value == target_value:
            return key
    return None  # or raise an exception if you prefer

def print_to_Serpent2_format(mix, suffix_for_temperature):
    """
    rewrite the info in Serpent2 format
    isotopes as keys and atomic density as values
    """
    print("$$$-- Output in Serpent2 format --$$$")
    association_dict = {
                        "1001" : "H1", "1002" : "H2", "1003" : "H3",
                        "2003" : "He3", "2004" : "He4",
                        "8016" : "O16", "8017" : "O17", 
                        "24050": "Cr50", "24052": "Cr52", "24053": "Cr53", "24054": "Cr54", 
                        "26054": "Fe54", "26056": "Fe56", "26057": "Fe57", "26058": "Fe58",
                        "28058": "Ni58", "28060": "Ni60", "28061": "Ni61", "28062": "Ni62", "28064": "Ni64",
                        "40090": "Zr90", "40091": "Zr91", "40092": "Zr92", "40094": "Zr94", "40096": "Zr96",
                        "50112": "Sn112", "50114": "Sn114", "50115": "Sn115", "50116": "Sn116","50117": "Sn117",  "50118": "Sn118", "50119": "Sn119", "50120": "Sn120", "50122": "Sn122", "50124": "Sn124",
                        "64154": "Gd154", "64155": "Gd155", "64156": "Gd156", "64157": "Gd157", "64158": "Gd158", "64160": "Gd160",
                        "92235": "U235", "92238": "U238", "92234": "U234", "92236": "U236"}
    for iso in mix.keys():
        iso_code = get_code_from_iso(association_dict, iso)
        print(f"{str(iso_code)}{suffix_for_temperature} {str((mix[iso]))}")


AT10_45Gd_pin = create_exact_BWR_pin("AT10_45Gd", 0.4435, 0.452, 0.514, isGd=True, height=1, pitch=1.295, path_to_composition_files="/home/loutre/RESEARCH/PolyBWR_project/PTT/compo_data/")

homogenized_45Gd = homogenize_pin(AT10_45Gd_pin, 1.295)

print(homogenized_45Gd)

print_to_Dragon5_format(homogenized_45Gd)
print_to_Serpent2_format(homogenized_45Gd, ".06c")
