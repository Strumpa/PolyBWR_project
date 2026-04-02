## Python3 file to group Serpent2 parsing functions
import os
import numpy as np
import serpentTools as st
import re


def parse_S2_pin_mat_det(assembly_model, assembly_id, XS_lib_S2, fission_isotopes, bu=0):
    """
    Serpent2 assembly detector post-treatment for Assembly reaction rates :
    use lattice detector definition
    assembly_model : (CartesianAssemblyModel) assembly model with geometry/material info
    XS_lib_S2 : (str) name of the Serpent2 XS library
    fission_isotopes : (list) list of fission isotopes
    ngamma_isotopes : (list) list of (n,gamma) isotopes
    bu : (int) burnup step
    returns : reaction rates for fission and (n,gamma) reactions for the specified isotopes
    """
    # Get lattice info from assembly model
    lattice_info = assembly_model.get_postprocessing_lattice_info()
    ordered_pin_indices = lattice_info['ordered_pin_indices']
    pin_idx_to_material = lattice_info['pin_idx_to_material_name']
    #print(f"pin_idx_to_material = {pin_idx_to_material}")
    pin_idx_to_composition = lattice_info['pin_idx_to_composition']
    #print(f"pin_idx_to_composition = {pin_idx_to_composition}")
    pin_idx_on_axis = lattice_info['pin_idx_on_symmetry_axis']
    n_unique_pins = lattice_info['n_unique_pins']

    MIXES_idx = [pidx for pidx in ordered_pin_indices]
    unique_mixes_on_diag = [pidx for pidx in pin_idx_on_axis]

    if assembly_id == "GE14_DOM":
        s2_case_name = f"GE14_DOM_assembly_{XS_lib_S2}_00.serp"
    elif assembly_id == "AT10":
        s2_case_name = f"AT10_assembly_starterDD_{XS_lib_S2}.serp"
    elif assembly_id == "GE14_DOM-C":
        s2_case_name = f"GE14_DOM-C_assembly_{XS_lib_S2}_00.serp"
    elif assembly_id == "AT10-C":
        s2_case_name = f"AT10_assembly_CTRL_00_{XS_lib_S2}.serp"
        
    detectorFile = st.read(f"{os.environ['SERPENT_RESULTS']}/{s2_case_name}_det{bu}.m")
    resultsFile = st.read(f"{os.environ['SERPENT_RESULTS']}/{s2_case_name}_res.m")
    keff = resultsFile.resdata["absKeff"].T[0]
    
    # try reading depletion 
    try:
        depletionFile = st.read(f"{os.environ['SERPENT_RESULTS']}/{assembly_id}_{XS_lib_S2}_BU_inp_dep.m")
    except:
        depletionFile = None


    #print(f"keff = {keff}")
    #print(detectorFile.detectors.keys())

    # Extracting the detector response
    if assembly_id == "GE14_DOM" or assembly_id == "GE14_DOM-C":
        ngroups = detectorFile.detectors["det_UOX16_51"].tallies.shape[0] # 2 = n_groups
    elif assembly_id == "AT10" or assembly_id == "AT10-C":
        ngroups = detectorFile.detectors["det_24UOX_1"].tallies.shape[0] # 2 = n_groups
    vol = np.pi * assembly_model.pin_geometry_dict.get("fuel_radius", None) ** 2
    if XS_lib_S2 == "endfb8r1":

        tally_index_to_react_pin_detectors = {
                                0: ("disappearance","U235"), 
                                1: ("disappearance","U238"), 
                                2: ("disappearance","Gd155"), 
                                3: ("disappearance","Gd157"), 
                                4: ("fission","U235"), 
                                5: ("fission","U238"),
                                6: ("n,gamma","U235"), 
                                7: ("n,gamma","U238"),  
                                8: ("n,gamma","Gd155"), 
                                9: ("n,gamma","Gd157"),
                                10: ("n,proton","U238"),
                                11: ("n,alpha","U235"),
                                12: ("n,alpha","U238"),
                                13: ('n,2n', "U235"),
                                14: ('n,2n', "U238"),
                                15: ('n,3n', "U235"),
                                16: ('n,3n', "U238"),
                                }
     
        tally_to_reaction_assembly = {
                                0: ("disappearance","U238"), 
                                1: ("fission","U238"),
                                2: ("n,gamma","U238"),  
                                3: ("n,proton","U238"),
                                4: ("n,alpha","U238"),
                                5: ('n,2n', "U238"),
                                6: ('n,3n', "U238"),
                            }
        
    npos = n_unique_pins
    disappearance = np.zeros((ngroups, npos))
    fission_rates = np.zeros((ngroups, npos))
    ngamma_rates = np.zeros((ngroups, npos))
    nproton_rates = np.zeros((ngroups, npos))
    nalpha_rates = np.zeros((ngroups, npos))
    n2n_rates = np.zeros((ngroups, npos))
    n3n_rates = np.zeros((ngroups, npos))
    n_isotopes_per_mat = {}
    groups = {}
    pattern = re.compile(r"^(\d{2}(?:UOX|Gd))_([A-Za-z])_(\d+)$")
    # Identify rings in groups of identical 'initial' materials ie same pin position, same initial enrichment.
    if depletionFile is not None:
        for id_nb, material_name in enumerate(depletionFile):
            #print(f"Depletion File Material name {material_name}")
            m = pattern.match(material_name)
            n_iso = 0
            if m:
                material = m.group(1)   # ie material code example : "24UOX"
                ring_number = m.group(2)   # ie "A", "B", "C" or "D" for UOX and + "E", "F" for Gd
                pos  = int(m.group(3))  # e.g. 1
                key = (material, pos)
                groups.setdefault(key, []).append(material_name)
    else:
        # NEED TO RECOVER COMPOSITION OF EACH MATERIAL MIXTURE DEFINED FROM THE ASSEMBLY MODEL, AND THEN GROUP THE MATERIAL MIXTURES BY PIN POSITION AND INITIAL ENRICHMENT TO RECOVER THE SAME GROUPS AS IN THE DEPLETION FILE
        print("Using the compositions recovered from the starterDD assembly models")
    
    # NEED TO EDIT THIS SO THAT IN THIS CASE USE THE COMPOSITIONS FROM THE ASSEMBLY MODEL
    if depletionFile is not None:
        for iso in fission_isotopes:
            n_iso = 0
            for key, material_names in groups.items():
                #print(f"Processing isotope {iso} for material {key}")
                material, pos = key
                n_iso = 0
                #print(f"Material {material}, Position {pos}, Material names : {material_names}")
                for material_name in material_names:
                    # Example of material_name : "24UOX_A_1" or "45Gd_D_31"
                    if key not in n_isotopes_per_mat.keys():
                        n_isotopes_per_mat[key] = {}
                    n_iso_ring = depletionFile[material_name]['adens'][iso_code_to_num[iso]][bu]*depletionFile[material_name]['volume'][bu]
                    #print(f"Material {material_name}, Isotope {iso}, Density {depletionFile[material_name]['adens'][iso_code_to_num[iso]][bu]}, Volume {depletionFile[material_name]['volume'][bu]}")
                    n_iso += n_iso_ring
                #print(f"Total number density of isotope {iso} in material {material} at position {pos} : {n_iso/vol} atoms/b-cm")
                #print(f"Comparing to pin definition : {N_iso[material][iso]} atoms/b-cm")
                n_isotopes_per_mat[key][iso] = n_iso/vol
                #print(f"Material {material_name}, Isotope {iso}, Density {n_iso}, pos {pos}")
            #print(n_isotopes_per_mat)
    
    # extract fision rates
    isotopes = ["U235", "U238", "Gd155", "Gd157"]
    for pos_idx in MIXES_idx:
        material = pin_idx_to_material[pos_idx]
        detector_name = f"det_{material}_{pos_idx}"
        n_groups = detectorFile.detectors[detector_name].tallies.shape[0]
        n_reactions = detectorFile.detectors[detector_name].tallies.shape[1]
        for g in range(n_groups):
            for r in range(n_reactions):
                rate = detectorFile.detectors[detector_name].tallies[g, r]
                reaction,isotope = tally_index_to_react_pin_detectors[r]
                if isotope in isotopes:
                    ### Sanity check : compare N_iso from pin definition and from depletion file
                    symmetry_factor = 1 if pos_idx in unique_mixes_on_diag else 2
                    try:
                        N_iso_pin = pin_idx_to_composition[pos_idx][isotope]
                    except KeyError:
                        print(f"Warning : isotope {isotope} not found in pin composition for material {material} at position {pos_idx}. Setting N_iso_pin to 0.")
                        N_iso_pin = 0.0
                        
                    if depletionFile is not None:
                        N_iso_dep = n_isotopes_per_mat[key][isotope] / symmetry_factor
                        if abs(N_iso_pin - N_iso_dep)/N_iso_pin > 0.1:
                            print(f"Warning : large discrepancy between N_iso from pin definition and from depletion file for isotope {isotope} in material {material} at position {pos_idx} : N_iso_pin = {N_iso_pin}, N_iso_dep = {N_iso_dep}")
                        #print(f"Material {material}, Position {pos_idx}, Isotope {isotope}, Rate {rate}, N_iso_pin = {N_iso_pin}, N_iso_dep = {N_iso_dep}")
                        N_iso = N_iso_dep
                    else: 
                        N_iso = N_iso_pin
                    ### Now fill the rates arrays
                    if reaction == "fission":
                        fission_rates[g, pos_idx-1] += rate / symmetry_factor * N_iso 
                    elif reaction == "n,gamma":
                        ngamma_rates[g, pos_idx-1] += rate / symmetry_factor * N_iso 
                    elif reaction == "n,proton":
                        nproton_rates[g, pos_idx-1] += rate / symmetry_factor * N_iso
                    elif reaction == "n,alpha":
                        nalpha_rates[g, pos_idx-1] += rate / symmetry_factor * N_iso
                    elif reaction == "n,2n":
                        n2n_rates[g, pos_idx-1] += rate / symmetry_factor * N_iso
                    elif reaction == "n,3n":
                        n3n_rates[g, pos_idx-1] += rate / symmetry_factor * N_iso
                    elif reaction == "disappearance":
                        disappearance[g, pos_idx-1] += rate / symmetry_factor * N_iso
                #print(f"Group {g+1}, Cell {pos_idx}, Reaction {reaction}, Rate: {rate}, N_iso = {N_iso}")
    #print(n_isotopes_per_mat)
    
    # recover 295g / 26g flux spectra
    flux_295g = detectorFile.detectors["flux_295g"].tallies
    flux_26g = detectorFile.detectors["flux_26g"].tallies
    print(f"Flux 295g : {flux_295g}")
    print(f"Flux 26g : {flux_26g}")
    fluxes = {"295g": flux_295g, "26g": flux_26g}
    
    # Recover U238 reaction rates on the dedicated 295g detector
    disappearance_U238_295g = np.zeros(295)
    fission_U238_295g = np.zeros(295)
    ngamma_U238_295g = np.zeros(295)
    nproton_U238_295g = np.zeros(295)
    nalpha_U238_295g = np.zeros(295)
    n2n_U238_295g = np.zeros(295)
    n3n_U238_295g = np.zeros(295)
    # recover absorption rate of U238 from the dedicated absorption detector
    n_reactions = detectorFile.detectors["det_assembly_295g"].tallies.shape[1]
    for r in range(n_reactions):
        rate = detectorFile.detectors["det_assembly_295g"].tallies[:, r]
        reaction,isotope = tally_to_reaction_assembly[r]
        if reaction == "n,gamma" and isotope == "U238":
            ngamma_U238_295g = detectorFile.detectors["det_assembly_295g"].tallies[:, r]
            #print(f"ngamma_U238_295g : {ngamma_U238_295g}")
        elif reaction == "fission" and isotope == "U238":
            fission_U238_295g = detectorFile.detectors["det_assembly_295g"].tallies[:, r]
            #print(f"fission_U238_295g : {fission_U238_295g}")
        elif reaction == "disappearance" and isotope == "U238":
            disappearance_U238_295g = detectorFile.detectors["det_assembly_295g"].tallies[:, r]
            #print(f"disappearance_U238_295g : {disappearance_U238_295g}")
        elif reaction == "n,proton" and isotope == "U238":
            nproton_U238_295g = detectorFile.detectors["det_assembly_295g"].tallies[:, r]
            #print(f"nproton_U238_295g : {nproton_U238_295g}")
        elif reaction == "n,alpha" and isotope == "U238":
            nalpha_U238_295g = detectorFile.detectors["det_assembly_295g"].tallies[:, r]
            #print(f"nalpha_U238_295g : {nalpha_U238_295g}")
        elif reaction == "n,2n" and isotope == "U238":
            n2n_U238_295g = detectorFile.detectors["det_assembly_295g"].tallies[:, r]
            #print(f"n2n_U238_295g : {n2n_U238_295g}")
        elif reaction == "n,3n" and isotope == "U238":
            n3n_U238_295g = detectorFile.detectors["det_assembly_295g"].tallies[:, r]
            #print(f"n3n_U238_295g : {n3n_U238_295g}")

    # Recover U238 reaction rates on the dedicated 296g detector
    disappearance_U238_26g = np.zeros(26)
    fission_U238_26g = np.zeros(26)
    ngamma_U238_26g = np.zeros(26)
    nproton_U238_26g = np.zeros(26)
    nalpha_U238_26g = np.zeros(26)
    n2n_U238_26g = np.zeros(26)
    n3n_U238_26g = np.zeros(26)
    # recover absorption rate of U238 from the dedicated absorption detector
    try:
        n_reactions = detectorFile.detectors["det_assembly_26g"].tallies.shape[1]
        for r in range(n_reactions):
            rate = detectorFile.detectors["det_assembly_26g"].tallies[:, r]
            reaction,isotope = tally_to_reaction_assembly[r]
            if reaction == "n,gamma" and isotope == "U238":
                ngamma_U238_26g = detectorFile.detectors["det_assembly_26g"].tallies[:, r]
                #print(f"ngamma_U238_26g : {ngamma_U238_26g}")
            elif reaction == "fission" and isotope == "U238":
                fission_U238_26g = detectorFile.detectors["det_assembly_26g"].tallies[:, r]
                #print(f"fission_U238_26g : {fission_U238_26g}")
            elif reaction == "disappearance" and isotope == "U238":
                disappearance_U238_26g = detectorFile.detectors["det_assembly_26g"].tallies[:, r]
                #print(f"disappearance_U238_26g : {disappearance_U238_26g}")
            elif reaction == "n,proton" and isotope == "U238":
                nproton_U238_26g = detectorFile.detectors["det_assembly_26g"].tallies[:, r]
                #print(f"nproton_U238_26g : {nproton_U238_26g}")
            elif reaction == "n,alpha" and isotope == "U238":
                nalpha_U238_26g = detectorFile.detectors["det_assembly_26g"].tallies[:, r]
                #print(f"nalpha_U238_26g : {nalpha_U238_26g}")
            elif reaction == "n,2n" and isotope == "U238":
                n2n_U238_26g = detectorFile.detectors["det_assembly_26g"].tallies[:, r]
                #print(f"n2n_U238_26g : {n2n_U238_26g}")
            elif reaction == "n,3n" and isotope == "U238":
                n3n_U238_26g = detectorFile.detectors["det_assembly_26g"].tallies[:, r]
                #print(f"n3n_U238_26g : {n3n_U238_26g}")
    except KeyError:
        print("Warning : det_assembly_26g not found in Serpent2 output. U238 absorption rates on 26g detector will be set to 0.")
        print(f"For file {s2_case_name}_det{bu}.m, available detectors are : {detectorFile.detectors.keys()}")
    
    
    # use Dragon's definitions of neutronic absorption : 
    # abs = fission + (n,gamma) + (n,proton) + (n,deutron) + (n,triton) + (n,alpha) + (n,2alpha) + (n,np) - (n,2n) - 2*(n,3n) - 3*(n,4n)
    print(f"Fission rates therm : {fission_rates[0]}")
    print(f"(n,gamma) rates therm : {ngamma_rates[0]}")
    print(f"(n,proton) rates therm : {nproton_rates[0]}")
    print(f"(n,alpha) rates therm : {nalpha_rates[0]}")
    print(f"(n,2n) rates therm : {n2n_rates[0]}")
    print(f"(n,3n) rates therm : {n3n_rates[0]}")
    total_absorptions = fission_rates + ngamma_rates + nproton_rates + nalpha_rates - n2n_rates -2*n3n_rates
    total_abs_U238_26g = fission_U238_26g + ngamma_U238_26g + nproton_U238_26g + nalpha_U238_26g - n2n_U238_26g - 2*n3n_U238_26g
    total_abs_U238_295g = fission_U238_295g + ngamma_U238_295g + nproton_U238_295g + nalpha_U238_295g - n2n_U238_295g - 2*n3n_U238_295g
    
    U238_abs = {"26g": total_abs_U238_26g, "295g": total_abs_U238_295g}    
    
    ### Compute fission over absoption on 2G detectors
    fiss_over_abs = np.zeros((ngroups, npos))
    for g in range(ngroups):
        for pos_idx in MIXES_idx:
            if total_absorptions[g, pos_idx-1] > 0:
                fiss_over_abs[g, pos_idx-1] = fission_rates[g, pos_idx-1] / total_absorptions[g, pos_idx-1]
            else:
                fiss_over_abs[g, pos_idx-1] = 0.0
    
    print(f"Fission rates : {fission_rates}")
    print(f"(n,gamma) rates : {ngamma_rates}")
    print(f"Fission over absorption ratio : {fiss_over_abs}")
    return keff, fission_rates, ngamma_rates, fiss_over_abs, fluxes, U238_abs


def sum_S2rates_over_iso(rates_dict):
    """
    rates_dict is a nested dictionary with the structure:
    {
        "mix1": {
            "iso1": [val1, val2, ...],
            "iso2": [...],
            ...
        },
        "mix2": {
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
    result = {}

    for mix, iso_dict in rates_dict.items():
        sum_vector = None
        for values in iso_dict.values():
            if sum_vector is None:
                sum_vector = [0.0] * len(values)
            for i, val in enumerate(values):
                sum_vector[i] += val
        result[mix] = sum_vector

    return result


