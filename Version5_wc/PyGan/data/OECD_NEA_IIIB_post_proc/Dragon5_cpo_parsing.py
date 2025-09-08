# Collection of function used to parse DRAGON5 COMPO output files for fission and neutron gamma rates post treatment.
import os
from re import U
import numpy as np
import lcm
from collections import defaultdict
from utilities import sum_rates_over_iso

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


def parse_DRAGON_DIAG_rates_regi_num(name_case, name_compo, composition_option, evaluation, ssh_method, correlation_option, fission_isotopes, n_gamma_isotopes, bu):
    """
    Parse DRAGON5 rates from the specified COMPO file.
    
    Parameters:
    - name_case (str): Name of the case.
    - name_compo (str): Name of the composition file.
    - composition_option (str): Isotopic compositon option.
    - evaluation (str): Nuclear data evaluation.
    - ssh_method (str): Self shielding method used.
    - correlation_option (str): Correlation option.
    - fission_isotopes (list): List of isotopes for fission rates.
    - n_gamma_isotopes (list): List of isotopes for neutron gamma rates.
    - bu (int): Burnup step.
    - unfold_symmetry (bool): Whether to unfold symmetry in the rates.
    
    Returns:
    - keff_D5 (float): The effective multiplication factor.
    - fiss_rates (dict): Dictionary of fission rates by isotope and mix.
    - n_gamma_rates (dict): Dictionary of neutron gamma rates by isotope and mix.
    - FLUX_295groups (float): The flux on the SHEM295 energy groups.
    
    """
    # Load the DRAGON rates
    path = os.getcwd()
    # AT10_void_0_J311_295_PT_NOCORR_region_num
    os.chdir(f"PYGAN_RESULTS/{name_case}_results/{composition_option}_{evaluation}_{ssh_method}_{correlation_option}_region_num")
    print(f"Loading {name_case} rates from {name_compo}")
    # Load the LCM file
    print(os.listdir())
    pyCOMPO = lcm.new('LCM_INP', name_compo, impx=0)
    os.chdir(path)
    # Retrieve the fission rates

    len_isotot = np.shape(pyCOMPO['EDIHOM_COND']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESDENS'])[0] - 1
    print(f"len_isotot = {len_isotot}")
    ########## CALCULATIONS ##########
    # Retrieve keff from pyCOMPO
    keff_D5 = pyCOMPO['EDIHOM_COND']['MIXTURES'][0]['CALCULATIONS'][0]['K-EFFECTIVE']
    print(f"keff_D5 = {keff_D5}")
    MIXES_idx = [0,1,2,3,4,5,6,7,8,9,
                 10,11,12,13,14,15,16,17,18,
                 19,20,21,22,23,24,25,26,
                 27,28,29,30,31,32,33,
                 # 35-1, 36-1, 37-1 are W1, WB, W2
                37, 38, 39,
                 # 41-1, 42-1 are  W0, WR
                42, 43, 44,
                # 46-5 is W4
                46, 47, 48,
                49, 50, 51, 
                52, 53, 
                54]
    MIX_unique_numbers = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                        11, 12, 13, 14, 15, 16, 17, 18, 19,
                        20, 21, 22, 23, 24, 25, 26, 27,
                        28, 29, 30, 31, 32, 33, 34,
                        38, 39, 40,
                        43, 44, 45,
                        47, 48, 49,
                        50, 51, 52,
                        53, 54,
                        55]
    unique_mixes_on_diag = [1, 11, 20, 28, 50, 53, 55]
    lattice_desc = [
    ["C1_1", "C2_2", "C3_3", "C5_4", "C6_5", "C5_6", "C4_7", "C3_8", "C2_9", "C1_10"],
    ["C2_2", "C4_11", "C7_12", "C6_13", "C7_14", "C6_15", "C6_16", "C7_17", "C4_18", "C2_19"],
    ["C3_3", "C7_12", "C6_20", "C6_21", "C6_22", "C7_23", "C6_24", "C6_25", "C7_26", "C3_27"],
    ["C5_4", "C6_13", "C6_21", "C6_28", "C6_29", "C6_30", "C7_31", "C6_32", "C5_33", "C4_34"],
    ["C6_5", "C7_14", "C6_22", "C6_29", "W1", "WB", "W2", "C4_38", "C6_39", "C4_40"],
    ["C5_6", "C6_15", "C7_23", "C6_30", "WL", "W0", "WR", "C3_43", "C7_44", "C4_45"],
    ["C4_7", "C6_16", "C6_24", "C7_31", "W3", "WT", "W4", "C4_47", "C4_48", "C4_49"],
    ["C3_8", "C7_17", "C6_25", "C6_32", "C4_38", "C3_43", "C4_47", "C4_50", "C8_51", "C3_52"],
    ["C2_9", "C4_18", "C7_26", "C5_33", "C6_39", "C7_44", "C4_48", "C8_51", "C4_53", "C2_54"],
    ["C1_10", "C2_19", "C3_27", "C4_34", "C4_40", "C4_45", "C4_49", "C3_52", "C2_54", "C1_55"],
    ]
    Iso_index_to_ALIAS = {}
    fiss_rates = {}
    n_gamma_rates = {}
    for iso in range(len_isotot):
        #print(f"iso index = {iso}, isotope = {Iso_index_to_ALIAS[iso]}")
        isotope = pyCOMPO['EDIHOM_COND']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][iso]['ALIAS'][0:5].strip()
        #print(f"isotope = {isotope}")
        if isotope in fission_isotopes:
            isotope_fission_rate = {}
            for mix in MIXES_idx:
                isotope_fission_rate[f"{mix+1}"] = {}
                NWT0 = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NWT0']
                N = pyCOMPO['H_EDI_REGI_1']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESDENS'][iso]
                vol = pyCOMPO['H_EDI_REGI_1']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESVOL'][iso]
                NFTOT = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NFTOT']
                if mix+1 in unique_mixes_on_diag:
                    sym_factor = 2
                else:
                    sym_factor = 1
                isotope_fission_rate[f"{mix+1}"] = np.array(NFTOT)*np.array(NWT0)*N*vol*sym_factor # multiply volume by 2 to account for diagonal symmetry of the assembly
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

    FLUX_295groups = pyCOMPO['EDIHOM_295']['MIXTURES'][0]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NWT0']
    print(f"FLUX_295groups = {FLUX_295groups}")
    #print(f"fiss_rates = {fiss_rates}")
    sum_summed_fiss_rates = sum_rates_over_iso(fiss_rates)

    sum_summed_n_gamma_rates = sum_rates_over_iso(n_gamma_rates)

    fiss_rates["TOT"] = sum_summed_fiss_rates
    #n_gamma_rates["TOT"] = sum_summed_n_gamma_rates

    return keff_D5, fiss_rates, n_gamma_rates, FLUX_295groups, lattice_desc, MIX_unique_numbers


def parse_Spectrum_D5(pyCOMPO):
    """
    Parse the spectrum from DRAGON5 results
    """
    print(pyCOMPO['EDI_HOM']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'])
    keff_D5 = pyCOMPO['EDI_HOM']['MIXTURES'][0]['CALCULATIONS'][0]['K-EFFECTIVE']
    energy_mesh = pyCOMPO['EDI_HOM']['MIXTURES'][0]['CALCULATIONS'][0]['ENERGY']
    energy_mesh = np.array(energy_mesh)
    print(f"energy_mesh = {energy_mesh}")
    FLUX_295groups = pyCOMPO['EDI_HOM']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][0]['NWT0']

    FLUX_295groups = np.array(FLUX_295groups)
    FLUX_295groups = FLUX_295groups / np.sum(FLUX_295groups)  # Normalize the flux
    print(f"FLUX_295groups = {FLUX_295groups}")
    return keff_D5[0], energy_mesh[::-1], FLUX_295groups[::-1]  # Reverse the order to match the increasing energy instead of decreasing energy

def parse_D5_rates(pyCOMPO, bu=0):
    """
    Parse fission rates for U235 & U238 from DRAGON5 results
    """
    fission_isotopes = ["U235", "U238"]
    len_isotot = np.shape(pyCOMPO['EDI_HOM']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESDENS'])[0] - 1
    print(f"len_isotot = {len_isotot}")
    ########## CALCULATIONS ##########
    # Retrieve keff from pyCOMPO
    keff_D5 = pyCOMPO['EDI_HOM']['MIXTURES'][0]['CALCULATIONS'][0]['K-EFFECTIVE']
    print(f"keff_D5 = {keff_D5}")
    ROD_ID_IDX = [0,1,2,3,4]
    POS_ID_IDX = [1,2,3,4,5,6,7,8,9]
    POS_ON_DIAG = [1, 5, 8]

    Iso_index_to_ALIAS = {}
    fiss_rates = {}
    n_groups = len(np.array(pyCOMPO['HOM2g']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][0]['NWT0']))
    n_scores = len(POS_ID_IDX)
    print(f"n_groups = {n_groups}, n_scores = {n_scores}")
    fiss_rates_D5 = np.zeros((n_groups, n_scores))
    fiss_iso_dens = {}

    for iso in range(len_isotot):
        isotope = pyCOMPO['EDI_HOM']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][iso]['ALIAS'][0:5].strip()
        #print(f"isotope = {isotope}, iso number = {iso}")
        Iso_index_to_ALIAS[iso] = isotope
    
        
    for mix in POS_ID_IDX:
        fiss_iso_dens[mix] = {}
        for iso in range(len_isotot):
            #print(f"iso index = {iso}, isotope = {Iso_index_to_ALIAS[iso]}")
            isotope = pyCOMPO['EDI_HOM']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][iso]['ALIAS'][0:5].strip()
            if isotope in fission_isotopes:
                #print(f"isotope = {isotope}")
                NWT0 = pyCOMPO['HOM2g']['MIXTURES'][mix-1]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NWT0']
                N = pyCOMPO['HOM2g']['MIXTURES'][mix-1]['CALCULATIONS'][bu]['ISOTOPESDENS'][iso]
                print(f"mix = {mix-1}, N = {N}, isotope = {isotope}")
                
                fiss_iso_dens[mix][isotope] = N 
                vol = pyCOMPO['HOM2g']['MIXTURES'][mix-1]['CALCULATIONS'][bu]['ISOTOPESVOL'][iso]
                NFTOT = pyCOMPO['HOM2g']['MIXTURES'][mix-1]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NFTOT']
                print(f"NFTOT[0] = {NFTOT[0]}, NWT0[0] = {NWT0[0]}, N = {N}, vol = {vol}")
                print(f"NFTOT[1] = {NFTOT[1]}, NWT0[1] = {NWT0[1]}")
                if mix in POS_ON_DIAG:
                    sym_factor = 2
                else:
                    sym_factor = 1
                #fiss_rates_D5[0, mix-1] += NFTOT[0]*NWT0[0]*N*vol*sym_factor  # multiply volume by 2 to account for diagonal symmetry of the assembly
                
                fiss_rates_D5[0, mix-1] = NFTOT[0]*NWT0[0]*N*vol*sym_factor
                fiss_rates_D5[1, mix-1] += NFTOT[1]*NWT0[1]*N*vol*sym_factor
    #            fiss_rates[isotope] = U235_fiss_rates_D5

    return fiss_rates_D5, fiss_iso_dens