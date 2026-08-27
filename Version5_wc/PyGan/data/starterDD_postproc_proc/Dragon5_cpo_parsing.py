# Collection of function used to parse DRAGON5 COMPO output files for fission and neutron gamma rates post treatment.
import os
import numpy as np
import lcm
from collections import defaultdict
import pandas as pd


def parse_COMPO(COMPO_name, assembly_model, fission_isotopes, path_to_compo=None):
    """
    Parse DRAGON5 rates from the specified COMPO file.
    
    Parameters:
    - name_case (str): Name of the case.
    - name_compo (str): Name of the composition file.
    - composition_option (str): Isotopic compositon option.
    - evaluation (str): Nuclear data evaluation.
    - ssh_method (str): Self shielding method used.
    - correlation_option (str): Correlation option.
    - geometry_refinement_option (str): identifictor for flux geometry refinement type.
    - fission_isotopes (list): List of isotopes for fission rates.
    - n_gamma_isotopes (list): List of isotopes for neutron gamma rates.
    - bu (int): Burnup step.
    - computational_scheme (str): Computational scheme used, default is "1L_MOC".
    
    Returns:
    - keff_D5 (float): The effective multiplication factor.
    - n_gamma_rates (dict): Dictionary of neutron gamma rates by isotope and mix.
    - FLUX_295groups (float): The flux on the SHEM295 energy groups.
    
    """
    bu=0 # Burnup step to extract from COMPO, should be 0 for initial condition
    # Load the DRAGON rates
    path = os.getcwd()
    if path_to_compo is not None:
        compo_path = path_to_compo
        os.chdir(compo_path)
    else:
        os.chdir(f"DRAGON_RESULTS")
    print(f"Loading {COMPO_name} rates from {COMPO_name}")
    print(f"Reading from directory : {os.getcwd()}")
    pyCOMPO = lcm.new('LCM_INP', COMPO_name, impx=0)
    os.chdir(path)

    # Get lattice info from assembly model
    lattice_info = assembly_model.get_postprocessing_lattice_info()
    ordered_pin_indices = lattice_info['ordered_pin_indices']
    pin_idx_on_axis = lattice_info['pin_idx_on_symmetry_axis']
    n_unique_pins = lattice_info['n_unique_pins']
    # For DRAGON, mixture indices are 1-based and should match pin_idx+1 if ordered
    MIXES_idx = [pidx-1 for pidx in ordered_pin_indices]
    unique_mixes_on_diag = [pidx-1 for pidx in pin_idx_on_axis]
    name_2g = "EDIR_2G" #  # "H_EDI_REGI_2"
    len_isotot = np.shape(pyCOMPO[name_2g]['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESDENS'])[0] - 1
    #print(f"len_isotot = {len_isotot}")
    ########## CALCULATIONS ##########
    # Retrieve keff from pyCOMPO
    keff_D5 = pyCOMPO[name_2g]['MIXTURES'][0]['CALCULATIONS'][0]['K-EFFECTIVE'][0]
    #print(f"keff_D5 = {keff_D5}")
    isotopes = ["U235", "U238", "Gd155", "Gd157"]
    
    fission_rates_1g = np.zeros(n_unique_pins)
    prod_over_abs_1g = np.zeros(n_unique_pins)
    
    fission_rates_2g = np.zeros((2, n_unique_pins))
    prod_over_abs_2g = np.zeros((2, n_unique_pins))
    
    n_gamma_rates = {}
    for iso in range(len_isotot):
        isotope = pyCOMPO[name_2g]['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][iso]['ALIAS'][0:5].strip()
        if isotope in isotopes:
            for idx, mix in enumerate(MIXES_idx):
                
                try:
                    NWT0_1g = pyCOMPO['H_EDI_REGI_1']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NWT0']
                except lcm.PyLcmError:
                    NWT0_1g = np.zeros(1)
                print(f"mix = {mix+1}, isotope = {isotope}, iso = {iso}")
                NWT0_2g = pyCOMPO[name_2g]['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NWT0']
                N = pyCOMPO[name_2g]['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESDENS'][iso]
                vol = pyCOMPO[name_2g]['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESVOL'][iso]
                try:
                    NFTOT_1g = pyCOMPO['H_EDI_REGI_1']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NFTOT']
                except lcm.PyLcmError:
                    NFTOT_1g = np.zeros(len(NWT0_1g))
                try:
                    NFTOT_2g = pyCOMPO[name_2g]['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NFTOT']
                except lcm.PyLcmError:
                    NFTOT_2g = np.zeros(len(NWT0_2g))
                    
                #print(f"mix = {mix+1}, NWT0 = {NWT0}, N = {N}, vol = {vol}, NFTOT = {NFTOT}")
                # absorption = sigma_g - sigma_sigs_g for multigroup keff estimate ?
                # = NTOT0 - SIGSS0
                try:
                    TOT_1g = pyCOMPO['H_EDI_REGI_1']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NTOT0']
                except lcm.PyLcmError:
                    TOT_1g = np.zeros(len(NWT0_1g))

                try:
                    TOT_2g = pyCOMPO[name_2g]['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NTOT0']
                except lcm.PyLcmError:
                    TOT_2g = np.zeros(len(NWT0_2g))

                try:
                    SIGS0_1g = pyCOMPO['H_EDI_REGI_1']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['SIGS00']
                except lcm.PyLcmError:
                    SIGS0_1g = np.zeros(len(NWT0_1g))

                try:
                    SIGS0_2g = pyCOMPO[name_2g]['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['SIGS00']
                except lcm.PyLcmError:
                    SIGS0_2g = np.zeros(len(NWT0_2g))

                ABS_1g = np.array(TOT_1g) - np.array(SIGS0_1g)
                ABS_2g = np.array(TOT_2g) - np.array(SIGS0_2g)

                try:
                    PROD_1g = pyCOMPO['H_EDI_REGI_1']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NUSIGF']
                except lcm.PyLcmError:
                    PROD_1g = np.zeros(len(NWT0_1g))
                    
                try:
                    PROD_2g = pyCOMPO[name_2g]['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NUSIGF']
                except lcm.PyLcmError:
                    PROD_2g = np.zeros(len(NWT0_2g))
                ### Recover different contributions to absorption to check consistency with TOT-SIGS0 definition of absorption and understand the importance of different reactions in the absorption rate.
                try:
                    NP_1g = pyCOMPO['H_EDI_REGI_1']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NP']
                except lcm.PyLcmError:
                    NP_1g = np.zeros(len(NWT0_1g))
                
                try:
                    NP_2g = pyCOMPO[name_2g]['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NP']
                except lcm.PyLcmError:
                    NP_2g = np.zeros(len(NWT0_2g))
                
                try: 
                    NG_1g = pyCOMPO['H_EDI_REGI_1']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NG']
                except lcm.PyLcmError:
                    NG_1g = np.zeros(len(NWT0_1g))
                try:
                    NG_2g = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NG']
                except lcm.PyLcmError:
                    NG_2g = np.zeros(len(NWT0_2g))
                try: 
                    ND_1g = pyCOMPO['H_EDI_REGI_1']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['ND']
                    ND_2g = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['ND']
                except lcm.PyLcmError:
                    ND_1g = np.zeros(len(NWT0_1g))
                    ND_2g = np.zeros(len(NWT0_2g))
                try: 
                    NT_1g = pyCOMPO['H_EDI_REGI_1']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NT']
                    NT_2g = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NT']
                except lcm.PyLcmError:
                    NT_1g = np.zeros(len(NWT0_1g))
                    NT_2g = np.zeros(len(NWT0_2g))
                try: 
                    NA_1g = pyCOMPO['H_EDI_REGI_1']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NA']
                    NA_2g = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NA']
                except lcm.PyLcmError:
                    NA_1g = np.zeros(len(NWT0_1g))
                    NA_2g = np.zeros(len(NWT0_2g))
                try: 
                    N2A_1g = pyCOMPO['H_EDI_REGI_1']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['N2A']
                    N2A_2g = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['N2A']
                except lcm.PyLcmError:
                    N2A_1g = np.zeros(len(NWT0_1g))
                    N2A_2g = np.zeros(len(NWT0_2g))
                try: 
                    NNP_1g = pyCOMPO['H_EDI_REGI_1']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NNP']
                    NNP_2g = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NNP']
                except lcm.PyLcmError:
                    NNP_1g = np.zeros(len(NWT0_1g))
                    NNP_2g = np.zeros(len(NWT0_2g))
                try:
                    N2N_1g = pyCOMPO['H_EDI_REGI_1']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['N2N']
                    N2N_2g = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['N2N']
                except lcm.PyLcmError:
                    N2N_1g = np.zeros(len(NWT0_1g))
                    N2N_2g = np.zeros(len(NWT0_2g))
                try:
                    N3N_1g = pyCOMPO['H_EDI_REGI_1']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['N3N']
                    N3N_2g = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['N3N']
                except lcm.PyLcmError:
                    N3N_1g = np.zeros(len(NWT0_1g))
                    N3N_2g = np.zeros(len(NWT0_2g))
                try:
                    N4N_1g = pyCOMPO['H_EDI_REGI_1']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['N4N']
                    N4N_2g = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['N4N']
                except lcm.PyLcmError:
                    N4N_1g = np.zeros(len(NWT0_1g))
                    N4N_2g = np.zeros(len(NWT0_2g))

                # Use DRAGON definitition of neutronic absorption
                # ABS = NFTOT + NG + NP + ND + NT + NA + N2A + NNP - N2N - 2*N3N - 3*N4N
                ABS_1g = np.array(NFTOT_1g) + np.array(NG_1g) + np.array(NP_1g) + np.array(NA_1g) - np.array(N2N_1g) - 2*np.array(N3N_1g) #- 3*N4N
                ABS_2g = np.array(NFTOT_2g) + np.array(NG_2g) + np.array(NP_2g) + np.array(NA_2g) - np.array(N2N_2g) - 2*np.array(N3N_2g) #- 3*N4N
                ABS_SIGS_1g = np.array(TOT_1g) - np.array(SIGS0_1g)
                ABS_SIGS_2g = np.array(TOT_2g) - np.array(SIGS0_2g)
                print(f"ABS from TOT-SIGS0 = {ABS_SIGS_1g}, ABS from sum over reactions = {ABS_1g}, difference = {np.array(ABS_1g)-np.array(ABS_SIGS_1g)}")
                print(f"ABS from TOT-SIGS0 = {ABS_SIGS_2g}, ABS from sum over reactions = {ABS_2g}, difference = {np.array(ABS_2g)-np.array(ABS_SIGS_2g)}")
                    
                if mix in unique_mixes_on_diag:
                    sym_factor = 2
                else:
                    sym_factor = 1
                
                fission_rates_1g[mix] += NFTOT_1g[0]*NWT0_1g[0]*N*vol*sym_factor
                #prod_over_abs_1g[mix] += (np.array(PROD_1g[0])/np.array(ABS_1g[0]))# *N*vol
                
                fission_rates_2g[0][mix] += NFTOT_2g[1]*NWT0_2g[1]*N*vol*sym_factor
                fission_rates_2g[1][mix] += NFTOT_2g[0]*NWT0_2g[0]*N*vol*sym_factor
                
                prod_over_abs_2g[0][mix] += (np.array(NFTOT_2g[1])/np.array(ABS_2g[1]))# *N*vol
                prod_over_abs_2g[1][mix] += (np.array(NFTOT_2g[0])/np.array(ABS_2g[0])) # *N*vol
    # Recover 295g absorption rates for U238
    FLUX_spectrum_ngroups = np.array(pyCOMPO['U238_295']['MIXTURES'][0]['CALCULATIONS'][bu]['ISOTOPESLIST'][0]['NWT0'])
    U238_TOT = pyCOMPO['U238_295']['MIXTURES'][0]['CALCULATIONS'][bu]['ISOTOPESLIST'][0]['NTOT0']
    U238_SIGS0 = pyCOMPO['U238_295']['MIXTURES'][0]['CALCULATIONS'][bu]['ISOTOPESLIST'][0]['SIGS00']
    U238_ABS_SIGS = (np.array(U238_TOT) - np.array(U238_SIGS0))*np.array(FLUX_spectrum_ngroups)
    try:
        NF_295 = np.array(pyCOMPO['U238_295']['MIXTURES'][0]['CALCULATIONS'][bu]['ISOTOPESLIST'][0]['NFTOT'])
    except lcm.PyLcmError:
        NF_295 = np.zeros(len(FLUX_spectrum_ngroups))
    try:
        NP_295 = np.array(pyCOMPO['U238_295']['MIXTURES'][0]['CALCULATIONS'][bu]['ISOTOPESLIST'][0]['NP'])
    except lcm.PyLcmError:
        NP_295 = np.zeros(len(FLUX_spectrum_ngroups))
    try:
        NG_295 = np.array(pyCOMPO['U238_295']['MIXTURES'][0]['CALCULATIONS'][bu]['ISOTOPESLIST'][0]['NG'])
    except lcm.PyLcmError:
        NG_295 = np.zeros(len(FLUX_spectrum_ngroups))
    try:
        ND_295 = np.array(pyCOMPO['U238_295']['MIXTURES'][0]['CALCULATIONS'][bu]['ISOTOPESLIST'][0]['ND'])
    except lcm.PyLcmError:
        ND_295 = np.zeros(len(FLUX_spectrum_ngroups))
    try:
        NT_295 = np.array(pyCOMPO['U238_295']['MIXTURES'][0]['CALCULATIONS'][bu]['ISOTOPESLIST'][0]['NT'])
    except lcm.PyLcmError:
        NT_295 = np.zeros(len(FLUX_spectrum_ngroups))
    try:
        NA_295 = np.array(pyCOMPO['U238_295']['MIXTURES'][0]['CALCULATIONS'][bu]['ISOTOPESLIST'][0]['NA'])
    except lcm.PyLcmError:
        NA_295 = np.zeros(len(FLUX_spectrum_ngroups))
    try:
        N2A_295 = np.array(pyCOMPO['U238_295']['MIXTURES'][0]['CALCULATIONS'][bu]['ISOTOPESLIST'][0]['N2A'])
    except lcm.PyLcmError:
        N2A_295 = np.zeros(len(FLUX_spectrum_ngroups))
    try:
        NNP_295 = np.array(pyCOMPO['U238_295']['MIXTURES'][0]['CALCULATIONS'][bu]['ISOTOPESLIST'][0]['NNP'])
    except lcm.PyLcmError:
        NNP_295 = np.zeros(len(FLUX_spectrum_ngroups))
    try:
        N2N_295 = np.array(pyCOMPO['U238_295']['MIXTURES'][0]['CALCULATIONS'][bu]['ISOTOPESLIST'][0]['N2N'])
    except lcm.PyLcmError:
        N2N_295 = np.zeros(len(FLUX_spectrum_ngroups))
    try:
        N3N_295 = np.array(pyCOMPO['U238_295']['MIXTURES'][0]['CALCULATIONS'][bu]['ISOTOPESLIST'][0]['N3N'])
    except lcm.PyLcmError:
        N3N_295 = np.zeros(len(FLUX_spectrum_ngroups))
    try:
        N4N_295 = np.array(pyCOMPO['U238_295']['MIXTURES'][0]['CALCULATIONS'][bu]['ISOTOPESLIST'][0]['N4N'])
    except lcm.PyLcmError:
        N4N_295 = np.zeros(len(FLUX_spectrum_ngroups))
    
    # Use DRAGON definitition of neutronic absorption
    #U238_ABS = (NF_295 + NP_295 + NG_295 + ND_295 + NT_295 + NA_295 + N2A_295 + NNP_295 - N2N_295 - 2*N3N_295 - 3*N4N_295)*np.array(FLUX_spectrum_ngroups)
    #U238_ABS = (NF_295 + NP_295 + NG_295 + ND_295 + NT_295 + NA_295 + N2A_295 + NNP_295 - N2N_295 - 2*N3N_295)*np.array(FLUX_spectrum_ngroups)
    U238_ABS = (NF_295 + NG_295 + NP_295 + NA_295 - N2N_295 - 2*N3N_295)*np.array(FLUX_spectrum_ngroups)
    #print(f"U238 absorption from TOT-SIGS0 = {U238_ABS_SIGS}, U238 absorption from sum over reactions = {U238_ABS}, difference = {U238_ABS - U238_ABS_SIGS}")

    FLUX_spectrum_ngroups = np.array(pyCOMPO['EDIHOM_295']['MIXTURES'][0]['CALCULATIONS'][bu]['ISOTOPESLIST'][0]['NWT0'])
    energy_mesh = np.array(pyCOMPO['EDIHOM_295']['MIXTURES'][0]['CALCULATIONS'][0]['ENERGY'])
    #print(f"D5 fiss over abs = {fiss_over_abs}")
    
    return keff_D5, fission_rates_1g, prod_over_abs_1g, fission_rates_2g, prod_over_abs_2g, U238_ABS[::-1], FLUX_spectrum_ngroups[::-1], energy_mesh[::-1]