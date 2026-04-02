# Collection of function used to parse DRAGON5 COMPO output files for fission and neutron gamma rates post treatment.
import os
import numpy as np
import lcm
from collections import defaultdict
import pandas as pd


def parse_COMPO(COMPO_name, assembly_model, fission_isotopes):
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

    len_isotot = np.shape(pyCOMPO['EDIHOM_COND']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESDENS'])[0] - 1
    #print(f"len_isotot = {len_isotot}")
    ########## CALCULATIONS ##########
    # Retrieve keff from pyCOMPO
    keff_D5 = pyCOMPO['EDIHOM_COND']['MIXTURES'][0]['CALCULATIONS'][0]['K-EFFECTIVE'][0]
    #print(f"keff_D5 = {keff_D5}")
    isotopes = ["U235", "U238", "Gd155", "Gd157"]
    fiss_over_abs = np.zeros((2, n_unique_pins))
    fission_rates = np.zeros((2, n_unique_pins))
    n_gamma_rates = {}
    for iso in range(len_isotot):
        isotope = pyCOMPO['EDIHOM_COND']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][iso]['ALIAS'][0:5].strip()
        if isotope in isotopes:
            for idx, mix in enumerate(MIXES_idx):
                NWT0 = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NWT0']
                N = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESDENS'][iso]
                vol = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESVOL'][iso]
                try:
                    NFTOT = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NFTOT']
                except lcm.PyLcmError:
                    NFTOT = np.zeros(len(NWT0))
                #print(f"mix = {mix+1}, NWT0 = {NWT0}, N = {N}, vol = {vol}, NFTOT = {NFTOT}")
                # absorption = sigma_g - sigma_sigs_g for multigroup keff estimate ?
                # = NTOT0 - SIGSS0
                TOT = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NTOT0']
                #print(f"mix = {mix+1}, TOT = {TOT}")
                SIGS0 = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['SIGS00']
                #print(f"mix = {mix+1}, SIGS0 = {SIGS0}")
                ABS = np.array(TOT) - np.array(SIGS0)
                #print(f"mix = {mix+1}, ABS = {ABS}")
                try:
                    PROD = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NUSIGF']
                except lcm.PyLcmError:
                    PROD = np.zeros(len(NWT0))
                ### Recover different contributions to absorption to check consistency with TOT-SIGS0 definition of absorption and understand the importance of different reactions in the absorption rate.
                try:
                    NP = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NP']
                except lcm.PyLcmError:
                    NP = np.zeros(len(NWT0))
                try: 
                    NG = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NG']
                except lcm.PyLcmError:
                    NG = np.zeros(len(NWT0))
                try: 
                    ND = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['ND']
                except lcm.PyLcmError:
                    ND = np.zeros(len(NWT0))
                try: 
                    NT = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NT']
                except lcm.PyLcmError:
                    NT = np.zeros(len(NWT0))
                try: 
                    NA = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NA']
                except lcm.PyLcmError:
                    NA = np.zeros(len(NWT0))
                try: 
                    N2A = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['N2A']
                except lcm.PyLcmError:
                    N2A = np.zeros(len(NWT0))
                try: 
                    NNP = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NNP']
                except lcm.PyLcmError:
                    NNP = np.zeros(len(NWT0))
                try:
                    N2N = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['N2N']
                except lcm.PyLcmError:
                    N2N = np.zeros(len(NWT0))
                try:
                    N3N = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['N3N']
                except lcm.PyLcmError:
                    N3N = np.zeros(len(NWT0))
                try:
                    N4N = pyCOMPO['H_EDI_REGI_2']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['N4N']
                except lcm.PyLcmError:
                    N4N = np.zeros(len(NWT0))
                
                # Use DRAGON definitition of neutronic absorption
                # ABS = NFTOT + NG + NP + ND + NT + NA + N2A + NNP - N2N - 2*N3N - 3*N4N
                ABS = np.array(NFTOT) + np.array(NG) + np.array(NP) + np.array(NA) - np.array(N2N) - 2*np.array(N3N) #- 3*N4N
                ABS_SIGS = np.array(TOT) - np.array(SIGS0)
                print(f"ABS from TOT-SIGS0 = {ABS_SIGS}, ABS from sum over reactions = {ABS}, difference = {np.array(ABS)-np.array(ABS_SIGS)}")
                    
                if mix in unique_mixes_on_diag:
                    sym_factor = 2
                else:
                    sym_factor = 1
                
                fission_rates[0][mix] += NFTOT[1]*NWT0[1]*N*vol*sym_factor
                fission_rates[1][mix] += NFTOT[0]*NWT0[0]*N*vol*sym_factor
                
                fiss_over_abs[0][mix] += (np.array(NFTOT[1])/np.array(ABS[1]))# *N*vol
                fiss_over_abs[1][mix] += (np.array(NFTOT[0])/np.array(ABS[0])) # *N*vol
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
    
    return keff_D5, fission_rates, fiss_over_abs, U238_ABS[::-1], FLUX_spectrum_ngroups[::-1], energy_mesh[::-1]