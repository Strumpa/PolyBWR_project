# Collection of function used to parse DRAGON5 COMPO output files for fission and neutron gamma rates post treatment.
import os
from re import U
import numpy as np
import lcm
from collections import defaultdict
from utilities import sum_rates_over_iso


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