## PyGan routine to handle flux calculations for the ATRIUM-10 BWR fuel assembly
# Author : R. Guasch
# Date : 2025-06-04
#


import lcm
import cle2000
import lifo 

def fluxCalculationMOC(track_lcm, track_binary, self_shielded_microlib, flux_lcm=None):
    """
    This function prepares a DRAGON5 MOC eigenvalue calculation.
    An option to include leakage model calculation is to be implemented.
    Parameters:
    ----------
    track_lcm : (lcm object)
        LCM object containing the tracking results.
    track_binary : (binary)
        Sequential binary tracking file containing the tracks lengths.
    self_shielded_microlib : (lcm object)
        LCM object containing the self-shielded cross sections library.
    
        
    Returns:
    ----------
    keff : (float)
        Effective multiplication factor from the flux calculation.
    flux_lcm : (lcm object)
        LCM object containing the flux calculation results.


    """

    # Run MOC_K_noL.c2m procedure
    ipLifo = lifo.new()
    if flux_lcm is None:
        ipLifo.pushEmpty("FLUX", "LCM")
    else:
        ipLifo.push(flux_lcm)
    ipLifo.push(track_lcm)
    ipLifo.push(track_binary)
    ipLifo.push(self_shielded_microlib)
    # Create a cle2000 object to handle the flux calculation
    moc_proc = cle2000.new('MOC_K_NOL', ipLifo, 1)
    # Execute the MOC calculation procedure
    moc_proc.exec()

    # Recover the results from the LIFO stack
    ipLifo.lib()
    flux_lcm = ipLifo.node('FLUX')
    keff = flux_lcm["K-EFFECTIVE"][0] 
    
    ## Clear stack
    while ipLifo.getMax() > 0:
        ipLifo.pop();
        

    return keff, flux_lcm


def fluxCalculationPIJ(track_lcm, track_binary, self_shielded_microlib):
    """
    This function prepares a DRAGON5 PIJ flux calculation.
    Parameters:
    ----------
    track_lcm : (lcm object)
        LCM object containing the tracking results.
    track_binary : (binary)
        Sequential binary tracking file containing the tracks lengths.
    self_shielded_microlib : (lcm object)
        LCM object containing the self-shielded cross sections library.

    Returns:
    ----------
    keff : (float)
        Effective multiplication factor from the flux calculation.
    flux_lcm : (lcm object)
        LCM object containing the flux calculation results.

    """

    # Run PIJ.c2m procedure
    ipLifo = lifo.new()
    ipLifo.pushEmpty("FLUX", "LCM")
    ipLifo.push(track_lcm)
    ipLifo.push(track_binary)
    ipLifo.push(self_shielded_microlib)
    # Create a cle2000 object to handle the flux calculation
    pij_proc = cle2000.new('PIJ_K_NOL', ipLifo, 1)
    # Execute the PIJ calculation procedure
    pij_proc.exec()
    
    # Recover the results from the LIFO stack
    ipLifo.lib()
    flux_lcm = ipLifo.node('FLUX')
    keff = flux_lcm["K-EFFECTIVE"][0]

    return keff, flux_lcm


def fluxCalculationIC(track_lcm, track_binary, self_shielded_microlib):
    """
    This function prepares a DRAGON5 IC flux calculation.
    Parameters:
    ----------
    track_lcm : (lcm object)
        LCM object containing the tracking results.
    track_binary : (binary)
        Sequential binary tracking file containing the tracks lengths.
    self_shielded_microlib : (lcm object)
        LCM object containing the self-shielded cross sections library.

    Returns:
    ----------
    keff : (float)
        Effective multiplication factor from the flux calculation.
    flux_lcm : (lcm object)
        LCM object containing the flux calculation results.

    """

    # Run IC.c2m procedure
    ipLifo = lifo.new()
    ipLifo.pushEmpty("FLUX", "LCM")
    ipLifo.push(track_lcm)
    ipLifo.push(track_binary)
    ipLifo.push(self_shielded_microlib)
    # Create a cle2000 object to handle the flux calculation
    ic_proc = cle2000.new('IC_K_NOL', ipLifo, 1)
    # Execute the IC calculation procedure
    ic_proc.exec()

    # Recover the results from the LIFO stack
    ipLifo.lib()
    flux_lcm = ipLifo.node('FLUX')
    keff = flux_lcm["K-EFFECTIVE"][0]

    return keff, flux_lcm


def fluxCalculationICMOC(track_MOC, track_bin_MOC, track_IC, track_bin_IC, self_shielded_microlib):
    """ 
    This function runs a DRAGON5 2-level flux calculation.
    This first level is an IC calculation on 295 energy groups,
    the second level is a MOC calculation on 295 energy groups, using the IC obtained flux as an initial guess.
    
    Parameters:
    ----------
    track_MOC : (lcm object)
        LCM object containing the tracking results for the MOC flux calculation.
    track_bin_MOC : (binary)
        Sequential binary tracking file containing the tracks lengths for the MOC flux calculation.
    track_IC : (lcm object)
        LCM object containing the tracking results for the IC flux calculation.
    track_bin_IC : (binary)
        Sequential binary tracking file containing the tracks lengths for the IC flux calculation.
    self_shielded_microlib : (lcm object)
        LCM object containing the self-shielded cross sections library.

    Returns:
    ----------
    keff : (float)
        Effective multiplication factor from the flux calculation.
    flux_lcm : (lcm object)
        LCM object containing the flux calculation results.
    """
    # Run 2L_FLUX.c2m procedure
    ipLifo = lifo.new()
    ipLifo.pushEmpty("FLUX", "LCM")
    ipLifo.pushEmpty("LIB26G", "LCM")
    ipLifo.push(track_MOC)
    ipLifo.push(track_bin_MOC)
    ipLifo.push(track_IC)
    ipLifo.push(track_bin_IC)
    ipLifo.push(self_shielded_microlib)
    # Create a cle2000 object to handle the flux calculation
    flux2l_proc = cle2000.new('2L_FLUX_IC_I', ipLifo, 1)
    # Execute the 2-level flux calculation procedure
    flux2l_proc.exec()
    # Recover the results from the LIFO stack
    ipLifo.lib()
    flux_lcm = ipLifo.node('FLUX')
    while ipLifo.getMax() > 0:
        ipLifo.pop()
    keff = flux_lcm["K-EFFECTIVE"][0]

    return keff, flux_lcm

def fluxCalculation2LScheme(track_lcm, track_binary, track_lvl1, track_binary_lvl1, self_shielded_microlib, first_level_solution_door, sph_grmax=22):
    """
    This function runs a DRAGON5 2-level flux calculation.
    A IC/PIJ calculation is performed at the 1st level on 295 energy groups, 
    microscopic cross sections are condensed to 26 groups
    A SPH: equivalence is performed,
    and a MOC calculation is performed at the 2nd level on 26 energy groups.
    
    Parameters:
    ----------
    track_lcm : (lcm object)
        LCM object containing the tracking results for the 2nd level flux calculation.
    track_binary : (binary)
        Sequential binary tracking file containing the tracks lengths for the 2nd level flux calculation.
    track_lvl1 : (lcm object)
        LCM object containing the tracking results for the 1st level flux calculation.
    track_binary_lvl1 : (binary)
        Sequential binary tracking file containing the tracks lengths for the 1st level flux calculation.
    self_shielded_microlib : (lcm object)
        LCM object containing the self-shielded cross sections library.
    first_level_solution_door : (int)
        Solution door to be used for the 1st level flux calculation "PIJ" or "IC".
    Returns:
    ----------
    flux_lcm : (lcm object)
        LCM object containing the flux calculation results.
        
    Note:
    2nd level flux calculation is performed with the "MOC".
    """
    """ 
    PARAMETER FLUX LIB26G TRACK TF_EXC TRACK_L1 TF_EXC_L1 LIBRARY2 ::
    ::: LINKED_LIST FLUX ;
    ::: LINKED_LIST LIB26G ;
    ::: LINKED_LIST TRACK ;
    ::: SEQ_BINARY TF_EXC ;
    ::: LINKED_LIST TRACK_L1 ;
    ::: SEQ_BINARY TF_EXC_L1 ;
    ::: LINKED_LIST LIBRARY2 ; ;
    + Level 1 sol door
    """
    # Run 2L_FLUX.c2m procedure
    ipLifo = lifo.new()
    ipLifo.pushEmpty("FLUX", "LCM")
    ipLifo.pushEmpty("LIB26G", "LCM")
    ipLifo.push(track_lcm)
    ipLifo.push(track_binary)
    ipLifo.push(track_lvl1)
    ipLifo.push(track_binary_lvl1)
    ipLifo.push(self_shielded_microlib)
    ipLifo.push(first_level_solution_door)
    ipLifo.push(sph_grmax)
    
    # Create a cle2000 object to handle the flux calculation
    flux2l_proc = cle2000.new('2L_FLUX', ipLifo, 1)
    # Execute the 2-level flux calculation procedure
    flux2l_proc.exec()
    # Recover the results from the LIFO stack
    ipLifo.lib()
    flux_lcm = ipLifo.node('FLUX')
    lib_26g = ipLifo.node('LIB26G')
    while ipLifo.getMax() > 0:
        ipLifo.pop()
    keff = flux_lcm["K-EFFECTIVE"][0]

    return keff, flux_lcm, lib_26g