## PyGan routine to handle flux calculations for the ATRIUM-10 BWR fuel assembly
# Author : R. Guasch
# Date : 2025-06-04
#


import lcm
import cle2000
import lifo 

def fluxCalculationMOC(track_lcm, track_binary, self_shielded_microlib, flux_lcm=None, sys_lcm=None, burnup_point=0.0):
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
    flux_lcm : (lcm object)
        LCM object containing the flux calculation results from a previous calculation, if None a new one will be created.
    burnup_point : (float)
        Current burnup point in MWd/kgU, used to check if BU>0 to use previous flux as initial guess.
    
        
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
    if sys_lcm is None:
        ipLifo.pushEmpty("SYS", "LCM")
    else:
        ipLifo.push(sys_lcm)
    ipLifo.push(track_lcm)
    ipLifo.push(track_binary)
    ipLifo.push(self_shielded_microlib)
    ipLifo.push(int(burnup_point))
    # Create a cle2000 object to handle the flux calculation
    moc_proc = cle2000.new('MOC_K_NOL', ipLifo, 1)
    # Execute the MOC calculation procedure
    moc_proc.exec()

    # Recover the results from the LIFO stack
    ipLifo.lib()
    flux_lcm = ipLifo.node('FLUX')
    sys_lcm = ipLifo.node('SYS')
    keff = flux_lcm["K-EFFECTIVE"][0] 
    ## Clear stack
    while ipLifo.getMax() > 0:
        ipLifo.pop();

    return keff, flux_lcm, sys_lcm


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
