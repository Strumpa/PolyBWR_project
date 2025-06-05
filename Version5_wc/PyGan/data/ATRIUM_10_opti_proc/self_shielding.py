## PyGan script for handling self-shielding in DRAGON5
# For now, only supports USS: calls for PT, RSE and SUBG methods.
# Due to the 3-level BWR assembly geometry, solution door is limited to SALT: CP or MOC.
# Due to time considerations, CP will be used since MOC is generally not adequate for self-shielding applications.
# Author: R. Guasch
# Date: 2025-06-03 

import lcm
import lifo
import cle2000

def selfShieldingUSS(mix_numbering_option, microlib, track_lcm, track_binary, name_geom, ssh_option):
    """
    Parameters:
    ----------
    
    mix_numbering_option : (str)
        Option for mix numbering, can be "number_mix_families_per_enrichment" or "number_mix_families_per_region" or "number_mix_families_per_selfshielding_region".
    microlib : (lcm object)
        LCM object containing the library microscopic cross sections for self-shielding calculations ie to be self-shielded.
    track_lcm : (lcm object)
        LCM object containing the tracking results.
    track_binary : (binary)
        Sequential binary tracking file containing the tracks lengths.
    name_geom : (str)
        Name of the geometry to be used for self-shielding calculations.
        Only "AT10_ASSBLY" is available for now, will need to add "AT10_ASSBLY_ctrl" and "AT10_ASSBLY_part" later
        to be able to describe assemblies with control cross and with part length rods (would change ssh parameters?)
    ssh_option : (str)
        Option to specify specific groupings of self-shielding regions, to be tested.
        "default" is default from USS: based on LIB: data,
        "RSE" is for Resonant Spectrum Expansion method, and groups all U8, U5 and Zr isotopes in a single self-shielding region --> to be tested.
        "VVER_REGI" is based of work done by the VVER team, to be tested.
    
    Returns:
    ----------
    libary_ssh : (lcm object)
        LCM object the self-shielded cross sections library.
    """

    if mix_numbering_option not in ["number_mix_families_per_enrichment", "number_mix_families_per_region", "number_mix_families_per_selfshielding_region"]:
        raise ValueError("Invalid mix numbering option selected.")
    if mix_numbering_option == "number_mix_families_per_enrichment":
        
        # Use mixAUTOP1.c2m procedure
        print("Creating self-shielded library with mix numbering per enrichment.")
        ipLifo = lifo.new() 
        ipLifo.pushEmpty("LIBRARY2", "LCM")
        ipLifo.push(microlib)
        ipLifo.push(track_lcm)
        ipLifo.push(track_binary)
        ipLifo.push(name_geom)
        ipLifo.push(ssh_option)
        

        # Create a cle2000 object to handle the library creation
        uss_proc = cle2000.new('mixAUTOP1', ipLifo, 1)
        # Execute the USS: procedure
        uss_proc.exec()

        # Recover the results from the LIFO stack
        ipLifo.lib()
        library_ssh = ipLifo.node('LIBRARY2')

        return library_ssh



