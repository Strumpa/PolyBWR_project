## PyGan script for handling self-shielding in DRAGON5
# For now, only supports USS: calls for PT, RSE and SUBG methods.
# Due to the 3-level BWR assembly geometry, solution door is limited to SALT: CP or MOC.
# Due to time considerations, CP will be used since MOC is generally not adequate for self-shielding applications.
# Author: R. Guasch
# Date: 2025-06-03 

import lcm
import lifo
import cle2000

def selfShieldingUSS(microlib, track_lcm, track_binary):
    """
    Parameters:
    ----------
    
    microlib : (lcm object)
        LCM object containing the library microscopic cross sections for self-shielding calculations ie to be self-shielded.
    track_lcm : (lcm object)
        LCM object containing the tracking results.
    track_binary : (binary)
        Sequential binary tracking file containing the tracks lengths.
    
    Returns:
    ----------
    libary_ssh : (lcm object)
        LCM object the self-shielded cross sections library.
    """
    ipLifo = lifo.new() 
    ipLifo.pushEmpty("LIBRARY2", "LCM")
    ipLifo.push(microlib)
    ipLifo.push(track_lcm)
    ipLifo.push(track_binary)
    
    # Create a cle2000 object to handle the library creation
    uss_proc = cle2000.new('USS_call', ipLifo, 1)
    # Execute the USS: procedure
    uss_proc.exec()

    # Recover the results from the LIFO stack
    ipLifo.lib()
    library_ssh = ipLifo.node('LIBRARY2')

    return library_ssh