## PyGan edi_compo.py procedure to handle calls to EDI: and COMPO: modues
# Author : R. Guasch
# Date : 2025-06-04

import lcm
import cle2000
import lifo


def ediCompo(mix_numbering_option, track_lcm, self_shielded_microlib, name_compo, save_option):
    """
    This function prepares a DRAGON5 EDI: and COMPO: procedure call.
    It handles the composition and cross-section data for the ATRIUM-10 BWR fuel assembly.
    
    Parameters:
    ----------
    mix_numbering_option : str
        Option for mix numbering (e.g., "number_mix_families_per_enrichment").
    track_lcm : lcm object
        LCM object containing the tracking results.
    self_shielded_microlib : lcm object
        LCM object containing the self-shielded cross sections library.
    name_compo : str
        Name of the composition to be processed.

    Returns:
    ----------
    compo_lcm : lcm object
        LCM object containing the composition results.
    """
    
    if mix_numbering_option == "number_mix_families_per_enrichment":
        # Initialize the LIFO stack for EDI: and COMPO: procedure
        ipLifo = lifo.new()
        
        # Push necessary data onto the LIFO stack
        ipLifo.pushEmpty("COMPO", "LCM")
        ipLifo.push(track_lcm)
        ipLifo.push(self_shielded_microlib)
        ipLifo.push(name_compo)
        ipLifo.push(save_option)  # Save option for the composition data

        # Create a cle2000 object to handle the EDI: and COMPO: procedure call
        edi_compo_proc = cle2000.new('EDICPO1', ipLifo, 1)
        
        # Execute the EDI_COMPO procedure
        edi_compo_proc.exec()

        # Recover the results from the LIFO stack
        compo_lcm = ipLifo.node('Compo')

    return compo_lcm
