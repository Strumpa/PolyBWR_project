## PyGan script for handling self-shielding in DRAGON5
# For now, only supports USS: calls for PT, RSE and SUBG methods.
# Due to the 3-level BWR assembly geometry, solution door is limited to SALT: CP or MOC.
# Due to time considerations, CP will be used since MOC is generally not adequate for self-shielding applications.
# Author: R. Guasch
# Date: 2025-06-03 

import lcm
import lifo
import cle2000
from mix_handling import cle2000_variable_declaration

def selfShieldingUSS(mix_numbering_option, microlib, track_lcm, track_binary, name_geom, ssh_option, solution_door, mix_connectivity_dict=None):
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
    solution_door : (str)
        Solution door to be used for self-shielding calculations, can be "PIJ" for full pij reconstruction or "IC" for interface current.
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
    elif mix_numbering_option == "number_mix_families_per_region":
        # Create the self sheilding procedure based on mix numbering definition
        print("Creating self-shielded library with mix numbering per region.")
        proc_name = "SSH_R_A.c2m"

        fill_ssh_proc(proc_name, mix_connectivity_dict, mix_numbering_option, ssh_option, solution_door)

        ipLifo = lifo.new() 
        ipLifo.pushEmpty("LIBRARY2", "LCM")
        ipLifo.push(microlib)
        ipLifo.push(track_lcm)
        ipLifo.push(track_binary)

        # Create a cle2000 object to handle the library creation
        uss_proc = cle2000.new(proc_name.split(".")[0], ipLifo, 1)
        # Execute the USS: procedure
        uss_proc.exec()

        # Recover the results from the LIFO stack
        ipLifo.lib()
        library_ssh = ipLifo.node('LIBRARY2')

        return library_ssh

def fill_ssh_proc(file_name, mix_connectivity_dict, mix_numbering_option, ssh_option, solution_door):
    """
    Generate CLE-2000 procedure for self-shielding calculations.
    Parameters:
    ----------
    file_name : (str)
        Name of the file to save the procedure.
    mix_connectivity_dict : (dict)
        Dictionary containing the connectivity of the mixes.
    mix_numbering_option : (str)
        Option for mix numbering, can be "number_mix_families_per_enrichment" or "number_mix_families_per_region" or "number_mix_families_per_selfshielding_region".
    ssh_option : (str)
        Self-shielding option to be used in the procedure, for now only default is available. but should be developed to include optimized use of CALC / REGI keywords.
    solution_door : (str)
        Solution door to be used for self-shielding calculations, can be "PIJ" for full pij reconstruction or "IC" for interface current.
    Returns:
    """
    header = (
            "* --------------------------------\n"
            "* Procedure generated by mix_handling.py\n"
            "* Author: R. Guasch\n"  
            "* --------------------------------\n"
            "*    INPUT & OUTPUT PARAMETERS\n"
            "* --------------------------------\n"
            "PARAMETER LIBRARY2 LIBRARY TRACK_SS TF_SS_EXC ::\n"
            "::: LINKED_LIST LIBRARY2 ;\n"
            "::: LINKED_LIST LIBRARY ;\n"
            "::: LINKED_LIST TRACK_SS ;\n"
            "::: SEQ_BINARY TF_SS_EXC ; ;\n"

            "\n"
            "* -------------------------------\n"
            "*    STRUCTURES AND MODULES\n"
            "* -------------------------------\n"
            "MODULE USS: END: ;\n"
            "\n"
            )
    # Generate variable declaration based on mix numbering option
    variable_declaration = cle2000_variable_declaration(mix_connectivity_dict)
    ssh_proc = (
    f"*PROCEDURE {file_name}\n"
    f"{header}\n"
    "* --------------------------------\n"
    "*    MIX NUMBERING DEFINITION\n"
    "* --------------------------------\n"
    f"{variable_declaration}\n"
    "* --------------------------------\n"
    "*    SELF-SHIELDING OPERATOR DEFINITION\n"
    "*    CALL TO USS: MODULE\n"
    "* --------------------------------\n"
    f"{fill_uss_call(ssh_option, solution_door)}\n"
    "* -----------------------------------------\n"
    "*     END OF USS: CALL\n"
    "* -----------------------------------------\n"
    "END: ;\n"
    "QUIT ."
    )
    # Write the procedure to a file
    with open(file_name, 'w') as file:
        file.write(ssh_proc)
    file.close()    
    
def fill_uss_call(ssh_option, solution_door, iter_pass = 3, max_st_iter = 50):
    """
    define call to USS: module for self-shielding calculations.
    """
    if solution_door not in ["PIJ", "IC"]:
        raise ValueError("Invalid solution door option selected.")
    if solution_door == "PIJ":
        solver = "PIJ"
    elif solution_door == "IC":
        solver = "ARM"
    if ssh_option == "default":
        uss_call = (
        "LIBRARY2 := USS: LIBRARY TRACK_SS TF_SS_EXC ::\n"
        f"    EDIT 1 PASS {iter_pass} MAXST {max_st_iter} {solver}\n"
        ";\n"
        )
    ### need to make more options available
        
    return uss_call






