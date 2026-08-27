#####################################################################
#                                                                   #
# Purpose   : LIBRARY definition for HOM_UOX_Gd157 benchmark cell   #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def Mix(draglib_name, mix_option, ssh_method, edpmode_flag, N_H, N_O):
    """
    draglib_name : str, choice of draglib to create the library with the corresponding evaluation. USe .access script to declare symbolic link to the draglib.
    mix_option : choose to include Gd157 or not.
    ssh_method : choice of self-shielding method to be set in LIB:
    edpmode_flag : choice to activate EDP0 flag which takes lumoed isotopes energy contributions out of H-FACTORs.
    N_H : float, number density of H1 in H2O moderator (in atoms/barn-cm)
    N_O : float, number density of O16 in H2O moderator (in atoms/barn-cm)
    """
    # Lifo
    namLIB = "LIBRARY"
    myLifo=lifo.new()
    myLifo.pushEmpty(namLIB, "LCM")
    myLifo.push(draglib_name)
    myLifo.push(ssh_method)
    myLifo.push(edpmode_flag)
    myLifo.push(N_H)
    myLifo.push(N_O)
    myLifo.lib()

    # Execution
    if mix_option == "UO2_BWR_PINCELL":
        Mix_proc = cle2000.new('MIX_UOX',myLifo,1)
    elif mix_option == "GD_BWR_PINCELL":
        Mix_proc = cle2000.new('MIX_GD',myLifo,1)
    else:
        raise ValueError(f"Invalid mix option parameter : {mix_option}")
    
    Mix_proc.exec()
    # Recover
    myLifo.lib()
    pyLIB = myLifo.node(namLIB)

    # Clear stack before next execution
    while myLifo.getMax() > 0: 
        myLifo.pop()
    
    return pyLIB
