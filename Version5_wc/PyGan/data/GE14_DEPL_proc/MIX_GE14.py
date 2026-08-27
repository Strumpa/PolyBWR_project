#####################################################################
#                                                                   #
# Purpose   : LIBRARY definition GE14 DOM and VAN lattices          #
# Author    : R.Guasch                                              #
# Date      : 11/08/2026                                            #
#####################################################################
#
import lifo
import cle2000

def Mix(draglib_name, mix_option, ssh_method, anisotropy_level, transport_correction, TFuel, TCool, N_H, N_O, edpmode_flag):
    """
    draglib_name : str, choice of draglib to create the library with the corresponding evaluation. USe .access script to declare symbolic link to the draglib.
    mix_option : choose to include Gd157 or not.
    ssh_method : choice of self-shielding method to be set in LIB:
    edpmode_flag : choice to activate EDP0 flag which takes lumoed isotopes energy contributions out of H-FACTORs.
    """
    # Lifo
    namLIB = "LIBRARY"
    myLifo=lifo.new()
    myLifo.pushEmpty(namLIB, "LCM")
    myLifo.push(draglib_name)
    myLifo.push(ssh_method)
    myLifo.push(anisotropy_level)
    myLifo.push(transport_correction)
    myLifo.push(TFuel)
    myLifo.push(TCool)
    myLifo.push(N_H)
    myLifo.push(N_O)
    myLifo.push(edpmode_flag)
    myLifo.lib()

    # Execution
    if mix_option == "DOM":
        Mix_proc = cle2000.new('MIX_DOM',myLifo,1)
    elif mix_option == "VAN":
        Mix_proc = cle2000.new('MIX_VAN',myLifo,1)
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
