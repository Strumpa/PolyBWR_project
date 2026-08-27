#####################################################################
#                                                                   #
# Purpose   : LIBRARY definition for HOM_UOX_Gd157 benchmark cell   #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def Mix(draglib_name, self_shielding_method, mix_option, edepmode_flag):
    """
    draglib_name : str, choice of draglib to create the library with the corresponding evaluation. USe .access script to declare symbolic link to the draglib.
    self_shielding_method : PT, RSE or SUBG to be set in LIB:
    mix_option : choose to include Gd157 or not.
    edepmode_flag : flag to test lumped isotopes energy in LIB:
    """
    # Lifo
    namLIB = "LIBRARY"
    myLifo=lifo.new()
    myLifo.pushEmpty(namLIB, "LCM")
    myLifo.push(draglib_name)
    myLifo.push(edepmode_flag)
    myLifo.push(self_shielding_method)
    myLifo.lib()

    # Execution
    if mix_option == "HOM_UOX_Gd157":
        Mix_proc = cle2000.new('MixGd157',myLifo,1)
    elif mix_option == "HOM_UOX":
        Mix_proc = cle2000.new('MixUOX',myLifo,1)
    elif mix_option == "HOM_U5":
        Mix_proc = cle2000.new('MixU5',myLifo,1)
    elif mix_option == "HOM_U8":
        Mix_proc = cle2000.new('MixU8',myLifo,1)
    elif mix_option == "HOM_U45":
        Mix_proc = cle2000.new('MixU45',myLifo,1)
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
