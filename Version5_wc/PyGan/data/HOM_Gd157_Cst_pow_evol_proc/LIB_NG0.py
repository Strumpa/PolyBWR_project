#####################################################################
#                                                                   #
# Purpose   : LIBRARY definition for HOM_UOX_Gd157 benchmark cell   #
# Specificity : modified (n,gamma) reaction Q-values                #
#####################################################################
#
import lifo
import cle2000

def LIB_NG0(draglib_name):
    """
    draglib_name : str, choice of draglib to create the library with the corresponding evaluation. USe .access script to declare symbolic link to the draglib.
    ## WARNING : Compatibility with the DEPLTION structure is not guaranteed : only endfb8r1 draglibs are supported here.
    """
    # Lifo
    namLIB = "LIBRARY"
    myLifo=lifo.new()
    myLifo.pushEmpty(namLIB, "LCM")
    myLifo.push(draglib_name)
    myLifo.lib()

    # Execution
    LIB_NG0_proc = cle2000.new('LIB_NG0',myLifo,1)
    LIB_NG0_proc.exec()

    # Recover
    myLifo.lib()
    pyLIB = myLifo.node(namLIB)

    # Clear stack before next execution
    while myLifo.getMax() > 0: 
        myLifo.pop()
    
    return pyLIB
