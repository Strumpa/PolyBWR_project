#####################################################################
#                                                                   #
# Purpose   : LIBRARY definition for HOM_UOX_Gd157 benchmark cell   #
#           -Edited burnup chain : all (n,gamma) Q-values set to 0  #
#####################################################################
#
import lifo
import cle2000

def RSEC_NG0(draglib):
    """
    RSEC_NG0(draglib) -> pyLIB, 
    draglib : (str) identifier to draglib used, (J311_295 or nJ311_295 supported here)
    Has to be based on Jeff3.1.1 for depletion chain consistency
    by default : RSE+CORR is set
    """
    # Lifo
    namLIB = "LIBRARY"
    myLifo=lifo.new()
    myLifo.pushEmpty(namLIB, "LCM")
    myLifo.push(draglib)
    myLifo.lib()

    # Execution
    RSEC_NG0_proc = cle2000.new('RSEC_NG0',myLifo,1)
    RSEC_NG0_proc.exec()

    # Recover
    myLifo.lib()
    pyLIB = myLifo.node(namLIB)

    # Clear stack before next execution
    while myLifo.getMax() > 0: 
        myLifo.pop()
    
    return pyLIB
