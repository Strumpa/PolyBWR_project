#####################################################################
#                                                                   #
# Purpose   : LIBRARY definition for HOM_UOX_Gd157 benchmark cell   #
#           -Edited burnup chain : all (n,gamma) Q-values set to 0  #
#           -Fission Q-values edited to match S2 edepmode 0         #
#####################################################################
#
import lifo
import cle2000

def RSEC_qfiss(draglib):
    """
    RSEC_qfiss() -> pyLIB, 
    draglib : (str) identifier to draglib used, (J311_295 or nJ311_295 supported here)
    Has to be based on Jeff3.1.1 for depletion chain consistency
    by default : RSE+CORR is set
    DEPL structure is edited to match S2 edepmode 0
    """
    # Lifo
    namLIB = "LIBRARY"
    myLifo=lifo.new()
    myLifo.pushEmpty(namLIB, "LCM")
    myLifo.push(draglib)
    myLifo.lib()

    # Execution
    RSEC_qfiss_proc = cle2000.new('RSEC_qfiss',myLifo,1)
    RSEC_qfiss_proc.exec()

    # Recover
    myLifo.lib()
    pyLIB = myLifo.node(namLIB)

    # Clear stack before next execution
    while myLifo.getMax() > 0: 
        myLifo.pop()
    
    return pyLIB
