#####################################################################
#                                                                   #
# Purpose   : LIBRARY definition for HOM_UOX_Gd157 benchmark cell   #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def MixGd157(draglib_name, ssh_method, correlation):
    """
    draglib_name : str, choice of draglib to create the library with the corresponding evaluation. USe .access script to declare symbolic link to the draglib.
    ssh_method : str, choice between 'RSE', 'PT' and 'SUBG' to create the library with the corresponding self-shielding method.
    correlation : str, choice between 'CORR' and 'NOCORR' to create the library with correlated resonant effects between Gd157 and U238 or not.
    """
    # Lifo
    namLIB = "LIBRARY"
    myLifo=lifo.new()
    myLifo.pushEmpty(namLIB, "LCM")
    myLifo.push(draglib_name)
    myLifo.push(ssh_method)
    myLifo.push(correlation)
    myLifo.lib()

    # Execution
    MixGd157_proc = cle2000.new('MixGd157',myLifo,1)
    MixGd157_proc.exec()

    # Recover
    myLifo.lib()
    pyLIB = myLifo.node(namLIB)

    # Clear stack before next execution
    while myLifo.getMax() > 0: 
        myLifo.pop()
    
    return pyLIB
