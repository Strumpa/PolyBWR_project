#####################################################################
#                                                                   #
# Purpose     : LIBRARY definition for gduo2_kec1 benchmark cell    #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def lib_gduo2(evaluation_name,  ssh_option, correlation, NG0):
    """
    evaluation_name : str, choice between 'J311_295' or 'endfb8r1_295' etc, defined in .access file
    ssh_option : str, choice between 'PT', 'RSE'
    correlation : str, choice between 'CORR' and 'NOCORR' to create the library with correlated resonant effects between Gd157 and U238 or not.
    NG0 : bool, if True : use modified DEPL structure with all NG0 Qvalues set to 0.  
    """
    # Lifo
    namLIB = "LIBRARY"
    myLifo=lifo.new()
    myLifo.pushEmpty(namLIB, "LCM")
    myLifo.push(evaluation_name)
    myLifo.lib()

    # Execution
    if NG0 and (evaluation_name == 'endfb8r1_295' or evaluation_name == 'endfb81295K' or evaluation_name == 'endfb81295K2'):
        if ssh_option == 'PT' and correlation == 'CORR':
            lib_gduo2_proc = cle2000.new('NG0_PTC',myLifo,1)
        elif ssh_option == 'RSE' and correlation == 'CORR':
            lib_gduo2_proc = cle2000.new('NG0_RSEC',myLifo,1)
        elif ssh_option == 'PT' and correlation == 'NOCORR':
            lib_gduo2_proc = cle2000.new('NG0_PT',myLifo,1)
        elif ssh_option == 'RSE' and correlation == 'NOCORR':
            lib_gduo2_proc = cle2000.new('NG0_RSE',myLifo,1)
    else:
        if ssh_option == 'PT' and correlation == 'CORR':
            lib_gduo2_proc = cle2000.new('lib_PTC',myLifo,1)
        elif ssh_option == 'RSE' and correlation == 'CORR':
            lib_gduo2_proc = cle2000.new('lib_RSEC',myLifo,1)
        elif ssh_option == 'PT' and correlation == 'NOCORR':
            lib_gduo2_proc = cle2000.new('lib_PT',myLifo,1)
        elif ssh_option == 'RSE' and correlation == 'NOCORR':
            lib_gduo2_proc = cle2000.new('lib_RSE',myLifo,1)
    
    lib_gduo2_proc.exec()
    # Recover
    myLifo.lib()
    pyLIB = myLifo.node(namLIB)
    # Clear stack before next execution
    while myLifo.getMax() > 0:
        myLifo.pop()

    return pyLIB
