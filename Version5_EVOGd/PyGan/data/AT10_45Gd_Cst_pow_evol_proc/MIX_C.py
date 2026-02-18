#############################################################################
#                                                                           #
# Description : Create material mixtures for single Gd cell test case       #
# Author      : R.Guasch                                                    #
# Date        : 2025                                                        #
# Purpose     : Compare different reaction rates normalization options      #
#                                                                           #
#############################################################################
#
import lifo
import cle2000

def MIX_C(Library, ssh_option, correlation, rates_norm_opt):
  
    # Lifo
    myLifo=lifo.new()
    myLifo.pushEmpty("LIBRARY", "LCM")
    if rates_norm_opt != "EDP0":
        myLifo.push(Library)
    myLifo.lib()
    print("Using library : ", Library)
    # Execution
    #mixBWR = cle2000.new('QJ311_PT',myLifo,1)
    
    if "endfb8r1" in Library:
        if rates_norm_opt == "EDP0":
            if ssh_option == "RSE" and correlation == "CORR":
                print("Using Q81_RSE_CORR with EDP0 rates normalization")
                mixBWR = cle2000.new('Q81_RSEC',myLifo,1)
            elif ssh_option == "RSE" and correlation == "NOCORR":
                print("Using Q81_RSE with EDP0 rates normalization")
                mixBWR = cle2000.new('Q81_RSE',myLifo,1)
            elif ssh_option == "PT" and correlation == "CORR":
                print("Using Q81_PT_CORR with EDP0 rates normalization")
                mixBWR = cle2000.new('Q81_PTC',myLifo,1)
            elif ssh_option == "PT" and correlation == "NOCORR":
                print("Using Q81_PT with EDP0 rates normalization")
                mixBWR = cle2000.new('Q81_PT',myLifo,1)
        else:
            if ssh_option == "RSE" and correlation == "CORR":
                mixBWR = cle2000.new('MIX_C_RSEC',myLifo,1)
            elif ssh_option == "RSE" and correlation == "NOCORR":
                mixBWR = cle2000.new('MIX_C_RSE',myLifo,1)
            elif ssh_option == "PT" and correlation == "CORR":
                mixBWR = cle2000.new('MIX_C_PTC',myLifo,1)
            elif ssh_option == "PT" and correlation == "NOCORR":
                mixBWR = cle2000.new('MIX_C_PT',myLifo,1)
    elif "J311" in Library:
        if rates_norm_opt == "EDP0":
            print("Using MIX_C_J3 with EDP0 rates normalization")
            if ssh_option == "RSE" and correlation == "CORR":
                mixBWR = cle2000.new('QJ311_RSEC',myLifo,1)
            elif ssh_option == "RSE" and correlation == "NOCORR":
                mixBWR = cle2000.new('QJ311_RSE',myLifo,1)
            elif ssh_option == "PT" and correlation == "CORR":
                mixBWR = cle2000.new('QJ311_PTC',myLifo,1)
            elif ssh_option == "PT" and correlation == "NOCORR":
                mixBWR = cle2000.new('QJ311_PT',myLifo,1)
        else:
            if ssh_option == "RSE" and correlation == "CORR":
                mixBWR = cle2000.new('MIX_C_RSEC',myLifo,1)
            elif ssh_option == "RSE" and correlation == "NOCORR":
                mixBWR = cle2000.new('MIX_C_RSE',myLifo,1)
            elif ssh_option == "PT" and correlation == "CORR":
                mixBWR = cle2000.new('MIX_C_PTC',myLifo,1)
            elif ssh_option == "PT" and correlation == "NOCORR":
                mixBWR = cle2000.new('MIX_C_PT',myLifo,1)
        
    """
    if Library == "endfb8r1_295" or Library == "endfb81295K" or Library == "endfb81295K2":
        
    elif Library == "J311_295" or Library == "J311_295K":
        if ssh_option == "PT" and correlation == "NOCORR":
            print("Using QJ311_PT")
            mixBWR = cle2000.new('QJ311_PT',myLifo,1)
      
        elif ssh_option == "RSE" and correlation == "NOCORR":
            mixBWR = cle2000.new('MIX_C_RSE_J3',myLifo,1)
        elif ssh_option == "PT" and correlation == "CORR":
            mixBWR = cle2000.new('MIX_C_PTC_J3',myLifo,1)
        elif ssh_option == "PT" and correlation == "NOCORR":
            mixBWR = cle2000.new('MIX_C_PT_J3',myLifo,1)
    """
    mixBWR.exec()

    # Recover
    myLifo.lib()
    pyMIX = myLifo.node("LIBRARY")

    # Clear stack before next execution
    while myLifo.getMax() > 0:
        myLifo.pop()

    return pyMIX
