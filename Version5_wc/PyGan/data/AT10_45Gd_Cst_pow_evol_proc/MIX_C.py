#####################################################################
#                                                                   #
# Description : PyGan scritp for BWR simulation with DRAGON5        #
# Author      : L. Fede, used by R.Guasch                           #
# Date        : 2023                                                #
# Purpose     : MIX definition for pin cells                        #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def MIX_C(Library, ssh_option, correlation):
  
    # Lifo
    myLifo=lifo.new()
    myLifo.pushEmpty("LIBRARY", "LCM")
    myLifo.push(Library)
    myLifo.lib()

    # Execution
    if ssh_option == "RSE" and correlation == "CORR":
        mixBWR = cle2000.new('MIX_C_RSEC',myLifo,1)
    elif ssh_option == "RSE" and correlation == "NOCORR":
        mixBWR = cle2000.new('MIX_C_RSE',myLifo,1)
    elif ssh_option == "PT" and correlation == "CORR":
        mixBWR = cle2000.new('MIX_C_PTC',myLifo,1)
    elif ssh_option == "PT" and correlation == "NOCORR":
        mixBWR = cle2000.new('MIX_C_PT',myLifo,1)
    
    mixBWR.exec()

    # Recover
    myLifo.lib()
    pyMIX = myLifo.node("LIBRARY")

    # Clear stack before next execution
    while myLifo.getMax() > 0:
        myLifo.pop()

    return pyMIX
