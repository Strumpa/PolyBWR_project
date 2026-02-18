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

def MIX_C(Library,ssh_option,rates_norm_opt):
    # Lifo
    myLifo=lifo.new()
    myLifo.pushEmpty("LIBRARY", "LCM")
    if rates_norm_opt == "EDP0":
        mixBWR = cle2000.new('Q81_PT',myLifo,1)
    else:
        myLifo.push(Library)
        myLifo.push(ssh_option)
        mixBWR = cle2000.new('MIX_C',myLifo,1)
    # Execution
    mixBWR.exec()
    # Recover
    myLifo.lib()
    pyMIX = myLifo.node("LIBRARY")

    # Clear stack before next execution
    while myLifo.getMax() > 0:
        myLifo.pop()

    return pyMIX
