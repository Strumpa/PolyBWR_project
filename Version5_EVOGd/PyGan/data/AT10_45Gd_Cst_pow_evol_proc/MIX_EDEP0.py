#############################################################################################
#                                                                                           #
# Description : PyGan script for BWR simulation with DRAGON5                                #
# Author      : R.Guasch                                                                    #
# Date        : 2025                                                                        #
# Purpose     : MIX definition for pin cells                                                #
# Modified DEPL structure to set fission Q-values to Serpent2 values : ER-ENU.              #
#############################################################################################
#
import lifo
import cle2000

def MIX_EDEP0(ssh_option, correlation):
    """
    Library : ENDF/B-VIII.1 (only supported here as DEPL structure is edited)
    ssh_option : self-shielding option (RSE or PT)
    correlation : correlation option (CORR or NOCORR)
    """
    if ssh_option == "RSE" and correlation == "CORR":
        # Lifo
        myLifo=lifo.new()
        myLifo.pushEmpty("LIBRARY", "LCM")
        myLifo.lib()

        # Execution
        mixNG0 = cle2000.new('EDEP0_RSEC',myLifo,1)
        mixNG0.exec()

        # Recover
        myLifo.lib()
        pyMIX = myLifo.node("LIBRARY")
    elif ssh_option == "RSE" and correlation == "NOCORR":
        # Lifo
        myLifo=lifo.new()
        myLifo.pushEmpty("LIBRARY", "LCM")
        myLifo.lib()

        # Execution
        mixNG0 = cle2000.new('EDEP0_RSE',myLifo,1)
        mixNG0.exec()

        # Recover
        myLifo.lib()
        pyMIX = myLifo.node("LIBRARY")

    elif ssh_option == "PT" and correlation == "CORR":
        # Lifo
        myLifo=lifo.new()
        myLifo.pushEmpty("LIBRARY", "LCM")
        myLifo.lib()

        # Execution
        mixNG0 = cle2000.new('EDEP0_PTC',myLifo,1)
        mixNG0.exec()

        # Recover
        myLifo.lib()
        pyMIX = myLifo.node("LIBRARY")

    elif ssh_option == "PT" and correlation == "NOCORR":
        # Lifo
        myLifo=lifo.new()
        myLifo.pushEmpty("LIBRARY", "LCM")
        myLifo.lib()

        # Execution
        mixNG0 = cle2000.new('EDEP0_PT',myLifo,1)
        mixNG0.exec()

        # Recover
        myLifo.lib()
        pyMIX = myLifo.node("LIBRARY")

    # Clear stack before next execution
    while myLifo.getMax() > 0:
        myLifo.pop()

    return pyMIX
