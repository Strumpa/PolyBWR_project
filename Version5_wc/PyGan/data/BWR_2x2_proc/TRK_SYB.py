#####################################################################
#                                                                   #
# Description : PyGan scritp for BWR simulation with DRAGON5        #
# Author      : R.Guasch                                            #
# Date        : 2024                                                #
# Purpose     : Tracking with SYBILT                                #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def TRK_SYB(namTRACK, pyGEOM, iqua2, nseg):
    # Lifo
    myLifo=lifo.new()


    myLifo.pushEmpty(namTRACK, "LCM")

    myLifo.push(pyGEOM)
    myLifo.push(iqua2)
    myLifo.push(nseg)

    # Execution
    trackBWR = cle2000.new('TRK_SYB',myLifo,1)
    trackBWR.exec()

    # Recover
    myLifo.lib()

    pyTRACK = myLifo.node(namTRACK)

    #pyTF_EXC_SS = myLifo.node(namTF_EXC_SS)

    return pyTRACK
