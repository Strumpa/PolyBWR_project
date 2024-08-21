#####################################################################
#                                                                   #
# Description : PyGan scritp for BWR simulation with DRAGON5        #
# Author      : R.Guasch                                            #
# Date        : 2024                                                #
# Purpose     : Tracking with NXT module + MCCGT if sol_door = MOC  #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def TRK_NXT(namTRACK,namTF_EXC,pyGEOM, sol_door, an2d, calc_type):
    # Lifo
    myLifo=lifo.new()


    myLifo.pushEmpty(namTRACK, "LCM")
    myLifo.pushEmpty(namTF_EXC, "BINARY")

    myLifo.push(pyGEOM)

    myLifo.push(sol_door)
    myLifo.push(an2d)
    myLifo.push(calc_type)

    # Execution
    trackBWR = cle2000.new('TRK_NXT',myLifo,1)
    trackBWR.exec()

    # Recover
    myLifo.lib()

    pyTRACK = myLifo.node(namTRACK)
    pyTF_EXC = myLifo.node(namTF_EXC)

    #pyTF_EXC_SS = myLifo.node(namTF_EXC_SS)

    return [pyTRACK,pyTF_EXC]



