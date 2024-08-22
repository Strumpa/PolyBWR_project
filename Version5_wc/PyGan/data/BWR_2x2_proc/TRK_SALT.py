#####################################################################
#                                                                   #
# Description : PyGan scritp for BWR simulation with DRAGON5        #
# Author      : R.Guasch                                            #
# Date        : 29 July 2024                                        #
# Purpose     : Tracking for AT10 2x2 pins cluster                  #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def TRK_SALT(namGEOM_G2S, namFIG, namTRACK, namTF_EXC, pyGEOM, sol_door, an2d, calc_type):
    # Lifo
    
    myLifo=lifo.new()
    myLifo.pushEmpty(namGEOM_G2S, "ASCII")
    myLifo.pushEmpty(namFIG, "ASCII")
    myLifo.pushEmpty(namTRACK, "LCM")
    myLifo.pushEmpty(namTF_EXC, "BINARY")

    myLifo.push(pyGEOM)
    myLifo.push(sol_door)
    myLifo.push(an2d)
    myLifo.push(calc_type)

    #myLifo.lib()
    # Execution
    trackBWR = cle2000.new('TRK_SALT',myLifo,1)
    trackBWR.exec()

    # Recover
    myLifo.lib()

    pyGEOM_G2S = myLifo.node(namGEOM_G2S)
    pyFIG = myLifo.node(namFIG)
    pyTRACK = myLifo.node(namTRACK)
    pyTF_EXC = myLifo.node(namTF_EXC)


    return [pyGEOM_G2S,pyFIG,pyTRACK,pyTF_EXC]







