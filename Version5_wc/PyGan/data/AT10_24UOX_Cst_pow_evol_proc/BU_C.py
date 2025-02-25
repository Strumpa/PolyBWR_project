#####################################################################
#                                                                   #
# Description : PyGan script for BWR simulation with DRAGON5        #
# Author      : R. Guasch                                           #
# Date        : 2025/02/25                                          #
# Purpose     : Depletion calculation for pin cell                  #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def BU_C(namCOMPO, pyMIX, pyTRACK, pyTF_EXC, pyTRACK_SS, pyTF_EXC_SS, StepList, name_compo, ssh_option, depl_sol_option, glob_option, sat_option, val_exp):
    # Lifo
    myLifo=lifo.new()
    myLifo.pushEmpty(namCOMPO,"LCM")
    myLifo.push(pyMIX)
    myLifo.push(pyTRACK)
    myLifo.push(pyTF_EXC)
    myLifo.push(pyTRACK_SS)
    myLifo.push(pyTF_EXC_SS)
    myLifo.push(StepList)
    myLifo.push(name_compo)
    myLifo.push(ssh_option)
    myLifo.push(depl_sol_option)
    myLifo.push(glob_option)
    myLifo.push(sat_option)
    myLifo.push(val_exp)

    myLifo.lib()

    # Execution 
    calcBWR = cle2000.new('BU_C',myLifo,1)
    calcBWR.exec()

    # Recover
    myLifo.lib()
    pyCOMPO = myLifo.node(namCOMPO)

    # Clear stack before next execution
    while myLifo.getMax() > 0: 
        myLifo.pop()

    return pyCOMPO
