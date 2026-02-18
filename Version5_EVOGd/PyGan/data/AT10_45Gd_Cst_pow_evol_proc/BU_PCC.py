### Python3 function wrapping BU_PCC CLE-2000 procedures
# Author : R. Guasch
# Date : 2025/04/23
# Purpose : Depletion calculation for pin cell with predictor-corrector method in CLE-2000

import lifo
import cle2000

def CECM(namCOMPO, pyMIX, pyTRACK, pyTF_EXC, pyTRACK_SS, pyTF_EXC_SS, StepList, name_compo, depl_sol_option, glob_option):
    """
    CECM method : Predictor-Corrector method for depletion calculation :
           - Constant extrpolation up until dt/2, 
           - Constant midpoint extrapolation on dt, constant midpoint corrector
    """
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
    myLifo.push(depl_sol_option)
    myLifo.push(glob_option)

    myLifo.lib()

    # Execution 
    calcBWR = cle2000.new('CECM',myLifo,1)
    calcBWR.exec()

    # Recover
    myLifo.lib()
    pyCOMPO = myLifo.node(namCOMPO)

    # Clear stack before next execution
    while myLifo.getMax() > 0: 
        myLifo.pop()

    return pyCOMPO


def CECE(namCOMPO, pyMIX, pyTRACK, pyTF_EXC, pyTRACK_SS, pyTF_EXC_SS, StepList, name_compo, depl_sol_option, glob_option):
    """
    CECE method : Predictor-Corrector method for depletion calculation :
              - Predictor step :Constant extrpolation until EoS,
              - Corrector step : Constant extrapolation on dt, using EoS rates
    """
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
    myLifo.push(depl_sol_option)
    myLifo.push(glob_option)

    myLifo.lib()

    # Execution 
    calcBWR = cle2000.new('CECE',myLifo,1)
    calcBWR.exec()

    # Recover
    myLifo.lib()
    pyCOMPO = myLifo.node(namCOMPO)

    # Clear stack before next execution
    while myLifo.getMax() > 0: 
        myLifo.pop()

    return pyCOMPO