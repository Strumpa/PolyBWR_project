### Python3 function wrapping BU_PCC.c2m CLE-2000 procedure
# Author : R. Guasch
# Date : 2025/04/23
# Purpose : Depletion calculation for pin cell with predictor-corrector method in CLE-2000

import lifo
import cle2000

def BU_PCC(namCOMPO, pyMIX, pyTRACK, pyTF_EXC, pyTRACK_SS, pyTF_EXC_SS, StepList, name_compo, depl_sol_option, glob_option, bos_rates_opt, eos_rates_opt):
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
    myLifo.push(bos_rates_opt)
    myLifo.push(eos_rates_opt)

    myLifo.lib()

    # Execution 
    calcBWR = cle2000.new('BU_PCC',myLifo,1)
    calcBWR.exec()

    # Recover
    myLifo.lib()
    pyCOMPO = myLifo.node(namCOMPO)

    # Clear stack before next execution
    while myLifo.getMax() > 0: 
        myLifo.pop()

    return pyCOMPO