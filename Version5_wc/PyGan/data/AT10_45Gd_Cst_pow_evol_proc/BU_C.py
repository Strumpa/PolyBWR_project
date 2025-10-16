#                                                                   
# Description : PyGan script for BWR simulation with DRAGON5        
# Author      : R. Guasch                                           
# Date        : 2025/02/25  (updated 2025/05/06)                                        
# Purpose     : Depletion calculation for pin cell   
#               with predictor only methods in CLE-2000               
#                                                                   

#
import lifo
import cle2000

def BU_C(namCOMPO, pyMIX, pyTRACK, pyTF_EXC, pyTRACK_SS, pyTF_EXC_SS, StepList, name_compo, ssh_option, depl_sol_option, glob_option, sat_option, rates_extr, edep_mode):
    """
    0th order predictor : Constant extrapolation on dt, using BoS rates : NOEX
    1st order predictor : Linear extrapolation on dt, using BoS rates : EXTR
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
    myLifo.push(ssh_option)
    myLifo.push(depl_sol_option)
    myLifo.push(glob_option)
    myLifo.push(sat_option)
    myLifo.push(rates_extr)

    myLifo.lib()
    # Execution 
    if edep_mode == "edep0":
        print(f"edep_mode = {edep_mode}, running EDP0 version of BU_C")
        calcBWR = cle2000.new('BU_C_E0',myLifo,1)
        calcBWR.exec()
    else:
        calcBWR = cle2000.new('BU_C',myLifo,1)
        calcBWR.exec()

    # Recover
    myLifo.lib()
    pyCOMPO = myLifo.node(namCOMPO)

    # Clear stack before next execution
    while myLifo.getMax() > 0: 
        myLifo.pop()

    return pyCOMPO


def BU_EXTR2(namCOMPO, pyMIX, pyTRACK, pyTF_EXC, pyTRACK_SS, pyTF_EXC_SS, StepList, name_compo, ssh_option, depl_sol_option, glob_option, sat_option):
    """
    Testing PARABOLIC extrapolation for reaction rates in Predictor only method
    2nd order predictor : Parabolic extrapolation on dt, using BoS rates : EXTR 2
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
    myLifo.push(ssh_option)
    myLifo.push(depl_sol_option)
    myLifo.push(glob_option)
    myLifo.push(sat_option)

    myLifo.lib()

    # Execution 
    calcBWR = cle2000.new('BU_EXTR2',myLifo,1)
    calcBWR.exec()

    # Recover
    myLifo.lib()
    pyCOMPO = myLifo.node(namCOMPO)

    # Clear stack before next execution
    while myLifo.getMax() > 0: 
        myLifo.pop()

    return pyCOMPO