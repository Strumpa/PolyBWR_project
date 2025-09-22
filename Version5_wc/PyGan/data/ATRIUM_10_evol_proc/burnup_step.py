#                                                                   
# Description : PyGan script for BWR simulation with DRAGON5        
# Author      : R. Guasch                                           
# Date        : 18/09/2025                                      
# Purpose     : Depletion calculation for BWR fuel assemblies
#               with predictor only methods in CLE-2000               
#                                                                   

#
import lifo
import cle2000

def BU_A(ssh_lib_lcm, track_lcm, burnup_lcm, flux_lcm, BU_beg, BU_end, depl_sol_option, glob_option, sat_option, rates_extr, edep_mode):
    """
    Parameters:
    -----------
    0th order predictor : Constant extrapolation on dt, using BoS rates : NOEX
    1st order predictor : Linear extrapolation on dt, using BoS rates : EXTR
    """
    """
    PARAMETER LIBRARY2 BURNUP FLUX TRACK ::
    ::: LINKED_LIST LIBRARY2 ; 
    ::: LINKED_LIST BURNUP ;
    ::: LINKED_LIST FLUX ;
    ::: LINKED_LIST TRACK ; ;

    DOUBLE DBUbeg DBUend LastBU ;
    :: >>DBUbeg<< >>DBUend<< >>LastBU<< ;

    STRING deplsolopt glob_opt sat_opt rates_extr ;
    :: >>deplsolopt<< >>glob_opt<< >>sat_opt<< >>rates_extr<< ; 
    """
    # Lifo
    myLifo=lifo.new()
    myLifo.push(ssh_lib_lcm)
    if burnup_lcm is None:
        myLifo.pushEmpty("BURNUP","LCM")
    else:
        myLifo.push(burnup_lcm)
    myLifo.push(flux_lcm)
    myLifo.push(track_lcm)
    myLifo.push(BU_beg)
    myLifo.push(BU_end)
    myLifo.push(depl_sol_option)
    myLifo.push(glob_option)
    myLifo.push(sat_option)
    myLifo.push(rates_extr)

    myLifo.lib()
    # Execution 
    if edep_mode == "edep0":
        calcBWR = cle2000.new('BU_A_E0',myLifo,1)
        calcBWR.exec()
    else:
        calcBWR = cle2000.new('BU_A',myLifo,1)
        calcBWR.exec()

    # Recover
    myLifo.lib()
    burnup_lcm = myLifo.node("BURNUP")
    ssh_lib_lcm = myLifo.node("LIBRARY2")

    # Clear stack before next execution
    while myLifo.getMax() > 0: 
        myLifo.pop()

    return burnup_lcm, ssh_lib_lcm

def normalizeRates(flux_lcm, lib_ssh, burnup_lcm, track_lcm, rates_extr, power_density, burnup_step):
    """
    Normalization of reaction rates at the begining of the depletion step
    """
    # Lifo
    myLifo=lifo.new()
    myLifo.push(lib_ssh)
    if burnup_lcm is None:
        myLifo.pushEmpty("BURNUP","LCM")
    else:
        myLifo.push(burnup_lcm)
    myLifo.push(flux_lcm)
    myLifo.push(track_lcm)
    myLifo.push(burnup_step)
    myLifo.push(power_density)
    myLifo.push(rates_extr)
    myLifo.lib()
    # Execution
    normalisation = cle2000.new('NORM_RATES',myLifo,1)
    normalisation.exec()
    # Recover
    myLifo.lib()
    burnup_lcm = myLifo.node("BURNUP")
    lib_ssh = myLifo.node("LIBRARY2")
    # Clear stack before next execution
    while myLifo.getMax() > 0:
        myLifo.pop();
    return burnup_lcm, lib_ssh

    
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