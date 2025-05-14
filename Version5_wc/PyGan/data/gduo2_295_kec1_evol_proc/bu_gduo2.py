#####################################################################
#                                                                   #
# Purpose     : Depletion calculation for gduo2_kec1 cell           #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def bu_gduo2(namCOMPO, pyLIB, pyTRACK, pyTF_EXC, pyTRACK_SS,  StepList, name_compo, time_integrator, glob_opt):

    # Lifo
    myLifo=lifo.new()
    myLifo.pushEmpty(namCOMPO,"LCM")
    myLifo.push(pyLIB)
    myLifo.push(pyTRACK)
    myLifo.push(pyTF_EXC)
    myLifo.push(pyTRACK_SS)
    myLifo.push(StepList)
    myLifo.push(name_compo)
    myLifo.push(glob_opt)
    myLifo.lib()

    if time_integrator == "EXTR":
        # Execution 
        calcBWR = cle2000.new('EXTR_gduo2',myLifo,1)
        calcBWR.exec()
    elif time_integrator == "NOEX":
        calcBWR = cle2000.new('NOEX_gduo2',myLifo,1)
        calcBWR.exec()
    elif time_integrator == "EXTR2":
        calcBWR = cle2000.new('EXTR2_gduo2',myLifo,1)
        calcBWR.exec()
    elif time_integrator == "CECE":
        calcBWR = cle2000.new('CECE_gduo2',myLifo,1)
        calcBWR.exec()

    # Recover
    myLifo.lib()
    pyCOMPO = myLifo.node(namCOMPO)
    # Clear stack before next execution
    while myLifo.getMax() > 0: 
        myLifo.pop()

    return pyCOMPO
