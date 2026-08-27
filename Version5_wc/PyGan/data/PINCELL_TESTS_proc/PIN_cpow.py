###########################################################################
#                                                                         #
# Description : PyGan script for BWR simulation with DRAGON5              #
# Author      : R. Guasch                                                 #
# Date        : 05/08/2026                                                #
# Purpose     : Depletion calculation for simplified PINCELLS             #
#                                                                         #
###########################################################################
#
import lifo
import cle2000

def cpow_EXTR(case_name, namCOMPO, pyLIB, pyTRACK, pyTF_EXC, pyTRACK_SS, pyTF_EXC_SS, StepList, name_compo):

    # Lifo
    myLifo=lifo.new()
    myLifo.pushEmpty(namCOMPO,"LCM")
    myLifo.push(pyLIB)
    myLifo.push(pyTRACK)
    myLifo.push(pyTF_EXC)
    myLifo.push(pyTRACK_SS)
    myLifo.push(pyTF_EXC_SS)
    myLifo.push(StepList)
    myLifo.push(name_compo)
    myLifo.lib()

    # Execution 
    if case_name == "UO2_BWR_PINCELL":
        calcBWR = cle2000.new('evolUOX',myLifo,1)
    elif case_name == "GD_BWR_PINCELL":
        calcBWR = cle2000.new('evolGD',myLifo,1)
    calcBWR.exec()

    # Recover
    myLifo.lib()
    pyCOMPO = myLifo.node(namCOMPO)
    StepList.close()

    # Clear stack before next execution
    while myLifo.getMax() > 0: 
        myLifo.pop()
        
    return pyCOMPO

def cpow_CELI(namCOMPO, pyLIB, StepList, name_compo, glob_opt):

  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty(namCOMPO,"LCM")
  myLifo.push(pyLIB)
  myLifo.push(StepList)
  myLifo.push(name_compo)
  myLifo.push(glob_opt)
  myLifo.lib()

  # Execution 
  calcBWR = cle2000.new('pcc1',myLifo,1)
  calcBWR.exec()

  # Recover
  myLifo.lib()
  pyCOMPO = myLifo.node(namCOMPO)
  StepList.close()

  # Clear stack before next execution
  while myLifo.getMax() > 0: 
      myLifo.pop()
      
  return pyCOMPO

def flux_295g(namCOMPO_295, pyLIB, name_compo):
    
    # Lifo
    myLifo=lifo.new()
    myLifo.pushEmpty(namCOMPO_295,"LCM")
    myLifo.push(pyLIB)
    myLifo.pushEmpty("LIBRARY2","LCM")
    myLifo.push(name_compo)
    myLifo.lib()

    # Execution 
    calcBWR = cle2000.new('flux',myLifo,1)
    calcBWR.exec()

    # Recover
    myLifo.lib()
    pyCOMPO_295 = myLifo.node(namCOMPO_295)
    pyLIB2_295 = myLifo.node("LIBRARY2")

    # Clear stack before next execution
    while myLifo.getMax() > 0: 
        myLifo.pop()
        
    return pyCOMPO_295, pyLIB2_295