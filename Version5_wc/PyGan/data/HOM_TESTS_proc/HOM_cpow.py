###########################################################################
#                                                                         #
# Description : PyGan script for BWR simulation with DRAGON5              #
# Author      : R. Guasch                                                 #
# Date        : 14/05/2026                                                #
# Purpose     : Depletion calculation for homogenized pincells            #
#                                                                         #
###########################################################################
#
import lifo
import cle2000

def HOM_cpow_EXTR(namCOMPO, pyLIB, StepList, name_compo, glob_opt):

  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty(namCOMPO,"LCM")
  myLifo.push(pyLIB)
  myLifo.push(StepList)
  myLifo.push(name_compo)
  myLifo.push(glob_opt)
  myLifo.lib()

  # Execution 
  calcBWR = cle2000.new('evol',myLifo,1)
  calcBWR.exec()

  # Recover
  myLifo.lib()
  pyCOMPO = myLifo.node(namCOMPO)
  StepList.close()

  # Clear stack before next execution
  while myLifo.getMax() > 0: 
      myLifo.pop()
      
  return pyCOMPO

def HOM_cpow_CELI(namCOMPO, pyLIB, StepList, name_compo, glob_opt):

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