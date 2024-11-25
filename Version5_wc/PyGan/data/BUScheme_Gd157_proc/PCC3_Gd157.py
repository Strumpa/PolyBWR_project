#####################################################################
#                                                                   #
# Description : PyGan script for BWR simulation with DRAGON5        #
# Author      : R. Guasch                                           #
# Date        : 2024-11-21                                          #
# Purpose     : Depletion calculation for homogeneous for pin cell  #
#               Implementing predictor-corrector method             #
#####################################################################
#
import lifo
import cle2000

def PCC3_Gd157(namCOMPO, pyMIX, StepList, name_compo, ssh_module):

  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty(namCOMPO,"LCM")
  myLifo.push(pyMIX)
  myLifo.push(StepList)
  myLifo.push(name_compo)
  myLifo.push(ssh_module)
  myLifo.lib()

  # Execution 
  calcBWR = cle2000.new('PCC3_Gd157',myLifo,1)
  calcBWR.exec()

  # Recover
  myLifo.lib()
  pyCOMPO = myLifo.node(namCOMPO)
  StepList.close()
  return pyCOMPO
