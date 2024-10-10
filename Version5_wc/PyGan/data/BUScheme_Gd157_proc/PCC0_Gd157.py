#####################################################################
#                                                                   #
# Description : PyGan scritp for BWR simulation with DRAGON5        #
# Author      : R. Guasch                                           #
# Date        : 2023/2024                                           #
# Purpose     : Depletion calculation for homogeneous for pin cell  #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def PCC0_Gd157(namCOMPO, StepList, name_compo, ssh_module, ssh_method, sat, depl_sol):

  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty(namCOMPO,"LCM")
  myLifo.push(StepList)
  myLifo.push(name_compo)
  myLifo.push(ssh_module)
  myLifo.push(ssh_method)
  myLifo.push(sat)
  myLifo.push(depl_sol)
  myLifo.lib()

  # Execution 
  calcBWR = cle2000.new('PCC0_Gd157',myLifo,1)
  calcBWR.exec()

  # Recover
  myLifo.lib()
  pyCOMPO = myLifo.node(namCOMPO)
  StepList.close()
  return pyCOMPO
