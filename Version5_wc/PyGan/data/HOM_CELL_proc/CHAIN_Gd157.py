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

def CHAIN_Gd157(namCOMPO, StepList, name_compo, ssh_module, sat, depl_sol, chain_modif):

  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty(namCOMPO,"LCM")
  myLifo.push(StepList)
  myLifo.push(name_compo)
  myLifo.push(ssh_module)
  myLifo.push(sat)
  myLifo.push(depl_sol)
  myLifo.push(chain_modif)
  myLifo.lib()

  # Execution 
  calcBWR = cle2000.new('CHAIN_Gd157',myLifo,1)
  calcBWR.exec()

  # Recover
  myLifo.lib()
  pyCOMPO = myLifo.node(namCOMPO)
  StepList.close()
  return pyCOMPO
