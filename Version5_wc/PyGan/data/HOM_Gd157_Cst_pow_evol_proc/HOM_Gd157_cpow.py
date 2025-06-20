###########################################################################
#                                                                         #
# Description : PyGan scritp for BWR simulation with DRAGON5              #
# Author      : R. Guasch                                                 #
# Date        : 2023/2024                                                 #
# Purpose     : Depletion calculation for homogeneous Gd loaded pin cell  #
#                                                                         #
###########################################################################
#
import lifo
import cle2000

def HOM_Gd157_cpow(namCOMPO, pyLIB, StepList, name_compo, ssh_module, sat, depl_sol, rates_extr, DIRA_opt):

  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty(namCOMPO,"LCM")
  myLifo.push(pyLIB)
  myLifo.push(StepList)
  myLifo.push(name_compo)
  myLifo.push(ssh_module)
  myLifo.push(sat)
  myLifo.push(depl_sol)
  myLifo.push(rates_extr)
  myLifo.push(DIRA_opt)
  myLifo.lib()

  # Execution 
  calcBWR = cle2000.new('Gd157_cpow',myLifo,1)
  calcBWR.exec()

  # Recover
  myLifo.lib()
  pyCOMPO = myLifo.node(namCOMPO)
  StepList.close()

  # Clear stack before next execution
  while myLifo.getMax() > 0: 
      myLifo.pop()
      
  return pyCOMPO
