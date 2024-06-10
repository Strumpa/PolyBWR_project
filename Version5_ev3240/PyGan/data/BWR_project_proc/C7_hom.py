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

def C7_hom(namCOMPO,StepList, name_compo):

  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty(namCOMPO,"LCM")
  myLifo.push(StepList)
  myLifo.push(name_compo)
  myLifo.lib()

  # Execution 
  calcBWR = cle2000.new('C7_hom',myLifo,1)
  calcBWR.exec()

  # Recover
  myLifo.lib()
  pyCOMPO = myLifo.node(namCOMPO)
  return pyCOMPO
