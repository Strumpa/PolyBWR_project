#####################################################################
#                                                                   #
# Description : PyGan scritp for BWR simulation with DRAGON5        #
# Author      : R. Guasch                                           #
# Date        : 2023/2024                                           #
# Purpose     : Depletion calculation for 2x2 AT10 Cluster          #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def BU_UOX(namCOMPO, pyMIX, pyTRACK, pyTF_EXC, pyTRACK_SS, pyTF_EXC_SS, StepList, Multicompo, name_compo, flx_sol, CALC_opt):

  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty(namCOMPO,"LCM")
  myLifo.push(pyMIX)
  myLifo.push(pyTRACK)
  myLifo.push(pyTF_EXC)
  myLifo.push(pyTRACK_SS)
  myLifo.push(pyTF_EXC_SS)
  myLifo.push(StepList)
  myLifo.push(Multicompo)
  myLifo.push(name_compo)
  myLifo.push(flx_sol)
  myLifo.push(CALC_opt)

  myLifo.lib()

  # Execution 
  calcBWR = cle2000.new('BU_UOX',myLifo,1)
  calcBWR.exec()

  # Recover
  myLifo.lib()
  pyCOMPO = myLifo.node(namCOMPO)
  return pyCOMPO
