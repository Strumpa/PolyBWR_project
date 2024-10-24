#####################################################################
#                                                                   #
# Description : PyGan scritp for BWR simulation with DRAGON5        #
# Author      : L. Fede, adapted by R. Guasch                       #
# Date        : 2023/2024                                           #
# Purpose     : Depletion calculus for pin cell                     #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def BU_C_SYBNXT(namCOMPO, pyMIX, pyTRACK, pyTF_EXC, pyTRACK_SS, StepList, Multicompo, name_geom, name_compo, ssh_option, ssh_sol, flx_sol):

  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty(namCOMPO,"LCM")
  myLifo.push(pyMIX)
  myLifo.push(pyTRACK)
  myLifo.push(pyTF_EXC)
  myLifo.push(pyTRACK_SS)
  myLifo.push(StepList)
  myLifo.push(Multicompo)
  myLifo.push(name_geom)
  myLifo.push(name_compo)
  myLifo.push(ssh_option)
  myLifo.push(ssh_sol)
  myLifo.push(flx_sol)
  myLifo.lib()

  # Execution 
  calcBWR = cle2000.new('BU_C_SYBNXT',myLifo,1)
  calcBWR.exec()

  # Recover
  myLifo.lib()
  pyCOMPO = myLifo.node(namCOMPO)
  return pyCOMPO
