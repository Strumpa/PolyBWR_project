#####################################################################
#                                                                   #
# Description : PyGan scritp for VVER simulation with DRAGON5       #
# Author      : L. Fede                                             #
# Date        : 2023                                                #
# Purpose     : Depletion calculus for NATIVE assemblies            #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def BU_A_N(namCOMPO,pyMIX,pyTRACK_1L,pyTRACK_SS,StepList,Multicompo,name_geom,name_compo):

  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty(namCOMPO,"LCM")
  myLifo.push(pyMIX)
  myLifo.push(pyTRACK_1L)
  myLifo.push(pyTRACK_SS)
  myLifo.push(StepList)
  myLifo.push(Multicompo)
  myLifo.push(name_geom)
  myLifo.push(name_compo)
  myLifo.lib()

  # Execution 
  calcVVER = cle2000.new('BU_A_N',myLifo,1)
  calcVVER.exec()

  # Recover
  myLifo.lib()
  pyCOMPO = myLifo.node(namCOMPO)
  return pyCOMPO
