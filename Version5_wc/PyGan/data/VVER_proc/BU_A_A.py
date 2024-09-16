#####################################################################
#                                                                   #
# Description : PyGan scritp for VVER simulation with DRAGON5       #
# Author      : L. Fede                                             #
# Date        : 2023                                                #
# Purpose     : Depletion calculus for ALAMOS assemblies            #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def BU_A_A(namCOMPO,pyMIX,pyTRACK,pyTF_EXC,pyTRACK_SS,pyTF_EXC_SS,StepList,Multicompo,name_geom,name_compo,AUTOP):

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
  myLifo.push(name_geom)
  myLifo.push(name_compo)
  myLifo.push(AUTOP)
  myLifo.lib()

  # Execution 
  calcVVER = cle2000.new('BU_A_A',myLifo,1)
  calcVVER.exec()

  # Recover
  myLifo.lib()
  pyCOMPO = myLifo.node(namCOMPO)
  return pyCOMPO
