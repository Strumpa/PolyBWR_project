#####################################################################
#                                                                   #
# Purpose     : Depletion calculation for gduo2_kec1 cell           #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def bu_gduo2(namCOMPO, pyLIB, pyTRACK, pyTF_EXC, pyTRACK_SS,  StepList, name_compo):

  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty(namCOMPO,"LCM")
  myLifo.push(pyLIB)
  myLifo.push(pyTRACK)
  myLifo.push(pyTF_EXC)
  myLifo.push(pyTRACK_SS)
  myLifo.push(StepList)
  myLifo.push(name_compo)

  myLifo.lib()

  # Execution 
  calcBWR = cle2000.new('bu_gduo2',myLifo,1)
  calcBWR.exec()

  # Recover
  myLifo.lib()
  pyCOMPO = myLifo.node(namCOMPO)
  return pyCOMPO
