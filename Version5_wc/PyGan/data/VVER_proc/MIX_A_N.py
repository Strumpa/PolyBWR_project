#####################################################################
#                                                                   #
# Description : PyGan scritp for VVER simulation with DRAGON5       #
# Author      : L. Fede                                             #
# Date        : 2023                                                #
# Purpose     : Geometry definition for NATIVE assemblies           #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def MIX_A_N(namMIX,Library,name_geom):
  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty(namMIX, "LCM")
  myLifo.push(Library)
  myLifo.push(name_geom)
  #myLifo.lib()

  # Execution
  mixVVER = cle2000.new('MIX_A_N',myLifo,1)
  mixVVER.exec()

  # Recover
  #myLifo.lib()
  pyMIX = myLifo.node(namMIX)
  return pyMIX
