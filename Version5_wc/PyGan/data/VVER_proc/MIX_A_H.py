#####################################################################
#                                                                   #
# Description : PyGan scritp for VVER simulation with DRAGON5       #
# Author      : L. Fede                                             #
# Date        : 2023                                                #
# Purpose     : MIX definition for REL2005 assemblies          #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def MIX_A_H(namMIX,Library,name_geom):
  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty(namMIX, "LCM")
  myLifo.push(Library)
  myLifo.push(name_geom)
  #myLifo.lib()

  # Execution
  mixVVER = cle2000.new('MIX_A_H',myLifo,1)
  mixVVER.exec()

  # Recover
  #myLifo.lib()
  pyMIX = myLifo.node(namMIX)
  return pyMIX
