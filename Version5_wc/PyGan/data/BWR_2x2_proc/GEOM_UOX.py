#####################################################################
#                                                                   #
# Description : PyGan scritp for BWR simulation with DRAGON5        #
# Author      : R.Guasch                                            #
# Date        : 29 July 2024                                        #
# Purpose     : Geometry definition for 2x2 Cells UOX cluster       #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def GEOM_UOX(namGEOM,namGEOM_SS):

  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty(namGEOM, "LCM")
  myLifo.pushEmpty(namGEOM_SS, "LCM")
  myLifo.lib()

  # Execution 
  geoBWR = cle2000.new('GEOM_UOX',myLifo,1)
  geoBWR.exec()

  # Recover
  myLifo.lib()
  pyGEOM = myLifo.node(namGEOM)
  pyGEOM_SS = myLifo.node(namGEOM_SS)
  return [pyGEOM,pyGEOM_SS]