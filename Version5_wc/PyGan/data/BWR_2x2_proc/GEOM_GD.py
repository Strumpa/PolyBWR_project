#####################################################################
#                                                                   #
# Description : PyGan scritp for BWR simulation with DRAGON5        #
# Author      : R.Guasch                                            #
# Date        : 29 July 2024                                        #
# Purpose     : Geometry definition for 2x2 Cells UOX+Gd cluster    #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def GEOM_GD(namGEOM,namGEOM_SS,name_geom):

  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty(namGEOM, "LCM")
  myLifo.pushEmpty(namGEOM_SS, "LCM")
  myLifo.push(name_geom)
  myLifo.lib()

  # Execution 
  geoBWR = cle2000.new('GEOM_GD',myLifo,1)
  geoBWR.exec()

  # Recover
  myLifo.lib()
  pyGEOM = myLifo.node(namGEOM)
  pyGEOM_SS = myLifo.node(namGEOM_SS)
  return [pyGEOM,pyGEOM_SS]