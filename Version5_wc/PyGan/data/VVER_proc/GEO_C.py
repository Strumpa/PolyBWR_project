#####################################################################
#                                                                   #
# Description : PyGan scritp for VVER simulation with DRAGON5       #
# Author      : L. Fede                                             #
# Date        : 2023                                                #
# Purpose     : Geometry definition for pin cells                   #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def GEO_C(namGEOM,namGEOM_SS,name_geom):

  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty(namGEOM, "LCM")
  myLifo.pushEmpty(namGEOM_SS, "LCM")
  myLifo.push(name_geom)
  myLifo.lib()

  # Execution 
  geoVVER = cle2000.new('GEO_C',myLifo,1)
  geoVVER.exec()

  # Recover
  myLifo.lib()
  pyGEOM = myLifo.node(namGEOM)
  pyGEOM_SS = myLifo.node(namGEOM_SS)
  return [pyGEOM,pyGEOM_SS]
