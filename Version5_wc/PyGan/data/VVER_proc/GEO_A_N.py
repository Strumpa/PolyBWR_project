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

def GEO_A_N(namGEOM_1L,namGEOM_SS,name_geom):

  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty(namGEOM_1L, "LCM")
  myLifo.pushEmpty(namGEOM_SS, "LCM")
  myLifo.push(name_geom)
  myLifo.lib()

  # Execution 
  geoVVER = cle2000.new('GEO_A_N',myLifo,1)
  geoVVER.exec()

  # Recover
  myLifo.lib()
  pyGEOM_1L = myLifo.node(namGEOM_1L)
  pyGEOM_SS = myLifo.node(namGEOM_SS)
  return [pyGEOM_1L,pyGEOM_SS]
