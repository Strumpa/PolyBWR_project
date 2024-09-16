#####################################################################
#                                                                   #
# Description : PyGan scritp for VVER simulation with DRAGON5       #
# Author      : L. Fede                                             #
# Date        : 2023                                                #
# Purpose     : Geometry definition for REL2005 assemblies          #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def GEO_A_H(namGEOM_1L,namGEOM_SS,namGEOM_1L_NXT,namGEOM_SS_NXT,name_geom):

  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty(namGEOM_1L, "LCM")
  myLifo.pushEmpty(namGEOM_SS, "LCM")
  myLifo.pushEmpty(namGEOM_1L_NXT, "LCM")
  myLifo.pushEmpty(namGEOM_SS_NXT, "LCM")
  myLifo.push(name_geom)
  myLifo.lib()

  # Execution 
  geoVVER = cle2000.new('GEO_A_H',myLifo,1)
  geoVVER.exec()

  # Recover
  myLifo.lib()
  pyGEOM_1L = myLifo.node(namGEOM_1L)
  pyGEOM_SS = myLifo.node(namGEOM_SS)
  pyGEOM_1L_NXT = myLifo.node(namGEOM_1L_NXT)
  pyGEOM_SS_NXT = myLifo.node(namGEOM_SS_NXT)
  return [pyGEOM_1L,pyGEOM_SS,pyGEOM_1L_NXT,pyGEOM_SS_NXT]
