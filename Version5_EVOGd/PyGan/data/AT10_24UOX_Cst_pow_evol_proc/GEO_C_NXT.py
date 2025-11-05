#####################################################################
#                                                                   #
# Description : PyGan scritp for BWR simulation with DRAGON5        #
# Author      : L. Fede, adapted by R.Guasch                        #
# Date        : 2024                                                #
# Purpose     : Geometry definition for pin cells                   #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def GEO_C_NXT():

  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty("GEOM", "LCM")
  myLifo.pushEmpty("GEOM_SS", "LCM")
  myLifo.lib()

  # Execution 
  geoBWR = cle2000.new('GEO_C_NXT',myLifo,1)
  geoBWR.exec()

  # Recover
  myLifo.lib()
  pyGEOM = myLifo.node("GEOM")
  pyGEOM_SS = myLifo.node("GEOM_SS")

  # Clear stack before next execution
  while myLifo.getMax() > 0: 
        myLifo.pop()
        
  return pyGEOM,pyGEOM_SS
