#####################################################################
#                                                                   #
# Description : PyGan scritp for BWR simulation with DRAGON5        #
# Author      : L. Fede, used by R.Guasch                           #
# Date        : 2023                                                #
# Purpose     : MIX definition for pin cells                        #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def MIX_C(Library,ssh_option,name_geom,setQFIS=False,setNGtozero=False):
  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty("LIBRARY", "LCM")
  myLifo.push(Library)
  myLifo.push(ssh_option)
  myLifo.lib()

  # Execution
  if name_geom == "AT10_24UOX":
    mixBWR = cle2000.new('MIX_24UOX',myLifo,1)
  elif name_geom == "AT10_45Gd":
    if setQFIS and setNGtozero:
      mixBWR = cle2000.new('MIX_45GdQNG0',myLifo,1)
    elif setQFIS:
      mixBWR = cle2000.new('MIX_45GdQ',myLifo,1)
    else:
      mixBWR = cle2000.new('MIX_45Gd',myLifo,1)
  mixBWR.exec()

  # Recover
  myLifo.lib()
  pyMIX = myLifo.node("LIBRARY")

  # Clear stack before next execution
  while myLifo.getMax() > 0:
    myLifo.pop()

  return pyMIX
