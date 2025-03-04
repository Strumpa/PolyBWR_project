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

def MIX_NG0(Library):
  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty("LIBRARY", "LCM")
  myLifo.push(Library)
  myLifo.lib()

  # Execution
  mixBWR = cle2000.new('MIX_NG0',myLifo,1)
  mixBWR.exec()

  # Recover
  myLifo.lib()
  pyMIX = myLifo.node("LIBRARY")

  # Clear stack before next execution
  while myLifo.getMax() > 0:
    myLifo.pop()

  return pyMIX
