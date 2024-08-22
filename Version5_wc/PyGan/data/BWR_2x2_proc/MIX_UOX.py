#####################################################################
#                                                                   #
# Description : PyGan scritp for BWR simulation with DRAGON5        #
# Author      : R.Guasch                                            #
# Date        : 29 July 2024                                        #
# Purpose     : MIX definition for 2x2 cluster of pins              #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def MIX_UOX(namMIX, Library, ssh_option):
  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty(namMIX, "LCM")
  myLifo.push(Library)
  myLifo.push(ssh_option)
  myLifo.lib()

  # Execution
  mixBWR = cle2000.new('MIX_UOX',myLifo,1)
  mixBWR.exec()

  # Recover
  myLifo.lib()
  pyMIX = myLifo.node(namMIX)
  return pyMIX
