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

def MIX_C(Library,ssh_option):
  # Lifo
  namLIB = "LIBRARY"
  myLifo=lifo.new()
  myLifo.pushEmpty(namLIB, "LCM")
  myLifo.push(Library)
  myLifo.push(ssh_option)
  myLifo.lib()

  # Execution
  mixBWR = cle2000.new('MIX_C',myLifo,1)
  mixBWR.exec()

  # Recover
  myLifo.lib()
  pyMIX = myLifo.node(namLIB)

  # Clear stack before next execution
  while myLifo.getMax() > 0: 
        myLifo.pop()
  
  return pyMIX
