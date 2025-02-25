#####################################################################
#                                                                   #
# Description : PyGan scritp for BWR simulation with DRAGON5        #
# Author      : R.Guasch                                            #
# Date        : 2024                                                #
# Purpose     : Tracking for pin cells with SYB+NXT modules         #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def TRK_C_SYBNXT(pyGEOM,pyGEOM_SS):
  # Lifo
  myLifo=lifo.new()

  myLifo.push(pyGEOM)
  myLifo.push(pyGEOM_SS)

  # Execution
  trackBWR = cle2000.new('TRK_C_SYBNXT',myLifo,1)
  trackBWR.exec()

  # Recover nodes
  myLifo.lib()
  pyTRACK = myLifo.node("TRACK")
  pyTF_EXC = myLifo.node("TF_EXC")
  pyTRACK_SS = myLifo.node("TRACK_SS")

  # Clear stack before next execution
  while myLifo.getMax() > 0: 
        myLifo.pop()

  return pyTRACK,pyTF_EXC,pyTRACK_SS



