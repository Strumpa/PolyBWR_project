#####################################################################
#                                                                   #
# Description : PyGan scritp for VVER simulation with DRAGON5       #
# Author      : L. Fede                                             #
# Date        : 2023                                                #
# Purpose     : Tracking for NATIVE assemblies                      #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def TRK_A_N(namTRACK_1L,namTRACK_SS,namFIG_1L,namFIG_SS,pyGEOM_SS,pyGEOM_1L,an2d):

  # Lifo
  myLifo=lifo.new()


  myLifo.pushEmpty(namTRACK_1L, "LCM")
#  myLifo.pushEmpty(namTF_EXC, "BINARY")
  myLifo.pushEmpty(namTRACK_SS, "LCM")
  myLifo.pushEmpty(namFIG_1L, "ASCII")
  myLifo.pushEmpty(namFIG_SS, "ASCII")

  myLifo.push(pyGEOM_SS)
  myLifo.push(pyGEOM_1L)

  myLifo.push(an2d)

  myLifo.lib()

  # Execution
  trackVVER = cle2000.new('TRK_A_N',myLifo,1)
  trackVVER.exec()

  # Recover
  myLifo.lib()

  pyTRACK_1L = myLifo.node(namTRACK_1L)
#  pyTF_EXC = myLifo.node(namTF_EXC)
  pyTRACK_SS = myLifo.node(namTRACK_SS)
  pyFIG_1L = myLifo.node(namFIG_1L)
  pyFIG_SS = myLifo.node(namFIG_SS)
  
  return [pyTRACK_1L,pyTRACK_SS,pyFIG_1L,pyFIG_SS]







