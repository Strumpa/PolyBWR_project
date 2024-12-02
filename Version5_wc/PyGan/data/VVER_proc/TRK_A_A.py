#####################################################################
#                                                                   #
# Description : PyGan scritp for VVER simulation with DRAGON5       #
# Author      : L. Fede                                             #
# Date        : 2023                                                #
# Purpose     : Tracking for ALAMOS assemblies                      #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def TRK_A_A(namTRACK,namTF_EXC,namTRACK_SS,namTF_EXC_SS,namFIG,namFIG_SS,name_alamos_MOC_ap3,name_alamos_MOC_za,name_alamos_SSH_ap3,name_alamos_SSH_za):

  # Lifo
  myLifo=lifo.new()


  myLifo.pushEmpty(namTRACK, "LCM")
  myLifo.pushEmpty(namTF_EXC, "BINARY")
  myLifo.pushEmpty(namTRACK_SS, "LCM")
  myLifo.pushEmpty(namTF_EXC_SS, "BINARY")
  myLifo.pushEmpty(namFIG, "ASCII")
  myLifo.pushEmpty(namFIG_SS, "ASCII")

  myLifo.push(name_alamos_MOC_ap3)
  myLifo.push(name_alamos_MOC_za)
  myLifo.push(name_alamos_SSH_ap3)
  myLifo.push(name_alamos_SSH_za)

  # Execution
  trackVVER = cle2000.new('TRK_A_A',myLifo,1)
  trackVVER.exec()

  # Recover
  myLifo.lib()

  pyTRACK = myLifo.node(namTRACK)
  pyTF_EXC = myLifo.node(namTF_EXC)
  pyTRACK_SS = myLifo.node(namTRACK_SS)
  pyTF_EXC_SS = myLifo.node(namTF_EXC_SS)
  pyFIG = myLifo.node(namFIG)
  pyFIG_SS = myLifo.node(namFIG_SS)
  
  return [pyTRACK,pyTF_EXC,pyTRACK_SS,pyTF_EXC_SS,pyFIG,pyFIG_SS]







