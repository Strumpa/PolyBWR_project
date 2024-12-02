#####################################################################
#                                                                   #
# Description : PyGan scritp for VVER simulation with DRAGON5       #
# Author      : L. Fede                                             #
# Date        : 2023                                                #
# Purpose     : Tracking for REL2005 assemblies                     #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def TRK_A_H(namTRACK_2L,namTF_EXC,namTRACK_SS,namFIG_2L,namFIG_SS,namTRACK_1L,namFIG_1L,pyGEOM_SS,pyGEOM_1L,pyGEOM_SS_NXT,pyGEOM_1L_NXT,namFIG_SS_NXT,namFIG_1L_NXT,name_alamos_SSH_ap3,name_alamos_MOC_ap3,name_alamos_SSH_za,name_alamos_MOC_za,an2d):

  # Lifo
  myLifo=lifo.new()


  myLifo.pushEmpty(namTRACK_2L, "LCM")
  myLifo.pushEmpty(namTF_EXC, "BINARY")
  myLifo.pushEmpty(namTRACK_SS, "LCM")
  myLifo.pushEmpty(namFIG_2L, "ASCII")
  myLifo.pushEmpty(namFIG_SS, "ASCII")
  myLifo.pushEmpty(namTRACK_1L, "LCM")
  myLifo.pushEmpty(namFIG_1L, "ASCII")

  myLifo.push(pyGEOM_SS)
  myLifo.push(pyGEOM_1L)

  myLifo.push(pyGEOM_SS_NXT)
  myLifo.push(pyGEOM_1L_NXT)
  myLifo.pushEmpty(namFIG_SS_NXT, "ASCII")
  myLifo.pushEmpty(namFIG_1L_NXT, "ASCII")

  myLifo.push(name_alamos_SSH_ap3)
  myLifo.push(name_alamos_MOC_ap3)
  myLifo.push(name_alamos_SSH_za)
  myLifo.push(name_alamos_MOC_za)

  myLifo.push(an2d)

  myLifo.lib()

  # Execution
  trackVVER = cle2000.new('TRK_A_H',myLifo,1)
  trackVVER.exec()

  # Recover
  myLifo.lib()

  pyTRACK_2L = myLifo.node(namTRACK_2L)
  pyTF_EXC = myLifo.node(namTF_EXC)
  pyTRACK_SS = myLifo.node(namTRACK_SS)
  pyFIG_2L = myLifo.node(namFIG_2L)
  pyFIG_SS = myLifo.node(namFIG_SS)
  pyTRACK_1L = myLifo.node(namTRACK_1L)
  pyFIG_1L = myLifo.node(namFIG_1L)
  pyFIG_1L_NXT = myLifo.node(namFIG_1L_NXT)
  pyFIG_SS_NXT = myLifo.node(namFIG_SS_NXT)
  
  return [pyTRACK_2L,pyTF_EXC,pyTRACK_SS,pyFIG_2L,pyFIG_SS,pyTRACK_1L,pyFIG_1L,pyFIG_1L_NXT,pyFIG_SS_NXT]







