#####################################################################
#                                                                   #
# Description : PyGan script for Gd157 evolution study              #
# Author      : R.Guasch                                            #
# Date        : 2024/11/21                                          #
# Purpose     : MIX definition for HOM_UOX_Gd157 case               #
# Options : Library : - J311_281                                    #
#                     - J311_295                                    #
#                     - J311_315                                    #
#                     - J311_172 (based on Salino's PyNjoy2016)     #
#                                                                   #
#           ssh_option : - RSE/PT/SUBG (for AUTO select SUBG)       #
#####################################################################
#

import lifo
import cle2000

def MixHGd157(namMIX,Library,ssh_option,correlation):
  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty(namMIX, "LCM")
  myLifo.push(Library)
  myLifo.push(ssh_option)
  myLifo.push(correlation)
  myLifo.lib()

  # Execution
  mixHOM = cle2000.new('MixHGd157',myLifo,1)
  mixHOM.exec()

  # Recover
  myLifo.lib()
  pyMIX = myLifo.node(namMIX)
  return pyMIX
