#############################################################################################
#                                                                                           #
# Description : PyGan scritp for BWR simulation with DRAGON5                                #
# Author      : L. Fede, used by R.Guasch                                                   #
# Date        : 2023                                                                        #
# Purpose     : MIX definition for pin cells                                                #
# Modified DEPL structure to exclude (n,gamma) reactions from constant power normalization  #
#############################################################################################
#
import lifo
import cle2000

def MIX_NG0(Library):
  """
  Due to length of DEPL structure : only RSE + CORR ssh_option is supported here.
  """
  # Lifo
  myLifo=lifo.new()
  myLifo.pushEmpty("LIBRARY", "LCM")
  myLifo.push(Library)
  myLifo.lib()

  # Execution
  mixNG0 = cle2000.new('MIX_NG0',myLifo,1)
  mixNG0.exec()

  # Recover
  myLifo.lib()
  pyMIX = myLifo.node("LIBRARY")

  # Clear stack before next execution
  while myLifo.getMax() > 0:
    myLifo.pop()

  return pyMIX
