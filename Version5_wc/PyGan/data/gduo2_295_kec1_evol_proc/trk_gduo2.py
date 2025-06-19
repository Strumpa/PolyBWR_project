#####################################################################
#                                                                   #
# Purpose     : Tracking for pin cells with SALT module             #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def trk_gduo2(pyGEOM):
  # Lifo
  myLifo=lifo.new()

  myLifo.pushEmpty("UOX_TBH", "ASCII")
  myLifo.pushEmpty("TRACK", "LCM")
  myLifo.pushEmpty("TF_EXC", "BINARY")

  myLifo.pushEmpty("TRACK_SS", "LCM")

  myLifo.push(pyGEOM)

  # Execution
  trk_gduo2_proc = cle2000.new('trk_gduo2',myLifo,1)
  trk_gduo2_proc.exec()

  # Recover
  myLifo.lib()

  pyUOX_TBH = myLifo.node("UOX_TBH")
  pyTRACK = myLifo.node("TRACK")
  pyTF_EXC = myLifo.node("TF_EXC")
  pyTRACK_SS = myLifo.node("TRACK_SS")

  return pyUOX_TBH, pyTRACK, pyTF_EXC, pyTRACK_SS







