#####################################################################
#                                                                   #
# Purpose     : LIBRARY definition for gduo2_kec1 benchmark cell    #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def lib_gduo2(evaluation_name, correlation):
  """
  evaluation_name : str, choice between 'J311_295' or 'ENDFb8r1_295' or 'E8R1295K' or 'J311295K' (same but with KERMA) to create the library with the corresponding evaluation.
  correlation : str, choice between 'CORR' and 'NOCORR' to create the library with correlated resonant effects between Gd157 and U238 or not.
  """
  # Lifo
  namLIB = "LIBRARY"
  myLifo=lifo.new()
  myLifo.pushEmpty(namLIB, "LCM")
  myLifo.push(evaluation_name)
  myLifo.push(correlation)
  myLifo.lib()

  # Execution
  lib_gduo2_proc = cle2000.new('lib_gduo2',myLifo,1)
  lib_gduo2_proc.exec()

  # Recover
  myLifo.lib()
  pyLIB = myLifo.node(namLIB)
  return pyLIB
