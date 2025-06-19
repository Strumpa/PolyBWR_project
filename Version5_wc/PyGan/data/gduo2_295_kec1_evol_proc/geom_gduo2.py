#####################################################################
#
# Description     : Geometry definition for gduo2_kec1 benchmark cell            
#                                                                   
#####################################################################
#
import lifo
import cle2000

def geom_gduo2():
  """
  Geometry definition for gduo2_kec1 benchmark cell
  """
  namGEOM = "GEOM"
  # Lifo  
  myLifo=lifo.new()
  myLifo.pushEmpty(namGEOM, "LCM")
  myLifo.lib()

  # Execution 
  geo_gduo2_proc = cle2000.new('geom_gduo2',myLifo,1)
  geo_gduo2_proc.exec()

  # Recover
  myLifo.lib()
  pyGEOM = myLifo.node(namGEOM)
  return pyGEOM
