#####################################################################
#                                                                   #
# Description : PyGan script for BWR simulation with DRAGON5        #
# Author      : R.Guasch                                            #
# Date        : 2025/02/25                                          #
# Purpose     : Geometry definition of BWR pincell                  #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def GEO_C_SALT():

    # Lifo
    myLifo=lifo.new()
    myLifo.pushEmpty("GEOM","LCM")
    myLifo.pushEmpty("GEOML1","LCM")
    myLifo.pushEmpty("GEOM_SS","LCM")
    myLifo.lib()

    # Execution 
    geoBWR = cle2000.new('GEO_C_SALT',myLifo,1)
    geoBWR.exec()

    # Recover
    myLifo.lib()
    pyGEOM = myLifo.node("GEOM")
    pyGEOML1 = myLifo.node("GEOML1")
    pyGEOM_SS = myLifo.node("GEOM_SS")
    
    # Clear stack before next execution
    while myLifo.getMax() > 0: 
        myLifo.pop()
    
    return pyGEOM, pyGEOML1, pyGEOM_SS
