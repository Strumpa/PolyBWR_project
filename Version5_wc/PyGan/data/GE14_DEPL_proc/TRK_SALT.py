#####################################################################
#                                                                   #
# Description : PyGan script for BWR simulation with DRAGON5        #
# Author      : R.Guasch                                            #
# Date        : 11/08/2026                                          #
# Purpose     : Tracking of BWR assemblies with the SALT module     #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def TRK(tdt_file_ssh, tdt_file_l1, tdt_file_l2):
    # Lifo
    myLifo=lifo.new()
    myLifo.pushEmpty("TRKSSH","LCM")
    myLifo.pushEmpty("TRKFILSSH","BINARY")
    myLifo.pushEmpty("TRKFLUXL1","LCM")
    myLifo.pushEmpty("TRKFILFLUXL1","BINARY")
    myLifo.pushEmpty("TRKFLUXL2","LCM")
    myLifo.pushEmpty("TRKFILFLUXL2","BINARY")
    myLifo.push(tdt_file_ssh)
    myLifo.push(tdt_file_l1)
    myLifo.push(tdt_file_l2)

    # Execution
    trackBWR = cle2000.new('TRK',myLifo,1)
    trackBWR.exec()

    # Recover
    myLifo.lib()

    trk_ssh = myLifo.node("TRKSSH")
    trkfil_ssh = myLifo.node("TRKFILSSH")
    trk_L1 = myLifo.node("TRKFLUXL1")
    trkfil_L1 = myLifo.node("TRKFILFLUXL1")
    trk_L2 = myLifo.node("TRKFLUXL2")
    trkfil_L2 = myLifo.node("TRKFILFLUXL2")

    # Clear stack before next execution
    while myLifo.getMax() > 0: 
            myLifo.pop()

    return trk_ssh, trkfil_ssh, trk_L1, trkfil_L1, trk_L2, trkfil_L2