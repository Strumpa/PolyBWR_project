#####################################################################
#                                                                   #
# Description : Python/PyGan script for treating 2x2 UOX cluster    #
# Author      : R. Guasch                                           #
# Date        : 2024                                                #
# Purpose     : Test tracking and paralleliation params             #
#                                                                   #
#####################################################################
#
import lifo
import cle2000

def UOX_2x2(batch_flx, batch_ssh, an2d, an2d_ss, densur, densur_ss):

    # Lifo
    myLifo=lifo.new()
    myLifo.push(batch_flx)
    myLifo.push(batch_ssh)
    myLifo.push(an2d)
    myLifo.push(an2d_ss)
    myLifo.push(densur)
    myLifo.push(densur_ss)
    myLifo.lib()

    # Execution 
    AT10_2x2_UOX = cle2000.new('UOX_2x2',myLifo,1)
    AT10_2x2_UOX.exec()

    # Recover
    myLifo.lib()
    pyKinf = myLifo.node('Kinf')
    return pyKinf
