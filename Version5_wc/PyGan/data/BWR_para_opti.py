##################################################################################
#                                                                                #
# Description : PyGan scritp for parallelization test on BWR bencharks           #
# Author      : R. Guasch					                                     #
# Date        : 2024                                                             #
# Purpose     : Find optimal parallelization options for BWR benchmark cases     #
#                                                                                #
##################################################################################
#
########################################################################################################################################################################################
#                                #
#             IMPORTS            #
#                                #
##################################
import sys
import numpy as np
import os, shutil
import lifo
import lcm
import cle2000
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from assertS import *
# POST-PROCESSING class

# --- CELLS

# --- CLUSTERS of CELLS
from UOX_2x2 import *

# --- ASSEMBLIES



########################################################################################################################################################################################
#                                                        #
#   PARAMETERS OF FLUX CALCULATION : TO BE MODIFIED      #
#                                                        #
##########################################################
#
# Library = name of the used Library
Library = 'J311_295'
# 
#
# name_geom = name of the geometry for the simulated case 
# - cells :
#           - 'AT10_24UOX'
#           - 'AT10_32UOX'
#           - 'AT10_42UOX'
#           - 'AT10_45UOX'
#           - 'AT10_48UOX'
#           - 'AT10_50UOX'
#           - 'AT10_45Gd'
#           - 'AT10_42Gd'
# - clusters :
# 		 	- 'AT10_2x2_UOX'
# 		 	- 'AT10_2x2_UOX_Gd'
# 		 	- 'AT10_3x3_UOX_Gd'
#
# - assemblies : 
#           - 'ATRIUM-10XM_NL24-2'
#           - 'ATRIUM-10XM_NL24-2_ctrl'
case = "AT10_2x2_UOX"


tracking_module = "SALT" # "NXT"
tracking_ssh_module = "SALT" # "NXT"
#
#
#


# main flux tracking parameters
an2d = [12, 14, 18, 20, 24, 30]  # number of azimuthal angles
densur = [30.0, 35.0, 40.0, 45.0, 50.0] # density of integration lines (line/cm)

# self sielding tracking parameters
an2d_ss = [6, 8, 12] # number of azimuthal angles
densur_ss = [10.0, 15.0, 25.0] # density of integration lines (line/cm)

# parallelisation parameters
batch_flx = [50, 100, 150, 200, 250, 300] # number of batches for flux tracking
batch_ssh = [50, 100, 150, 200, 250, 300] # number of batches for self shielding tracking



########################################################################################################################################################################################
#                                                 #
#             EXECUTION OF PROCEDURES             #
#                                                 #
###################################################

if case == 'AT10_2x2_UOX':
	pyKinf = UOX_2x2(batch_flx[1], batch_ssh[1], an2d[0], an2d_ss[0], densur[0], densur_ss[0])
    
	

