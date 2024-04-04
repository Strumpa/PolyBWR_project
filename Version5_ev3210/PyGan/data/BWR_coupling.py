##################################################################################
#                                                                                #
# Description : PyGan scritp for BWR simulation with DRAGON5                     #
# Author      : R. Guasch, C. Huet                                               #
# Date        : 2024                                                             #
# Purpose     : Perform coupled Dragon/Donjon:THM calculations for a BWR pincell #
#                                                                                #
##################################################################################
#
########################################################################################################################################################################################
#                                #
#             IMPORTS            #
#                                #
##################################
import os, shutil, sys
import lifo
import lcm
import cle2000
import numpy as np
import matplotlib
matplotlib.use('Agg')
#from matplotlib import pyplot
import matplotlib.pyplot as plt
from assertS import *
# POST-PROCESSING class
#from POSTPROC import POSTPROC

# --- CELL
# GEOMETRY class
from GEO_C_SALT import *
from GEO_C_NXT import *
# MIXTURES class
from MIX_C_THM import *
# TRACKING class
from TRK_C_SALT import *
from TRK_C_SYBNXT import *
# BU DEPLETION class
from FLX_C import *

# --- OTHERS
from getLists import *


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
# TYPE = type of calculus
#           - 'ALAMOS'
#           - 'HYBRIDE'
#           - 'NATIVE'
TYPE = 'NATIVE'
#
# name_geom = 'AT10_UOX' or 'AT10_UOX_Gd2O3' of the simulated case, different discretization for UOX and UOX+Gd2O3 cells.
# - cells = name_mix:
#           - 'AT10_24UOX'
#           - 'AT10_32UOX'
#           - 'AT10_42UOX'
#           - 'AT10_45UOX'
#           - 'AT10_48UOX'
#           - 'AT10_50UOX'
#           - 'AT10_45Gd'
#           - 'AT10_42Gd'
# - assemblies : 
#           - 'ATRIUM-10XM_NL24-2'
#           - 'ATRIUM-10XM NL24-2_ctrl'
name_geom = 'AT10_UOX' # 'AT10_UOX_Gd2O3'
name_mix = "AT10_24UOX"
tracking_module = "SALT"

case="CELL"
#
# Multicompo = 1 if you want to generate a MULTICOMPO object 
#            = 0 either
Multicompo = 1
# Initialize pin (normalization factor) power :
pin_power = 38.6 # MW/t
# Initialize moderator relative density, Fuel and coolant temperature :
mode_rel_dens = 1.0 # = current moderator density / nominal moderator density.
Tfuel = 750.0 # K, it is assumed in MIX_C.c2m that Tgap = Tfuel
Tcool = 559.0 # K, it is assumed in MIX_C.c2m that Tclad = Tcool, we sould check the validity of this claim 
# --> could use Python protoype from ENE6107A project with modified correlation for drift flux model in moderator before modifying THM.
# suffix = suffix added to name_geom for creation of figures, MULTICOMPO and BU vector
suffix = "_"+tracking_module+"_THM"
#
# names for exportation
name_compo = './_COMPO_'+name_geom+'_'+name_mix+suffix
name_BUvector = name_geom+'_'+name_mix+suffix
name_fig = name_geom+'_'+name_mix+suffix+'_2L.ps'
name_fig_SS = name_geom+'_'+name_mix+suffix+'_SS.ps'
name_fig_1L = name_geom+'_'+name_mix+suffix+'_1L.ps'
namFIG_SS_NXT=name_geom+'_'+name_mix+suffix+'_SS_NXT.ps'
namFIG_1L_NXT=name_geom+'_'+name_mix+suffix+'_1L_NXT.ps'

#
# Creation of results directory
path=os.getcwd()

a=os.path.exists('BWR_Coupling_PyGan_'+name_geom)
if a==False:
	os.mkdir('BWR_Coupling_PyGan_'+name_geom)

SAVE_DIR='BWR_Coupling_PyGan_'+name_geom+'/'+name_mix+'_'+suffix+'_postprocess'
a=os.path.exists(SAVE_DIR)
if a==False:
	os.mkdir(SAVE_DIR)


###################################################
#                                                 #
#             EXECUTION OF PROCEDURES             #
#                                                 #
###################################################
#
# --------------------------------------
#                 CELL                 
# --------------------------------------
#
if case == 'CELL' :
    if tracking_module == 'SALT':

            # ---------- GEOMETRY
            [pyGEOM,pyGEOM_SS] = GEO_C_SALT("GEOM","GEOM_SS",name_geom)

            # ---------- MIXTURES
            pyMIX = MIX_C_THM("LIBRARY",Library,name_mix,mode_rel_dens,Tfuel,Tcool) # we want to find a way to recover 

            # ---------- TRACKING
            [pyUOX_TBH,pyFIG,pyTRACK,pyTF_EXC,pyUOX_TBH_SS,pyFIG_SS,pyTRACK_SS,pyTF_EXC_SS] = TRK_C_SALT("UOX_TBH",name_fig,"TRACK","TF_EXC","UOX_TBH_SS",name_fig_SS,"TRACK_SS","TF_EXC_SS",pyGEOM,pyGEOM_SS)

            # ---------- FLUX CALCULATION FOR 
            #pyCOMPO = BU_C("COMPO",pyMIX,pyTRACK,pyTF_EXC,pyTRACK_SS,pyTF_EXC_SS,StepList,Multicompo,name_geom,name_compo)

    elif tracking_module == "SYBNXT":

                        # ---------- GEOMETRY
            [pyGEOM,pyGEOM_SS] = GEO_C_NXT("GEOM","GEOM_SS",name_geom)

            # ---------- MIXTURES
            pyMIX = MIX_C_THM("LIBRARY",Library,name_mix,mode_rel_dens,Tfuel,Tcool) # we want to find a way to recover

            # ---------- TRACKING
            [pyUOX_TBH,pyFIG,pyTRACK,pyTF_EXC,pyUOX_TBH_SS,pyFIG_SS,pyTRACK_SS,pyTF_EXC_SS] = TRK_C_SYBNXT("UOX_TBH",name_fig,"TRACK","TF_EXC","UOX_TBH_SS",name_fig_SS,"TRACK_SS","TF_EXC_SS",pyGEOM,pyGEOM_SS)

            # ---------- BU DEPLETION
            pyCOMPO = BU_C("COMPO",pyMIX,pyTRACK,pyTF_EXC,pyTRACK_SS,pyTF_EXC_SS,StepList,Multicompo,name_geom,name_compo)
