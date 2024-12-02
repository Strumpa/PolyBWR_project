##################################################################################
#                                                                                #
# Description : PyGan scritp for BWR simulation with DRAGON5                     #
# Author      : R. Guasch                                                        #
# Date        : 2024                                                             #
# Purpose     : Perform lattice calculations for BWR 2x2 UOX/UOX+Gd cluster      #
#                                                                                #
##################################################################################

# Python3 imports
import sys
import numpy as np
import os, shutil
import lifo
import lcm
import cle2000
import matplotlib


# --- DRAGON5 procedures
from getLists import *
from GEOM_UOX import *
from GEOM_GD import *
from MIX_UOX import *
from MIX_GD import *
from TRK_SALT import *
# from TRK_NXT import *
# from TRK_SYBILT import *
#from FLU_2x2 import *
from BU_UOX import *
from BU_GD import *

# --- Post-Processing 

from POSTPROC_2x2 import *

#
# Visualisation options for post-processing functions
# 	= 1 if you want to plot and save the results
# 	= 0 either
# NB : if visu_XXX_Keff==1 or visu_XXX_BU==1 or visu_XXX_ISOTOPOESDENS==1 then visu_XXX must be equal to 1
visu_DRAGON=0
visu_SERPENT=0
visu_COMP=1
visu_DELTA=1

# ---------- GENERAL OPTIONS
# Multicompo = 1 if you want to generate a MULTICOMPO object 
#            = 0 either
Multicompo = 1

# Choice of parameters 
Library = 'J311_295'

case_name = "2x2_UOX" # "2x2_UOX_Gd" or "2x2_UOX"

# SSH method
ssh_method = "PT" # "PT", "RSE", "SUBG"
CALC_opt = "noCalc" # "CALC1", "CALC2", "CALC3", "CALC4", "noCalc"

# BU options 
#
# Burnup points
burnup_points = "UOx6_autop5"

# ---------- TRACKING
# -- Self-shielding options :

# ssh_trk = "SYBILT" # allows for SYBIL tracking of self-shielding geometry, allows for CP and IC solutons
# ssh_trk = "NXT" # allows for NXT: tracking of self-shielding geometry, allows for CP and MOC techniques (if MCCGT called)
ssh_trk = "SALT" # allows for SALT: tracking of self-shielding gometry, allows for CP and MOC techniques (if MCCGT called)
ssh_sol = "CP" # "CP" or "MOC" or "IC" 

# -- ssh tracking parameters :

if ssh_trk == "SALT" or ssh_trk == "NXT":
    an2d_ssh = 8 # number of azimuthal angles for self-shielding tracking
    densur_ssh = 15.0 # track density for self-shielding tracking
elif ssh_trk == "SYBILT":
    iqua2_ssh = 12
    nseg_ssh = 5

# -- Flux calculation options :

flx_trk = "SALT" # "NXT" or "SALT" or "SYBILT" : all allow for CP flux calculation, NXT and SALT allow for MOC flux calculation if MCCGT is used too
flx_sol = "MOC" # "CP" or "MOC" or "IC" allows for CP flux calculation

if flx_trk == "SALT" or flx_trk == "NXT":
    an2d = 18 # number of azimuthal angles for self-shielding tracking
    densur = 35.0 # track density for self-shielding tracking
elif ssh_trk == "SYBILT":
    iqua2 = 24
    nseg = 10

if case_name == "2x2_UOX":
    name_fig = "BWR_AT10_2x2_UOX.ps"
    name_fig_ss = "BWR_AT10_2x2_UOX_SS.ps"
elif case_name == "2x2_UOX_Gd":
    name_fig = "BWR_AT10_2x2_UOX_Gd.ps"
    name_fig_ss = "BWR_AT10_2x2_UOX_Gd_SS.ps"


    
name_BUvector = 'BUvector_'+case_name+'_'+burnup_points
name_compo = './_COMPO_AT10_'+case_name+'_'+burnup_points+"_"+ssh_method+"_"+flx_sol+"_"+CALC_opt
    
## Post treatment options
# Visualisation parameters stored in VISU_param
VISU_param=[visu_DRAGON,visu_SERPENT,visu_COMP,visu_DELTA]



# Recovering ListBU ListAUTOP ListCOMPO
[ListBU,ListSSH,ListCOMPO]=getLists(burnup_points)
#
# Create Steplist for BU - SELFSHIELDING - COMPO save 
StepList = lcm.new('LCM','burnup_steps')
StepList['ListBU']    = np.array(ListBU, dtype='f')
StepList['ListAutop'] = np.array(ListSSH, dtype='f')
StepList['ListCompo'] = np.array(ListCOMPO, dtype='f')
StepList.close() # close without erasing
#



# ---------- DRAGON5 PROCEDURES : BEGIN CALCULATION ROUTE
if case_name == "2x2_UOX":
    # ---------- GEOMETRY, 2x2 cluster with 4 UOX pins, 
    [pyGEOM,pyGEOM_SS] = GEOM_UOX("GEOM_2x2", "GEOM_2x2_SS")

    # ---------- MIXTURES
    # -- UOX with 3 different compositions (C1: 1, C2: 2 and C4: 1)
    pyMIX = MIX_UOX("MIX_2x2", Library, ssh_method)

elif case_name == "2x2_UOX_Gd":
    # ---------- GEOMETRY : 2x2 cluster with 3 UOX and 1 Gd pins
    [pyGEOM,pyGEOM_SS] = GEOM_GD("GEOM_2x2", "GEOM_2x2_SS")
    # ---------- MIXTURES
    # -- 3 UOX and 1 Gd with 3 different compositions

# ---------- SELF-SHIELDING TRACKING

if ssh_trk == "SALT":
    [pyGEOM_SS_G2S,pyFIG_SS,pyTRACK_SS,pyTF_EXC_SS] = TRK_SALT("GEOM_SS_G2S", name_fig_ss, "TRACK_SS", "TF_EXC_SS", pyGEOM_SS, ssh_sol, an2d_ssh, "SSH")

# ---------- FLUX CALCULATION TRACKING
if flx_trk == "SALT":
    [pyGEOM_G2S,pyFIG,pyTRACK,pyTF_EXC] = TRK_SALT("GEOM_G2S", name_fig, "TRACK", "TF_EXC",pyGEOM, flx_sol, an2d, "FLUX")
elif flx_trk == "NXT":
    [pyTRACK,pyTF_EXC,pyTRACK_SS] = TRK_NXT("TRACK", "TF_EXC", "TRACK_SS", pyGEOM, flx_sol, an2d, "FLUX")


if case_name == "2x2_UOX":
    # SSH + FLX + BU CALCULATION (+SSH at given burnup points)
    pyCOMPO = BU_UOX("COMPO", pyMIX, pyTRACK, pyTF_EXC, pyTRACK_SS, pyTF_EXC_SS, StepList, Multicompo, name_compo, flx_sol, CALC_opt)
elif case_name == "2x2_UOX_Gd":
    # SSH + FLX + BU CALCULATION (+SSH at given burnup points)
    pyCOMPO = BU_GD("COMPO", pyMIX, pyTRACK, pyTF_EXC, pyTRACK_SS, pyTF_EXC_SS, StepList, Multicompo, name_compo, flx_sol, CALC_opt)


# ---------- POST-PROCESSING
POSTPROC_2x2(pyCOMPO, ListCOMPO, ListSSH, burnup_points, case_name, ssh_method, ssh_sol, flx_sol, CALC_opt)

