# PyGan script for gduo2_295_kec1 evolution study
# Authors : R.Guasch, A. Hébert
# Date : 2025-01-09
# Hypothesis : RSE method with CORR should elimiate mutual self shielding effect between U8 and Gd157 in group 100.

# Importing libraries
import sys
import numpy as np
import os, shutil
import lifo
import lcm
import cle2000
import matplotlib

# Importing Procedures
# --- Geometry
from geom_gduo2 import *
# --- Materials
from lib_gduo2 import *
# --- Traking
from trk_gduo2 import *
# --- Burnup
from bu_gduo2 import *

# --- OTHERS
from getLists import *

# --- Post-treatment

from Serpent2_edep_pcc import Serpent2_case 
from Serpent2_edep_pcc import SERPENT2_comparisons
from POSTPROC_gduo2_evol import DRAGON_case
from POSTPROC_gduo2_evol import D5S2_comparisons

# Input data
burnup_points = 'Gd_autop3' # 'Gd_autop4' or 'Gd_autop3', refine if necessary
# Recovering BUlist SSHlist COMPOlist
[BUlist, SSHlist, COMPOlist]=getLists(burnup_points)

post_process_D5S2 = True # True or False

# Evaluation name
evaluation_name = 'J311_295' # 'J311_295' or 'ENDFb8r1_295'
correlation = 'CORR' # 'CORR' or 'NOCORR'

# StepList
StepList = lcm.new('LCM','burnup_steps')
StepList['BUlist']    = np.array(BUlist, dtype='f')
StepList['SSHlist'] = np.array(SSHlist, dtype='f')
StepList['COMPOlist'] = np.array(COMPOlist, dtype='f')
StepList.close() # close without erasing

# Creation of results directory
path=os.getcwd()

a=os.path.exists(f'{path}/gduo2_295_kec1_evol_results_PyGan')
if a==False:
	os.mkdir(f'{path}/gduo2_295_kec1_evol_results_PyGan')

# Begin calculation procedure for gduo2_kec1 evolution study

# Geometry
pyGEOM = geom_gduo2()

# Materials library
pyLIB = lib_gduo2(evaluation_name, correlation)

# Tracking
pyUOX_TBH, pyTRACK, pyTF_EXC, pyTRACK_SS = trk_gduo2(pyGEOM)

# Burnup evolution
pyCOMPO = bu_gduo2("COMPO", pyLIB, pyTRACK, pyTF_EXC, pyTRACK_SS, StepList, name_compo=f'_COMPO_gduo2_295_kec1_{evaluation_name}_{correlation}')


# End of calculation procedure for gduo2_kec1 evolution study
# Post treat results once Serpent2 / OpenMC has run.

if post_process_D5S2:
	D5_savedir = "gduo2_295_kec1_evol_results_PyGan/D5_results"
	S2_savedir = "gduo2_295_kec1_evol_results_PyGan/S2_results"
	COMP_savedir = "gduo2_295_kec1_evol_results_PyGan/D5_S2_comparison"
	# Check if the directories exist and create them if necesary
	if not os.path.exists(D5_savedir):
		os.makedirs(D5_savedir)
	if not os.path.exists(S2_savedir):
		os.makedirs(S2_savedir)
	if not os.path.exists(COMP_savedir):
		os.makedirs(COMP_savedir)
	edep_id = 2
	pcc_id = 1
	BUlists = {}
	draglib_name = "J311_295"
	D5_PCC_ids = ["1"]
	bu_points = burnup_points
	lists = getLists(bu_points)
	BUlists["BU"] = lists[0]
	BUlists["AUTOP"] = lists[1]
	BUlists["COMPO"] = lists[2]
	D5_case = DRAGON_case(pyCOMPO, draglib_name, pcc_id, bu_points, "USS", "RSE", correlation, ["U235", "U238", "Gd155", "Gd157"], BUlists, save_dir=D5_savedir)
	S2_case = Serpent2_case("gduo2_295_kec1", "PyNjoy2016", edep_id = 2, pcc_id = 1, specific_power = 38.6, tracked_nuclides=["U235","U238", "Gd155", "Gd157"], save_dir=S2_savedir)

	D5S2_case = D5S2_comparisons("gduo2_295_kec1_evol", D5_case, [S2_case], ["U235", "U238", "Gd155", "Gd157"], save_dir=COMP_savedir)
	D5S2_case.compare_keffs()
	D5S2_case.compare_Ni()
	D5S2_case.plot_delta_Keff()
	D5S2_case.plot_delta_Ni()