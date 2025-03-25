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
import matplotlib.pyplot as plt

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
evaluation_name = 'J311_295' # 'J311_295' or 'J311295K' or 'ENDFb8r1_295' or 'E8R1295K'
evaluation_name_kerma = 'J311295K'
correlation = 'NOCORR' # 'CORR' or 'NOCORR'

tracked_nuclides = ["U235", "U234", "U238", "Pu239", "Pu240", "Pu241", "Pu242", "Gd155", "Gd157", "Xe135", "Sm149"]

# StepList
StepList = lcm.new('LCM','burnup_steps')
StepList['BUlist']    = np.array(BUlist, dtype='f')
StepList['SSHlist'] = np.array(SSHlist, dtype='f')
StepList['COMPOlist'] = np.array(COMPOlist, dtype='f')
StepList.close() # close without erasing

# StepList
StepList2 = lcm.new('LCM','burnup_steps')
StepList2['BUlist']    = np.array(BUlist, dtype='f')
StepList2['SSHlist'] = np.array(SSHlist, dtype='f')
StepList2['COMPOlist'] = np.array(COMPOlist, dtype='f')
StepList2.close() # close without erasing

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
puLIB_KERMA = lib_gduo2(evaluation_name_kerma, correlation)

# Tracking
pyUOX_TBH, pyTRACK, pyTF_EXC, pyTRACK_SS = trk_gduo2(pyGEOM)

# Burnup evolution
pyCOMPO = bu_gduo2("COMPO", pyLIB, pyTRACK, pyTF_EXC, pyTRACK_SS, StepList, name_compo=f'_COMPO_gduo2_295_kec1_{evaluation_name}_{correlation}')
pyCOMPO_KERMA = bu_gduo2("COMPO", puLIB_KERMA, pyTRACK, pyTF_EXC, pyTRACK_SS, StepList2, name_compo=f'_COMPO_gduo2_295_kec1_{evaluation_name_kerma}_{correlation}')


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
	draglib_name = evaluation_name
	draglib_name_kerma = evaluation_name_kerma
	bu_points = burnup_points
	lists = getLists(bu_points)
	BUlists["BU"] = lists[0]
	BUlists["AUTOP"] = lists[1]
	BUlists["COMPO"] = lists[2]
	
	D5_case = DRAGON_case(pyCOMPO, draglib_name, pcc_id, bu_points, "USS", "RSE", correlation, tracked_nuclides, BUlists, save_dir=D5_savedir)
	D5_case.plot_keffs()
	for iso in tracked_nuclides:
		D5_case.plot_Ni(iso)
	
	D5_case_kerma = DRAGON_case(pyCOMPO_KERMA, draglib_name_kerma, pcc_id, bu_points, "USS", "RSE", correlation, tracked_nuclides, BUlists, save_dir=D5_savedir)
	D5_case_kerma.plot_keffs()
	for iso in tracked_nuclides:
		D5_case_kerma.plot_Ni(iso)
	

	### Compare KERMA vs non-KERMA
	## Keffs :
	diff_keffs = (D5_case_kerma.DRAGON_Keff - D5_case.DRAGON_Keff)*1e5
	plt.figure()
	plt.plot(D5_case.DRAGON_BU, diff_keffs, label = f"KERMA - non-KERMA", marker = "x", linestyle = "--")
	plt.xlabel("Burnup [MWd/kgU]")
	plt.ylabel("Keff difference [pcm]")
	plt.title(f"Keff difference between KERMA and non-KERMA cases for {draglib_name} and {draglib_name_kerma} cases")
	plt.legend()
	plt.grid()
	plt.savefig(f"{D5_savedir}/DRAGON_Keff_diff_{draglib_name}_{draglib_name_kerma}_{correlation}_PCC{pcc_id}.png")
	plt.close()

	## Ni :
	delta_Niso = {}
	for iso in tracked_nuclides:
		delta_Niso = [(D5_case_kerma.DRAGON_ISOTOPESDENS[iso][idx] - D5_case.DRAGON_ISOTOPESDENS[iso][idx]) * 100 / D5_case.DRAGON_ISOTOPESDENS[iso][idx]
			if D5_case.DRAGON_ISOTOPESDENS[iso][idx] != 0 else 0
			for idx in range(len(D5_case.DRAGON_ISOTOPESDENS[iso]))]
		plt.figure()
		plt.plot(D5_case.DRAGON_BU, delta_Niso, label = f"KERMA - non-KERMA", marker = "x", linestyle = "--")
		plt.xlabel("Burnup [MWd/kgU]")
		plt.ylabel(f"{iso} density difference [%]")
		plt.title(f"{iso} density difference between KERMA and non-KERMA cases for {draglib_name} and {draglib_name_kerma} cases")
		plt.legend()
		plt.grid()
		plt.savefig(f"{D5_savedir}/DRAGON_{iso}_diff_{draglib_name}_{draglib_name_kerma}_{correlation}_PCC{pcc_id}.png")
		plt.close()

	S2_case = Serpent2_case("gduo2_295_kec1", "J311_pynjoy2016", edep_id = 2, pcc_id = 1, specific_power = 38.6, tracked_nuclides=tracked_nuclides, save_dir=S2_savedir)

	# Compare D5 and S2
	
	D5S2_case = D5S2_comparisons("gduo2_295_kec1_evol", D5_case, [S2_case], tracked_nuclides, save_dir=COMP_savedir)
	D5S2_case.compare_keffs()
	D5S2_case.compare_Ni()
	D5S2_case.plot_delta_Keff()
	D5S2_case.plot_delta_Ni()
	

	D5S2_case_kerma = D5S2_comparisons("gduo2_295_kec1_evol_KERMA", D5_case_kerma, [S2_case], tracked_nuclides, save_dir=COMP_savedir)
	D5S2_case_kerma.compare_keffs()
	D5S2_case_kerma.compare_Ni()
	D5S2_case_kerma.plot_delta_Keff()
	D5S2_case_kerma.plot_delta_Ni()