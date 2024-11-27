# PyGan procedure to import compos from the BUScheme_Gd157 case and post-process / compare them to Serpent2 BU evolution files

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import lifo
import lcm
import cle2000

from Serpent2_edep_pcc import Serpent2_case 
from Serpent2_edep_pcc import SERPENT2_comparisons
from POSTPROC_BUScheme import DRAGON_case
from POSTPROC_BUScheme import D5S2_comparisons
from getLists import *

post_treat_S2 = True
post_treat_D5 = True
post_treat_D5S2_comparison = True

D5_savedir = "PyGan_results_BUScheme/D5_HOM_UOX_Gd157_edep_pcc"
S2_savedir = "PyGan_results_BUScheme/S2_HOM_UOX_Gd157_edep_pcc"
COMP_savedir = "PyGan_results_BUScheme/COMP_HOM_UOX_Gd157_edep_pcc"
# Check if the directories exist and create them if necesary
if not os.path.exists(D5_savedir):
	os.makedirs(D5_savedir)
if not os.path.exists(S2_savedir):
	os.makedirs(S2_savedir)
if not os.path.exists(COMP_savedir):
	os.makedirs(COMP_savedir)

original_specific_power = 38.6 # W/gU

if post_treat_D5:
	BUlists = {}
	draglib_names = ["J311_172", "J311_281", "J311_295"] #, "J311_315"]
	D5_PCC_ids = ["0","1","2","3","3b"]
	bu_points = "BOC_fine_autop5"
	lists = getLists(bu_points)
	BUlists["BU"] = lists[0]
	BUlists["AUTOP"] = lists[1]
	BUlists["COMPO"] = lists[2]
	ssh_module = "USS"
	ssh_method = "RSE"
	correlation = "noCORR"
	D5_CASES = []
	for draglib_name in draglib_names:
		for pcc_id in D5_PCC_ids:
			D5_CASES.append(DRAGON_case(draglib_name, pcc_id, bu_points, ssh_module, ssh_method, correlation, ["U235", "U238", "Gd157", "Gd158"], BUlists, D5_savedir))


if post_treat_S2:
	S2_CASES = []
	PyNjoy_cases = []
	oldlib_cases = []
	S2_libs = ["PyNjoy2016","oldlib"]
	for lib_name in S2_libs:
		if lib_name == "PyNjoy2016":
			pcc_ids = [1,2]
			edep_ids = [0,1,2]
		elif lib_name == "oldlib":
			pcc_ids = [1]
			edep_ids = [0,1]
		for pcc in pcc_ids:
			for edep in edep_ids:
				S2_CASES.append(Serpent2_case("HOM_UOX_Gd157", lib_name, edep, pcc, original_specific_power, ["U235", "U238", "Gd157", "Gd158"], S2_savedir))
				if lib_name == "PyNjoy2016":
					PyNjoy_cases.append(S2_CASES[-1])
				elif lib_name == "oldlib":
					oldlib_cases.append(S2_CASES[-1])


	for S2_case in S2_CASES:
		S2_case.plot_keff()
		S2_case.plot_concentrations(["Gd157", "Gd158"])
		#S2_case.rescale_BU("MWd/tU")	

	S2_comp_PyNjoy = SERPENT2_comparisons("Gd157 evolution schemes Pynjoy", PyNjoy_cases, COMP_savedir)
	S2_comp_PyNjoy.set_accepted_params_for_lib("PyNjoy2016", [0,1,2], [1,2,3,4,6])
	S2_comp_PyNjoy.set_ref_case(ref_lib="PyNjoy2016", ref_edep=0, ref_pcc=1)
	S2_comp_PyNjoy.plot_keffs()
	S2_comp_PyNjoy.compute_delta_keffs()
	S2_comp_PyNjoy.compute_delta_Niso(["U235", "U238", "Gd157", "Gd158"])
	S2_comp_PyNjoy.plot_delta_keffs()
	S2_comp_PyNjoy.plot_delta_Niso(["U235", "U238", "Gd157", "Gd158"])

	S2_comp_oldlib = SERPENT2_comparisons("Gd157 evolution schemes", oldlib_cases, COMP_savedir)
	S2_comp_oldlib.set_ref_case(ref_lib="oldlib", ref_edep=0, ref_pcc=1)
	S2_comp_oldlib.plot_keffs()
	S2_comp_oldlib.compute_delta_keffs()
	S2_comp_oldlib.compute_delta_Niso(["U235", "U238", "Gd157", "Gd158"])
	S2_comp_oldlib.plot_delta_keffs()
	S2_comp_oldlib.plot_delta_Niso(["U235", "U238", "Gd157", "Gd158"])


	print("Serpent2 comparisons done")

if post_treat_D5S2_comparison:
	D5_S2_comparisons = []
	for D5_case in D5_CASES:
		D5_S2_comparisons.append(D5S2_comparisons(f"{D5_case.draglib_name} {D5_case.ssh_module} {D5_case.ssh_method} PCC{D5_case.pcc_id}", D5_case, S2_CASES, ["U235", "U238", "Gd157", "Gd158"], COMP_savedir))
		D5_S2_comparisons[-1].compare_keffs() 
		D5_S2_comparisons[-1].compare_Ni()
		D5_S2_comparisons[-1].plot_delta_Keff()
		D5_S2_comparisons[-1].plot_delta_Ni()
	print("D5-S2 comparisons done")
print("All post-treatments done")
