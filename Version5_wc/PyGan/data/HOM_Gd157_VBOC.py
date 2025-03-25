### post treat HOM Gd157 Very Beginning of Cycle results
# open compos from Dragon/Linux_*/
# open Serpent2 results

import serpentTools as st
import lcm
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

from postproc_cst_pow_evol import Serpent2_case as S2_case
from postproc_cst_pow_evol import DRAGON_case as D5_case


#### Running the post-treatment : options

eval_name = "J311" # or "J311" "endfb8r1"
origin_of_S2_data = "pynjoy2016" #"pynjoy2012_kerma" # or "pynjoy2016" or "NJOY2016"

# Serpent2 results
post_process_S2 = False # True or False
post_process_D5 = True # True or False

# Creation of results directory
path=os.getcwd()

a=os.path.exists(f'{path}/HOM_Gd157_VBOC_results_{eval_name}')
if a==False:
	os.mkdir(f'{path}/HOM_Gd157_VBOC_results_{eval_name}')

save_dir = f'{path}/HOM_Gd157_VBOC_results_{eval_name}'

save_dir_Serpent2 = f'{save_dir}/Serpent2'
a=os.path.exists(save_dir_Serpent2)
if a==False:
    os.mkdir(save_dir_Serpent2)

save_dir_DRAGON = f'{save_dir}/DRAGON'
a=os.path.exists(save_dir_DRAGON)
if a==False:
    os.mkdir(save_dir_DRAGON)



tracked_nuclides = ["U235", "U234", "U238", "Pu239", "Pu240", "Pu241", "Pu242", "Gd157", "Gd158", "Xe135", "Sm149"]

if post_process_S2:
    # Serpent2 results
    # ----------------
    Serpent2_case = S2_case(case_name = "HOM_Gd157_VBOC", lib_name=f"{eval_name}_{origin_of_S2_data}", 
                            edep_id=2, isEcaptSet=False, pcc_id=2, specific_power=38.6, 
                            tracked_nuclides=tracked_nuclides, save_dir=save_dir_Serpent2)	


    Serpent2_case.plot_keff()
    for iso in tracked_nuclides:
        Serpent2_case.plot_concentrations([iso])

if post_process_D5:
    # DRAGON results
    # --------------
        
    ## Compo without KERMAs
    name_compo = f"COMPO_HOM_Gd157_VBOC_{eval_name}_295"
    ## Compo with KERMAs
    name_compo_Kermas = f"COMPO_HOM_Gd157_VBOC_{eval_name}_295K"
    print(name_compo)
    print(name_compo_Kermas)
    os.chdir("Linux_aarch64/")
    print(os.listdir())
    #print(name_compo)
    pyCOMPO = lcm.new('LCM_INP', name_compo, impx=0)
    pyCOMPO_KERMA = lcm.new('LCM_INP', name_compo_Kermas, impx=0)
    print(f"pyCOMPO['EDIBU']['GLOBAL']['pval000000001'] are = {pyCOMPO['EDIBU']['GLOBAL']['pval00000001']}")
    os.chdir(path)
    BU_lists = {"BU":[0.0, 25.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0,
                    750.0, 800.0, 850.0, 900.0, 950.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0,    
                    5500.0, 6000.0, 6500.0, 7000.0, 7500.0, 8000.0, 8500.0, 9000.0, 9500.0, 10000.0 ],
                "AUTOP":[25.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0,
                    750.0, 800.0, 850.0, 900.0, 950.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0,    
                    5500.0, 6000.0, 6500.0, 7000.0, 7500.0, 8000.0, 8500.0, 9000.0, 9500.0, 10000.0 ],
                "COMPO":[0.0, 25.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0,
                    750.0, 800.0, 850.0, 900.0, 950.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0,    
                    5500.0, 6000.0, 6500.0, 7000.0, 7500.0, 8000.0, 8500.0, 9000.0, 9500.0, 10000.0 ]}
    # Load the composition
    DRAGON_case = D5_case(pyCOMPO=pyCOMPO, dlib_name=f"{eval_name}_295", bu_points="VBOC", 
                        ssh_opt="PT", correlation="CORR", sat="def", depl_sol="def", 
                        tracked_nuclides=tracked_nuclides, BU_lists=BU_lists, save_dir=save_dir_DRAGON)

    DRAGON_caseK = D5_case(pyCOMPO=pyCOMPO_KERMA, dlib_name=f"{eval_name}_295K", bu_points="VBOC",
                            ssh_opt="PT", correlation="CORR", sat="def", depl_sol="def",
                            tracked_nuclides=tracked_nuclides, BU_lists=BU_lists, save_dir=save_dir_DRAGON)


    ## Comparison between D5 results :

    delta_keff_D5 = (DRAGON_caseK.DRAGON_Keff - DRAGON_case.DRAGON_Keff)*1e5 # pcm

    plt.figure()
    plt.plot(DRAGON_case.DRAGON_BU, delta_keff_D5, label = f"KERMA - non-KERMA", marker = "x", linestyle = "--")
    plt.xlabel("Burnup [MWd/kgU]")
    plt.ylabel("Keff difference [pcm]")
    plt.title(f"$\\Delta$ Keff for KERMA and non-KERMA HOM_Gd157_VBOC")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir_DRAGON}/DRAGON_Keff_diff_{DRAGON_case.draglib_name}_{DRAGON_caseK.draglib_name}.png")

    for iso in tracked_nuclides:
        delta_Niso = [(DRAGON_caseK.DRAGON_ISOTOPESDENS[iso][idx] - DRAGON_case.DRAGON_ISOTOPESDENS[iso][idx]) * 100 / DRAGON_case.DRAGON_ISOTOPESDENS[iso][idx]
                        if DRAGON_case.DRAGON_ISOTOPESDENS[iso][idx] != 0 else 0
                        for idx in range(len(DRAGON_case.DRAGON_ISOTOPESDENS[iso]))]
        print(f"For iso = {iso}, deltaNiso = {delta_Niso} %")
        fig, ax1 = plt.subplots()
        ax1.plot(DRAGON_case.DRAGON_BU, DRAGON_case.DRAGON_ISOTOPESDENS[iso], label = f"non-KERMA", marker = "x", linestyle = "--")
        ax1.plot(DRAGON_caseK.DRAGON_BU, DRAGON_caseK.DRAGON_ISOTOPESDENS[iso], label = f"KERMA", marker = "x", linestyle = "--")
        ax1.set_xlabel("Burnup [MWd/kgU]")
        ax1.set_ylabel(f"{iso} density [iso/cm*b]", color='tab:blue')

        ax2 = ax1.twinx()
        ax2.plot(DRAGON_case.DRAGON_BU, delta_Niso, label = f"KERMA - non-KERMA", marker = "x", linestyle = "--", color='tab:red')
        ax2.set_ylabel(f"{iso} density difference [%]", color='tab:red')
        plt.title(f"$\\Delta$ N{iso} beteween KERMA and non-KERMA HOM_Gd157_VBOC, {eval_name}")
        plt.legend()
        plt.grid()
        plt.savefig(f"{save_dir_DRAGON}/DRAGON_{iso}_diff_{DRAGON_case.draglib_name}_{DRAGON_caseK.draglib_name}.png")
        plt.close()

        delta_Niso_abs = [(DRAGON_caseK.DRAGON_ISOTOPESDENS[iso][idx] - DRAGON_case.DRAGON_ISOTOPESDENS[iso][idx])
                        for idx in range(len(DRAGON_case.DRAGON_ISOTOPESDENS[iso]))]
        print(f"For iso = {iso}, absolute deltaNiso = {delta_Niso_abs}")
        fig, ax1 = plt.subplots()
        ax1.plot(DRAGON_case.DRAGON_BU, DRAGON_case.DRAGON_ISOTOPESDENS[iso], label = f"non-KERMA", marker = "x", linestyle = "--")
        ax1.plot(DRAGON_caseK.DRAGON_BU, DRAGON_caseK.DRAGON_ISOTOPESDENS[iso], label = f"KERMA", marker = "x", linestyle = "--")
        ax1.set_xlabel("Burnup [MWd/kgU]")
        ax1.set_ylabel(f"{iso} density [iso/cm*b]", color='tab:blue')

        ax2 = ax1.twinx()
        ax2.plot(DRAGON_case.DRAGON_BU, delta_Niso_abs, label = f"KERMA - non-KERMA", marker = "x", linestyle = "--", color='tab:red')
        ax2.set_ylabel(f"{iso} density difference [iso/cm*b]", color='tab:red')
        plt.title(f"$\\Delta$ N{iso} beteween KERMA and non-KERMA HOM_Gd157_VBOC, {eval_name}")
        plt.legend()
        plt.grid()
        plt.savefig(f"{save_dir_DRAGON}/DRAGON_{iso}_diff_abs_{DRAGON_case.draglib_name}_{DRAGON_caseK.draglib_name}.png")
        plt.close()



if post_process_S2 and post_process_D5:
    ## Plot Serpent2 and DRAGON results
    for iso in tracked_nuclides:
        fig, ax1 = plt.subplots()
        ax1.plot(Serpent2_case.BU, Serpent2_case.Ni[iso], label = f"Serpent2", marker = "x", linestyle = "--")
        ax1.plot(DRAGON_case.DRAGON_BU, DRAGON_case.DRAGON_ISOTOPESDENS[iso], label = f"DRAGON", marker = "x", linestyle = "--")
        ax1.plot(DRAGON_caseK.DRAGON_BU, DRAGON_caseK.DRAGON_ISOTOPESDENS[iso], label = f"DRAGON KERMA", marker = "x", linestyle = "--")
        ax1.set_xlabel("Burnup [MWd/kgU]")
        ax1.set_ylabel(f"{iso} density [iso/cm*b]", color='tab:blue')
        plt.title(f"N{iso} for Serpent2 and DRAGON cases for HOM_Gd157_VBOC, {eval_name}")
        plt.legend()
        plt.grid()
        plt.savefig(f"{save_dir}/Serpent2_DRAGON_{iso}_{Serpent2_case.case_name}_{DRAGON_case.draglib_name}_{DRAGON_caseK.draglib_name}.png")
        plt.close()

    # Comparison between Serpent2 and DRAGON results
    # ----------------------------------------------
    delta_keff = (DRAGON_case.DRAGON_Keff - Serpent2_case.keffs)*1e5 # pcm
    delta_keffK = (DRAGON_caseK.DRAGON_Keff - Serpent2_case.keffs)*1e5 # pcm

    plt.figure()
    plt.plot(Serpent2_case.BU, delta_keff, label = f"DRAGON - Serpent2", marker = "x", linestyle = "--")
    plt.plot(Serpent2_case.BU, delta_keffK, label = f"DRAGON KERMA - Serpent2", marker = "x", linestyle = "--")
    plt.xlabel("Burnup [MWd/kgU]")
    plt.ylabel("Keff difference [pcm]")
    plt.title(f"$\\Delta$ Keff between Serpent2 and DRAGON5 for HOM_Gd157_VBOC, {eval_name}")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir}/Serpent2_DRAGON_Keff_diff_{Serpent2_case.case_name}_{DRAGON_case.draglib_name}_{DRAGON_caseK.draglib_name}.png")
    plt.close()

    for iso in tracked_nuclides:
        
        delta_Niso = [(DRAGON_case.DRAGON_ISOTOPESDENS[iso][idx] - Serpent2_case.Ni[iso][idx]) for idx in range(len(Serpent2_case.Ni[iso]))]
        delta_NisoK = [(DRAGON_caseK.DRAGON_ISOTOPESDENS[iso][idx] - Serpent2_case.Ni[iso][idx]) for idx in range(len(Serpent2_case.Ni[iso]))]
        fig, ax1 = plt.subplots()

        ax1.plot(Serpent2_case.BU, delta_Niso, label = f"DRAGON - Serpent2", marker = "D", linestyle = "--", color='tab:red')
        ax1.plot(Serpent2_case.BU, delta_NisoK, label = f"DRAGON KERMA - Serpent2", marker = "D", linestyle = "--", color='tab:green')
        ax1.set_ylabel(f"{iso} density difference [#/b*cm]", color='tab:red')
        plt.title(f"$\\Delta$ N{iso} between Serpent2 and DRAGON5 for HOM_Gd157_VBOC, {eval_name}")
        plt.legend()
        plt.grid()
        plt.savefig(f"{save_dir}/Serpent2_DRAGON_{iso}_diff_{Serpent2_case.case_name}_{DRAGON_case.draglib_name}_{DRAGON_caseK.draglib_name}.png")
        plt.close()



print(Serpent2_case.BU)
print(DRAGON_case.DRAGON_BU)
print(DRAGON_caseK.DRAGON_BU)

# Relative diferences between Serpent2 and DRAGON5 BU points

delta_BUs = [(DRAGON_case.DRAGON_BU[idx] - Serpent2_case.BU[idx]) * 100 / Serpent2_case.BU[idx]  if Serpent2_case.BU[idx] != 0 else 0 for idx in range(len(Serpent2_case.BU))]
    
plt.figure()
plt.plot(Serpent2_case.BU, delta_BUs, label = f"DRAGON - Serpent2", marker = "x", linestyle = "--")
plt.xlabel("Burnup [MWd/kgU]")
plt.ylabel("Burnup difference [%]")
plt.title(f"$\\Delta$ BU between Serpent2 and DRAGON5 for HOM_Gd157_VBOC, {eval_name}")
plt.legend()
plt.grid()
plt.savefig(f"{save_dir}/Serpent2_DRAGON_BU_diff_{Serpent2_case.case_name}_{DRAGON_case.draglib_name}_{DRAGON_caseK.draglib_name}.png")
plt.close()

# Consider interpolating between Serpent2 and DRAGON5 BU points : calculated BU in S2 edepmode 2 is not the same as in D5 --> to be checked
# Interpolating between Serpent2 and DRAGON5 BU points

# keff

val_keff_interp = np.interp(Serpent2_case.BU, DRAGON_case.DRAGON_BU, DRAGON_case.DRAGON_Keff)
val_keff_interp_Kerma = np.interp(Serpent2_case.BU, DRAGON_caseK.DRAGON_BU, DRAGON_caseK.DRAGON_Keff)

# NGd157
val_NGd157_interp = np.interp(Serpent2_case.BU, DRAGON_case.DRAGON_BU, DRAGON_case.DRAGON_ISOTOPESDENS["Gd157"])
val_NGd157_interp_Kerma = np.interp(Serpent2_case.BU, DRAGON_caseK.DRAGON_BU, DRAGON_caseK.DRAGON_ISOTOPESDENS["Gd157"])

# compute new differences

delta_keff_interp = (val_keff_interp - Serpent2_case.keffs)*1e5 # pcm
delta_keff_interp_Kerma = (val_keff_interp_Kerma - Serpent2_case.keffs)*1e5 # pcm

delta_NGd157_interp = (val_NGd157_interp - Serpent2_case.Ni["Gd157"])*100 / Serpent2_case.Ni["Gd157"]
delta_NGd157_interp_Kerma = (val_NGd157_interp_Kerma - Serpent2_case.Ni["Gd157"])*100 / Serpent2_case.Ni["Gd157"]

plt.figure()
plt.plot(Serpent2_case.BU, delta_keff_interp, label = f"DRAGON interp - Serpent2", marker = "x", linestyle = "--")
plt.plot(Serpent2_case.BU, delta_keff_interp_Kerma, label = f"DRAGON KERMA interp - Serpent2", marker = "x", linestyle = "--")
plt.xlabel("Burnup [MWd/kgU]")
plt.ylabel("Keff difference [pcm]")
plt.title(f"$\\Delta$ Keff between Serpent2 and interpolated DRAGON5 for HOM_Gd157_VBOC, {eval_name}")
plt.legend()
plt.grid()
plt.savefig(f"{save_dir}/Serpent2_DRAGON_interp_Keff_diff_{Serpent2_case.case_name}_{DRAGON_case.draglib_name}_{DRAGON_caseK.draglib_name}.png")
plt.close()

plt.figure()
plt.plot(Serpent2_case.BU, delta_NGd157_interp, label = f"DRAGON interp - Serpent2", marker = "x", linestyle = "--")
plt.plot(Serpent2_case.BU, delta_NGd157_interp_Kerma, label = f"DRAGON KERMA interp - Serpent2", marker = "x", linestyle = "--")
plt.xlabel("Burnup [MWd/kgU]")
plt.ylabel("Gd157 density difference [%]")
plt.title(f"$\\Delta$ NGd157 between Serpent2 and interpolated DRAGON5 for HOM_Gd157_VBOC, {eval_name}")
plt.legend()
plt.grid()
plt.savefig(f"{save_dir}/Serpent2_DRAGON_interp_Gd157_diff_{Serpent2_case.case_name}_{DRAGON_case.draglib_name}_{DRAGON_caseK.draglib_name}.png")
plt.close()
