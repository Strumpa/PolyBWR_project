### Main comparison script used to run S2 vsOpenMC comaprisons


import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import serpentTools as st
from post_process import Serpent2_case as S2_case
from post_process import OpenMC_case as OMC_case


def plot_delta_keffs(S2_edep0, S2_edep0_fissh, S2_edep2, OMC_fissq, OMC_fissq_set, OMC_energy_deposition):
    integrator = OMC_fissq.integrator
    # delta keff (pcm) comparison
    ## Might need to interpolate to compare with edepmode 2 results ?
    # Check delta in days and delta in reconstructed BU [GWj/tU] for S2 and OpenMC
    delta_Days_S2_edep0_OpenMC_fissq = S2_edep0.BUdays - OMC_fissq.BUdays
    delta_Days_S2_edep0_fissh_OpenMC_fissq_set = S2_edep0_fissh.BUdays - OMC_fissq_set.BUdays
    delta_Days_S2_edep2_OpenMC_energy_deposition = S2_edep2.BUdays - OMC_energy_deposition.BUdays

    print(f"BEGIN BU POINTS COMPARISONS using time integrator {integrator}")
    print(f"BU in S2 edep0 : {S2_edep0.BU} GWd/tU")
    print(f"BU in OpenMC fissq: {OMC_fissq.BUdays*specific_power} GWd/tU")
    print(f"delta BU S2 edep0 - OpenMC fissq: {S2_edep0.BU - OMC_fissq.BUdays*specific_power} GWd/tU")
    print(f"delta days S2 edep0 - OpenMC fissq: {delta_Days_S2_edep0_OpenMC_fissq} days")

    print(f"BU in S2 edep0 set fissh: {S2_edep0_fissh.BU} GWd/tU")
    print(f"BU in OpenMC set fissq: {OMC_fissq_set.BUdays*specific_power} GWd/tU")
    print(f"delta BU S2 edep0 set fissh - OpenMC set fissq: {S2_edep0_fissh.BU - OMC_fissq_set.BUdays*specific_power} GWd/tU")
    print(f"delta days S2 edep0 - OpenMC fissq set: {delta_Days_S2_edep0_fissh_OpenMC_fissq_set} days")

    print(f"BU in S2 edep2: {S2_edep2.BU} GWd/tU")
    print(f"BU in OpenMC energy_deposition: {OMC_energy_deposition.BUdays*specific_power} GWd/tU")
    print(f"delta BU S2 edep2 - OpenMC fissq: {S2_edep2.BU - OMC_energy_deposition.BUdays*specific_power} GWd/tU")
    print(f"delta days S2 edep2 - OpenMC fissq: {delta_Days_S2_edep2_OpenMC_energy_deposition} days")
    print("END BU POINTS COMPARISONS")


    delta_keff_edep0 = (S2_edep0.keffs - OMC_fissq.keffs)*1e5
    delta_keff_edep0_set = (S2_edep0_fissh.keffs - OMC_fissq_set.keffs)*1e5
    delta_keff_edep2 = (S2_edep2.keffs - OMC_energy_deposition.keffs)*1e5
    time_days = S2_edep0.BUdays
    plt.figure()
    plt.plot(time_days, delta_keff_edep0, label="Serpent2 - OpenMC, edep0", marker="x", linestyle="--")
    plt.plot(time_days, delta_keff_edep0_set, label="Serpent2 - OpenMC, edep0 \n set fission q-values", marker="x", linestyle="--")
    plt.plot(time_days, delta_keff_edep2, label="Serpent2 - OpenMC, edep2", marker="x", linestyle="--")
    plt.plot(time_days, 200*np.ones_like(time_days), label="200 pcm", linestyle="--", color="red")
    plt.plot(time_days, -200*np.ones_like(time_days), label="-200 pcm", linestyle="--", color="red")
    plt.xlabel("Time [days]")
    plt.ylabel("Delta Keff [pcm]")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir}/delta_keff_energy_deposition_comparison_S2_OMC_{integrator}.png")

    print(f"PRINTING RESULTS FOR {integrator}")
    print(f"Printing max, RMS and avg delta keff values for S2 edep0 vs OpenMC")
    print(f"max delta keff S2 edep0 vs OpenMC = {np.max(np.abs(delta_keff_edep0))} pcm")
    print(f"RMS delta keff S2 edep0 vs OpenMC = {np.sqrt(np.mean(delta_keff_edep0**2))} pcm")
    print(f"avg delta keff S2 edep0 vs OpenMC = {np.mean(delta_keff_edep0)} pcm")
    print(f"delta keff S2 edep0 vs OpenMC = {delta_keff_edep0} pcm")

    print("Printing max, RMS and avg delta keff values for S2 edep0 set fissh0 vs OpenMC with set q-fiss values")
    print(f"max delta keff S2 edep0 set vs OpenMC set qfiss = {np.max(np.abs(delta_keff_edep0_set))} pcm")
    print(f"RMS delta keff S2 edep0 set vs OpenMC set qfiss = {np.sqrt(np.mean(delta_keff_edep0_set**2))} pcm")
    print(f"avg delta keff S2 edep0 set vs OpenMC set qfiss = {np.mean(delta_keff_edep0_set)} pcm")
    print(f"delta keff S2 edep0 set vs OpenMC set qfiss = {delta_keff_edep0_set} pcm")

    print("Prining max, RMS and avg delta keff values for S2 edep2 vs OpenMC")
    print(f"max delta keff S2 edep2 vs OpenMC = {np.max(np.abs(delta_keff_edep2))} pcm")
    print(f"RMS delta keff S2 edep2 vs OpenMC = {np.sqrt(np.mean(delta_keff_edep2**2))} pcm")
    print(f"avg delta keff S2 edep2 vs OpenMC = {np.mean(delta_keff_edep2)} pcm")
    print(f"delta keff S2 edep2 vs OpenMC = {delta_keff_edep2} pcm")


def plot_keffs(S2_edep0, S2_edep0_fissh, S2_edep2, OMC_fissq, OMC_fissq_set, OMC_energy_deposition):
    integrator = OMC_fissq.integrator
    ## Prepare comparison
    # compare keffs for S2 edepmode 0 and OpenMC fissq
    plt.figure()
    plt.errorbar(OMC_fissq.BUdays, # time axis in days
                OMC_fissq.keffs, # keff values
                yerr=OMC_fissq.sigmas_keff, # statistical uncertainties on keff
                label=f"OpenMC, fissq {OMC_fissq.integrator}", # labels and plot customization options
                marker=".", linestyle="--")
    plt.errorbar(S2_edep0.BUdays, 
                S2_edep0.keffs, 
                yerr=S2_edep0.sigmas_keff, 
                label=f"Serpent2, edep0 pcc{S2_edep0.pcc_id}", 
                marker=".", linestyle="--")
    plt.xlabel("Time [days]")
    plt.ylabel("Keff")
    plt.title("Comparison of keffs for S2 and OpenMC, \n fission Q-values normalization with default values")
    plt.legend()
    plt.grid()
    plt.show()
    plt.savefig(f"{save_dir}/keff_comparison_S2_edep0_OMC_fissq_{integrator}.png")

    # compare keffs for S2 edepmode 0 set fissh and OpenMC fissq set
    plt.figure()
    plt.errorbar(OMC_fissq_set.BUdays, # time axis in days
                OMC_fissq_set.keffs, # keff values
                yerr=OMC_fissq_set.sigmas_keff, # statistical uncertainties on keff
                label=f"OpenMC, set fissq {OMC_fissq_set.integrator}", # labels and plot customization options
                marker=".", linestyle="--")
    plt.errorbar(S2_edep0_fissh.BUdays, 
                S2_edep0_fissh.keffs, 
                yerr=S2_edep0_fissh.sigmas_keff, 
                label=f"Serpent2, edep0, fissh pcc{S2_edep0_fissh.pcc_id}", 
                marker=".", linestyle="--")
    plt.xlabel("Time [days]")
    plt.ylabel("Keff")
    plt.title("Comparison of keffs for S2 and OpenMC, \n fission Q-values normalization with set values")
    plt.legend()
    plt.grid()
    plt.show()
    plt.savefig(f"{save_dir}/keff_comparison_S2_edep0_setfissh_OMC_setfissq_{integrator}.png")

    # compare keffs for S2 edepmode 2 and OpenMC energy_deposition
    plt.figure()
    plt.errorbar(OMC_energy_deposition.BUdays, 
                OMC_energy_deposition.keffs, 
                yerr=OMC_energy_deposition.sigmas_keff, 
                label=f"OpenMC, energy_deposition {OMC_energy_deposition.integrator}", 
                marker=".", linestyle="--")
    plt.errorbar(S2_edep2.BUdays, 
                S2_edep2.keffs, 
                yerr=S2_edep2.sigmas_keff, 
                label=f"Serpent2, edep2 pcc{S2_edep2.pcc_id}", 
                marker=".", linestyle="--")
    plt.xlabel("Time [days]")
    plt.ylabel("Keff")
    plt.title("Comparison of keffs for S2 and OpenMC, \n energy deposition normalization with MT301-MT318 + MT458 data")
    plt.legend()
    plt.grid()
    plt.show()
    plt.savefig(f"{save_dir}/keff_comparison_S2_edep2_OMC_energy_deposition_{integrator}.png")

    return

def plot_delta_Ni(S2_edep0, S2_edep0_fissh, S2_edep2, OMC_fissq, OMC_fissq_set, OMC_energy_deposition):
    # Compare isotopic densities for tracked nuclides
    integrator = OMC_fissq.integrator
    time_days = S2_edep0.BUdays
    for iso in tracked_nuclides:
        ## Compare Niso :
        delta_Ni_edep0 = [(S2_edep0.Ni[iso][i] - OMC_fissq.Ni[iso][i])*100.0/OMC_fissq.Ni[iso][i] 
                        if OMC_fissq.Ni[iso][i] != 0 else 0 for i in range(len(OMC_fissq.Ni[iso]))]
        delta_Ni_edep0_set = [(S2_edep0_fissh.Ni[iso][i] - OMC_fissq_set.Ni[iso][i])*100.0/OMC_fissq_set.Ni[iso][i]
                            if OMC_fissq_set.Ni[iso][i] != 0 else 0 for i in range(len(OMC_fissq_set.Ni[iso]))]
        delta_Ni_edep2 = [(S2_edep2.Ni[iso][i] - OMC_energy_deposition.Ni[iso][i])*100.0/OMC_energy_deposition.Ni[iso][i] 
                        if OMC_energy_deposition.Ni[iso][i] != 0 else 0 for i in range(len(OMC_energy_deposition.Ni[iso]))]

        plt.figure()
        plt.plot(time_days, delta_Ni_edep0, label="Serpent2, edep0 - OpenMC fissq", marker="x", linestyle="--")
        plt.plot(time_days, delta_Ni_edep0_set, label="Serpent2, edep0 - OpenMC fissq \n set fisson q-values", marker="x", linestyle="--")
        plt.plot(time_days, delta_Ni_edep2, label="Serpent2, edep2 - OpenMC edep", marker="x", linestyle="--")
        plt.plot(time_days, 2*np.ones_like(time_days), label="2%", linestyle="--", color="red")
        plt.plot(time_days, -2*np.ones_like(time_days), label="-2%", linestyle="--", color="red")
        plt.xlabel("Time [days]")
        plt.ylabel(f"Delta N {iso} [%]")
        plt.title(f"Delta N {iso} comparison for S2 and OpenMC, {integrator}")
        plt.legend()
        plt.grid()
        plt.savefig(f"{save_dir}/delta_{iso}_comparison_{integrator}.png")

    return

def checkGdsum_correlation(S2_edep0, S2_edep0_fissh, S2_edep2, OMC_fissq, OMC_fissq_set, OMC_energy_deposition):
    # Check correlation between Gd157 and Gd158
    time_days = S2_edep0.BUdays
    integrator = OMC_fissq.integrator
    ### Analyse correlation between Gd157 depletion and Gd158 production
    

    plt.figure()
    plt.plot(time_days, S2_edep2.Ni["Gd157"], label="Serpent2, NGd157", marker="x", linestyle="--")
    plt.plot(time_days, S2_edep2.Ni["Gd158"], label="Serpent2, NGd158", marker="x", linestyle="--")
    plt.plot(time_days, S2_edep2.Ni["Gd157"]+S2_edep2.Ni["Gd158"], label="Serpent2, sum Gd", marker="x", linestyle="--")
    plt.plot(time_days, OMC_energy_deposition.Ni["Gd157"], label="OpenMC, NGd157", marker="o", linestyle="--")
    plt.plot(time_days, OMC_energy_deposition.Ni["Gd158"], label="OpenMC, Gd158", marker="o", linestyle="--")
    plt.plot(time_days, OMC_energy_deposition.Ni["Gd157"]+OMC_energy_deposition.Ni["Gd158"], label="OpenMC, sum Gd", marker="o", linestyle="--")
    plt.xlabel("Time [days]")
    plt.ylabel("Number density [atom/b-cm]")
    #plt.yscale("log")
    plt.title(f"Evolution of Gd157 and Gd158 densities, \n edepmode 2, {integrator}")
    plt.legend()
    plt.grid()
    plt.show()
    plt.savefig(f"{save_dir}/Gd157_Gd158_correlation_edep2_{integrator}.png")

    plt.figure()
    plt.plot(time_days, S2_edep0.Ni["Gd157"], label="Serpent2, NGd157", marker="x", linestyle="--")
    plt.plot(time_days, S2_edep0.Ni["Gd158"], label="Serpent2, NGd158", marker="x", linestyle="--")
    plt.plot(time_days, S2_edep0.Ni["Gd157"]+S2_edep0.Ni["Gd158"], label="Serpent2, sum Gd", marker="x", linestyle="--")
    plt.plot(time_days, OMC_fissq.Ni["Gd157"], label="OpenMC, NGd157", marker="o", linestyle="--")
    plt.plot(time_days, OMC_fissq.Ni["Gd158"], label="OpenMC, Gd158", marker="o", linestyle="--")
    plt.plot(time_days, OMC_fissq.Ni["Gd157"]+OMC_fissq.Ni["Gd158"], label="OpenMC, sum Gd", marker="o", linestyle="--")
    plt.xlabel("Time [days]")
    plt.ylabel("Number density [atom/b-cm]")
    #plt.yscale("log")
    plt.title(f"Evolution of Gd157 and Gd158 densities, \n edep0 default fission q-values, {integrator}")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir}/Gd157_Gd158_correlation_edep0_{integrator}.png")

    plt.figure()
    plt.plot(time_days, S2_edep0_fissh.Ni["Gd157"], label="Serpent2, NGd157", marker="x", linestyle="--")
    plt.plot(time_days, S2_edep0_fissh.Ni["Gd158"], label="Serpent2, NGd158", marker="x", linestyle="--")
    plt.plot(time_days, S2_edep0_fissh.Ni["Gd157"]+S2_edep0_fissh.Ni["Gd158"], label="Serpent2, sum Gd", marker="x", linestyle="--")
    plt.plot(time_days, OMC_fissq_set.Ni["Gd157"], label="OpenMC, NGd157", marker="o", linestyle="--")
    plt.plot(time_days, OMC_fissq_set.Ni["Gd158"], label="OpenMC, Gd158", marker="o", linestyle="--")
    plt.plot(time_days, OMC_fissq_set.Ni["Gd157"]+OMC_fissq_set.Ni["Gd158"], label="OpenMC, sum Gd", marker="o", linestyle="--")
    plt.xlabel("Time [days]")
    plt.ylabel("Number density [atom/b-cm]")
    #plt.yscale("log")
    plt.title(f"Evolution of Gd157 and Gd158 densities, \n edep0 with set fission q-values, {integrator}")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir}/Gd157_Gd158_correlation_edep0_fissh_{integrator}.png")

    ### Differences on sumGd for S2 and OpenMC for edep0
    sumGd_S2_edep0 = S2_edep0.Ni["Gd157"] + S2_edep0.Ni["Gd158"]
    sumGd_OMC_fissq = OMC_fissq.Ni["Gd157"]+OMC_fissq.Ni["Gd158"]

    delta_sumGd = (sumGd_S2_edep0 - sumGd_OMC_fissq)

    plt.figure()
    plt.plot(time_days, delta_sumGd, label="Serpent2, edep0 - OpenMC", marker="x", linestyle="--")
    plt.xlabel("Time [days]")
    plt.ylabel("Delta sumGd [atom/b-cm]")
    plt.title(f"Delta sumGd, edep0 {integrator}")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir}/delta_sumGd_comparison_edep0_{integrator}.png")

    print(f"delta sum Gd max: {np.max(delta_sumGd)}")
    print(f"delta sum Gd: {delta_sumGd}")

    ### Differences on sumGd for S2 and OpenMC for edep0 with set fission q-values
    sumGd_S2_edep0_fissh = S2_edep0_fissh.Ni["Gd157"] + S2_edep0_fissh.Ni["Gd158"]
    sumGd_OMC_set_fissq = OMC_fissq_set.Ni["Gd157"]+OMC_fissq_set.Ni["Gd158"]

    delta_sumGd = (sumGd_S2_edep0_fissh - sumGd_OMC_set_fissq)

    plt.figure()
    plt.plot(time_days, delta_sumGd, label="Serpent2 - OpenMC \n normlization with set fission q-values", marker="x", linestyle="--")
    plt.xlabel("Time [days]")
    plt.ylabel("Delta sumGd [atom/b-cm]")
    plt.title(f"Delta sumGd, edep0 set fission q-values, {integrator}")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir}/delta_sumGd_comparison_edep0_set_fission_qvalues_{integrator}.png")

    print(f"delta sum Gd max: {np.max(delta_sumGd)}")
    print(f"delta sum Gd: {delta_sumGd}")

    ### Differences on sumGd for S2 and OpenMC for edep2
    sumGd_S2_edep2 = S2_edep2.Ni["Gd157"] + S2_edep2.Ni["Gd158"]
    sumGd_OMC = S2_edep2.Ni["Gd157"]+S2_edep2.Ni["Gd158"]

    delta_sumGd = (sumGd_S2_edep2 - sumGd_OMC)

    plt.figure()
    plt.plot(time_days, delta_sumGd, label="Serpent2 - OpenMC, edep2", marker="x", linestyle="--")
    plt.xlabel("Time [days]")
    plt.ylabel("Delta sumGd [atom/b-cm]")
    plt.title(f"Delta sumGd, edep2, {integrator}")
    plt.legend()
    plt.grid()
    plt.savefig(f"{save_dir}/delta_sumGd_comparison_edep2_{integrator}.png")

    print(f"delta sum Gd max: {np.max(delta_sumGd)}")
    print(f"delta sum Gd: {delta_sumGd}")

    return


def compare_energy_deposition(case_name, lib_name, OMC_integrator,tracked_nuclides,specific_power,save_dir,openMC_save_dir,S2_save_dir):
    integrator_to_pccid = {"Predictor": 0, "CELI": 1, "LE": 2, "LELI": 3, "LEQI": 4, "CECE": 6}
    S2_pcc_id = integrator_to_pccid[OMC_integrator]

    # Constructor input for S2_case obj : case_name, lib_name, edep_id, areQfissSet, isEcaptSet, pcc_id, specific_power, tracked_nuclides, save_dir
    S2_edep0 = S2_case(case_name, lib_name, 
                            edep_id=0, 
                            areQfissSet=False,
                            isEcaptSet=False, 
                            pcc_id=S2_pcc_id, 
                            specific_power=specific_power, 
                            tracked_nuclides=tracked_nuclides, 
                            save_dir=S2_save_dir)

    S2_edep0_fissh = S2_case(case_name, lib_name, 
                            edep_id=0, 
                            areQfissSet=True,
                            isEcaptSet=False, 
                            pcc_id=S2_pcc_id, 
                            specific_power=specific_power, 
                            tracked_nuclides=tracked_nuclides, 
                            save_dir=S2_save_dir)

    S2_edep2 = S2_case(case_name, lib_name, 
                            edep_id=2, 
                            areQfissSet=False,
                            isEcaptSet=False,
                            pcc_id=S2_pcc_id, 
                            specific_power=specific_power, 
                            tracked_nuclides=tracked_nuclides, 
                            save_dir=S2_save_dir)

    # OpenMC results
    OMC_fissq = OMC_case(case_name, lib_name, 
                                edep_id="fissq", 
                                areQfissSet=False, 
                                integrator=OMC_integrator, 
                                specific_power=specific_power, 
                                tracked_nuclides=tracked_nuclides, 
                                save_dir=openMC_save_dir)

    OMC_fissq_set = OMC_case(case_name, lib_name, 
                                edep_id="fissq", 
                                areQfissSet=True, 
                                integrator=OMC_integrator, 
                                specific_power=specific_power, 
                                tracked_nuclides=tracked_nuclides, 
                                save_dir=openMC_save_dir)

    OMC_energy_deposition = OMC_case(case_name, lib_name, 
                                edep_id="energy_deposition", 
                                areQfissSet=False, 
                                integrator=OMC_integrator, 
                                specific_power=specific_power, 
                                tracked_nuclides=tracked_nuclides, 
                                save_dir=openMC_save_dir)
    
    plot_keffs(S2_edep0, S2_edep0_fissh, S2_edep2, OMC_fissq, OMC_fissq_set, OMC_energy_deposition)
    plot_delta_keffs(S2_edep0, S2_edep0_fissh, S2_edep2, OMC_fissq, OMC_fissq_set, OMC_energy_deposition)
    plot_delta_Ni(S2_edep0, S2_edep0_fissh, S2_edep2, OMC_fissq, OMC_fissq_set, OMC_energy_deposition)
    checkGdsum_correlation(S2_edep0, S2_edep0_fissh, S2_edep2, OMC_fissq, OMC_fissq_set, OMC_energy_deposition)



    return


if __name__ == "__main__":
    # Set up directories for saving results
    save_dir = "OMC_S2_comparison_results"
    OpenMC_save_dir = "OMC_S2_comparison_results/OpenMC_results/"
    if not os.path.exists(OpenMC_save_dir):
        os.makedirs(OpenMC_save_dir)

    S2_save_dir = "OMC_S2_comparison_results/S2_results/"
    if not os.path.exists(S2_save_dir):
        os.makedirs(S2_save_dir)

    ### Open S2 results : HOM_Gd157_VBOC_Cst_pow_depl_edep0_pcc0 and HOM_Gd157_VBOC_Cst_pow_depl_edep2_pcc0
        # Set time integrator to simple Euler predictor and in both S2 and OpenMC, 
        # study effect of energy deposition mode in Serpent2-OpenMC comaprisons :
        # Expected results : 
        # - S2 edep0 and OpenMC fissq should be similar : fission q-values are used for constant power evolution normalization of reaction rates
        # - S2 edep2 and OpenMC energy_deposition should be similar : MT301-MT318 + (MT458 data) used for constant power evolution normalization of reaction rates
    case_name = "HOM_Gd157_VBOC_OMC"
    lib_name = "endfb8r1_pynjoy2012_kerma"
    specific_power = 38.6
    tracked_nuclides = ["Gd157", "Gd158", "U235", "U238", "Pu239", "Xe135", "Sm149"]

    compare_energy_deposition(case_name, lib_name, "Predictor", tracked_nuclides, specific_power, save_dir, OpenMC_save_dir, S2_save_dir)
    compare_energy_deposition(case_name, lib_name, "CELI", tracked_nuclides, specific_power, save_dir, OpenMC_save_dir, S2_save_dir)






