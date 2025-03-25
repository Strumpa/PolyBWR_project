### Main comparison script used to run S2 vsOpenMC comaprisons


import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import serpentTools as st
from post_process import Serpent2_case as S2_case


### For now just open txt files for openMC results

path_to_openmc_res = f"{os.environ['OPENMC_RESULTS']}/HOM_Gd157_VBOC/energy_deposition_mode/HOM_Gd157_VBOC_Predictor"
openMC_save_dir = "OMC_S2_comparison_results/OpenMC_results/"

time_days = np.loadtxt(f"{path_to_openmc_res}/HOM_Gd157_Cst_pow_depl_time_days.txt")
keffs = np.loadtxt(f"{path_to_openmc_res}/HOM_Gd157_Cst_pow_depl_keff.txt")
# keffs = n_BU_steps x 2 numpy array with keffs in first column and sigmas in second column
# extract keffs and sigmas
OpenMC_keffs = keffs[:,0]
OpenMC_sigmas_keffs = keffs[:,1]
NGd157 = np.loadtxt(f"{path_to_openmc_res}/HOM_Gd157_Cst_pow_depl_NGd157.txt") # atoms
NGd158 = np.loadtxt(f"{path_to_openmc_res}/HOM_Gd157_Cst_pow_depl_NGd158.txt") # atoms
NU235 = np.loadtxt(f"{path_to_openmc_res}/HOM_Gd157_Cst_pow_depl_NU235.txt") # atoms
vol = 1*1.295**2 # cm3
NGd157 = NGd157*1e-24/vol # atom/b-cm
NGd158 = NGd158*1e-24/vol # atom/b-cm
NU235 = NU235*1e-24/vol # atom/b-cm
#NU238 = np.loadtxt(f"{path_to_openmc_res}/HOM_Gd157_Cst_pow_depl_NU238.txt")
#NU238 = np.loadtxt(f"{path_to_openmc_res}/HOM_Gd157_VBOC_Cst_pow_depl_NU238.txt")


# Check shapes of numpy arrays
print(time_days.shape)
print(keffs.shape)
print(NGd157.shape)
print(NGd158.shape)
print(NU235.shape)
#print(NU238.shape)

### Open S2 results : HOM_Gd157_VBOC_Cst_pow_depl_edep0_pcc2

case_name = "HOM_Gd157_VBOC"
lib_name = "endfb8r1_pynjoy2012_kerma"
specific_power = 38.6
tracked_nuclides = ["Gd157", "Gd158", "U235"]
S2_save_dir = "OMC_S2_comparison_results/S2_results" 
if not os.path.exists(S2_save_dir):
    os.makedirs(S2_save_dir)

S2_edep0_pcc0 = S2_case(case_name, lib_name, 0, False, 0, specific_power, tracked_nuclides, S2_save_dir)
S2_edep0_pcc2 = S2_case(case_name, lib_name, 0, False, 2, specific_power, tracked_nuclides, S2_save_dir)
S2_edep2_pcc0 = S2_case(case_name, lib_name, 2, False, 0, specific_power, tracked_nuclides, S2_save_dir)
S2_edep2_pcc2 = S2_case(case_name, lib_name, 2, False, 2, specific_power, tracked_nuclides, S2_save_dir)
print(S2_edep0_pcc2.BU.shape)
print(f"S2 BU in days: {S2_edep0_pcc2.BU/S2_edep0_pcc2.specific_power}")

## Prepare comparison

save_dir = "OMC_S2_comparison_results"


# compare keffs 

plt.figure()
plt.errorbar(time_days, OpenMC_keffs, yerr=OpenMC_sigmas_keffs, label="OpenMC", marker="o", linestyle="--")
plt.errorbar(S2_edep0_pcc2.BUdays, S2_edep0_pcc2.keffs, yerr=S2_edep0_pcc2.sigmas_keff, label="Serpent2, edep0 pcc2", marker="x", linestyle="--")
plt.errorbar(S2_edep2_pcc2.BUdays, S2_edep2_pcc2.keffs, yerr=S2_edep2_pcc2.sigmas_keff, label="Serpent2, edep2 pcc2", marker="x", linestyle="--")
plt.xlabel("Time [days]")
plt.ylabel("Keff")
plt.legend()
plt.grid()
plt.show()
plt.savefig(f"{save_dir}/keff_comparison.png")

# delta keff (pcm) comparison
## Might need to interpolate to compare with edepmode 2 results ?

delta_keff_edep0_pcc2 = (S2_edep0_pcc2.keffs - OpenMC_keffs)*1e5
delta_keff_edep2_pcc2 = (S2_edep2_pcc2.keffs - OpenMC_keffs)*1e5
delta_keff_edep0_pcc0 = (S2_edep0_pcc0.keffs - OpenMC_keffs)*1e5
delta_keff_edep2_pcc0 = (S2_edep2_pcc0.keffs - OpenMC_keffs)*1e5
plt.figure()
plt.plot(time_days, delta_keff_edep0_pcc2, label="Serpent2, edep0 pcc2 - OpenMC", marker="x", linestyle="--")
plt.plot(time_days, delta_keff_edep2_pcc2, label="Serpent2, edep2 pcc2 - OpenMC", marker="x", linestyle="--")
plt.plot(time_days, delta_keff_edep0_pcc0, label="Serpent2, edep0 pcc0 - OpenMC", marker="x", linestyle="--")
plt.plot(time_days, delta_keff_edep2_pcc0, label="Serpent2, edep2 pcc0 - OpenMC", marker="x", linestyle="--")
plt.plot(time_days, 200*np.ones_like(time_days), label="200 pcm", linestyle="--", color="red")
plt.plot(time_days, -200*np.ones_like(time_days), label="-200 pcm", linestyle="--", color="red")
plt.xlabel("Time [days]")
plt.ylabel("Delta Keff [pcm]")
plt.legend()
plt.grid()
plt.show()
plt.savefig(f"{save_dir}/delta_keff_comparison.png")

print(f"max delta keff S2 edep2 pcc0 vs OpenMC = {np.max(np.abs(delta_keff_edep2_pcc0))} pcm")
print(f"RMS delta keff S2 edep2 pcc0 vs OpenMC = {np.sqrt(np.mean(delta_keff_edep2_pcc0**2))} pcm")
print(f"avg delta keff S2 edep2 pcc0 vs OpenMC = {np.mean(delta_keff_edep2_pcc0)} pcm")
print(f"delta keff S2 edep2 pcc0 vs OpenMC = {delta_keff_edep2_pcc0} pcm")



## Compare NGd157 : 
delta_NGd157_edep0_pcc2 = (S2_edep0_pcc2.Ni["Gd157"] - NGd157)*100.0/NGd157
delta_NGd157_edep0_pcc0 = (S2_edep0_pcc0.Ni["Gd157"] - NGd157)*100.0/NGd157
delta_NGd157_edep2_pcc2 = (S2_edep2_pcc2.Ni["Gd157"] - NGd157)*100.0/NGd157
delta_NGd157_edep2_pcc0 = (S2_edep2_pcc0.Ni["Gd157"] - NGd157)*100.0/NGd157

plt.figure()
plt.plot(time_days, delta_NGd157_edep0_pcc2, label="Serpent2, edep0 pcc2 - OpenMC", marker="x", linestyle="--")
plt.plot(time_days, delta_NGd157_edep2_pcc2, label="Serpent2, edep2 pcc2 - OpenMC", marker="x", linestyle="--")
plt.plot(time_days, delta_NGd157_edep0_pcc0, label="Serpent2, edep0 pcc0 - OpenMC", marker="x", linestyle="--")
plt.plot(time_days, delta_NGd157_edep2_pcc0, label="Serpent2, edep2 pcc0 - OpenMC", marker="x", linestyle="--")
plt.plot(time_days, 2*np.ones_like(time_days), label="2%", linestyle="--", color="red")
plt.plot(time_days, -2*np.ones_like(time_days), label="-2%", linestyle="--", color="red")
plt.xlabel("Time [days]")
plt.ylabel("Delta NGd157 [%]")
plt.legend()
plt.grid()
plt.show()
plt.savefig(f"{save_dir}/delta_NGd157_comparison.png")

# plot NGd157 for S2 edep2 pcc0 and OpenMC
plt.figure()
plt.plot(time_days, S2_edep2_pcc0.Ni["Gd157"], label="Serpent2, edep2 pcc0", marker="x", linestyle="--")
plt.plot(time_days, NGd157, label="OpenMC", marker="o", linestyle="--")
plt.xlabel("Time [days]")
plt.ylabel("NGd157 [atom/b-cm]")
plt.yscale("log")
plt.legend()
plt.grid()
plt.show()
plt.savefig(f"{save_dir}/NGd157_comparison.png")


## Compare NGd158 :
delta_NGd158_edep0_pcc2 = [(S2_edep0_pcc2.Ni["Gd158"][i] - NGd158[i])*100.0/NGd158[i] if NGd158[i] != 0 else 0 for i in range(len(NGd158))]
delta_NGd158_edep0_pcc0 = [(S2_edep0_pcc0.Ni["Gd158"][i] - NGd158[i])*100.0/NGd158[i] if NGd158[i] != 0 else 0 for i in range(len(NGd158))]
delta_NGd158_edep2_pcc2 = [(S2_edep2_pcc2.Ni["Gd158"][i] - NGd158[i])*100.0/NGd158[i] if NGd158[i] != 0 else 0 for i in range(len(NGd158))]
delta_NGd158_edep2_pcc0 = [(S2_edep2_pcc0.Ni["Gd158"][i] - NGd158[i])*100.0/NGd158[i] if NGd158[i] != 0 else 0 for i in range(len(NGd158))]

plt.figure()
plt.plot(time_days, delta_NGd158_edep0_pcc2, label="Serpent2, edep0 pcc2 - OpenMC", marker="x", linestyle="--")
plt.plot(time_days, delta_NGd158_edep2_pcc2, label="Serpent2, edep2 pcc2 - OpenMC", marker="x", linestyle="--")
plt.plot(time_days, delta_NGd158_edep0_pcc0, label="Serpent2, edep0 pcc0 - OpenMC", marker="x", linestyle="--")
plt.plot(time_days, delta_NGd158_edep2_pcc0, label="Serpent2, edep2 pcc0 - OpenMC", marker="x", linestyle="--")
plt.plot(time_days, 2*np.ones_like(time_days), label="2%", linestyle="--", color="red")
plt.plot(time_days, -2*np.ones_like(time_days), label="-2%", linestyle="--", color="red")
plt.xlabel("Time [days]")
plt.ylabel("Delta NGd158 [%]")
plt.legend()
plt.grid()
plt.show()
plt.savefig(f"{save_dir}/delta_NGd158_comparison.png")

# plot NGd158 for S2 edep2 pcc0 and OpenMC
plt.figure()
plt.plot(time_days, S2_edep2_pcc0.Ni["Gd158"], label="Serpent2, edep2 pcc0", marker="x", linestyle="--")
plt.plot(time_days, NGd158, label="OpenMC", marker="o", linestyle="--")
plt.xlabel("Time [days]")
plt.ylabel("NGd158 [atom/b-cm]")
plt.yscale("log")
plt.legend()
plt.grid()
plt.show()
plt.savefig(f"{save_dir}/NGd158_comparison.png")



## Compare NU235 :
delta_NU235_edep0_pcc2 = (S2_edep0_pcc2.Ni["U235"] - NU235)*100.0/NU235
delta_NU235_edep0_pcc0 = (S2_edep0_pcc0.Ni["U235"] - NU235)*100.0/NU235
delta_NU235_edep2_pcc2 = (S2_edep2_pcc2.Ni["U235"] - NU235)*100.0/NU235
delta_NU235_edep2_pcc0 = (S2_edep2_pcc0.Ni["U235"] - NU235)*100.0/NU235

plt.figure()
plt.plot(time_days, delta_NU235_edep0_pcc2, label="Serpent2, edep0 pcc2 - OpenMC", marker="x", linestyle="--")
plt.plot(time_days, delta_NU235_edep2_pcc2, label="Serpent2, edep2 pcc2 - OpenMC", marker="x", linestyle="--")
plt.plot(time_days, delta_NU235_edep0_pcc0, label="Serpent2, edep0 pcc0 - OpenMC", marker="x", linestyle="--")
plt.plot(time_days, delta_NU235_edep2_pcc0, label="Serpent2, edep2 pcc0 - OpenMC", marker="x", linestyle="--")
plt.plot(time_days, 2*np.ones_like(time_days), label="2%", linestyle="--", color="red")
plt.plot(time_days, -2*np.ones_like(time_days), label="-2%", linestyle="--", color="red")
plt.xlabel("Time [days]")
plt.ylabel("Delta NU235 [%]")
plt.legend()
plt.grid()
plt.show()
plt.savefig(f"{save_dir}/delta_NU235_comparison.png")

### Analyse correlation between Gd157 depletion and Gd158 production

plt.figure()
plt.plot(time_days, S2_edep2_pcc0.Ni["Gd157"], label="Serpent2, edep0 pcc2", marker="x", linestyle="--")
plt.plot(time_days, S2_edep2_pcc0.Ni["Gd158"], label="Serpent2, edep0 pcc2", marker="x", linestyle="--")
plt.plot(time_days, S2_edep2_pcc0.Ni["Gd157"]+S2_edep2_pcc0.Ni["Gd158"], label="Serpent2, sum Gd", marker="x", linestyle="--")
plt.plot(time_days, NGd157, label="OpenMC, Gd157", marker="o", linestyle="--")
plt.plot(time_days, NGd158, label="OpenMC, Gd158", marker="o", linestyle="--")
plt.plot(time_days, NGd157+NGd158, label="OpenMC, sum Gd", marker="o", linestyle="--")
plt.xlabel("Time [days]")
plt.ylabel("Number density [atom/b-cm]")
#plt.yscale("log")
plt.legend()
plt.grid()
plt.show()
plt.savefig(f"{save_dir}/Gd157_Gd158_correlation.png")

### Differences on sumGd for S2 and OpenMC

sumGd_S2_edep2_pcc0 = S2_edep2_pcc0.Ni["Gd157"] + S2_edep2_pcc0.Ni["Gd158"]
sumGd_OMC = NGd157 + NGd158

delta_sumGd = (sumGd_S2_edep2_pcc0 - sumGd_OMC)*100.0/sumGd_OMC

plt.figure()
plt.plot(time_days, delta_sumGd, label="Serpent2, edep2 pcc0 - OpenMC", marker="x", linestyle="--")
plt.plot(time_days, 2*np.ones_like(time_days), label="2%", linestyle="--", color="red")
plt.plot(time_days, -2*np.ones_like(time_days), label="-2%", linestyle="--", color="red")
plt.xlabel("Time [days]")
plt.ylabel("Delta sumGd [%]")
plt.legend()
plt.grid()
plt.show()
plt.savefig(f"{save_dir}/delta_sumGd_comparison.png")

print(f"delta sum Gd max: {np.max(delta_sumGd)} %")
print(f"delta sum Gd: {delta_sumGd} %")
