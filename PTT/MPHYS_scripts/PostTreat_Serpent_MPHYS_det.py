# Python3 script to post treat serpent2 rates output
# Author : R. Guasch
# Date : 11/10/2024
# Purpose : Post treat serpent2 rates output to validate power distribution in Multi-physics simulations

import numpy as np
import matplotlib.pyplot as plt
import serpentTools as st
import os
import sys


def parse_detector_rates(number_axial_slices, det):
    #print(det.detectors["_FUEL_1G_1"].tallies.shape)
    fiss_rates = {"U235":[], "U238":[], "U234":[]}
    n_gamma_rates = {"U235":[], "U238":[], "U234":[]}
    heat_production = {"U235":[], "U238":[], "U234":[]}
    for i in range(number_axial_slices):
        #detectors.append(det.detectors[f"_FUEL_1G_{i+1}"])
        n_gamma_rates["U235"].append(det.detectors[f"_FUEL_1G_{i+1}"].tallies[0])
        n_gamma_rates["U238"].append(det.detectors[f"_FUEL_1G_{i+1}"].tallies[1])
        n_gamma_rates["U234"].append(det.detectors[f"_FUEL_1G_{i+1}"].tallies[2])
        fiss_rates["U235"].append(det.detectors[f"_FUEL_1G_{i+1}"].tallies[3])
        fiss_rates["U238"].append(det.detectors[f"_FUEL_1G_{i+1}"].tallies[4])
        fiss_rates["U234"].append(det.detectors[f"_FUEL_1G_{i+1}"].tallies[5])
        heat_production["U235"].append(det.detectors[f"_FUEL_1G_{i+1}"].tallies[6])
        heat_production["U238"].append(det.detectors[f"_FUEL_1G_{i+1}"].tallies[7])
        heat_production["U234"].append(det.detectors[f"_FUEL_1G_{i+1}"].tallies[8])
    print(fiss_rates["U235"])           

    sum_fiss_rates = np.zeros(number_axial_slices)
    sum_n_gamma_rates = np.zeros(number_axial_slices)
    sum_heat_production = np.zeros(number_axial_slices)
    for i in range(number_axial_slices):
        for iso in fiss_rates.keys():
            sum_fiss_rates[i] += fiss_rates[iso][i]
            sum_n_gamma_rates[i] += n_gamma_rates[iso][i]
            sum_heat_production[i] += heat_production[iso][i]
    # normalize the fission rates
    fiss_rates = sum_fiss_rates/sum(sum_fiss_rates)
    n_gamma_rates = sum_n_gamma_rates/sum(sum_n_gamma_rates)
    heat_production = sum_heat_production/sum(sum_heat_production)

    return fiss_rates, n_gamma_rates, heat_production

def parse_detector_fluxes(number_axial_slices, det):

    flux_G1 = []
    flux_G2 = []

    for i in range(number_axial_slices):
        flux_G1.append(det.detectors[f"_FLUX_2G_{i+1}"].tallies[0])
        flux_G2.append(det.detectors[f"_FLUX_2G_{i+1}"].tallies[1])

    return flux_G1, flux_G2

def parse_DONJON_power(mesh_size, pow_relax, power_scaling_factor, th_relax):

    # retrieve DONJON5 power distribution
    donjon5_power = np.loadtxt(f"/home/p117902/working_dir/PolyBWR_project/Version5_wc/PyGan/Linux_aarch64/multiPhysics_PyGan_24UOX_cell/BiCG/EPRIvoidModel_Churchill_HEM1/mesh{mesh_size}_{power_scaling_factor}/Data/Power_Distrib_24UOX_mesh{mesh_size}_BiCG_EPRIvoidModel_relaxedPOW_{pow_relax}_relaxedTH_{th_relax}.txt")
    # normalize the donjon5 power distribution
    donjon5_power = donjon5_power/sum(donjon5_power)
    
    return donjon5_power

def parse_DONJON_fluxes(mesh_size, power_scaling_factor):
    path_to_DISTR = f"/home/p117902/working_dir/PolyBWR_project/Version5_wc/PyGan/Linux_aarch64/multiPhysics_PyGan_24UOX_cell/BiCG/EPRIvoidModel_Churchill_HEM1/mesh{mesh_size}_{power_scaling_factor}/Data"
    # retrieve DONJON5 flux distribution
    donjon5_flux_G1 = []
    donjon5_flux_G2 = []
    with open(f"{path_to_DISTR}/Flux01.res", "r") as f:
        donjon5_flux_file1 = f.readlines()
        for line in donjon5_flux_file1:
            if "PLANE-Z" in line:
                plane_number = int(line.split("#")[-1])
                line_index = donjon5_flux_file1.index(line)
                donjon5_flux_G1.append(float(donjon5_flux_file1[line_index+2].strip()))
    with open(f"{path_to_DISTR}/Flux02.res", "r") as f:
        donjon5_flux_file2 = f.readlines()
        for line in donjon5_flux_file2:
            if "PLANE-Z" in line:
                plane_number = int(line.split("#")[-1])
                line_index = donjon5_flux_file2.index(line)
                donjon5_flux_G2.append(float(donjon5_flux_file2[line_index+2]))
    return donjon5_flux_G1, donjon5_flux_G2

def normalize_fluxes(flux_G1, flux_G2):
    # normalize the fluxes
    flux_G1 = np.array(flux_G1)/sum(np.array(flux_G1))
    flux_G2 = np.array(flux_G2)/sum(np.array(flux_G2))
    return flux_G1, flux_G2                

def compare_results(donjon5_powers, sum_fiss_rates):
    # compare the results
    fig3, ax3 = plt.subplots()
    ax3.plot(range(number_axial_slices), sum_fiss_rates, label="S2 Fission rates")
    for i in range(len(donjon5_powers)):
        ax3.plot(range(number_axial_slices), donjon5_powers[i], label=f"Donjon5 power {i}")
    fig3.legend()
    fig3.savefig(f"AT10_24UOX_MPHYS_comparison.png")

def compute_relative_error(donjon5_rates, S2_rates):
    # calculate the relative error
    relative_error = np.abs(donjon5_rates-S2_rates)*100/S2_rates
    return relative_error

def plot_relative_error(relative_error_dict, height, type):
    fig4, ax4 = plt.subplots()
    for key in relative_error_dict.keys():
        if "10" in key:
            # Compute the boundaries and midpoints of each control volume
            z_boundaries = np.linspace(0, height, 10 + 1)
        elif "20" in key:
            z_boundaries = np.linspace(0, height, 20 + 1)
        elif "40" in key:
            z_boundaries = np.linspace(0, height, 40 + 1)
        z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes in z
    
        ax4.plot(z_values, relative_error_dict[key], label=f"Relative error {key}", marker="x", linestyle="--", linewidth=0.5)
    ax4.set_ylabel("Relative error (%)")
    ax4.set_xlabel("Height m")
    fig4.legend()
    fig4.savefig(f"AT10_24UOX_MPHYS_relative_errors_{type}.png")
    return



path_to_S2results = "/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/" # path to the serpent2 results

number_axial_slices = [10,20,40]
power_scaling_factor = [1,2] # run 4 and 8
pow_relax = 0.9
th_relax = 0.1
ERRORS_r = {}
ERRORS_G1 = {}
ERRORS_G2 = {}
for n in number_axial_slices:
    for p in power_scaling_factor:
        S2_res = f"AT10_24UOX_3D_MPHYS_mesh{n}_{p}_mc"
        det = st.read(f'{path_to_S2results}/{S2_res}_det0.m') # read the serpent2 output file
        # Parse Reaction rates (fission, n_gamma, heat production) for S2 and power for DONJON5
        sum_fiss_rates, sum_n_gamma_rates, sum_heat_prod = parse_detector_rates(n, det)
        donjon5_power = parse_DONJON_power(n, pow_relax, p, th_relax)
        relative_error = compute_relative_error(donjon5_power, sum_heat_prod)
        ERRORS_r[f"Power : Mesh {n}, power {p}"] = relative_error

        # Parse 2G fluxes for S2 and DONJON5
        S2_flux_G1, S2_flux_G2 = parse_detector_fluxes(n, det)
        donjon5_flux_G1, donjon5_flux_G2 = parse_DONJON_fluxes(n, p)
        S2_flux_G1, S2_flux_G2 = normalize_fluxes(S2_flux_G1, S2_flux_G2)
        donjon5_flux_G1, donjon5_flux_G2 = normalize_fluxes(donjon5_flux_G1, donjon5_flux_G2)
        relative_error_G1 = compute_relative_error(donjon5_flux_G2, S2_flux_G1)
        relative_error_G2 = compute_relative_error(donjon5_flux_G1, S2_flux_G2)
        ERRORS_G1[f"G1 Flux mesh {n}, power {p}"] = relative_error_G1
        ERRORS_G2[f"G2 Flux mesh {n}, power {p}"] = relative_error_G2


plot_relative_error(ERRORS_r, height=3.8, type="power")

plot_relative_error(ERRORS_G1, height=3.8, type="fluxG1")
plot_relative_error(ERRORS_G2, height=3.8, type="fluxG2")

