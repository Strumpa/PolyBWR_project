# Python3 script to post treat serpent2 rates output
# Author : R. Guasch
# Date : 11/10/2024
# Purpose : Post treat serpent2 rates output to validate power distribution in Multi-physics simulations

import numpy as np
import matplotlib.pyplot as plt
import serpentTools as st
import os
import sys

### Parsing functions
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
    #print(fiss_rates["U235"])           

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

def parse_DONJON_power(mesh_size, pow_relax, power_scaling_factor, th_relax, Ptot=40000):
    if pow_relax == False:
        power_relax_id = "non_relaxedPOW"
    else:
        power_relax_id = f"relaxedPow_{pow_relax}"
    if th_relax == False:
        th_relax_id = "non_relaxedTH"
    else:
        th_relax_id = f"relaxedTh_{th_relax}"
    # retrieve DONJON5 power distribution
    donjon5_power = np.loadtxt(f"/home/p117902/working_dir/PolyBWR_project/Version5_wc/PyGan/Linux_aarch64/multiPhysics_PyGan_24UOX_cell/BiCG/EPRIvoidModel_Churchill_HEM1/mesh{mesh_size}_{power_scaling_factor}/Data/Power_Distrib_24UOX_mesh{mesh_size}_BiCG_EPRIvoidModel_{power_relax_id}_{th_relax_id}.txt")
    # normalize the donjon5 power distribution
    donjon5_power = np.array(donjon5_power)/np.sum(donjon5_power)*Ptot # Renormalize to Ptot
    
    return donjon5_power

def parse_MPHYS_output_field(mesh_size, power_scaling_factor, field, pow_relax, th_relax):
    if pow_relax == False:
        power_relax_id = "non_relaxedPOW"
    else:
        power_relax_id = f"relaxedPOW_{pow_relax}"
    if th_relax == False:
        th_relax_id = "non_relaxedTH"
    else:
        th_relax_id = f"relaxedTH_{th_relax}"
    # retrieve DONJON5 power distribution
    field_values = np.loadtxt(f"/home/p117902/working_dir/PolyBWR_project/Version5_wc/PyGan/Linux_aarch64/multiPhysics_PyGan_24UOX_cell/BiCG/EPRIvoidModel_Churchill_HEM1/mesh{mesh_size}_{power_scaling_factor}/Data/{field}_24UOX_mesh{mesh_size}_BiCG_EPRIvoidModel_{power_relax_id}_{th_relax_id}.txt")
    return field_values

def parse_DONJON_Keffs(mesh_size, power_scaling_factor, pow_relax,th_relax):
    if pow_relax == False:
        power_relax_id = "non_relaxedPOW"
    else:
        power_relax_id = f"relaxedPow_{pow_relax}"
    if th_relax == False:
        th_relax_id = "non_relaxedTH"
    else:
        th_relax_id = f"relaxedTh_{th_relax}"
    # retrieve DONJON5 Keffs
    donjon5_keffs = np.loadtxt(f"/home/p117902/working_dir/PolyBWR_project/Version5_wc/PyGan/Linux_aarch64/multiPhysics_PyGan_24UOX_cell/BiCG/EPRIvoidModel_Churchill_HEM1/mesh{mesh_size}_{power_scaling_factor}/Data/Keffs_24UOX_mesh{mesh_size}_BiCG_EPRIvoidModel_{power_relax_id}_{th_relax_id}.txt")
    
    return donjon5_keffs

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


### Normalization functions
def normalize_fluxes(flux_G1, flux_G2):
    # normalize the fluxes to Ptot
    flux_G1 = np.array(flux_G1)/sum(np.array(flux_G1))
    flux_G2 = np.array(flux_G2)/sum(np.array(flux_G2))
    return flux_G1, flux_G2             

def normalize_S2_power(Serpent2_power, Ptot):
    Serpent2_power = np.array(Serpent2_power)/sum(np.array(Serpent2_power))
    Serpent2_power = Serpent2_power*Ptot
    return Serpent2_power   


### Homogenization/interpolation functions

def sum_neighbors(arr):
    # Ensure the input array has an even length
    if len(arr) % 2 != 0:
        raise ValueError("Input array length must be even.")
    
    # Sum every two neighboring elements
    result = np.add(arr[::2], arr[1::2])
    
    return result


# write interpolation function : take a list from larger mesh size and interpolate to get the values for the smaller mesh size


### Error computing functions
def compute_relative_error(donjon5_rates, S2_rates):
    # calculate the relative error
    relative_error = np.abs(donjon5_rates-S2_rates)*100/S2_rates
    return relative_error

def compute_quadratic_error(donjon5_rates, S2_rates):


    # calculate the quadratic error in units of field
    quadratic_error = np.sqrt(np.sum(((donjon5_rates-S2_rates))**2)/len(S2_rates))
    return quadratic_error


def error_Keffs_pcm(D5_keff, S2_keff):
    # compute delta_keff in pcm for converged Keff value
    delta_keff = (D5_keff - S2_keff)*1e5

    return delta_keff

### Plotting functions
def compare_results(field, unit, donjon5_fields, S2_field, n_slices, name, save_name):
    if " " in field:
        save_field = field.replace(" ", "_")
    else:
        save_field = field
    # compute z_mesh
    z_boundaries = np.linspace(0, 3.8, n_slices + 1)
    z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes

    # compare the results
    fig3, ax3 = plt.subplots()
    if S2_field is not None:
        ax3.plot(z_values, S2_field, label=f"S2 {field}", marker="x", linestyle="--", linewidth=0.5)
        S2_id = "S2"
    else:
        S2_id = ""
    for i in range(len(donjon5_fields)):
        ax3.plot(z_values, donjon5_fields[i], label=f"Donjon5 {field} {i}")
    fig3.legend(loc="best")
    ax3.set_title(f"Comparison of power distributions {name}")
    ax3.set_ylabel(f"{field} ({unit})")
    ax3.set_xlabel("Height (m)")
    ax3.grid()

    fig3.savefig(f"Figures_PT/AT10_24UOX_MPHYS_D5{S2_id}_{save_field}_comparison_{save_name}.png")

def plot_relative_error(relative_error_dict, height, type, power_scaling_factors_to_plot):
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
        power_scaling_factor = int(key.split(" ")[5])
        if power_scaling_factor in power_scaling_factors_to_plot:
            ax4.plot(z_values, relative_error_dict[key], label=f"Relative error {key}", marker="x", linestyle="--", linewidth=0.5)
    ax4.plot(z_values, 3.00*np.ones(len(z_values)), color="red", linestyle="--")
    ax4.set_ylabel("Relative error (%)")
    ax4.set_xlabel("Height m")
    fig4.legend(loc="best")
    fig4.savefig(f"Figures_PT/AT10_24UOX_MPHYS_relative_errors_{type}_4_8.png")
    plt.close(fig4)
    return

def plot_quadratic_errors(quadratic_error_dict, type, max_slices):
    """
    quadratic_error_dict : dictionary containing the quadratic errors for each mesh size
    type : string, type of the plot, flux or power
    max_slices : int, maximum number of axial slices
    """
    fig5, ax5 = plt.subplots()
    z_slices = np.linspace(0, max_slices, 10)
    #print(f"z_slices : {z_slices}")
    for key in quadratic_error_dict.keys():
        
        num_axial_slices = key.split(" ")[3].split(",")[0]
        #print(f"num_axial_slices : {num_axial_slices}")
        print(f"key = {key}, quadratic_error_dict[key] : {quadratic_error_dict[key]}")
        ax5.plot(float(num_axial_slices), quadratic_error_dict[key], label=f"{key}", marker="x", linestyle="--", linewidth=0.5)
    #ax5.plot(z_slices, 3.00*np.ones(len(z_slices)), color="red", linestyle="--")
    ax5.set_ylabel("Quadratic error")
    ax5.set_xlabel("Number of axial slices")
    #fig5.legend()
    fig5.savefig(f"Figures_PT/AT10_24UOX_MPHYS_quadratic_errors_{type}.png")
    plt.close(fig5)
    return

def plot_single_distr(title, field, unit, n_slices, power_distr, name, save_name):
    # compute z_mesh
    z_boundaries = np.linspace(0, 3.8, n_slices + 1)
    z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes
    if " " in field:
        save_field = field.replace(" ", "_")
    else:
        save_field = field
    # plot distribution
    fig3, ax3 = plt.subplots()
    ax3.plot(z_values, power_distr, label=f"{field} distribution {name}", marker="x", linestyle="--", linewidth=0.5)
    ax3.set_ylabel(f"{field} ({unit})")
    ax3.set_xlabel("Height (m)")
    ax3.set_title(title)
    ax3.grid()
    fig3.legend(loc="best")
    fig3.savefig(f"Figures_PT/AT10_24UOX_MPHYS_{save_field}_distribution_{save_name}.png")
    plt.close(fig3)
    return

def plot_Keff_convergence(keffs, name, save_name):
    # plot keffs
    fig6, ax6 = plt.subplots()
    ax6.plot(keffs, label=f"Keffs {name}", marker="x", linestyle="--", linewidth=0.5)
    ax6.set_ylabel("Keff")
    ax6.set_xlabel("Number of iterations")
    ax6.grid()
    fig6.legend()
    fig6.savefig(f"Figures_PT/AT10_24UOX_MPHYS_keffs_conv_{save_name}.png")
    plt.close(fig6)
    return


def plot_keff_error_spatial_convergence(delta_keffs_dict, name, save_name):
    """
    idea is to see how delta_Keff evolves with the number of axial slices increasing
    """
    fig7, ax7 = plt.subplots()
    for key in delta_keffs_dict.keys():
        if "10" in key:
            n = 10
        elif "20" in key:
            n = 20
        elif "40" in key:
            n = 40
        elif "70" in key:
            n = 70
        elif "80" in key:
            n = 80
        elif "160" in key:
            n = 160
        ax7.plot(n, delta_keffs_dict[key], marker="x", linestyle="--", linewidth=0.5)
    ax7.set_ylabel("Delta Keff (pcm)")
    ax7.set_xlabel("Number of axial slices")
    ax7.set_title(f"Delta Keffs {name}")
    ax7.grid()
    fig7.legend()
    fig7.savefig(f"Figures_PT/AT10_24UOX_MPHYS_delta_keffs_{save_name}.png")
    plt.close(fig7)

    return

def plot_keff_D5S2_spatial_convergence(n_slices, D5_keffs, S2_keffs, name, save_name):
    """
    idea is to see how Keff evolves with the number of axial slices increasing : should reach spatial convergence !
    """

    fig7, ax7 = plt.subplots()
    ax7.plot(n_slices, D5_keffs, label=f"Keffs D5 {name}", marker="x", linestyle="--", linewidth=0.5)
    ax7.plot(n_slices, S2_keffs, label=f"Keffs S2 {name}", marker="*", linestyle="--", linewidth=0.5)
    ax7.set_ylabel("Keff")
    ax7.set_xlabel("Number of axial slices")
    ax7.set_title(f"Keffs {name}")
    ax7.grid()
    fig7.legend()
    fig7.savefig(f"Figures_PT/AT10_24UOX_MPHYS_D5S2_keffs_{save_name}.png")
    plt.close(fig7)

    return

def plot_nodal_errors(n_slices_list, nodal_errors_dict, field, unit):
    """
    plot the nodal errors for each number of axial slices, error dict built from comparison with S2
    """
    z_meshes = []
    for n in n_slices_list:
        z_boundaries = np.linspace(0, 3.8, n + 1)
        z_values = (z_boundaries[:-1] + z_boundaries[1:]) / 2  # Midpoints of control volumes
        z_meshes.append(z_values)
    fig8, ax8 = plt.subplots()
    for i in range(len(n_slices_list)):
        key = f"Nodal error on power : {n_slices_list[i]} axial slices"
        ax8.plot(z_meshes[i], nodal_errors_dict[key], label=f"{key}", marker="x", linestyle="--", linewidth=0.5)
    ax8.set_ylabel(f"Nodal error on {field} {unit}")
    ax8.set_xlabel("Height (m)")
    ax8.grid()
    fig8.legend()
    fig8.savefig(f"Figures_PT/AT10_24UOX_MPHYS_nodal_errors_{field}.png")
    plt.close(fig8)
    return

def plot_quadratic_errors(n_slices_list, quadratic_error_dict):
    """
    plot the quadratic errors for each number of axial slices
    """
    fig9, ax9 = plt.subplots()
    for i in range(len(n_slices_list)):
        key = f"RMS error on power : {n_slices_list[i]} axial slices"
        ax9.plot(n_slices_list[i], quadratic_error_dict[key], label=f"{key}", marker="x", linestyle="--", linewidth=0.5)
    ax9.set_ylabel("RMS error on power")
    ax9.set_xlabel("Number of axial slices")
    ax9.grid()
    fig9.legend()
    fig9.savefig(f"Figures_PT/AT10_24UOX_MPHYS_quadratic_errors_power.png")
    plt.close(fig9)
    return

def compare_all_numbers_slices(field, unit, donjon_field_list, serpent_field_list, n_slices_list):
    """
    Compare the distributions for all the number of axial slices
    """
    if " " in field:
        save_field = field.replace(" ", "_")
    else:
        save_field = field
    fig10, ax10 = plt.subplots()
    for i in range(len(n_slices_list)):
        ax10.plot(donjon_field_list[i], label=f"Donjon5 {field} {n_slices_list[i]} axial slices")
        ax10.plot(serpent_field_list[i], label=f"S2 {field} {n_slices_list[i]} axial slices")
    ax10.set_ylabel(f"{field} ({unit})")
    ax10.set_xlabel("Height (m)")
    ax10.grid()
    fig10.legend()
    fig10.savefig(f"Figures_PT/AT10_24UOX_MPHYS_D5S2_comparison_{save_field}_all_slices.png")
    plt.close(fig10)
    return


### LaTeX table printing functions time !!!
def print_latex_table(n_slices_list, quadratic_error_dict):
    """
    print a latex table with the quadratic errors for each number of axial slices
    """
    print("\\begin{table}[H]")
    print("\\centering")
    print("\\begin{tabular}{|c|c|}")
    print("\\hline")
    print("Number of axial slices & RMS error on power \\\\")
    print("\\hline")
    for i in range(len(n_slices_list)):
        key = f"RMS error on power : {n_slices_list[i]} axial slices"
        print(f"{n_slices_list[i]} & {quadratic_error_dict[key]} \\\\")
    print("\\hline")
    print("\\end{tabular}")
    print("\\end{table}")
    return

# Post treatment and checker functions

def post_treat_D5vsS2(power_scaling_factor, number_axial_slices, pow_relax, th_relax):

    NODAL_errors = {}
    quadratic_errors_Power = {}

    delta_keffs = {}

    donjon_powers = []
    serpent_powers = []

    D5_keffs = []
    #D5_keffs_dict = {}
    S2_keffs = []


    for p in power_scaling_factor:
        for n in number_axial_slices:
            if p == 1:
                power = "40 kW"
                ptot = 40000
            elif p == 2:
                power = "20 kW"
                ptot = 20000
            elif p == 4:
                power = "10 kW"
                ptot = 10000
            elif p == 8:
                power = "5 kW"
                ptot = 5000
            # Recover S2 data
            path_to_S2results = f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64"#/MPHYS/mesh{n}" # path to the serpent2 results
            S2_res = f"AT10_24UOX_MPHYS_mesh{n}_{p}_mc"
            # Recover S2 Keff
            res = st.read(f'{path_to_S2results}/{S2_res}_res.m') # read the serpent2 output file

            serpent_keff=res.resdata["absKeff"]
            S2_keffs.append(serpent_keff[0])
            #print(f"S2 Keff : {serpent_keff}")

            # Recover S2 detector data
            det = st.read(f'{path_to_S2results}/{S2_res}_det0.m') # read the serpent2 output file
            # Parse Reaction rates (fission, n_gamma, heat production) for S2 and power for DONJON5
            sum_fiss_rates, sum_n_gamma_rates, sum_heat_prod = parse_detector_rates(n, det)
            # Normalize the S2 power distribution
            Serpent2_power = normalize_S2_power(sum_heat_prod, ptot)
            serpent_powers.append(Serpent2_power)
            # plot the power distribution
            plot_single_distr(f"Axial Serpent2 power distribution, normalized to {power}", "Power", "W", n, Serpent2_power, f"S2 : {n} axial slices", f"S2_{n}_axial_slices_{p}")

            # Recover DONJON5 data
            # Parse DONJON5 Keffs
            DONJON_keffs = parse_DONJON_Keffs(n,  p, pow_relax, th_relax)
            #print(f"DONJON5 Keffs : {DONJON_keffs}")
            D5_keffs.append(DONJON_keffs[-1])
            #D5_keffs_dict[f"Keff : {n} axial slices, power = {power}"] = DONJON_keffs[-1]

            # Compare Keffs
            delta_keffs[f"Delta Keffs : {n} axial slices"] = error_Keffs_pcm(DONJON_keffs[-1], serpent_keff[0])

            # multiphysics convergence of Keffs :
            plot_Keff_convergence(DONJON_keffs, f"D5 : {n} axial slices", f"D5_{n}_axial_slices")


            # Parse DONJON5 power distribution
            donjon5_power = parse_DONJON_power(n, pow_relax, p, th_relax)
            donjon_powers.append(donjon5_power)
            # plot the power distribution
            plot_single_distr(f"Axial DONJON5 power distribution, normalized to {power}", "Power", "W", n, donjon5_power, f"D5 : {n} axial slices", f"D5_{n}_axial_slices_{p}")


            # Compare D5-S2 results        
            #relative_error = compute_relative_error(donjon5_power, sum_heat_prod)
            compare_results("Power", "unit", [donjon5_power], Serpent2_power, n, f"{n} axial slices", f"{n}_axial_slices")
            nodal_error = donjon5_power - Serpent2_power
            NODAL_errors[f"Nodal error on power : {n} axial slices"] = nodal_error
            quadratic_error = compute_quadratic_error(donjon5_power, Serpent2_power)
            quadratic_errors_Power[f"RMS error on power : {n} axial slices"] = quadratic_error
            

            # Compare axial power distributions
            #ERRORS_r[f"Power : Mesh {n}, power {p}"] = relative_error
            #quadratic_error_r[f"Power : Mesh {n}, power {p}"] = compute_quadratic_error(donjon5_power, sum_heat_prod)

        plot_keff_D5S2_spatial_convergence(number_axial_slices, D5_keffs, S2_keffs, f"spatial convergence, p = {p}", f"spatial_convergence_p{p}")
    # plot keff spatial convergence
    plot_keff_error_spatial_convergence(delta_keffs, "spatial convergence", "spatial_convergence")

    # plot nodal errors on power for each node
    plot_nodal_errors(number_axial_slices, NODAL_errors, field="power", unit="W")

    # plot quadratic errors on power, see how it evolves with the number of axial slices
    plot_quadratic_errors(number_axial_slices, quadratic_errors_Power)

    # plot keff spatial convergence

    
    return

def consistency_check_Qfiss(n,p):
    Qfiss_distribs = []
    # Recover MPHYS outputs for n=10, p=1
    relax_options = [(False,False), (False, 0.9), (0.5, 0.9), (0.5, False), (0.9, 0.9), (0.5, 0.8), (0.5,0.7), (0.5,0.5), (0.5,0.2), (0.2,0.2)]
    for relax in relax_options:
        Qfiss = parse_MPHYS_output_field(n, p, "Qfiss", relax[0], relax[1])
        Qfiss_distribs.append(Qfiss)
        # plot the power density distribution
        #plot_single_distr("Power density", "W/m3", n, Qfiss, f"Qfiss : {n} axial slices", f"Qfiss_{n}_axial_slices_{relax[0]}_{relax[1]}")

    compare_results("Power density", "W/m3", Qfiss_distribs, None, n, f"{n} axial slices", f"{n}_axial_slices")

    diffQfiss = []
    for i in range(1,len(Qfiss_distribs)):
        diffQfiss.append((Qfiss_distribs[i] - Qfiss_distribs[0])*100/Qfiss_distribs[0])

    fig,ax = plt.subplots()
    for i in range(len(diffQfiss)):
        ax.plot(diffQfiss[i], label=f"diffQfiss_{i}")
    ax.set_ylabel("Relative difference (%)")
    fig.legend()
    fig.savefig(f"Figures_PT/MPHYS_diffs_Qdiff_mesh{n}_{p}_RELAXStudy.png")
    plt.close(fig)

def consistency_check_TH_params(n,p):
    TF_list = []
    TC_list = []
    DC_list = []

    # Recover MPHYS outputs for n=10, p=1
    relax_options = [(False,False), (False, 0.9), (0.9, False), (0.9, 0.9), (0.5, 0.8), (0.5,0.7), (0.5,0.5), (0.5,0.2), (0.2,0.2)]
    for relax in relax_options:
        TF = parse_MPHYS_output_field(n, p, "TeffFuel", relax[0], relax[1])
        TC = parse_MPHYS_output_field(n, p, "Twater", relax[0], relax[1])
        DC = parse_MPHYS_output_field(n, p, "rho", relax[0], relax[1])
        TF_list.append(TF)
        TC_list.append(TC)
        DC_list.append(DC)
        # plot the power density distribution
        #plot_single_distr("TFuel", "K", n, TF, f"TF : {n} axial slices", f"TF_{n}_axial_slices_{relax[0]}_{relax[1]}")
        #plot_single_distr("TCoolant", "K", n, TC, f"TC : {n} axial slices", f"TC_{n}_axial_slices_{relax[0]}_{relax[1]}")
        #plot_single_distr("DCoolant", "kg/m3", n, DC, f"DC : {n} axial slices", f"DC_{n}_axial_slices_{relax[0]}_{relax[1]}")

    compare_results("TFuel", "K", TF_list, None, n, f"{n} axial slices", f"{n}_axial_slices")
    compare_results("TCoolant", "K", TC_list, None, n, f"{n} axial slices", f"{n}_axial_slices")
    compare_results("DCoolant", "kg/m3", DC_list, None, n, f"{n} axial slices", f"{n}_axial_slices")

    diffTF = []
    diffTC = []
    diffDC = []
    for i in range(1,len(TF_list)):
        
        diffTF.append((TF_list[i] - TF_list[0])*100/TF_list[0])
        diffTC.append((TC_list[i] - TC_list[0])*100/TC_list[0])
        diffDC.append((DC_list[i] - DC_list[0])*100/DC_list[0])

    fig,ax = plt.subplots()
    for i in range(len(diffTF)):
        ax.plot(diffTF[i], label=f"diffTF_{i}")
    ax.set_ylabel("Relative difference (%)")
    fig.legend()
    fig.savefig(f"Figures_PT/MPHYS_diffs_TF_mesh{n}_{p}_RELAXStudy.png")
    plt.close(fig)

    fig,ax = plt.subplots()
    for i in range(len(diffTC)):
        ax.plot(diffTC[i], label=f"diffTC_{i+1}")
    ax.set_ylabel("Relative difference (%)")
    fig.legend()
    fig.savefig(f"Figures_PT/MPHYS_diffs_TC_mesh{n}_{p}_RELAXStudy.png")
    plt.close(fig)

    fig,ax = plt.subplots()
    for i in range(len(diffDC)):
        ax.plot(diffDC[i], label=f"diffDC_{i+1}")
    ax.set_ylabel("Relative difference (%)")
    fig.legend()
    fig.savefig(f"Figures_PT/MPHYS_diffs_DC_mesh{n}_{p}_RELAXStudy.png")
    plt.close(fig)

    return


if __name__ == "__main__" :

    number_axial_slices = [10, 20, 40, 70, 80, 160]
    power_scaling_factor = [1] #, 2, 4, 8] # run 4 and 8
    
    
    pow_relax = False #pow_relax = [0.2,0.5,0.8, False]
    th_relax = False #th_relax = [0.2,0.5,0.8, False]

    post_treat_D5vsS2(power_scaling_factor, number_axial_slices, pow_relax, th_relax)
    consistency_check_Qfiss(n=10, p=1)
    consistency_check_TH_params(n=10,p=1)
    consistency_check_TH_params(n=20,p=1)




    




