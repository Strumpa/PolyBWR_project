# Python3 post treatment script for Serpent2 output files
# makes use of SerpentTools to parse the output files
# Author : R. Guasch
# Date: 2024-12-10
# Purpose : Study spectrum, reaction rates and XS for HOM_UOX_Gd157 benchmark
# Investigate Serpent2 set ures 0/1 option.

# Importing libraries
import os
import numpy as np
import serpentTools as st
import matplotlib.pyplot as plt


# Define the path to the output files
path_to_S2_output = "/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/HOM_CELL_study/HOM_UOX_Gd157/URES_study"
energy_mesh_name = "SHEM295"
ures_opt = [0, 1]
S2_lib = ["oldlib", "PyNjoy2016"]

Gd157_102rates = {"oldlib": {"ures0":None,"ures1":None}, "PyNjoy2016": {"ures0":None,"ures1":None}}
Gd157_102XS = {"oldlib": {"ures0":None,"ures1":None}, "PyNjoy2016": {"ures0":None,"ures1":None}}
U238_102rates = {"oldlib": {"ures0":None,"ures1":None}, "PyNjoy2016": {"ures0":None,"ures1":None}}
U238_102XS = {"oldlib": {"ures0":None,"ures1":None}, "PyNjoy2016": {"ures0":None,"ures1":None}}
spectra = {"oldlib": {"ures0":None,"ures1":None}, "PyNjoy2016": {"ures0":None,"ures1":None}}


for lib in S2_lib:
    for ures in ures_opt:
        # Parse S2/PyNjoy2016 XS lib results
        # parse detector output file
        rates_spectrum_det = st.read(f"{path_to_S2_output}/HOM_UOX_Gd157_XS_{energy_mesh_name}_{lib}_ures{ures}_mc_det0.m")
        # parse XS output file
        XS_det = st.read(f"{path_to_S2_output}/HOM_UOX_Gd157_XS_{energy_mesh_name}_{lib}_ures{ures}_mc_mdx0.m")
        XS_values = XS_det.xsVal

        # recover XS values
        for key in XS_values['1'].keys():
            Gd157_XS_102 = XS_values['1'][key]
        for key in XS_values['2'].keys():
            Gd157_XS_101 = XS_values['2'][key]
        for key in XS_values['3'].keys():
            U238_XS_102 = XS_values['3'][key]
        for key in XS_values['4'].keys():
            U238_XS_101 = XS_values['4'][key]
        
        # recover reaction rates
        Gd_det = rates_spectrum_det.detectors["Gd_det"].tallies
        n_groups = Gd_det.shape[0]
        n_reactions = Gd_det.shape[1]
        Gd157_nGamma_rates = []
        Gd157_abs_rates = []
        U238_nGamma_rates = []
        U238_abs_rates = []
        for group in range(n_groups):
            Gd157_nGamma_rates.append(Gd_det[group, 0])
            Gd157_abs_rates.append(Gd_det[group, 1])
            U238_nGamma_rates.append(Gd_det[group, 2])
            U238_abs_rates.append(Gd_det[group, 3])

        # recover spectrum
        spectrum_det = rates_spectrum_det.detectors["spectrum_295G"] # spectrum_295G is the detector name in the detector file
        n_groups = spectrum_det.tallies.shape[0]
        spectrum_scores = spectrum_det.tallies

        # reconstruct the energy grid
        energy_grid = []
        spectrum_to_plot = []
        for energy_bin in spectrum_det.grids['E']: # spectrum_det.grids['E'] : list of list of values forming the energy grid
        # in each entry, the first value is the lower bound of the energy bin, the second value is the upper bound of the energy bin and the third value is the average
            energy_grid.extend([energy_bin[0], energy_bin[1]])
        for i in range(n_groups):
            spectrum_to_plot.extend([spectrum_scores[i], spectrum_scores[i]])
        energy_grid = np.array(energy_grid)
        spectrum_to_plot = np.array(spectrum_to_plot)        

        # Organize the results in dictionaries
        Gd157_102rates[lib][f"ures{ures}"] = Gd157_nGamma_rates
        U238_102rates[lib][f"ures{ures}"] = U238_nGamma_rates
        Gd157_102XS[lib][f"ures{ures}"] = Gd157_XS_102[1:]
        U238_102XS[lib][f"ures{ures}"] = U238_XS_102[1:]
        N_Gd157 = Gd157_XS_102[0]
        N_U238 = U238_XS_102[0]
        spectra[lib][f"ures{ures}"] = spectrum_to_plot

if energy_mesh_name == "SHEM295":
    n_groups = 295
# plot the spectrum using piecewise constant distribution
plt.figure()
for lib in S2_lib:
    for ures in ures_opt:
        plt.step(energy_grid, spectra[lib][f"ures{ures}"], where='post', label=f"{lib} ures{ures}")
plt.xlabel("Energy [MeV]")
plt.ylabel("Flux")
plt.title("Spectrum")
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.savefig("HOM_UOX_Gd157_spectrum_lib_ures_study.png")
plt.show()
plt.close()

# plot the reaction rates
# Gd157 MT=102 (n,gamma)
plt.figure()
for lib in S2_lib:
    for ures in ures_opt:
        rates_to_plot = []
        for i in range(n_groups):
            rates_to_plot.extend([Gd157_102rates[lib][f"ures{ures}"][i], Gd157_102rates[lib][f"ures{ures}"][i]])
        plt.step(energy_grid, rates_to_plot, where='post', label=f"{lib} ures{ures}")
plt.xlabel("Energy [MeV]")
plt.ylabel("Reaction rates")
plt.title("Gd157 $\\tau (n,\\gamma)$")
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.savefig("HOM_UOX_Gd157_nGamma_rates_Gd157_lib_ures_study.png")
plt.show()
plt.close()

# U238 MT=102 (n,gamma)
plt.figure()
for lib in S2_lib:
    for ures in ures_opt:
        rates_to_plot = []
        for i in range(n_groups):
            rates_to_plot.extend([U238_102rates[lib][f"ures{ures}"][i], U238_102rates[lib][f"ures{ures}"][i]])
        plt.step(energy_grid, rates_to_plot, where='post', label=f"{lib} ures{ures}")
plt.xlabel("Energy [MeV]")
plt.ylabel("Reaction rates")
plt.title("U238 $\\tau (n,\\gamma)$")
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.savefig("HOM_UOX_Gd157_nGamma_rates_U238_lib_ures_study.png")
plt.show()
plt.close()

# plot the XS
# Gd157 MT=102 (n,gamma)
plt.figure()
for lib in S2_lib:
    for ures in ures_opt:
        XS_to_plot = []
        for i in range(n_groups):
            XS_to_plot.extend([Gd157_102XS[lib][f"ures{ures}"][i], Gd157_102XS[lib][f"ures{ures}"][i]])
        plt.step(energy_grid, XS_to_plot, where='post', label=f"{lib} ures{ures}")
plt.xlabel("Energy [MeV]")
plt.ylabel("Cross section (barns)")
plt.title("Gd157 $\\sigma (n,\\gamma)$")
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.savefig("HOM_UOX_Gd157_nGamma_XS_Gd157_lib_ures_study.png")
plt.show() 
plt.close()

# U238 MT=102 (n,gamma)
plt.figure()
for lib in S2_lib:
    for ures in ures_opt:
        XS_to_plot = []
        for i in range(n_groups):
            XS_to_plot.extend([U238_102XS[lib][f"ures{ures}"][i], U238_102XS[lib][f"ures{ures}"][i]])
        plt.step(energy_grid, XS_to_plot, where='post', label=f"{lib} ures{ures}")
plt.xlabel("Energy [MeV]")
plt.ylabel("Cross section (barns)")
plt.title("U238 $\\sigma (n,\\gamma)$")
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.savefig("HOM_UOX_Gd157_nGamma_XS_U238_lib_ures_study.png")
plt.show()
plt.close()


# Compute errors on XS and reaction rates : for each energy bin compute the relative error between ures0 and ures1 for a given lib
Gd157_102rates_error = {"oldlib:ures1-ures0":None, "PyNjoy2016:ures1-ures0":None, "ures0:PyNjoy2016-oldlib":None, "ures1:PyNjoy2016-oldlib":None}
U238_102rates_error = {"oldlib:ures1-ures0":None, "PyNjoy2016:ures1-ures0":None, "ures0:PyNjoy2016-oldlib":None, "ures1:PyNjoy2016-oldlib":None}
Gd157_102XS_error = {"oldlib:ures1-ures0":None, "PyNjoy2016:ures1-ures0":None, "ures0:PyNjoy2016-oldlib":None, "ures1:PyNjoy2016-oldlib":None}
U238_102XS_error = {"oldlib:ures1-ures0":None, "PyNjoy2016:ures1-ures0":None, "ures0:PyNjoy2016-oldlib":None, "ures1:PyNjoy2016-oldlib":None}

# compute the relative errors for Gd157 MT=102 (n,gamma)
Gd157_102rates_error["oldlib:ures1-ures0"] = 100*(np.array(Gd157_102rates["oldlib"]["ures1"]) - np.array(Gd157_102rates["oldlib"]["ures0"])) / np.array(Gd157_102rates["oldlib"]["ures0"])
Gd157_102rates_error["PyNjoy2016:ures1-ures0"] = 100*(np.array(Gd157_102rates["PyNjoy2016"]["ures1"]) - np.array(Gd157_102rates["PyNjoy2016"]["ures0"])) / np.array(Gd157_102rates["PyNjoy2016"]["ures0"])
Gd157_102rates_error["ures0:PyNjoy2016-oldlib"] = 100*(np.array(Gd157_102rates["PyNjoy2016"]["ures0"]) - np.array(Gd157_102rates["oldlib"]["ures0"])) / np.array(Gd157_102rates["oldlib"]["ures0"])
Gd157_102rates_error["ures1:PyNjoy2016-oldlib"] = 100*(np.array(Gd157_102rates["PyNjoy2016"]["ures1"]) - np.array(Gd157_102rates["oldlib"]["ures1"])) / np.array(Gd157_102rates["oldlib"]["ures1"])

Gd157_102XS_error["oldlib:ures1-ures0"] = 100*(np.array(Gd157_102XS["oldlib"]["ures1"]) - np.array(Gd157_102XS["oldlib"]["ures0"])) / np.array(Gd157_102XS["oldlib"]["ures0"])
Gd157_102XS_error["PyNjoy2016:ures1-ures0"] = 100*(np.array(Gd157_102XS["PyNjoy2016"]["ures1"]) - np.array(Gd157_102XS["PyNjoy2016"]["ures0"])) / np.array(Gd157_102XS["PyNjoy2016"]["ures0"])
Gd157_102XS_error["ures0:PyNjoy2016-oldlib"] = 100*(np.array(Gd157_102XS["PyNjoy2016"]["ures0"]) - np.array(Gd157_102XS["oldlib"]["ures0"])) / np.array(Gd157_102XS["oldlib"]["ures0"])
Gd157_102XS_error["ures1:PyNjoy2016-oldlib"] = 100*(np.array(Gd157_102XS["PyNjoy2016"]["ures1"]) - np.array(Gd157_102XS["oldlib"]["ures1"])) / np.array(Gd157_102XS["oldlib"]["ures1"])

# compute the relative errors for U238 MT=102 (n,gamma)
U238_102rates_error["oldlib:ures1-ures0"] = 100*(np.array(U238_102rates["oldlib"]["ures1"]) - np.array(U238_102rates["oldlib"]["ures0"])) / np.array(U238_102rates["oldlib"]["ures0"])
U238_102rates_error["PyNjoy2016:ures1-ures0"] = 100*(np.array(U238_102rates["PyNjoy2016"]["ures1"]) - np.array(U238_102rates["PyNjoy2016"]["ures0"])) / np.array(U238_102rates["PyNjoy2016"]["ures0"])
U238_102rates_error["ures0:PyNjoy2016-oldlib"] = 100*(np.array(U238_102rates["PyNjoy2016"]["ures0"]) - np.array(U238_102rates["oldlib"]["ures0"])) / np.array(U238_102rates["oldlib"]["ures0"])
U238_102rates_error["ures1:PyNjoy2016-oldlib"] = 100*(np.array(U238_102rates["PyNjoy2016"]["ures1"]) - np.array(U238_102rates["oldlib"]["ures1"])) / np.array(U238_102rates["oldlib"]["ures1"])

U238_102XS_error["oldlib:ures1-ures0"] = 100*(np.array(U238_102XS["oldlib"]["ures1"]) - np.array(U238_102XS["oldlib"]["ures0"])) / np.array(U238_102XS["oldlib"]["ures0"])
U238_102XS_error["PyNjoy2016:ures1-ures0"] = 100*(np.array(U238_102XS["PyNjoy2016"]["ures1"]) - np.array(U238_102XS["PyNjoy2016"]["ures0"])) / np.array(U238_102XS["PyNjoy2016"]["ures0"])
U238_102XS_error["ures0:PyNjoy2016-oldlib"] = 100*(np.array(U238_102XS["PyNjoy2016"]["ures0"]) - np.array(U238_102XS["oldlib"]["ures0"])) / np.array(U238_102XS["oldlib"]["ures0"])
U238_102XS_error["ures1:PyNjoy2016-oldlib"] = 100*(np.array(U238_102XS["PyNjoy2016"]["ures1"]) - np.array(U238_102XS["oldlib"]["ures1"])) / np.array(U238_102XS["oldlib"]["ures1"])


# plot the relative errors
# Gd157 MT=102 (n,gamma) rates
plt.figure()
for key in Gd157_102rates_error.keys():
    errors = Gd157_102rates_error[key]
    errors_to_plot = []
    for i in range(n_groups):
        errors_to_plot.extend([errors[i], errors[i]])
    plt.step(energy_grid, errors_to_plot, where='post', label=f"{key}")
plt.xlabel("Energy [MeV]")
plt.ylabel("Relative error (%)")
plt.title("Gd157 $\\tau (n,\\gamma)$ relative errors")
plt.grid()
plt.xscale('log')
plt.legend()
plt.savefig("HOM_UOX_Gd157_nGamma_rates_Gd157_relative_errors_lib_ures_study.png")
plt.show()
plt.close()

# plot relative errors for Gd157 MT=102 (n,gamma) XS

plt.figure()
for key in Gd157_102XS_error.keys():
    errors = Gd157_102XS_error[key]
    errors_to_plot = []
    for i in range(n_groups):
        errors_to_plot.extend([errors[i], errors[i]])
    plt.step(energy_grid, errors_to_plot, where='post', label=f"{key}")
plt.xlabel("Energy [MeV]")
plt.ylabel("Relative error (%)")
plt.title("Gd157 $\\sigma (n,\\gamma)$ relative errors")
plt.grid()
plt.xscale('log')
plt.legend()
plt.savefig("HOM_UOX_Gd157_nGamma_XS_Gd157_relative_errors_lib_ures_study.png")
plt.show()
plt.close()


# U238 MT=102 (n,gamma) rates
plt.figure()
for key in U238_102rates_error.keys():
    errors = U238_102rates_error[key]
    errors_to_plot = []
    for i in range(n_groups):
        errors_to_plot.extend([errors[i], errors[i]])
    plt.step(energy_grid, errors_to_plot, where='post', label=f"{key}")
plt.xlabel("Energy [MeV]")
plt.ylabel("Relative error (%)")
plt.title("U238 $\\tau (n,\\gamma)$ relative errors")
plt.grid()
plt.xscale('log')
plt.legend()
plt.savefig("HOM_UOX_Gd157_nGamma_rates_U238_relative_errors_lib_ures_study.png")
plt.show()
plt.close()

# plot relative errors for U238 MT=102 (n,gamma) XS
plt.figure()
for key in U238_102XS_error.keys():
    errors = U238_102XS_error[key]
    errors_to_plot = []
    for i in range(n_groups):
        errors_to_plot.extend([errors[i], errors[i]])
    plt.step(energy_grid, errors_to_plot, where='post', label=f"{key}")
plt.xlabel("Energy [MeV]")
plt.ylabel("Relative error (%)")
plt.title("U238 $\\sigma (n,\\gamma)$ relative errors")
plt.grid()
plt.xscale('log')
plt.legend()
plt.savefig("HOM_UOX_Gd157_nGamma_XS_U238_relative_errors_lib_ures_study.png")
plt.show()
plt.close()

# Focus on relative errors for Gd157 MT=102 (n,gamma) rates, 
# plot individual relative errors for each lib and ures option
plt.figure()
errors = Gd157_102rates_error["oldlib:ures1-ures0"]
errors_to_plot = []
for i in range(n_groups):
    errors_to_plot.extend([errors[i], errors[i]])
plt.step(energy_grid, errors_to_plot, where='post', label="oldlib ures1-ures0")
plt.xlabel("Energy [MeV]")
plt.ylabel("Relative error (%)")
plt.title("Gd157 $\\tau (n,\\gamma)$ relative errors")
plt.grid()
plt.xscale('log')
plt.legend()
plt.savefig("HOM_UOX_Gd157_nGamma_rates_Gd157_relative_errors_oldlib_ures10_study.png")
plt.show()
plt.close()

plt.figure()
errors = Gd157_102rates_error["PyNjoy2016:ures1-ures0"]
errors_to_plot = []
for i in range(n_groups):
    errors_to_plot.extend([errors[i], errors[i]])
plt.step(energy_grid, errors_to_plot, where='post', label="PyNjoy2016 ures1-ures0")
plt.xlabel("Energy [MeV]")
plt.ylabel("Relative error (%)")
plt.title("Gd157 $\\tau (n,\\gamma)$ relative errors")
plt.grid()
plt.xscale('log')
plt.legend()
plt.savefig("HOM_UOX_Gd157_nGamma_rates_Gd157_relative_errors_PyNjoy2016_ures10_study.png")
plt.show()
plt.close()


plt.figure()
errors = Gd157_102rates_error["ures0:PyNjoy2016-oldlib"]
errors_to_plot = []
for i in range(n_groups):
    errors_to_plot.extend([errors[i], errors[i]])
plt.step(energy_grid, errors_to_plot, where='post', label="ures0 PyNjoy2016-oldlib")
plt.xlabel("Energy [MeV]")
plt.ylabel("Relative error (%)")
plt.title("Gd157 $\\tau (n,\\gamma)$ relative errors")
plt.grid()
plt.xscale('log')
plt.legend()
plt.savefig("HOM_UOX_Gd157_nGamma_rates_Gd157_relative_errors_ures0_PyNjoy2016-oldlib_study.png")
plt.show()
plt.close()

plt.figure()
errors = Gd157_102rates_error["ures1:PyNjoy2016-oldlib"]
errors_to_plot = []
for i in range(n_groups):
    errors_to_plot.extend([errors[i], errors[i]])
plt.step(energy_grid, errors_to_plot, where='post', label="ures1 PyNjoy2016-oldlib")
plt.xlabel("Energy [MeV]")
plt.ylabel("Relative error (%)")
plt.title("Gd157 $\\tau (n,\\gamma)$ relative errors")
plt.grid()
plt.xscale('log')
plt.legend()
plt.savefig("HOM_UOX_Gd157_nGamma_rates_Gd157_relative_errors_ures1_PyNjoy2016-oldlib_study.png")
plt.show()
plt.close()






