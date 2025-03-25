## Python3 post treatment script for the outputs of the Serpent2 code
# Purpose : Read and compare outputs of different Serpent2 calculations, based on different input ACE files.


import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import serpentTools as st

# Define the Serpent2_case class

class Serpent2_case:
    def __init__(self, case_name, lib_name, edep_id, isEcaptSet, pcc_id, specific_power, tracked_nuclides, save_dir):
        self.case_name = case_name
        self.lib_name = lib_name        
        self.isEcaptSet = isEcaptSet # Option to be activated when energy deposition mode is set to 1 in Serpent2 and Ecapt has been manually set in Serpent2.
        if self.isEcaptSet:
            self.edep_id = f"{edep_id}_Ecapt"
        else:
            self.edep_id = edep_id
        self.pcc_id = pcc_id
        self.specific_power = specific_power
        self.tracked_nuclides = tracked_nuclides
        self.save_dir = save_dir
        
        self.unitsBU = "MWd/kgU" # Default Serpent2 units for BU

        if self.edep_id == 0:
            self.edep_mode = "Constant energy deposition per fission"
        elif self.edep_id == 1:
            self.edep_mode = "Local energy deposition based on ENDF MT458 data"
        elif self.edep_id == 2:
            self.edep_mode = "Local photon energy deposition, using KERMA coefficients"
        
        if self.pcc_id == 0:
            self.BUScheme = "CE"
        elif self.pcc_id == 1:
            self.BUScheme = "CELI"
        elif self.pcc_id == 2:
            self.BUScheme = "LE"
        elif self.pcc_id == 3:
            self.BUScheme = "LELI"
        elif self.pcc_id == 4:
            self.BUScheme = "LEQI"
        elif self.pcc_id == "6":
            self.BUScheme = "CECE"

        self.read_S2_outputs()
        self.rescale_BU("MWd/tU")
        

    
    def read_S2_outputs(self):
        # Read the results and depletion files
        self.res = st.read(f"{os.environ['SERPENT_RESULTS']}/{self.case_name}/{self.case_name}_{self.lib_name}_edep{self.edep_id}_pcc{self.pcc_id}_mc_res.m")
        self.depl = st.read(f"{os.environ['SERPENT_RESULTS']}/{self.case_name}/{self.case_name}_{self.lib_name}_edep{self.edep_id}_pcc{self.pcc_id}_mc_dep.m")
        
        self.keffs = self.res.resdata["absKeff"].T[0]
        self.sigmas_keff = self.res.resdata["absKeff"].T[1]
        self.mat = self.depl.materials['total']
        self.BU = self.mat.burnup
        data_frame = self.mat.toDataFrame("adens", names = self.tracked_nuclides, time = "days")
        self.Ni = {}
        for isotope in data_frame.columns:
            self.BUdays = data_frame.index.to_numpy()
            self.Ni[isotope] = data_frame[isotope].to_numpy()
        return 
    
    def rescale_BU(self, units):
        if units == "MWd/tU" and self.unitsBU == "MWd/kgU":
            self.BU = self.BU*1000
            self.unitsBU = "MWd/tU"


# Define the Serpent2_comparison class
            
class Serpent2_comparison:
    def __init__(self, comparison_name, cases_to_test, case_ref, tracked_nulides, save_dir):
        self.comparison_name = comparison_name
        self.cases_to_test = cases_to_test
        self.case_ref = case_ref
        self.tracked_nulides = tracked_nulides
        self.save_dir = save_dir

        if not os.path.exists(f"{self.save_dir}/{self.comparison_name}"):
            print(f"Creating directory {self.save_dir}/{self.comparison_name}")
            os.makedirs(f"{self.save_dir}/{self.comparison_name}")


        self.plot_keff()
        self.plot_keff_with_sigma()
        for iso in tracked_nulides:
            self.plot_Ni(iso)


        return
    
    def plot_BU(self):
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        for case in self.cases_to_test:
            ax.plot(case.BU, case.BU, label = f"{case.case_name} - {case.edep_mode} - {case.lib_name}", marker = "x", linestyle = "--")
        ax.plot(self.case_ref.BU, self.case_ref.BU, label = f"{self.case_ref.case_name} - {self.case_ref.edep_mode} - {self.case_ref.lib_name}", marker = "x", linestyle = "--")
        ax.set_xlabel("Burnup [MWd/tU]")
        ax.set_ylabel("Burnup [MWd/tU]")
        ax.legend()
        ax.grid()
        plt.show()
        plt.savefig(f"{self.save_dir}/{self.comparison_name}/{self.comparison_name}_BU.png")
        plt.close()
        return

    
    def plot_keff(self):
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        for case in self.cases_to_test:
            ax.plot(case.BU, case.keffs, label = f"{case.case_name} - {case.edep_mode} - {case.lib_name}", marker = "x", linestyle = "--")
        ax.plot(self.case_ref.BU, self.case_ref.keffs, label = f"{self.case_ref.case_name} - {self.case_ref.edep_mode} - {self.case_ref.lib_name}", marker = "x", linestyle = "--")
        ax.set_xlabel("Burnup [MWd/tU]")
        ax.set_ylabel("Keff")
        ax.legend()
        ax.grid()
        plt.show()
        plt.savefig(f"{self.save_dir}/{self.comparison_name}/{self.comparison_name}_keff.png")
        plt.close()
        return
    
    def plot_keff_with_sigma(self):
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        for case in self.cases_to_test:
            ax.errorbar(case.BU, case.keffs, yerr = case.sigmas_keff, label = f"{case.case_name} - {case.edep_mode} - {case.lib_name}", marker = "x", linestyle = "--")
        ax.errorbar(self.case_ref.BU, self.case_ref.keffs, yerr = self.case_ref.sigmas_keff, label = f"{self.case_ref.case_name} - {self.case_ref.edep_mode} - {self.case_ref.lib_name}", marker = "x", linestyle = "--")
        ax.set_xlabel("Burnup [MWd/tU]")
        ax.set_ylabel("Keff")
        ax.legend()
        ax.grid()
        plt.show()
        plt.savefig(f"{self.save_dir}/{self.comparison_name}/{self.comparison_name}_keff_with_sigma.png")
        plt.close()
        return
    
    def plot_Ni(self, iso):
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        for case in self.cases_to_test:
            ax.plot(case.BU, case.Ni[iso], label = f"{case.case_name} - {iso}, {case.lib_name}", marker = "x", linestyle = "--")
        ax.plot(self.case_ref.BU, self.case_ref.Ni[iso], label = f"{self.case_ref.case_name} - {iso}, {case.lib_name}", marker = "x", linestyle = "--")
        ax.set_xlabel("Burnup [MWd/tU]")
        ax.set_ylabel("Atom density [#/b*cm]")
        ax.legend()
        ax.grid()
        plt.savefig(f"{self.save_dir}/{self.comparison_name}/{self.comparison_name}_Ni_{iso}.png")
        plt.show()
        plt.close()
        return
    
    def compute_and_plot_diff_Ni(self, isotope):
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        for case in self.cases_to_test:
            diff_Ni = [(case.Ni[isotope][idx] - self.case_ref.Ni[isotope][idx])*100/self.case_ref.Ni[isotope][idx] if self.case_ref.Ni[isotope][idx] != 0 else 0 for idx in range(len(case.BU))]
            ax.plot(case.BU, diff_Ni, label = f"{case.case_name} - {isotope} : {case.lib_name} vs {self.case_ref.lib_name}", marker = "x", linestyle = "--")
        ax.set_xlabel("Burnup [MWd/tU]")
        ax.set_ylabel("Atom density difference [%]")
        ax.legend()
        ax.grid()
        plt.savefig(f"{self.save_dir}/{self.comparison_name}/{self.comparison_name}_diff_Ni_{isotope}.png")
        plt.show()
        plt.close()
        return

    def compute_and_plot_keff_diff(self):
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        for case in self.cases_to_test:
            diff_keff = (case.keffs - self.case_ref.keffs)*1e5
            sig_err = np.sqrt(case.sigmas_keff**2 + self.case_ref.sigmas_keff**2)*1e5
            ax.errorbar(case.BU, diff_keff, yerr = sig_err, label = f"{case.case_name} - {case.lib_name} vs {self.case_ref.lib_name}", marker = "x", linestyle = "--")
        ax.set_xlabel("Burnup [MWd/tU]")
        ax.set_ylabel("Keff difference [pcm]")
        ax.legend()
        ax.grid()
        plt.savefig(f"{self.save_dir}/{self.comparison_name}/{self.comparison_name}_diff_keff.png")
        plt.show()
        plt.close()
        return
    
if __name__ == "__main__":
    save_dir_comparison = "ACE_files_comparison"
    if not os.path.exists(save_dir_comparison):
        os.makedirs(save_dir_comparison)
    # Define the Serpent2_case instances
    tracked_nulides = ["U235", "U238", "Pu239", "Gd157", "Gd158", "Xe135", "Sm149"]
    endfb8r1_2012_edep0 = Serpent2_case("HOM_Gd157_VBOC", "endfb8r1_pynjoy2012", 0, False, 2, 38.6, tracked_nulides, "HOM_Gd157_VBOC")
    
    endfb8r1_2012_kerma_edep0 = Serpent2_case("HOM_Gd157_VBOC", "endfb8r1_pynjoy2012_kerma", 0, False, 2, 38.6, tracked_nulides, "HOM_Gd157_VBOC") 
    endfb8r1_2012_kerma_edep2 = Serpent2_case("HOM_Gd157_VBOC", "endfb8r1_pynjoy2012_kerma", 2, False, 2, 38.6, tracked_nulides, "HOM_Gd157_VBOC")

    endfb8r1_2016_edep0 = Serpent2_case("HOM_Gd157_VBOC", "endfb8r1_NJOY2016", 0, False, 2, 38.6, tracked_nulides, "HOM_Gd157_VBOC")
    print(endfb8r1_2012_edep0.BU)
    print(endfb8r1_2012_kerma_edep0.BU)
    print(endfb8r1_2012_kerma_edep2.BU)
    print(endfb8r1_2016_edep0.BU)


    # Define the Serpent2_comparison instance

    comparison = Serpent2_comparison(comparison_name="HOM_Gd157_VBOC_pynjoy_kerma_vs_nokerma", 
                                    cases_to_test=[endfb8r1_2012_kerma_edep0], 
                                    case_ref=endfb8r1_2012_edep0,
                                    tracked_nulides=tracked_nulides, save_dir=save_dir_comparison)
    for iso in tracked_nulides:
        comparison.compute_and_plot_diff_Ni(iso)
    comparison.compute_and_plot_keff_diff()


    comparison = Serpent2_comparison(comparison_name="HOM_Gd157_VBOC_NJOY2016_vs_pynjoy2012",
                                    cases_to_test=[endfb8r1_2016_edep0],
                                    case_ref=endfb8r1_2012_edep0,
                                    tracked_nulides=tracked_nulides, save_dir=save_dir_comparison)
    for iso in tracked_nulides:
        comparison.compute_and_plot_diff_Ni(iso)
    comparison.compute_and_plot_keff_diff()

    comparison = Serpent2_comparison(comparison_name="HOM_Gd157_VBOC_pynjoy2012_kerma_edep0_vs_edep2",
                                    cases_to_test=[endfb8r1_2012_kerma_edep2],
                                    case_ref=endfb8r1_2012_kerma_edep0,
                                    tracked_nulides=tracked_nulides, save_dir=save_dir_comparison)
    
    for iso in tracked_nulides:
        comparison.compute_and_plot_diff_Ni(iso)
    comparison.compute_and_plot_keff_diff()

    values_BU_edep2 = endfb8r1_2012_kerma_edep2.BU
    values_NGd157_edep2 = endfb8r1_2012_kerma_edep2.Ni["Gd157"]
    
    values_BU_edep0 = endfb8r1_2012_kerma_edep0.BU
    values_NGd157_edep0 = endfb8r1_2012_kerma_edep0.Ni["Gd157"]

    # interpolate Gd157 atom density values for edep2 on the BU values of edep0
    values_NGd157_edep2_interp = np.interp(values_BU_edep0, values_BU_edep2, values_NGd157_edep2)
    # difference between the interpolated values and the edep0 values
    diff_NGd157 = (values_NGd157_edep2_interp-values_NGd157_edep0)*100/values_NGd157_edep0

    diff_no_interp = (values_NGd157_edep2 - values_NGd157_edep0)*100/values_NGd157_edep0

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.plot(values_BU_edep0, diff_NGd157, label = f"Diff Gd157 atom density [%]", marker = "x", linestyle = "--")
    ax.plot(values_BU_edep0, diff_no_interp, label = f"Diff Gd157 atom density [%] - no interpolation", marker = "x", linestyle = "--")
    ax.set_xlabel("Burnup [MWd/tU]")
    ax.set_ylabel("Atom density difference [%]")
    ax.legend()
    ax.grid()
    plt.savefig(f"{save_dir_comparison}/HOM_Gd157_VBOC_pynjoy2012_kerma_edep0_vs_edep2/interp_HOM_Gd157_VBOC_pynjoy_kerma_edep0_vs_edep2_diff_NGd157.png")
    plt.close()


    val_keff_edep2 = endfb8r1_2012_kerma_edep2.keffs
    val_keff_edep0 = endfb8r1_2012_kerma_edep0.keffs

    # interpolate keff values for edep2 on the BU values of edep0
    val_keff_edep2_interp = np.interp(values_BU_edep0, values_BU_edep2, val_keff_edep2)

    # difference between the interpolated values and the edep0 values
    diff_keff_interp = (val_keff_edep2_interp - val_keff_edep0)*1e5
    diff_keff = (val_keff_edep2 - val_keff_edep0)*1e5

    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    ax.plot(values_BU_edep0, diff_keff_interp, label = f"Diff keff [pcm]", marker = "x", linestyle = "--")
    ax.plot(values_BU_edep0, diff_keff, label = f"Diff keff [pcm] - no interpolation", marker = "x", linestyle = "--")
    ax.set_xlabel("Burnup [MWd/tU]")
    ax.set_ylabel("Keff difference [pcm]")
    ax.legend()
    ax.grid()
    plt.savefig(f"{save_dir_comparison}/HOM_Gd157_VBOC_pynjoy2012_kerma_edep0_vs_edep2/interp_HOM_Gd157_VBOC_pynjoy_kerma_edep0_vs_edep2_diff_keff.png")
    plt.close()
    

    print(f"number of BU points = {len(endfb8r1_2012_edep0.BU)}")
    


