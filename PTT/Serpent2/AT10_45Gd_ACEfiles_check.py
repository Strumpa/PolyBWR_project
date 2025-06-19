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

# S2 PT class, updated 16/04/2025 to include option to post treat cases with manually  set fission heating values

class Serpent2_case: 
    def __init__(self, case_name, lib_name, edep_id, areQfissSet, isEcaptSet, pcc_id, specific_power, tracked_nuclides, save_dir):
        self.case_name = case_name
        self.lib_name = lib_name        
        self.isEcaptSet = isEcaptSet # Option to be activated when energy deposition mode is set to 1 in Serpent2 and Ecapt has been manually set in Serpent2.
        if self.isEcaptSet and edep_id == 1:
            self.edep_id = f"{edep_id}_Ecapt"
        else:
            self.edep_id = edep_id
        if self.edep_id == 2 or self.edep_id == 1 or self.edep_id == "1_Ecapt":
            self.areQfissSet = False
        else:
            self.areQfissSet = areQfissSet
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
        if self.areQfissSet:
            self.res = st.read(f"{os.environ['SERPENT_RESULTS']}/{self.case_name}/{self.case_name}_{self.lib_name}_edep{self.edep_id}_qfiss_pcc{self.pcc_id}_mc_res.m")
            self.depl = st.read(f"{os.environ['SERPENT_RESULTS']}/{self.case_name}/{self.case_name}_{self.lib_name}_edep{self.edep_id}_qfiss_pcc{self.pcc_id}_mc_dep.m")
        else:
            self.res = st.read(f"{os.environ['SERPENT_RESULTS']}/{self.case_name}/{self.case_name}_{self.lib_name}_edep{self.edep_id}_pcc{self.pcc_id}_mc_res.m")
            self.depl = st.read(f"{os.environ['SERPENT_RESULTS']}/{self.case_name}/{self.case_name}_{self.lib_name}_edep{self.edep_id}_pcc{self.pcc_id}_mc_dep.m")
        self.keff = self.res.resdata["absKeff"].T[0] # Retrieve the keffs
        self.sigmas_keff = self.res.resdata["absKeff"].T[1] # Retrieve the uncertainties (standard deviation) on the keffs
        self.mat = self.depl.materials['total'] # Retrieve the total, homogenized material
        self.BU = self.mat.burnup # Burnup in MWd/kgU
        
        data_frame = self.mat.toDataFrame("adens", names = self.tracked_nuclides, time = "days") # Convert the material to a pandas dataframe
        self.BUdays = data_frame.index.to_numpy() # retrieve the burnup points in days
        
        self.Ni = {} # Dictionary to store the isotopic densities
        for isotope in data_frame.columns:
            self.Ni[isotope] = data_frame[isotope].to_numpy() # Convert the pandas series to numpy array
        
        return 
    
    def rescale_BU(self, units):
        if units == "MWd/tU" and self.unitsBU == "MWd/kgU":
            self.BU = self.BU*1000
            self.unitsBU = "MWd/tU"

    def plot_keff(self):
        plt.figure()
        plt.plot(self.BU, self.keff, label = f"edepmode {self.edep_id} pcc {self.pcc_id}", marker = "x", linestyle = "--")
        plt.xlabel(f"Burnup [{self.unitsBU}]")
        plt.ylabel("Keff")
        plt.title(f"Keff evolution for {self.case_name} case")
        plt.legend()
        plt.grid()
        plt.savefig(f"{self.save_dir}/S2_Keff_{self.case_name}_edep{self.edep_id}_pcc{self.pcc_id}.png")
        plt.close()
        return
    
    def plot_concentrations(self, isotopes_list):
        plt.figure()
        for nuclide in isotopes_list:
            plt.plot(self.BU, self.Ni[nuclide], label = nuclide)
        plt.xlabel(f"Burnup [{self.unitsBU}]")
        plt.ylabel("Concentration [atom/b-cm]")
        plt.title(f"Concentration evolution for {self.case_name} case")
        plt.legend()
        plt.grid()
        plt.savefig(f"{self.save_dir}/S2_Concentrations_{self.case_name}_edep{self.edep_id}_pcc{self.pcc_id}.png")
        plt.close()
        return
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
            ax.plot(case.BU, case.keff, label = f"{case.case_name} - {case.edep_mode} - {case.lib_name}", marker = "x", linestyle = "--")
        ax.plot(self.case_ref.BU, self.case_ref.keff, label = f"{self.case_ref.case_name} - {self.case_ref.edep_mode} - {self.case_ref.lib_name}", marker = "x", linestyle = "--")
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
            ax.errorbar(case.BU, case.keff, yerr = case.sigmas_keff, label = f"{case.case_name} - {case.edep_mode} - {case.lib_name}", marker = "x", linestyle = "--")
        ax.errorbar(self.case_ref.BU, self.case_ref.keff, yerr = self.case_ref.sigmas_keff, label = f"{self.case_ref.case_name} - {self.case_ref.edep_mode} - {self.case_ref.lib_name}", marker = "x", linestyle = "--")
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
            diff_keff = (case.keff - self.case_ref.keff)*1e5
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
    

# --- Main script
# Define the input parameters
if __name__ == "__main__":
    """
    Post treatment script for the outputs of the Serpent2 code, AT10_45Gd cell

    goal is to assess consistency between different ACE files, generated by NJOY2016, pynjoy2012 and pynjoy2012_kerma
    evaluation used in implied : "endfb8r1"
    """

    save_dir_comparison = "ACE_files_comparison"
    if not os.path.exists(save_dir_comparison):
        os.makedirs(save_dir_comparison)

    tracked_nulides = ["U235", "U238", "Pu239", "Gd157", "Xe135", "Sm149"]
    # Define the input parameters
    # --- Serpent2 case
    # --- ATRIUM-10 : 45Gd cell

    # constructor of S2_case :
    # (self, case_name, lib_name, edep_id, areQfissSet, isEcaptSet, pcc_id, specific_power, tracked_nuclides, save_dir)

    S2_edep0_setqfiss_pcc0_pynjoy2012 = Serpent2_case("AT10_45Gd", "endfb8r1_pynjoy2012", edep_id=0, areQfissSet=True, isEcaptSet=False, pcc_id=0, 
                                                      specific_power=26.5, tracked_nuclides=tracked_nulides, save_dir="AT10_45Gd_Serpent2_results")
    
    S2_edep0_setqfiss_pcc0_pynjoy2012_kerma = Serpent2_case("AT10_45Gd", "endfb8r1_pynjoy2012_kerma", edep_id=0, areQfissSet=True, isEcaptSet=False, pcc_id=0, 
                                                      specific_power=26.5, tracked_nuclides=tracked_nulides, save_dir="AT10_45Gd_Serpent2_results")
    
    S2_edep0_setqfiss_pcc0_NJOY2016 = Serpent2_case("AT10_45Gd", "endfb8r1_NJOY2016", edep_id=0, areQfissSet=True, isEcaptSet=False, pcc_id=0, 
                                                      specific_power=26.5, tracked_nuclides=tracked_nulides, save_dir="AT10_45Gd_Serpent2_results")
    
    S2_edep0_setqfiss_pcc1_pynjoy2012 = Serpent2_case("AT10_45Gd", "endfb8r1_pynjoy2012", edep_id=0, areQfissSet=True, isEcaptSet=False, pcc_id=1, 
                                                      specific_power=26.5, tracked_nuclides=tracked_nulides, save_dir="AT10_45Gd_Serpent2_results")
    
    S2_edep0_setqfiss_pcc1_pynjoy2012_kerma = Serpent2_case("AT10_45Gd", "endfb8r1_pynjoy2012_kerma", edep_id=0, areQfissSet=True, isEcaptSet=False, pcc_id=1, 
                                                      specific_power=26.5, tracked_nuclides=tracked_nulides, save_dir="AT10_45Gd_Serpent2_results")
    
    S2_edep0_setqfiss_pcc1_NJOY2016 = Serpent2_case("AT10_45Gd", "endfb8r1_NJOY2016", edep_id=0, areQfissSet=True, isEcaptSet=False, pcc_id=1, 
                                                      specific_power=26.5, tracked_nuclides=tracked_nulides, save_dir="AT10_45Gd_Serpent2_results")
    
    # pynjoy2012_kerma vs pynjoy2012 acefiles, no kerma data used
    comparison = Serpent2_comparison(comparison_name="AT10_45Gd_pynjoy_kerma_vs_nokerma", 
                                    cases_to_test=[S2_edep0_setqfiss_pcc0_pynjoy2012_kerma], 
                                    case_ref=S2_edep0_setqfiss_pcc0_pynjoy2012,
                                    tracked_nulides=tracked_nulides, save_dir=save_dir_comparison)
    
    for iso in tracked_nulides:
        comparison.compute_and_plot_diff_Ni(iso)
    comparison.compute_and_plot_keff_diff()

    # pynjoy2012 vs NJOY2016 acefiles, no kerma data used
    comparison = Serpent2_comparison(comparison_name="AT10_45Gd_pynjoy_vs_NJOY2016",
                                    cases_to_test=[S2_edep0_setqfiss_pcc0_pynjoy2012], 
                                    case_ref=S2_edep0_setqfiss_pcc0_NJOY2016,
                                    tracked_nulides=tracked_nulides, save_dir=save_dir_comparison)
    
    for iso in tracked_nulides:
        comparison.compute_and_plot_diff_Ni(iso)
    comparison.compute_and_plot_keff_diff()

    # pynjoy2012_kerma vs NJOY2016 acefiles, no kerma data used
    comparison = Serpent2_comparison(comparison_name="AT10_45Gd_pynjoy_kerma_vs_NJOY2016",
                                    cases_to_test=[S2_edep0_setqfiss_pcc0_pynjoy2012_kerma], 
                                    case_ref=S2_edep0_setqfiss_pcc0_NJOY2016,
                                    tracked_nulides=tracked_nulides, save_dir=save_dir_comparison)
    
    for iso in tracked_nulides:
        comparison.compute_and_plot_diff_Ni(iso)
    comparison.compute_and_plot_keff_diff()

    ## Do the same for pcc1 case :
    # pynjoy2012_kerma vs pynjoy2012 acefiles, no kerma data used
    comparison = Serpent2_comparison(comparison_name="AT10_45Gd_pynjoy_kerma_vs_nokerma_pcc1", 
                                    cases_to_test=[S2_edep0_setqfiss_pcc1_pynjoy2012_kerma], 
                                    case_ref=S2_edep0_setqfiss_pcc1_pynjoy2012,
                                    tracked_nulides=tracked_nulides, save_dir=save_dir_comparison)
    for iso in tracked_nulides:
        comparison.compute_and_plot_diff_Ni(iso)
    comparison.compute_and_plot_keff_diff()
    # pynjoy2012 vs NJOY2016 acefiles, no kerma data used
    comparison = Serpent2_comparison(comparison_name="AT10_45Gd_pynjoy_vs_NJOY2016_pcc1",
                                    cases_to_test=[S2_edep0_setqfiss_pcc1_pynjoy2012], 
                                    case_ref=S2_edep0_setqfiss_pcc1_NJOY2016,
                                    tracked_nulides=tracked_nulides, save_dir=save_dir_comparison)
    for iso in tracked_nulides:
        comparison.compute_and_plot_diff_Ni(iso)
    comparison.compute_and_plot_keff_diff()
    # pynjoy2012_kerma vs NJOY2016 acefiles, no kerma data used
    comparison = Serpent2_comparison(comparison_name="AT10_45Gd_pynjoy_kerma_vs_NJOY2016_pcc1",
                                    cases_to_test=[S2_edep0_setqfiss_pcc1_pynjoy2012_kerma], 
                                    case_ref=S2_edep0_setqfiss_pcc1_NJOY2016,
                                    tracked_nulides=tracked_nulides, save_dir=save_dir_comparison)
    
    for iso in tracked_nulides:
        comparison.compute_and_plot_diff_Ni(iso)
    comparison.compute_and_plot_keff_diff()


    ## Make a simple comparison for different aec files, same pcc and edep mode
    # pcc1 : pynjoy2012, pynjoy2012_kerma vs NJOY2016

    comparision = Serpent2_comparison(comparison_name="AT10_45Gd_pynjoy_kerma_pynjoy_vs_NJOY2016_pcc1",
                                    cases_to_test=[S2_edep0_setqfiss_pcc1_pynjoy2012_kerma, S2_edep0_setqfiss_pcc1_pynjoy2012], 
                                    case_ref=S2_edep0_setqfiss_pcc1_NJOY2016,
                                    tracked_nulides=tracked_nulides, save_dir=save_dir_comparison)
    
    for iso in tracked_nulides:
        comparision.compute_and_plot_diff_Ni(iso)
    comparision.compute_and_plot_keff_diff()
    
