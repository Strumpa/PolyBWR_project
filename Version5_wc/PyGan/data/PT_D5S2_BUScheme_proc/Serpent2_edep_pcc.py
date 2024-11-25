# Python3 class used to post treat Serpent2 BU evolution files
# Author: R. Guasch
# Date: 2024/11/19
# Purpose : Post treat Serpent2 BU evolution files for HOM_UOX_Gd157 case, 
# Two main parameters are considered : - energy deposition mode : set edep keyword in Serpent2
#                                      - bu evolution scheme predictor corrector mode : set pcc keyword in Serpent2
import serpentTools as st
import numpy as np
import matplotlib.pyplot as plt
import os

class Serpent2_case:
    def __init__(self, case_name, lib_name, edep_id, pcc_id, tracked_nuclides):
        self.case_name = case_name
        self.lib_name = lib_name
        self.edep_id = edep_id
        self.pcc_id = pcc_id
        self.tracked_nuclides = tracked_nuclides
        self.unitsBU = "MWd/kgU" # Default Serpent2 units for BU
        self.read_S2_outputs()
        path=os.getcwd()
        a=os.path.exists('BWRresults_PyGan_Serpent2_edep_pcc')
        if a==False:
            os.mkdir('BWRresults_PyGan_Serpent2_edep_pcc')
        self.save_path = "BWRresults_PyGan_Serpent2_edep_pcc"
    
    def read_S2_outputs(self):
        # Read the results
        self.res = st.read(f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/HOM_CELL_study/{self.case_name}/ENERGY_DEPOSITION_study/{self.case_name}_{self.lib_name}_edep{self.edep_id}_pcc{self.pcc_id}_mc_res.m")
        self.keffs = self.res.resdata["absKeff"].T[0]
        self.sigmas_keff = self.res.resdata["absKeff"].T[1]
        self.depl = st.read(f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/HOM_CELL_study/{self.case_name}/ENERGY_DEPOSITION_study/{self.case_name}_{self.lib_name}_edep{self.edep_id}_pcc{self.pcc_id}_mc_dep.m")
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

    def plot_keff(self):
        plt.figure()
        #print(self.keffs)
        plt.plot(self.BU, self.keffs, label = f"edepmode {self.edep_id} pcc {self.pcc_id}", marker = "x", linestyle = "--")
        plt.xlabel(f"Burnup [{self.unitsBU}]")
        plt.ylabel("Keff")
        plt.title(f"Keff evolution for {self.case_name} case")
        plt.legend()
        plt.grid()
        plt.savefig(f"{self.save_path}/S2_Keff_{self.case_name}_edep{self.edep_id}_pcc{self.pcc_id}.png")
        plt.close()
        return
    
    def plot_concentrations(self, isotopes_list):
        plt.figure()
        for nuclide in isotopes_list:
            plt.plot(self.BU, self.Ni[nuclide], label = nuclide)
        plt.xlabel("Burnup [MWd/kgU]")
        plt.ylabel("Concentration [atom/b-cm]")
        plt.title(f"Concentration evolution for {self.case_name} case")
        plt.legend()
        plt.grid()
        plt.savefig(f"{self.save_path}/S2_Concentrations_{self.case_name}_edep{self.edep_id}_pcc{self.pcc_id}.png")
        plt.close()
        return
 
    
class SERPENT2_comparisons:
    def __init__(self, comparison_name, S2_cases_list):
        self.comparison_name = comparison_name
        self.S2_cases_list = S2_cases_list
        self.accepted_libs = []
        self.accepted_edep = {}
        self.accepted_pcc = {}
        self.save_path = "BWRresults_PyGan_Serpent2_edep_pcc/S2Comparisons"
        path=os.getcwd()
        a=os.path.exists('BWRresults_PyGan_Serpent2_edep_pcc/S2Comparisons')
        if a==False:
            os.mkdir('BWRresults_PyGan_Serpent2_edep_pcc/S2Comparisons')
        
        return
    
    def set_accepted_params_for_lib(self, library, accepted_edep, accepted_pcc):
        self.accepted_libs.append(library)
        self.accepted_edep[library] = accepted_edep
        self.accepted_pcc[library] = accepted_pcc
        return
    
    def plot_keffs(self):
        plt.figure()
        for S2_case in self.S2_cases_list:
            plt.plot(S2_case.BU, S2_case.keffs, label = f"edepmode {S2_case.edep_id} pcc {S2_case.pcc_id}", marker = "x", linestyle = "--")
        plt.xlabel("Burnup [MWd/kgU]")
        plt.ylabel("Keff")
        plt.title(f"Keff evolution for {self.comparison_name} case")
        plt.legend()
        plt.grid()
        plt.savefig(f"{self.save_path}/S2_Keff_{self.comparison_name}.png")
        plt.close()
        return
    
    def set_ref_case(self, ref_lib, ref_edep, ref_pcc):
        self.ref_lib = ref_lib
        self.ref_edep = ref_edep
        self.ref_pcc = ref_pcc
        for S2_case in self.S2_cases_list:
            if S2_case.edep_id == ref_edep and S2_case.pcc_id == ref_pcc and S2_case.lib_name == ref_lib:
                print(f"Reference case set to {S2_case.case_name} for library {S2_case.lib_name} with edepmode {S2_case.edep_id} and pcc {S2_case.pcc_id}")
                self.ref_case = S2_case
                break
        return
    def get_case(self, lib, edep, pcc):
        for S2_case in self.S2_cases_list:
            if S2_case.edep_id == edep and S2_case.pcc_id == pcc and S2_case.lib_name == lib:
                print(f"Extracted case : {S2_case.case_name} with library {S2_case.lib_name}, edepmode {S2_case.edep_id} and pcc {S2_case.pcc_id}")
                return S2_case
            else:
                print(f"Case not found for library {lib}, edepmode {edep} and pcc {pcc}")
        return None
    
    def compute_delta_keffs(self):
        self.delta_keffs = {}
        for S2_case in self.S2_cases_list:
            if S2_case != self.ref_case:
                delta_keff = (S2_case.keffs - self.ref_case.keffs)*1e5
                self.delta_keffs[f"{S2_case.lib_name} edep{S2_case.edep_id}_pcc{S2_case.pcc_id}"] = delta_keff
        return
    
    def compute_delta_Niso(self, isotopes_list):
        self.delta_Niso = {}
        for iso in isotopes_list:
            self.delta_Niso[iso] = {}
            for S2_case in self.S2_cases_list:
                if S2_case != self.ref_case:
                    delta_Niso = [(S2_case.Ni[iso][idx] - self.ref_case.Ni[iso][idx]) * 100 / self.ref_case.Ni[iso][idx]
                    if self.ref_case.Ni[iso][idx] != 0 else 0
                    for idx in range(len(self.ref_case.Ni[iso]))]
                    self.delta_Niso[iso][f"{S2_case.lib_name} edep{S2_case.edep_id} pcc{S2_case.pcc_id}"] = delta_Niso
        return

    def make_pcc_comparison(self, edep_id, lib_name, isotopes_list):
        """
        For a given energy deposition mode, compare the pcc schemes for a given library
        edep_id : int : energy deposition mode,
        lib_name : str : library name,
        isotopes_list : list : list of isotopes to compare
        additionally this funtion compares the keffs as well as the isotopes concentrations
        """
        comparison_ref_case = self.get_case(lib_name, edep_id, pcc=1) # by default the reference case is the one with pcc=1

    
    def plot_delta_keffs(self):
        plt.figure()
        for case, delta_keff in self.delta_keffs.items():
            plt.plot(self.ref_case.BU, delta_keff, label = case, marker = "x", linestyle = "--")
        plt.xlabel("Burnup [MWd/kgU]")
        plt.ylabel("Delta Keff")
        plt.axhline(y = 300.0, color = 'r', linestyle = '-')
        plt.axhline(y = -300.0, color = 'r', linestyle = '-')  
        plt.title(f"Delta Keff evolution for {self.comparison_name} case")
        plt.legend()
        plt.grid()
        plt.savefig(f"{self.save_path}/S2_DeltaKeff_{self.comparison_name}.png")
        plt.close()
        return
    
    def plot_delta_Niso(self, isotopes_list):
        for iso in isotopes_list:
            plt.figure()
            for case, delta_Niso in self.delta_Niso[iso].items():
                print(delta_Niso)
                plt.plot(self.ref_case.BU, delta_Niso, label = case, marker = "x", linestyle = "--")
            plt.xlabel("Burnup [MWd/kgU]")
            plt.ylabel(f"Delta {iso} [%]")
            plt.axhline(y = 2.0, color = 'r', linestyle = '-')
            plt.axhline(y = -2.0, color = 'r', linestyle = '-')  
            plt.title(f"Delta {iso} evolution for {self.comparison_name} case")
            plt.legend()
            plt.grid()
            plt.savefig(f"{self.save_path}/S2_Delta{iso}_{self.comparison_name}.png")
            plt.close()
        return