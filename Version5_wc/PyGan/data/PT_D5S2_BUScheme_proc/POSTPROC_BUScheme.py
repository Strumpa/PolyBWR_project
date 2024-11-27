# Python3 class for post treating DRAGON5-Serpent2 comparison of Burnup Schemes
# Author : R. Guasch
# Date : 2024/11/22
# Purpose : Post treat comaprison of DRAGON5-Serpent2 BU schemes

import sys
import numpy as np
import os, shutil
import lifo
import lcm
import cle2000
import matplotlib.pyplot as plt

class DRAGON_case:
    def __init__(self, dlib_name, pcc_id, bu_points, ssh_module, ssh_method, correlation, tracked_nuclides, BU_lists, save_dir):
        path = os.getcwd()
        self.DIR = "EDIBU"
        self.draglib_name = dlib_name 
        self.pcc_id = pcc_id
        self.bu_points = bu_points
        self.ssh_module = ssh_module
        self.ssh_method = ssh_method
        self.correlation = correlation
        self.specific_power = 38.6 # W/g, test with lower values to see if system is more numerically stable : maybe less important varaitions in N_Gd157 ?
        self.tracked_nuclides = tracked_nuclides
        self.save_dir = save_dir

        self.BU_list = BU_lists["BU"]
        self.AUTOP_list = BU_lists["AUTOP"] # only relevant for PCC0 and PCC2 schemes as ssh is performed at all preditor-corrector steps in PCC1, PCC3 and PCC3b 
        self.COMPO_list = BU_lists["COMPO"]
    
        os.chdir("PyGan_results_HOM_UOX_Gd157_PCC")
        # list all files in the directory
        files = os.listdir()
        print(f"Files in the directory : {files}")
        name_compo = f"COMPO_HOM_UOX_Gd157_PCC{pcc_id}_{dlib_name}_{bu_points}_{ssh_module}_{ssh_method}_{correlation}"
        print(f"Name of the COMPO file : {name_compo}") 
        self.pyCOMPO = lcm.new('LCM_INP',name_compo,impx=0)
        os.chdir(path)

        self.DRAGON_BU = None
        self.DRAGON_Keff = None
        self.DRAGON_ISOTOPESDENS = {}

        self.parse_compo()
        if self.pcc_id == 0:
            self.BUScheme = "CE"
        elif self.pcc_id == 1:
            self.BUScheme = "CELI"
        elif self.pcc_id == 2:
            self.BUScheme = "LE"
        elif self.pcc_id == 3:
            self.BUScheme = "LELI"

        return
        

    def parse_compo(self):
        lenBU_DRAGON=np.shape(self.COMPO_list)[0]
        ISOTOPES=self.pyCOMPO[self.DIR]['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESDENS']
        print(f"Dragon isotopes = {ISOTOPES}")
        #lenISOT_DRAGON=np.shape(ISOTOPES)[0]-1
        lenISOT_DRAGON=np.shape(ISOTOPES)[0]
        DRAGON_BU=self.COMPO_list
        DRAGON_ISOTOPESDENS=np.zeros((lenISOT_DRAGON,lenBU_DRAGON))
        DRAGON_Keff=np.zeros(lenBU_DRAGON)

        #print("$$$ ---------------- DRAGON_ISOTOPESDENS shape = ",lenISOT_DRAGON,lenBU_DRAGON)

        for k in range(lenBU_DRAGON):
            DRAGON_Keff[k]=self.pyCOMPO[self.DIR]['MIXTURES'][0]['CALCULATIONS'][k]['K-EFFECTIVE']
            #print("$$$ ---------------- ISOTOPES BU step ",k,"/",lenBU_DRAGON," = ",COMPO_py['EDIBU_HOM']['MIXTURES'][0]['CALCULATIONS'][k]['ISOTOPESDENS'])    
            for j in range(lenISOT_DRAGON):
                #print("$$$ ---------------- ISOTOPES ",j,"/",lenISOT_DRAGON," = ",COMPO_py['EDIBU']['MIXTURES'][0]['CALCULATIONS'][k]['ISOTOPESDENS'][j])
                DRAGON_ISOTOPESDENS[j][k]=self.pyCOMPO[self.DIR]['MIXTURES'][0]['CALCULATIONS'][k]['ISOTOPESDENS'][j]


        # --------- List of isotopes from the DRAGON Multicompo results
        isotopes2=[]
        isotopes=[]
        for k in range(lenISOT_DRAGON):
            isotopes2=isotopes2+[self.pyCOMPO[self.DIR]['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][k]['ALIAS']]
        for k in range(lenISOT_DRAGON):
            if isotopes2[k][0]=='U':
                isotopes=isotopes+[isotopes2[k][0:4]]
            else:
                isotopes=isotopes+[isotopes2[k][0:5]]

        indices=np.zeros(len(self.tracked_nuclides))
        for n in range(len(self.tracked_nuclides)):
            for m in range(len(isotopes)):
                if self.tracked_nuclides[n]==isotopes[m]:
                    indices[n]=m

        for k in range(len(self.tracked_nuclides)):
            self.DRAGON_ISOTOPESDENS[self.tracked_nuclides[k]] = DRAGON_ISOTOPESDENS[int(indices[k])]
        self.DRAGON_BU = DRAGON_BU
        self.DRAGON_Keff = DRAGON_Keff

        return

    def plot_keffs(self):
        plt.figure()
        plt.plot(self.DRAGON_BU, self.DRAGON_Keff, label = f"DRAGON {self.draglib_name} pcc {self.pcc_id}", marker = "x", linestyle = "--")
        plt.xlabel("Burnup [MWd/kgU]")
        plt.ylabel("Keff")
        plt.title(f"Keff evolution for {self.draglib_name} case")
        plt.legend()
        plt.grid()
        plt.savefig(f"DRAGON_Keff_{self.draglib_name}_{self.ssh_module}_{self.ssh_method}_PCC{self.pcc_id}.png")
        plt.close()
        return
    
    def plot_Ni(self, isotope):
        plt.figure()
        plt.plot(self.DRAGON_BU, self.DRAGON_ISOTOPESDENS[isotope], label = f"D5 : {isotope} atomic density", marker = "x", linestyle = "--")
        plt.xlabel("Burnup [MWd/kgU]")
        plt.ylabel("Isotope density")
        plt.title(f"{isotope} evolution for {self.draglib_name} case")
        plt.legend()
        plt.grid()
        plt.savefig(f"DRAGON_{isotope}_{self.draglib_name}_{self.ssh_module}_{self.ssh_method}_PCC{self.pcc_id}.png")
        plt.close()
        return
    
class D5S2_comparisons:
    def __init__(self, comparison_name, D5_case, S2_cases, tracked_nuclides, save_dir):
        self.comparison_name = comparison_name
        self.D5_case = D5_case
        self.S2_cases = S2_cases
        self.delta_keffs = {}
        self.delta_Niso = {}
        self.tracked_nuclides = tracked_nuclides
        self.save_dir = save_dir
        # check if the save directory exists and create it if not
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        
        return

    def compare_keffs(self):
        for case in self.S2_cases:
            delta_keff = (self.D5_case.DRAGON_Keff - case.keffs)*1e5 # error on Keff in pcm
            self.delta_keffs[f"{case.lib_name}_edep_{case.edep_id}_pcc_{case.pcc_id}"] = delta_keff
        return
    def compare_Ni(self):
        for iso in self.tracked_nuclides:
            delta_Niso_case = {}
            for case in self.S2_cases:
                delta_Niso = [(self.D5_case.DRAGON_ISOTOPESDENS[iso][idx] - case.Ni[iso][idx]) * 100 / case.Ni[iso][idx]
                    if case.Ni[iso][idx] != 0 else 0
                    for idx in range(len(case.Ni[iso]))]
                delta_Niso_case[f"{case.lib_name}_edep{case.edep_id}_PCC{case.pcc_id}"] = delta_Niso
            self.delta_Niso[f"{iso}"] = delta_Niso_case
        return

    def plot_delta_Keff(self):
        """
        Plot the delta Keff for all cases : 1 D5 case compared to several S2 cases
        """
        plt.figure()
        for comparison_case in self.delta_keffs.keys():
            plt.plot(self.D5_case.DRAGON_BU, self.delta_keffs[comparison_case], label = f"{comparison_case}".replace("_"," "), marker = "x", linestyle = "--")
        plt.xlabel("Burnup [MWd/tU]")
        plt.ylabel("$\\Delta$ Keff [pcm]")
        plt.axhline(y = 300.0, color = 'r', linestyle = '-')
        plt.axhline(y = -300.0, color = 'r', linestyle = '-')  
        plt.title(f"$\\Delta$ Keff evolution for {self.comparison_name} case")
        plt.legend()
        plt.grid()
        plt.savefig(f"{self.save_dir}/Delta_Keff_{self.comparison_name}.png")
        plt.close()
        return
    
    def plot_delta_Ni(self):
        """
        Plot the delta Ni for all cases : 1 D5 case compared to several S2 cases
        """
        for iso in self.tracked_nuclides:
            delta_case = self.delta_Niso[iso]
            plt.figure()
            for comparison_case in delta_case.keys():
                plt.plot(self.D5_case.DRAGON_BU, self.delta_Niso[iso][comparison_case], label = f"{comparison_case}".replace("_"," "), marker = "x", linestyle = "--")
            plt.xlabel("Burnup [MWd/tU]")
            plt.ylabel(f"$\\Delta$ N{iso} [%]")
            plt.axhline(y = 2.0, color = 'r', linestyle = '-')
            plt.axhline(y = -2.0, color = 'r', linestyle = '-')  
            plt.title(f"$\\Delta$ N{iso} evolution for {self.comparison_name} case")
            plt.legend()
            plt.grid()
            plt.savefig(f"{self.save_dir}/Delta_{iso}_{self.comparison_name}.png")
            plt.close()
        return