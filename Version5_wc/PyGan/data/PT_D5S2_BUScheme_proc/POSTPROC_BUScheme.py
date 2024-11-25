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
    def __init__(self, lib_name, pcc_id, bu_points, ssh_module, ssh_method, correlation, tracked_nuclides, BU_lists, save_dir):
        path = os.getcwd()
        self.DIR = "EDIBU"
        self.lib_name = lib_name 
        self.pcc_id = pcc_id
        self.bu_points = bu_points
        self.ssh_module = ssh_module
        self.ssh_method = ssh_method
        self.correlation = correlation
        self.tracked_nuclides = tracked_nuclides
        self.save_dir = save_dir
        
        self.BU_list = BU_lists["BU"]
        self.AUTOP_list = BU_lists["AUTOP"] # only relevant for PCC0 and PCC2 schemes as ssh is performed at all preditor-corrector steps in PCC1, PCC3 and PCC3b 
        self.COMPO_list = BU_lists["COMPO"]
    
        os.chdir("PyGan_results_HOM_UOX_Gd157_PCC")
        name_compo = f"_COMPO_HOM_UOX_Gd157_PCC{pcc_id}_{lib_name}_{bu_points}_{ssh_module}_{ssh_method}_{correlation}"
        os.chdir(path)
        self.pyCOMPO = lcm.new('LCM_INP',name_compo,impx=0)

        self.DRAGON_BU = None
        self.DRAGON_Keff = None
        self.DRAGON_ISOTOPESDENS = {}

        self.parse_compo()
        

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
        plt.plot(self.DRAGON_BU, self.DRAGON_Keff, label = f"DRAGON {self.lib_name} pcc {self.pcc_id}", marker = "x", linestyle = "--")
        plt.xlabel("Burnup [MWd/kgU]")
        plt.ylabel("Keff")
        plt.title(f"Keff evolution for {self.lib_name} case")
        plt.legend()
        plt.grid()
        plt.savefig(f"DRAGON_Keff_{self.lib_name}_{self.ssh_module}_{self.ssh_method}_PCC{self.pcc_id}.png")
        plt.close()
        return
    
    def plot_Ni(self, isotope):
        plt.figure()
        plt.plot(self.DRAGON_BU, self.DRAGON_ISOTOPESDENS[isotope], label = f"D5 : {isotope} atomic density", marker = "x", linestyle = "--")
        plt.xlabel("Burnup [MWd/kgU]")
        plt.ylabel("Isotope density")
        plt.title(f"{isotope} evolution for {self.lib_name} case")
        plt.legend()
        plt.grid()
        plt.savefig(f"DRAGON_{isotope}_{self.lib_name}_{self.ssh_module}_{self.ssh_method}_PCC{self.pcc_id}.png")
        plt.close()
        return
    
class D5S2_comparisons:
    def __init__(self, comparison_name, D5_case, S2_cases, tracked_nuclides):
        self.comparison_name = comparison_name
        self.D5_case = D5_case
        self.S2_cases = S2_cases
        self.delta_keffs = {}
        self.delta_Niso = {}
        self.tracked_nuclides = tracked_nuclides

    def compare_keffs(self):
        for case in self.S2_cases:
            delta_keff = (self.D5_case.DRAGON_Keff - case.S2_Keff)*1e5 # error on Keff in pcm
            self.delta_keffs[f"{case.lib_name}_{case.edep_id}_{case.pcc_id}"] = delta_keff
        return
    def compare_Ni(self):
        for case in self.S2_cases:
            delta_Niso = {}
            for iso in self.tracked_nuclides:
                delta_Niso = [(self.D5_case.ISOTOPESDENS[iso][idx] - case[iso][idx]) * 100 / case[iso][idx]
                    if self.case[iso][idx] != 0 else 0
                    for idx in range(len(case[iso]))]
            self.delta_Niso[f"{iso}_{case.lib_name}_{case.edep_id}_{case.pcc_id}"] = delta_Niso
        return

    def plot_delta_Keff(self):
        """
        Plot the delta Keff for all cases : 1 D5 case compared to several S2 cases

        """
        plt.figure()
        for comparison_case in self.delta_keffs.keys():
            plt.plot(self.D5_case.DRAGON_BU, self.delta_keffs[comparison_case], label = f"{comparison_case}".replace("_"," "), marker = "x", linestyle = "--")
        plt.xlabel("Burnup [MWd/tU]")
        plt.ylabel("Delta Keff [pcm]")
        plt.title(f"Delta Keff evolution for {self.comparison_name} case")
        plt.legend()
        plt.grid()
        plt.savefig(f"Delta_Keff_{self.comparison_name}.png")
        plt.close()
        return