# Python3 post treatment script for constant flux evolution case
# Author : R. Guasch
# Date : 2025/02/04

import numpy as np
import matplotlib.pyplot as plt
import serpentTools as st
import os

# S2 PT class :

class Serpent2_case: # Adapted frorm edep_pcc class.
    def __init__(self, case_name, lib_name, pcc_id, tracked_nuclides, save_dir):
        self.case_name = case_name
        self.lib_name = lib_name # PyNjoy2016 Serpent2 lib library
        self.pcc_id = pcc_id # Default is 1 for predictor corrector CE-LI
        self.tracked_nuclides = tracked_nuclides
        self.save_dir = save_dir
        self.unitsBU = "days" # Default Serpent2 units for BU
        self.read_S2_outputs()
        
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
    
    def read_S2_outputs(self):
        # Read the results
        self.res = st.read(f"{os.environ['SERPENT_RESULTS']}/HOM_CELL_study/{self.case_name}/BUScheme_EDEP_PCC_study/{self.case_name}_{self.lib_name}_Cst_flx_time_steps_mc_res.m")
        self.keffs = self.res.resdata["absKeff"].T[0]
        self.sigmas_keff = self.res.resdata["absKeff"].T[1]
        self.depl = st.read(f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/HOM_CELL_study/{self.case_name}/BUScheme_EDEP_PCC_study/{self.case_name}_{self.lib_name}_Cst_flx_time_steps_mc_dep.m")
        self.mat = self.depl.materials['total']
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
        plt.plot(self.BUdays, self.keffs, label = f"Constant flux pcc {self.pcc_id}", marker = "x", linestyle = "--")
        plt.xlabel(f"Burnup [{self.unitsBU}]")
        plt.ylabel("Keff")
        plt.title(f"Keff evolution for {self.case_name} case")
        plt.legend()
        plt.grid()
        plt.savefig(f"{self.save_dir}/S2_Keff_{self.case_name}_cst_flx_pcc{self.pcc_id}.png")
        plt.close()
        return
    
    def plot_concentrations(self, isotopes_list):
        plt.figure()
        for nuclide in isotopes_list:
            plt.plot(self.BUdays, self.Ni[nuclide], label = nuclide)
        plt.xlabel(f"Burnup [{self.unitsBU}]")
        plt.ylabel("Concentration [atom/b-cm]")
        plt.title(f"Concentration evolution for {self.case_name} case")
        plt.legend()
        plt.grid()
        plt.savefig(f"{self.save_dir}/S2_Concentrations_{self.case_name}_cst_flx_pcc{self.pcc_id}.png")
        plt.close()
        return


# Dragon case post treatment class :
    
class DRAGON_case:
    def __init__(self, pyCOMPO, dlib_name, bu_points, ssh_module, ssh_method, correlation, sat, depl_sol, tracked_nuclides, BU_lists, save_dir):
        """
        DRAGON5 output related variables
        DIR = : (str) name of the directory storing DRAGON5 outputs in COMPO object
        pyCOMPO = : (lcm object) COMPO object containing DRAGON5 outputs
        """
        self.DIR = "EDIBU"
        self.pyCOMPO = pyCOMPO

        """
        Calculation options
        """
        self.draglib_name = dlib_name # (str) name of the DRAGON library
        self.bu_points = bu_points # (str) identifier for the burnup points
        self.ssh_module = ssh_module
        self.ssh_method = ssh_method
        self.correlation = correlation
        self.sat = sat
        self.depl_sol = depl_sol

        """
        time step to BU step normalization factor
        """
        self.specific_power = 38.6 # W/gU
        
        """
        Output options
        """
        self.tracked_nuclides = tracked_nuclides
        self.save_dir = save_dir

        self.BU_list = BU_lists["BU"]
        self.AUTOP_list = BU_lists["AUTOP"] # only relevant for PCC0 and PCC2 schemes as ssh is performed at all preditor-corrector steps in PCC1, PCC3 and PCC3b 
        self.COMPO_list = BU_lists["COMPO"]

        # convert BU points to days using specific power.
        self.BU_days = np.array(self.BU_list)/self.specific_power
        self.AUTOP_days = np.array(self.AUTOP_list)/self.specific_power
        self.COMPO_days = np.array(self.COMPO_list)/self.specific_power


        self.DRAGON_BU = None
        self.DRAGON_Keff = None
        self.DRAGON_ISOTOPESDENS = {}

        self.parse_compo()

        self.BUScheme = "LE" # Default DRAGON burnup scheme is LE (Linearl Extrapolation of reaction rates over time steps)

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
        plt.plot(self.DRAGON_BU, self.DRAGON_Keff, label = f"D5 {self.draglib_name} : sat {self.sat}, depl sol {self.depl_sol}", marker = "x", linestyle = "--")
        plt.xlabel("Burnup [days]")
        plt.ylabel("Keff")
        plt.title(f"Keff evolution for {self.draglib_name} case")
        plt.legend()
        plt.grid()
        plt.savefig(f"{self.save_dir}/DRAGON_Keff_{self.draglib_name}_{self.ssh_module}_{self.ssh_method}_sat_{self.sat}_sol{self.depl_sol}.png")
        plt.close()
        return
    
    def plot_Ni(self, isotope):
        plt.figure()
        plt.plot(self.DRAGON_BU, self.DRAGON_ISOTOPESDENS[isotope], label = f"D5 : {isotope} atomic density", marker = "x", linestyle = "--")
        plt.xlabel("Burnup [days]")
        plt.ylabel("Isotope density")
        plt.title(f"{isotope} evolution for {self.draglib_name} case")
        plt.legend()
        plt.grid()
        plt.savefig(f"{self.save_dir}/DRAGON_{isotope}_{self.draglib_name}_{self.ssh_module}_{self.ssh_method}_sat_{self.sat}_sol{self.depl_sol}.png")
        plt.close()
        return
    

# Comparison between Serpent2 and DRAGON cases :
    
class D5S2_comparisons:
    def __init__(self, comparison_name, D5_cases, S2_case, tracked_nuclides, save_dir):
        self.comparison_name = comparison_name
        self.D5_cases = D5_cases
        self.S2_case = S2_case
        self.delta_keffs = {}
        self.delta_Niso = {}
        self.tracked_nuclides = tracked_nuclides
        self.save_dir = save_dir
        # check if the save directory exists and create it if not
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        
        self.compare_keffs()
        self.compare_Ni()
        
        return

    def compare_keffs(self):
        for case in self.D5_cases:
            delta_keff = (case.DRAGON_Keff - self.S2_case.keffs)*1e5 # error on Keff in pcm
            self.delta_keffs[f"{case.ssh_method}_{case.correlation}_{case.sat}_{case.depl_sol}"] = delta_keff
        return
    def compare_Ni(self):
        for iso in self.tracked_nuclides:
            delta_Niso_case = {}
            for case in self.D5_cases:
                delta_Niso = [(case.DRAGON_ISOTOPESDENS[iso][idx] - self.S2_case.Ni[iso][idx]) * 100 /  self.S2_case.Ni[iso][idx]
                    if self.S2_case.Ni[iso][idx] != 0 else 0
                    for idx in range(len(self.S2_case.Ni[iso]))]
                delta_Niso_case[f"{case.ssh_method}_{case.correlation}_{case.sat}_{case.depl_sol}"] = delta_Niso
            self.delta_Niso[f"{iso}"] = delta_Niso_case
        return

    def plot_delta_Keff(self):
        """
        Plot the delta Keff for all cases : all D5 cases compared to S2 case
        """
        plt.figure()
        for comparison_case in self.delta_keffs.keys():
            plt.plot(self.S2_case.BUdays, self.delta_keffs[comparison_case], label = f"{comparison_case}".replace("_"," "), marker = "x", linestyle = "--")
        plt.xlabel("Burnup [days]")
        plt.ylabel("$\\Delta$ Keff [pcm]")
        plt.axhline(y = 300.0, color = 'r', linestyle = '-')
        plt.axhline(y = -300.0, color = 'r', linestyle = '-')  
        plt.title(f"$\\Delta$ Keff evolution for {self.comparison_name} case")
        plt.legend()
        plt.grid()
        plt.savefig(f"{self.save_dir}/Delta_Keff_{self.comparison_name}_cst_flx.png")
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
                plt.plot(self.S2_case.BUdays, self.delta_Niso[iso][comparison_case], label = f"{comparison_case}".replace("_"," "), marker = "x", linestyle = "--")
            plt.xlabel("Burnup [days]")
            plt.ylabel(f"$\\Delta$ N{iso} [%]")
            plt.axhline(y = 2.0, color = 'r', linestyle = '-')
            plt.axhline(y = -2.0, color = 'r', linestyle = '-')  
            plt.title(f"$\\Delta$ N{iso} evolution for {self.comparison_name} case")
            plt.legend()
            plt.grid()
            plt.savefig(f"{self.save_dir}/Delta_{iso}_{self.comparison_name}_cst_flx.png")
            plt.close()
        return