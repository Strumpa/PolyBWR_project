# Python3 post treatment script for constant power evolution case
# Author : R. Guasch
# Date : 2025/02/10
import numpy as np
import matplotlib.pyplot as plt
import serpentTools as st
import os

# S2 PT class :

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

# Dragon case post treatment class :
    
class DRAGON_case:
    def __init__(self, pyCOMPO, dlib_name, bu_points, ssh_opt, correlation, sat, depl_sol, tracked_nuclides, BU_lists, save_dir):
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
        self.ssh_opt = ssh_opt
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

        self.EVO_BU_steps = BU_lists["BU"]
        self.AUTOP_BU_steps = BU_lists["AUTOP"] # only relevant for PCC0 and PCC2 schemes as ssh is performed at all preditor-corrector steps in PCC1, PCC3 and PCC3b 
        self.COMPO_BU_steps = BU_lists["COMPO"]

        self.BU = None
        self.keff = None
        self.DRAGON_ISOTOPESDENS = {}

        self.parse_compo()

        self.BUScheme = "LE" # Default DRAGON burnup scheme is LE (Linear Extrapolation of reaction rates over time steps)

        return
        

    def parse_compo(self):
        lenBU_DRAGON=np.shape(self.COMPO_BU_steps)[0]
        ISOTOPES=self.pyCOMPO[self.DIR]['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESDENS']
        print(f"Dragon isotopes = {ISOTOPES}")
        #lenISOT_DRAGON=np.shape(ISOTOPES)[0]-1, check this for consistency : the last isotope is the total density "MAC"
        lenISOT_DRAGON=np.shape(ISOTOPES)[0]
        DRAGON_BU=self.COMPO_BU_steps
        DRAGON_ISOTOPESDENS=np.zeros((lenISOT_DRAGON,lenBU_DRAGON))
        DRAGON_Keff=np.zeros(lenBU_DRAGON)
        print(f"DRAGON_ISOTOPESDENS shape = {lenISOT_DRAGON} {lenBU_DRAGON}")
        for k in range(lenBU_DRAGON):
            DRAGON_Keff[k]=self.pyCOMPO[self.DIR]['MIXTURES'][0]['CALCULATIONS'][k]['K-EFFECTIVE']
            for j in range(lenISOT_DRAGON):
                DRAGON_ISOTOPESDENS[j][k]=self.pyCOMPO[self.DIR]['MIXTURES'][0]['CALCULATIONS'][k]['ISOTOPESDENS'][j]

        # --------- List of isotopes from the DRAGON Multicompo results
        isotopes2=[]
        isotopes=[]
        for k in range(lenISOT_DRAGON):
            isotopes2.append(self.pyCOMPO[self.DIR]['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][k]['ALIAS'])
        for k in range(lenISOT_DRAGON):
            if isotopes2[k][0]=='U':
                isotopes=isotopes+[isotopes2[k][0:4]]
            else:
                isotopes=isotopes+[isotopes2[k][0:5]]

        indices=np.zeros(len(self.tracked_nuclides))
        for n in range(len(self.tracked_nuclides)):
            for m in range(len(isotopes)):
                if self.tracked_nuclides[n]==isotopes[m]:
                    #print(f"tracked_nuclides[{n}] = {self.tracked_nuclides[n]} is isotopes[{m}] = {isotopes[m]}")
                    indices[n]=m

        for k in range(len(self.tracked_nuclides)):
            self.DRAGON_ISOTOPESDENS[self.tracked_nuclides[k]] = DRAGON_ISOTOPESDENS[int(indices[k])]
        self.BU = DRAGON_BU
        self.keff = DRAGON_Keff

        print(f"DRAGON_BU = {self.BU} with shape = {np.shape(self.BU)}")
        print(f"DRAGON_Keff = {self.keff} with shape = {np.shape(self.keff)}")

        return
    
    def set_BUscheme(self, BUscheme, extrapolation_type):
        """
        Set the burnup scheme for the DRAGON case
        BUscheme : (str) "predictor" or "precitor-corrector"
        extrapolation_type : (str) "EXTR" or "NOEX" if predictor or ("EXTR" or "NOEX", "EXTR" no "EXTR") if predictor-corrector
        """
        self.BUScheme = BUscheme
        self.extrapolation_type = extrapolation_type
        return

    def plot_keffs(self):
        plt.figure()
        plt.plot(self.BU, self.keff, label = f"D5 {self.draglib_name} : sat {self.sat}, depl sol {self.depl_sol}", marker = "x", linestyle = "--")
        plt.xlabel("Burnup [MWd/tU]")
        plt.ylabel("Keff")
        plt.title(f"Keff evolution for {self.draglib_name} case")
        plt.legend()
        plt.grid()
        plt.savefig(f"{self.save_dir}/DRAGON_Keff_{self.draglib_name}_{self.ssh_opt}_sat_{self.sat}_sol{self.depl_sol}.png")
        plt.close()
        return
    
    def plot_Ni(self, isotope):
        plt.figure()
        plt.plot(self.BU, self.DRAGON_ISOTOPESDENS[isotope], label = f"D5 : {isotope} atomic density", marker = "x", linestyle = "--")
        plt.xlabel("Burnup [MWd/tU]")
        plt.ylabel("Isotope density")
        plt.title(f"{isotope} evolution for {self.draglib_name} case")
        plt.legend()
        plt.grid()
        plt.savefig(f"{self.save_dir}/DRAGON_{isotope}_{self.draglib_name}_{self.ssh_opt}_sat_{self.sat}_sol{self.depl_sol}.png")
        plt.close()
        return
    
class OpenMC_case:
    def __init__(self, case_name, lib_name, edep_id, areQfissSet, integrator, specific_power, tracked_nuclides, save_dir):
        self.case_name = case_name
        self.lib_name = lib_name
        self.edep_id = edep_id
        self.areQfissSet = areQfissSet
        if edep_id == "energy_deposition":
            self.areQfissSet = False # If energy deposition mode is set, Qfiss values are not set as KERMA coefficients are used instead
        self.integrator = integrator
        self.specific_power = specific_power
        self.tracked_nuclides = tracked_nuclides
        self.save_dir = save_dir
        self.unitsBU = "days"
        self.read_OpenMC_outputs()

    def read_OpenMC_outputs(self):
        # Read the results and depletion
        if self.areQfissSet == True and self.edep_id == "fissq":
            path_to_results = f"{os.environ['OPENMC_RESULTS']}/{self.case_name}/{self.edep_id}_{self.integrator}/results_set_qfiss"
        elif self.areQfissSet == False and self.edep_id == "fissq":
            path_to_results = f"{os.environ['OPENMC_RESULTS']}/{self.case_name}/{self.edep_id}_{self.integrator}/results_default_Q_values"
        elif self.areQfissSet == False and self.edep_id == "energy_deposition":
            path_to_results = f"{os.environ['OPENMC_RESULTS']}/{self.case_name}/{self.edep_id}_{self.integrator}/"
        self.keff = np.loadtxt(f"{path_to_results}/{self.case_name}_depl_keff.txt")
        self.sigmas_keff = self.keff[:,1]
        self.keff = self.keff[:,0]
        self.BUdays = np.loadtxt(f"{path_to_results}/{self.case_name}_depl_time.txt")
        self.BU = self.BUdays * self.specific_power # Convert the time in days to burnup in MWd/tU
        
        print(f"Reading OpenMC results for {self.case_name} with {self.edep_id} and {self.integrator}")
        # Sanity check
        print(f"OpenMC keffs shape : {self.keff.shape}")
        print(f"OpenMC BU_days shape : {self.BUdays.shape}")

        self.Ni = {}
        for isotope in self.tracked_nuclides:
            self.Ni[isotope] = np.loadtxt(f"{path_to_results}/{self.case_name}_depl_N{isotope}.txt")

        return
    

# Comparison between several DRAGON cases and 1 reference Serpent2 case:
class multiD5S2_comparisons:
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
        
        return

    def compare_keffs(self):
        for case in self.D5_cases:
            delta_keff = (case.DRAGON_Keff - self.S2_case.keffs)*1e5 # error on Keff in pcm
            self.delta_keffs[f"{case.draglib_name}_{case.ssh_opt}_{case.correlation}_to_S2_edep{self.S2_case.edep_id}"] = delta_keff
        return
    def compare_Ni(self):
        for iso in self.tracked_nuclides:
            delta_Niso_case = {}
            for case in self.D5_cases:
                delta_Niso = [(case.DRAGON_ISOTOPESDENS[iso][idx] - self.S2_case.Ni[iso][idx]) * 100 / self.S2_case.Ni[iso][idx]
                    if self.S2_case.Ni[iso][idx] != 0 else 0
                    for idx in range(len(self.S2_case.Ni[iso]))]
                delta_Niso_case[f"{case.draglib_name}_{case.ssh_opt}_{case.correlation}_to_S2_edep{self.S2_case.edep_id}"] = delta_Niso
            self.delta_Niso[f"{iso}"] = delta_Niso_case
        return

    def plot_delta_Keff(self):
        """
        Plot the delta Keff for all cases : 1 D5 case compared to several S2 cases
        """
        plt.figure()
        for comparison_case in self.delta_keffs.keys():
            plt.plot(self.S2_case.BU, self.delta_keffs[comparison_case], label = f"{comparison_case}".replace("_"," "), marker = "x", linestyle = "--")
        plt.xlabel(f"Burnup [{self.S2_case.unitsBU}]")
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
        Plot the delta Ni for all cases : several D5 cases compared to one S2 case
        """
        for iso in self.tracked_nuclides:
            plt.figure()
            for comparison_case in self.delta_Niso[iso].keys():
                plt.plot(self.S2_case.BU, self.delta_Niso[iso][comparison_case], label = f"{comparison_case}".replace("_"," "), marker = "x", linestyle = "--")
            plt.xlabel(f"Burnup [{self.S2_case.unitsBU}]")
            plt.ylabel(f"$\\Delta$ N{iso} [%]")
            plt.axhline(y = 2.0, color = 'r', linestyle = '-')
            plt.axhline(y = -2.0, color = 'r', linestyle = '-')  
            plt.title(f"$\\Delta$ N{iso} evolution for {self.comparison_name} case")
            plt.legend()
            plt.grid()
            plt.savefig(f"{self.save_dir}/Delta_{iso}_{self.comparison_name}.png")
            plt.close()
        return
    
# Comparison between 1 DRAGON case and several Serpent2 cases:
class D5multiS2_comparisons:
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
            self.delta_keffs[f"{self.D5_case.draglib_name}_{self.D5_case.ssh_opt}_{self.D5_case.correlation}_to_S2_edep{case.edep_id}"] = delta_keff
        return
    def compare_Ni(self):
        for iso in self.tracked_nuclides:
            delta_Niso_case = {}
            for case in self.S2_cases:
                delta_Niso = [(self.D5_case.DRAGON_ISOTOPESDENS[iso][idx] - case.Ni[iso][idx]) * 100 / case.Ni[iso][idx]
                    if case.Ni[iso][idx] != 0 else 0
                    for idx in range(len(case.Ni[iso]))]
                delta_Niso_case[f"{self.D5_case.draglib_name}_{self.D5_case.ssh_opt}_{self.D5_case.correlation}_to_S2_edep{case.edep_id}"] = delta_Niso
            self.delta_Niso[f"{iso}"] = delta_Niso_case
        return

    def plot_delta_Keff(self):
        """
        Plot the delta Keff for all cases : 1 D5 case compared to several S2 cases
        """
        plt.figure()
        for comparison_case in self.delta_keffs.keys():
            plt.plot(self.D5_case.DRAGON_BU, self.delta_keffs[comparison_case], label = f"{comparison_case}".replace("_"," "), marker = "x", linestyle = "--")
        plt.xlabel(f"Burnup [MWd/tU]")
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
        Plot the delta Ni for all cases : several D5 cases compared to one S2 case
        """
        for iso in self.tracked_nuclides:
            plt.figure()
            for comparison_case in self.delta_Niso[iso].keys():
                plt.plot(self.D5_case.DRAGON_BU, self.delta_Niso[iso][comparison_case], label = f"{comparison_case}".replace("_"," "), marker = "x", linestyle = "--")
            plt.xlabel(f"Burnup [MWd/tU]")
            plt.ylabel(f"$\\Delta$ N{iso} [%]")
            plt.axhline(y = 2.0, color = 'r', linestyle = '-')
            plt.axhline(y = -2.0, color = 'r', linestyle = '-')  
            plt.title(f"$\\Delta$ N{iso} evolution for {self.comparison_name} case")
            plt.legend()
            plt.grid()
            plt.savefig(f"{self.save_dir}/Delta_{iso}_{self.comparison_name}.png")
            plt.close()
        return
    
def interp_on_BU(BU_points_to_interp, BU_points_ref, values_ref):
    """
    Interpolate the values_ref on the BU_points_to_interp using linear interpolation
    """
    interp_values = np.interp(BU_points_to_interp, BU_points_ref, values_ref)
    return interp_values