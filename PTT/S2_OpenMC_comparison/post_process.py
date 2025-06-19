### Python3 classes for S2 / OpenMC results handling and plotting
## Author : R. Guasch
## Date : 2025-03-24


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
import serpentTools as st

# Define the Serpent2_case class

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
            self.res = st.read(f"{os.environ['SERPENT_RESULTS']}/{self.case_name}/{self.case_name}_qfiss_{self.lib_name}_edep{self.edep_id}_pcc{self.pcc_id}_mc_res.m")
            self.depl = st.read(f"{os.environ['SERPENT_RESULTS']}/{self.case_name}/{self.case_name}_qfiss_{self.lib_name}_edep{self.edep_id}_pcc{self.pcc_id}_mc_dep.m")
        else:
            self.res = st.read(f"{os.environ['SERPENT_RESULTS']}/{self.case_name}/{self.case_name}_{self.lib_name}_edep{self.edep_id}_pcc{self.pcc_id}_mc_res.m")
            self.depl = st.read(f"{os.environ['SERPENT_RESULTS']}/{self.case_name}/{self.case_name}_{self.lib_name}_edep{self.edep_id}_pcc{self.pcc_id}_mc_dep.m")
        self.keffs = self.res.resdata["absKeff"].T[0] # Retrieve the keffs
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
            path_to_results = f"{os.environ['OPENMC_RESULTS']}/{self.case_name}/{self.edep_id}_{self.integrator}"
        self.keffs = np.loadtxt(f"{path_to_results}/{self.case_name}_depl_keff.txt")
        self.sigmas_keff = self.keffs[:,1]
        self.keffs = self.keffs[:,0]
        self.BUdays = np.loadtxt(f"{path_to_results}/{self.case_name}_depl_time.txt")
        
        print(f"Reading OpenMC results for {self.case_name} with {self.edep_id} and {self.integrator}")
        # Sanity check
        print(f"OpenMC keffs shape : {self.keffs.shape}")
        print(f"OpenMC BU_days shape : {self.BUdays.shape}")

        self.Ni = {}
        for isotope in self.tracked_nuclides:
            self.Ni[isotope] = np.loadtxt(f"{path_to_results}/{self.case_name}_depl_N{isotope}.txt")

        return
    

