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
        self.res = st.read(f"{os.environ['SERPENT_RESULTS']}/OpenMC_comparisons/{self.case_name}_OMC_{self.lib_name}_edep{self.edep_id}_pcc{self.pcc_id}_mc_res.m")
        self.depl = st.read(f"{os.environ['SERPENT_RESULTS']}/OpenMC_comparisons/{self.case_name}_OMC_{self.lib_name}_edep{self.edep_id}_pcc{self.pcc_id}_mc_dep.m")
        
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



class OpenMC_case:
    def __init__(self, case_name, lib_name, edep_id, integrator, specific_power, tracked_nuclides, save_dir):
        self.case_name = case_name
        self.lib_name = lib_name
        self.edep_id = edep_id
        self.integrator = integrator
        self.specific_power = specific_power
        self.tracked_nuclides = tracked_nuclides
        self.save_dir = save_dir
        self.unitsBU = "days"



    def read_OpenMC_outputs(self):

        path_to_results = f"{os.environ['OPENMC_RESULTS']}/{self.case_name}/"
        