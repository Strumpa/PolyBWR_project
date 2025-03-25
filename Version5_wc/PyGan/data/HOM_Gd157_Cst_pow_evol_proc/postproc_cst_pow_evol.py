# Python3 post treatment script for constant power evolution case
# Author : R. Guasch
# Date : 2025/02/10
import numpy as np
import matplotlib.pyplot as plt
import serpentTools as st
import os
import warnings
from typing import Dict, List

# S2 PT class :

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
        elif self.edep_id == 1 or self.edep_id == "1_Ecapt":
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
        self.res = st.read(f"{os.environ['SERPENT_RESULTS']}/HOM_CELL_study/{self.case_name}/BUScheme_EDEP_PCC_study/{self.case_name}_{self.lib_name}_edep{self.edep_id}_pcc{self.pcc_id}_mc_res.m")
        self.depl = st.read(f"{os.environ['SERPENT_RESULTS']}/HOM_CELL_study/{self.case_name}/BUScheme_EDEP_PCC_study/{self.case_name}_{self.lib_name}_edep{self.edep_id}_pcc{self.pcc_id}_mc_dep.m")
        
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

    def plot_keff(self):
        plt.figure()
        plt.plot(self.BU, self.keffs, label = f"edepmode {self.edep_id} pcc {self.pcc_id}", marker = "x", linestyle = "--")
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

        self.EVO_BU_steps = BU_lists["BU"]
        self.AUTOP_BU_steps = BU_lists["AUTOP"] # only relevant for PCC0 and PCC2 schemes as ssh is performed at all preditor-corrector steps in PCC1, PCC3 and PCC3b 
        self.COMPO_BU_steps = BU_lists["COMPO"]

        self.DRAGON_BU = None
        self.DRAGON_Keff = None
        self.DRAGON_ISOTOPESDENS = {}

        self.parse_compo()

        self.BUScheme = "LE" # Default DRAGON burnup scheme is LE (Linear Extrapolation of reaction rates over time steps)

        return
        

    def parse_compo(self):
        """
        Parse and extract pyCOMPO object obtained from DRAGON5.

        Extracts burnup steps, keff values, and isotopic densities from the pyCOMPO object.
        """

        # Ensure required keys exist before proceeding
        try:
            BU_DRAGON = np.array(self.pyCOMPO[self.DIR]['GLOBAL']['pval00000001'], dtype=float)
        except KeyError as e:
            raise KeyError(f"Missing key in pyCOMPO structure: {e}")

        lenBU_DRAGON = BU_DRAGON.shape[0]

        # Check burnup step consistency
        if lenBU_DRAGON != len(self.COMPO_BU_steps):
            warnings.warn("Mismatch between BU steps from ListCOMPO and pyCOMPO.", RuntimeWarning)

        try:
            ISOTOPES = self.pyCOMPO[self.DIR]['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESDENS']
            ISOTOPE_LIST = []
            for i in range(len(ISOTOPES)):
                ISOTOPE_LIST.append(self.pyCOMPO[self.DIR]['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][i]['ALIAS'].split(" ")[0].strip())
        except KeyError as e:
            raise KeyError(f"Missing isotope data in pyCOMPO: {e}")

        lenISOT_DRAGON = len(ISOTOPES)

        # Initialize storage arrays
        DRAGON_BU = np.array(self.COMPO_BU_steps, dtype=float)
        DRAGON_ISOTOPESDENS = np.zeros((lenISOT_DRAGON, lenBU_DRAGON), dtype=float)
        DRAGON_Keff = np.zeros(lenBU_DRAGON, dtype=float)

        # Extract keff and isotope densities
        for k in range(lenBU_DRAGON):
            try:
                DRAGON_Keff[k] = self.pyCOMPO[self.DIR]['MIXTURES'][0]['CALCULATIONS'][k]['K-EFFECTIVE']
                for j in range(lenISOT_DRAGON):
                    DRAGON_ISOTOPESDENS[j, k] = self.pyCOMPO[self.DIR]['MIXTURES'][0]['CALCULATIONS'][k]['ISOTOPESDENS'][j]
            except KeyError as e:
                warnings.warn(f"Missing key at BU step {k}: {e}. Skipping this step.", RuntimeWarning)
                continue

        # Extract isotope names
        print(ISOTOPE_LIST)

        isotopes = ISOTOPE_LIST

        # Map tracked nuclides to their indices in the extracted isotope list
        indices = np.full(len(self.tracked_nuclides), -1, dtype=int)
        for n, nuclide in enumerate(self.tracked_nuclides):
            if nuclide in isotopes:
                indices[n] = isotopes.index(nuclide)

        # Check for missing nuclides
        missing_nuclides = [self.tracked_nuclides[n] for n, idx in enumerate(indices) if idx == -1]
        if missing_nuclides:
            warnings.warn(f"Some tracked nuclides not found in DRAGON data: {missing_nuclides}", RuntimeWarning)

        # Assign isotopic densities
        for k, idx in enumerate(indices):
            if idx != -1:  # Only assign if index was found
                self.DRAGON_ISOTOPESDENS[self.tracked_nuclides[k]] = DRAGON_ISOTOPESDENS[idx]

        # Store parsed data
        self.DRAGON_BU = DRAGON_BU
        self.DRAGON_Keff = DRAGON_Keff

    def plot_keffs(self):
        plt.figure()
        plt.plot(self.DRAGON_BU, self.DRAGON_Keff, label = f"D5 {self.draglib_name} : sat {self.sat}, depl sol {self.depl_sol}", marker = "x", linestyle = "--")
        plt.xlabel("Burnup [MWd/tU]")
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
        plt.xlabel("Burnup [MWd/tU]")
        plt.ylabel("Isotope density")
        plt.title(f"{isotope} evolution for {self.draglib_name} case")
        plt.legend()
        plt.grid()
        plt.savefig(f"{self.save_dir}/DRAGON_{isotope}_{self.draglib_name}_{self.ssh_module}_{self.ssh_method}_sat_{self.sat}_sol{self.depl_sol}.png")
        plt.close()
        return
    

# Comparison between several DRAGON cases and 1 reference Serpent2 case:
class multiD5S2_comparisons:
    def __init__(self, comparison_name, D5_cases, S2_case, tracked_nuclides, save_dir):
        self.comparison_name = comparison_name
        self.D5_cases = D5_cases
        self.S2_case = S2_case
        self.delta_keffs = {}
        self.delta_Niso = {} # absolute differences in atom/b-cm, consider plotting on log scale
        self.relative_delta_Niso = {} # relative differences in %
        self.tracked_nuclides = tracked_nuclides
        self.save_dir = f"{save_dir}/multiD5S2"
        # check if the save directory exists and create it if not
        if not os.path.exists(self.save_dir):
            os.makedirs(self.save_dir)
        self.D5_BU_lists = {}

        return
    
    def interpolate_values(self, x_source, y_source, x_target):
        """Interpolates values of y_source at x_target using x_source as reference."""
        interpolated_values = np.interp(x_target, x_source, y_source)
        return interpolated_values

    def compare_keffs(self, BU_interp: bool) -> None:
        """
        Compare keff values (in pcm) between the DRAGON cases and the S2 case.
        
        BU_interp : (bool) 
            If True, the DRAGON keffs are interpolated to the S2 burnup steps.
        """
        self.BU_interp = BU_interp

        # Ensure dictionary exists before modifying it
        if not hasattr(self, "delta_keffs"):
            self.delta_keffs = {}

        if not hasattr(self.S2_case, "BU") or not hasattr(self.S2_case, "keffs"):
            warnings.warn("S2_case is missing BU or keffs data. Aborting keff comparison.")
            return

        S2_BU = np.asarray(self.S2_case.BU)
        S2_keffs = np.asarray(self.S2_case.keffs)

        if len(S2_BU) == 0 or len(S2_keffs) == 0:
            warnings.warn("S2 burnup or keff data is empty. Aborting keff comparison.")
            return

        for case in self.D5_cases:
            case_key = f"{case.draglib_name}_{case.ssh_module}_{case.ssh_method}_{case.correlation}_to_S2_edep{self.S2_case.edep_id}"

            # Ensure required attributes exist in the case
            if not hasattr(case, "DRAGON_BU") or not hasattr(case, "DRAGON_Keff"):
                warnings.warn(f"Missing burnup or keff data for case {case_key}. Skipping...")
                continue

            dragon_BU = np.asarray(case.DRAGON_BU)
            self.D5_BU_lists[case_key] = dragon_BU
            dragon_Keff = np.asarray(case.DRAGON_Keff)

            if BU_interp:
                # It has been noticed that BU S2 goes over the last BU of the D5 cases, 
                # interpolate S2 keffs on D5 BU points
                S2_interp = self.interpolate_values(S2_BU, S2_keffs, dragon_BU)
            else:
                if len(dragon_Keff) != len(S2_keffs):
                    warnings.warn(f"Mismatch in burnup step sizes for case {case_key}. Skipping...")
                    continue
                S2_interp = S2_keffs

            # Compute difference in pcm
            delta_keff = (dragon_Keff - S2_interp) * 1e5

            # Store result
            self.delta_keffs[case_key] = delta_keff

            self.plot_delta_Keff()

    def compare_Ni(self, BU_interp: bool, abs: bool, rel: bool) -> None:
        """
        Compare the isotopic densities of the DRAGON cases to the S2 case for a given isotope.

        Depending on the energy deposition mode selected in Serpent2.1.32, the burnup steps may differ.
        If `BU_interp` is True, the DRAGON isotopic densities are interpolated to the S2 burnup steps.

        Results are stored in:
        - `self.relative_delta_Niso[iso][case_name]`: Relative differences (%) at each BU step.
        - `self.delta_Niso[iso][case_name]`: Absolute differences (atom/b-cm) at each BU step.

        Parameters:
        BU_interp (bool): If True, interpolate DRAGON isotopic densities to S2 burnup steps.
        abs (bool): If True, plot absolute differences.
        rel (bool): If True, plot relative differences.
        """
        self.BU_interp = BU_interp

        for iso in self.tracked_nuclides:
            print(f"treating {iso} in compare_Ni mutliD5S2")
            delta_Niso_case: Dict[str, np.ndarray] = {}
            relative_delta_Niso_case: Dict[str, np.ndarray] = {}
            for case in self.D5_cases:
                case_key = f"{case.draglib_name}_{case.ssh_module}_{case.ssh_method}_{case.correlation}_to_S2_edep{self.S2_case.edep_id}"

                # Ensure required attributes exist
                if not hasattr(case, "DRAGON_BU") or not hasattr(case, "DRAGON_ISOTOPESDENS") or iso not in case.DRAGON_ISOTOPESDENS:
                    warnings.warn(f"Missing burnup or isotope data for case {case_key}. Skipping...")
                    continue

                # Convert to numpy arrays for safety
                dragon_BU = np.asarray(case.DRAGON_BU)
                dragon_Niso = np.asarray(case.DRAGON_ISOTOPESDENS[iso])
                S2_BU = np.asarray(self.S2_case.BU)
                print(f"S2_BU : {S2_BU}")
                S2_Ni = np.asarray(self.S2_case.Ni.get(iso, []))

                if len(S2_BU) == 0 or len(S2_Ni) == 0:
                    warnings.warn(f"S2 burnup data is missing for isotope {iso}. Skipping...")
                    continue

                if BU_interp:
                    # Interpolating S2 values between D5 burnup points
                    S2_interp = self.interpolate_values(S2_BU, S2_Ni, dragon_BU)
                else:
                    if len(dragon_Niso) != len(S2_Ni):
                        warnings.warn(f"Mismatch in burnup step sizes for isotope {iso} in case {case_key}. Skipping...")
                        continue
                    S2_interp = S2_Ni # no interpolation

                # Compute absolute differences
                delta_Niso = dragon_Niso - S2_interp
                # Compute relative differences safely
                relative_delta_Niso = [delta_Niso[i] * 100 / S2_interp[i] if S2_interp[i] != 0 else 0 for i in range(len(S2_interp))]

                # Store results
                delta_Niso_case[case_key] = delta_Niso
                relative_delta_Niso_case[case_key] = relative_delta_Niso

                self.delta_Niso[iso] = delta_Niso_case
                self.relative_delta_Niso[iso] = relative_delta_Niso_case
        print(f"self.delta_Niso[iso] : {self.delta_Niso[iso]}")
        self.plot_delta_Ni(abs, rel)
        

    def plot_delta_Keff(self):
        """
        Plot the delta Keff for all cases : D5 cases compared to 1 "reference" S2 case

        if BU interpolation is activated, the keff differences plotted are based on D5 keffs interpolated on the S2 BU points
        """
        # plot absolute differences
        plt.figure()
        for comparison_case in self.delta_keffs.keys():
            plt.plot(self.D5_BU_lists[comparison_case], self.delta_keffs[comparison_case], label = f"{comparison_case}".replace("_"," "), marker = "x", linestyle = "--")
        plt.xlabel(f"Burnup [MWd/tU]")
        plt.ylabel("$\\Delta$ Keff [pcm]")
        plt.axhline(y = 300.0, color = 'r', linestyle = '-')
        plt.axhline(y = -300.0, color = 'r', linestyle = '-')  
        plt.title(f"$\\Delta$ Keff evolution for {self.comparison_name} case")
        plt.legend()
        plt.grid()
        if self.BU_interp:
            plt.savefig(f"{self.save_dir}/Delta_Keff_{self.comparison_name}_interp.png")  
        else:
            plt.savefig(f"{self.save_dir}/Delta_Keff_{self.comparison_name}.png")
        plt.close()

        return
    
    def plot_delta_Ni(self, plot_abs: bool = True, plot_rel: bool = True):
        """
        Plot the delta Ni for all cases : several D5 cases compared to one S2 case
        """
        for iso in self.tracked_nuclides:
            print(f"plotting {iso}")
            if plot_abs:
                # plot absolute differences
                plt.figure()
                for comparison_case in self.delta_Niso[iso].keys():
                    plt.plot(self.S2_case.BU, self.delta_Niso[iso][comparison_case], label = f"{comparison_case}".replace("_"," "), marker = "x", linestyle = "--")
                plt.xlabel(f"Burnup [{self.S2_case.unitsBU}]")
                plt.ylabel(f"$\\Delta$ N{iso} [atom/b-cm]") 
                plt.title(f"$\\Delta$ N{iso} evolution for {self.comparison_name} case")
                plt.legend()
                plt.grid()
                if self.BU_interp:
                    plt.savefig(f"{self.save_dir}/abs_Delta_{iso}_{self.comparison_name}_interp.png")  
                else: 
                    plt.savefig(f"{self.save_dir}/abs_Delta_{iso}_{self.comparison_name}.png")
                plt.close()

            if plot_rel:
                # plot absolute differences
                plt.figure()
                for comparison_case in self.delta_Niso[iso].keys():
                    plt.plot(self.D5_BU_lists[comparison_case], self.relative_delta_Niso[iso][comparison_case], label = f"{comparison_case}".replace("_"," "), marker = "x", linestyle = "--")
                plt.xlabel(f"Burnup [{self.S2_case.unitsBU}]")
                plt.ylabel(f"$\\Delta$ N{iso} [atom/b-cm]")
                plt.title(f"$\\Delta$ N{iso} evolution for {self.comparison_name} case")
                plt.ylabel(f"$\\Delta$ N{iso} [%]")
                plt.axhline(y = 2.0, color = 'r', linestyle = '-')
                plt.axhline(y = -2.0, color = 'r', linestyle = '-') 
                plt.legend()
                plt.grid()
                if self.BU_interp:
                    plt.savefig(f"{self.save_dir}/rel_Delta_{iso}_{self.comparison_name}_interp.png")  
                else:
                    plt.savefig(f"{self.save_dir}/rel_Delta_{iso}_{self.comparison_name}.png")
                plt.close()
        return


class D5multiS2_comparisons:
    def __init__(self, comparison_name, D5_case, S2_cases, tracked_nuclides, save_dir):
        self.comparison_name = comparison_name
        self.D5_case = D5_case
        self.S2_cases = S2_cases
        self.tracked_nuclides = tracked_nuclides
        self.save_dir = f"{save_dir}/D5multiS2"
        self.delta_keffs = {}
        self.delta_Niso_abs = {}  # Absolute differences
        self.delta_Niso_rel = {}  # Relative differences

        # Ensure save directory exists
        os.makedirs(self.save_dir, exist_ok=True)

    def interpolate_values(self, x_source, y_source, x_target):
        """Interpolates values of y_source at x_target using x_source as reference."""
        interpolated_values = np.interp(x_target, x_source, y_source)
        print(f"Interpolated values: {interpolated_values}")
        return interpolated_values

    def compare_keffs(self, BU_interp: bool = True):
        """
        Compare keff values between DRAGON5 and multiple Serpent2 cases with option to interpolate.
        
        BU_interp (bool): If True, interpolate DRAGON keffs to S2 burnup steps.
        """
        self.BU_interp = BU_interp
        burnup_D5 = self.D5_case.DRAGON_BU  # Reference burnup points

        for case in self.S2_cases:
            burnup_S2 = case.BU
            if BU_interp:
                keffs_S2_interp = self.interpolate_values(burnup_S2, case.keffs, burnup_D5)
            else:
                if len(burnup_D5) != len(case.keffs):
                    warnings.warn("Mismatch in burnup step sizes. Skipping comparison.")
                    continue
                keffs_S2_interp = case.keffs

            delta_keff = (self.D5_case.DRAGON_Keff - keffs_S2_interp) * 1e5  # Error in pcm
            key = f"{self.D5_case.draglib_name}_{self.D5_case.ssh_module}_{self.D5_case.ssh_method}_{self.D5_case.correlation}_to_S2_edep{case.edep_id}"
            self.delta_keffs[key] = delta_keff

        self.plot_delta_Keff()

    def compare_Ni(self, BU_interp: bool = True, abs: bool = True, rel: bool = True):
        """
        Compare isotope densities between DRAGON5 and multiple Serpent2 cases with interpolation.
        
        Parameters:
        BU_interp (bool): If True, interpolate DRAGON isotope densities to S2 burnup steps.
        abs (bool): If True, plot absolute differences.
        rel (bool): If True, plot relative differences.
        """
        self.BU_interp = BU_interp
        burnup_D5 = self.D5_case.DRAGON_BU  # Reference burnup points

        for iso in self.tracked_nuclides:
            delta_Niso_abs_case = {}
            delta_Niso_rel_case = {}
            for case in self.S2_cases:
                burnup_S2 = case.BU
                key = f"{self.D5_case.draglib_name}_{self.D5_case.ssh_module}_{self.D5_case.ssh_method}_{self.D5_case.correlation}_to_S2_edep{case.edep_id}"
                if iso not in self.D5_case.DRAGON_ISOTOPESDENS or iso not in case.Ni:
                    print(f"Warning: Isotope {iso} missing in one of the cases.")
                    delta_Niso_abs_case[key] = np.zeros(len(burnup_D5))
                    delta_Niso_rel_case[key] = np.zeros(len(burnup_D5))
                    continue

                # Interpolate S2 isotope densities at D5 burnup points
                if BU_interp:
                    Ni_S2_interp = self.interpolate_values(burnup_S2, case.Ni[iso], burnup_D5)
                else:
                    if len(burnup_D5) != len(case.Ni[iso]):
                        warnings.warn("Mismatch in burnup step sizes. Skipping comparison.")
                        continue
                    Ni_S2_interp = case.Ni[iso]

                # Absolute difference: ΔN = N_D5 - N_S2
                delta_N_abs = self.D5_case.DRAGON_ISOTOPESDENS[iso] - Ni_S2_interp

                # Relative difference: ΔN% = 100 × (N_D5 - N_S2) / N_S2 (avoiding division by zero)
                delta_N_rel = np.array([delta_N_abs[i] * 100 / Ni_S2_interp[i] 
                                        if Ni_S2_interp[i] != 0 else 0 for i in range(len(Ni_S2_interp))])

                key = f"{self.D5_case.draglib_name}_{self.D5_case.ssh_module}_{self.D5_case.ssh_method}_{self.D5_case.correlation}_to_S2_edep{case.edep_id}"
                delta_Niso_abs_case[key] = delta_N_abs
                delta_Niso_rel_case[key] = delta_N_rel

            self.delta_Niso_abs[iso] = delta_Niso_abs_case
            self.delta_Niso_rel[iso] = delta_Niso_rel_case

        self.plot_delta_Ni(abs,rel)


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
        if self.BU_interp:
            plt.savefig(f"{self.save_dir}/Delta_Keff_{self.comparison_name}_interp.png")
        else:
            plt.savefig(f"{self.save_dir}/Delta_Keff_{self.comparison_name}.png")
        plt.close()
        return
    
    def plot_delta_Ni(self, plot_abs: bool = True, plot_rel: bool = True):
        """
        Plot the delta Ni for all cases : several D5 cases compared to one S2 case
        """
        for iso in self.tracked_nuclides:
            print(f"plotting {iso}")
            if plot_rel:
                plt.figure()
                for comparison_case in self.delta_Niso_rel[iso].keys():
                    plt.plot(self.D5_case.DRAGON_BU, self.delta_Niso_rel[iso][comparison_case], label = f"{comparison_case}".replace("_"," "), marker = "x", linestyle = "--")
                plt.xlabel(f"Burnup [MWd/tU]")
                plt.ylabel(f"$\\Delta$ N{iso} [%]")
                plt.axhline(y = 2.0, color = 'r', linestyle = '-')
                plt.axhline(y = -2.0, color = 'r', linestyle = '-')  
                plt.title(f"$\\Delta$ N{iso} evolution for {self.comparison_name} case")
                plt.legend()
                plt.grid()
                if self.BU_interp:
                    plt.savefig(f"{self.save_dir}/rel_Delta_{iso}_{self.comparison_name}_interp.png")
                else:
                    plt.savefig(f"{self.save_dir}/rel_Delta_{iso}_{self.comparison_name}.png")
                plt.close()

            if plot_abs:
                plt.figure()
                for comparison_case in self.delta_Niso_abs[iso].keys():
                    plt.plot(self.D5_case.DRAGON_BU, self.delta_Niso_abs[iso][comparison_case], label = f"{comparison_case}".replace("_"," "), marker = "x", linestyle = "--")
                plt.xlabel(f"Burnup [MWd/tU]")
                plt.ylabel(f"$\\Delta$ N{iso} [atom/b-cm]")
                plt.title(f"$\\Delta$ N{iso} evolution for {self.comparison_name} case")
                plt.legend()
                plt.grid()
                if self.BU_interp:
                    plt.savefig(f"{self.save_dir}/abs_Delta_{iso}_{self.comparison_name}_interp.png")
                else:
                    plt.savefig(f"{self.save_dir}/abs_Delta_{iso}_{self.comparison_name}_abs.png")
                plt.close()
        return