# Python3 script for post treating Serpent2 test cases for JEFF311-PyNjoy2016 study
# Author : R. Guasch
# test 1 : removed all metastable isotopes that caused error ReadACEFile function
# test 2 : changed all metastable isotopes to fundamental states
# test 3 : changed all fundamentals to metastables
#####################################################
# Update : 21/08/2024
# Error found in acefiles for isotopes with metastable states : wrong name format, Serpent2 was not able to read them
# Purpose : Post-treatment of Serpent2 test cases for JEFF311-PyNjoy2016 study


import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import serpentTools
from serpentTools.settings import rc
import numpy as np
import os
import getBULists



class Serpent2_case:
    def __init__(self, path, file_name, case_name, temperature, library, isotopes, save_dir):
        self.path_to_data = path
        self.file_name = file_name
        self.case_name = case_name
        self.Temp = temperature
        self.library_used = library
        self.isotopes = isotopes
        self.savedir = save_dir

        self.case_is_ref = False # By default, the case is not a reference case

        self.load_Serpent_BU() # initialize SERPENT_BU attribute
        self.load_serpent2_keffs() # initialize SERPENT_Keff attribute
        self.load_Serpent_IsoDens() # initialize isotopic densities in evolutions attribute, stored in SERPENT_isotopes_dens dictionnary with key = isotope, values = list of atomic densities

    def setisRef(self, isRef): # use this method to set the case as a reference case
        self.case_is_ref = isRef
        return

    def load_Serpent_BU(self):   
        """
        Initialize SEREPENT_BU attribute = BU steps used for depletion calculations
        """
        depFile = self.path_to_data+self.file_name+'_dep.m'
        dep = serpentTools.read(depFile)
        fuel=dep.materials['total']
        serpent_BU=fuel.burnup
        np.savetxt(self.savedir+f'/serpent_BU_{self.file_name+"_"+str(self.Temp)}.txt',serpent_BU)
        self.SERPENT_BU=np.loadtxt(self.savedir+f'/serpent_BU_{self.file_name+"_"+str(self.Temp)}.txt',dtype=float)
        Ls1=np.shape(self.SERPENT_BU)
        lenBU_SERPENT=Ls1[0]
        for k in range(lenBU_SERPENT):
            self.SERPENT_BU[k]=1000*self.SERPENT_BU[k]
        return

    def load_serpent2_keffs(self):
        """
        Initialize SERPENT_Keff attribute
        """
        res=serpentTools.read(self.path_to_data+self.file_name+"_res.m")
        serpent_keff=res.resdata["absKeff"]
        np.savetxt(self.savedir+f'serpent_keff_{self.file_name+"_"+str(self.Temp)}.txt',serpent_keff)
        SERPENT_keff=np.loadtxt(self.savedir+f'serpent_keff_{self.file_name+"_"+str(self.Temp)}.txt',dtype=float)
        Ls2=np.shape(SERPENT_keff)
        lenBU_SERPENT=Ls2[0]
        self.SERPENT_Keff=np.zeros(lenBU_SERPENT)    
        for k in range(lenBU_SERPENT):
            self.SERPENT_Keff[k]=SERPENT_keff[k][0] # reshaping array to a 1D list
        return

    def load_Serpent_IsoDens(self):
        """
        Initialize isotopic evolution data
        """
        depFile = self.path_to_data+self.file_name+'_dep.m'
        dep = serpentTools.read(depFile)
        fuel=dep.materials['total']
        self.SERPENT_isotopes_dens = {}
        for isotope in self.isotopes:
            serpent_ISOTOPESDENS=fuel.toDataFrame("adens",names=isotope)
            np.savetxt(self.savedir+f'serpent_ISOTOPESDENS_{self.file_name+"_"+str(self.Temp)}_{isotope}.txt',serpent_ISOTOPESDENS)
            SERPENT_ISOTOPESDENS=np.loadtxt(self.savedir+f'serpent_ISOTOPESDENS_{self.file_name+"_"+str(self.Temp)}_{isotope}.txt',dtype=float)
            SERPENT_ISOTOPESDENS=np.transpose(SERPENT_ISOTOPESDENS)
            self.SERPENT_isotopes_dens[isotope] = np.array(SERPENT_ISOTOPESDENS)
        return



class Serpent2_post_treatment:
    def __init__(self, comparision_name, comparision_list, save_dir):
        """
        comparison_name = (str) to describe the comparison to plot
        comparision_list = (list) with a sequence of (Serpent2_case objects) to compare, 
        save_dir = (str) directory where to save figures

        Fit error on isotopic densities options: 
        BU_split_for_fit_dict = (dictionary) with (keys :) isotopes, (values :) lists of indices to split the BU list for fitting the error on sub-intervals
        error_models_to_fit_dict = (dictionary) with (keys :) isotopes, (values :) list of error models to fit as values, if empty, no fit is performed on the corresponding isotope/BU list
        """
        self.print_name = comparision_name.split("_")[0]
        self.ref_cases = [case for case in comparision_list if case.case_is_ref]
        self.test_cases = [case for case in comparision_list if not case.case_is_ref]           
        self.savedir = save_dir
        for ref in self.ref_cases:
            self.listBU = ref.SERPENT_BU
        self.save_name = comparision_name.replace(" ","_")
        # Compare test cases to reference case according to their temperature attribute
        # Generalize this to allow for the comparison of a test case with a reference case based on the .Temp attribute

        # Keff plots
        self.plot_comparisons_keff() # plot Keffs for reference case + all test cases given as input
        self.error_keffs_list = []
        if self.test_cases:
            for i in range(len(self.ref_cases)):
                self.compute_errors_Keff_to_ref(self.ref_cases[i], self.test_cases[i]) # creation of the error_keffs_list attribute
            self.plot_errorK_to_ref() # plot errors on Keff between reference case and test cases considered


        # Isotopic density plots
        self.plot_comparisons_isotopeDens() # plot isotopic densities for reference case + all test cases given as input
        self.error_iso_dens_for_tests = []
        if self.test_cases:
            for i in range(len(self.ref_cases)):
                self.error_iso_dens_for_tests.append(self.compute_error_isoDens(self.ref_cases[i], self.test_cases[i])) # creation of the error_iso_dens_for_tests attribute
            self.plot_error_isoDens()


    def fit_error_isoDens(self, BU_split_for_fit_dict, error_models_to_fit_dict):
        # Fit error on isotopic densities
        self.BUsplit_indx_dict = BU_split_for_fit_dict # subdivide the BU list for fitting the error on isotopic densities
        self.error_models_to_fit = error_models_to_fit_dict # associate to each iso a list of error models to fit to the error, can choose not to fit on specific BU sub-intervals, provided the BU_split_for_fit_dict is correctly set
        if self.test_cases:
            for iso in self.ref_cases[0].SERPENT_isotopes_dens.keys():
                self.fit_error_isoDens(iso)

    def plot_comparisons_keff(self): 
        print('$$$ -------- POSTPROC.py : Serpent2 figures (Keff) ')
        fig,ax = plt.subplots(dpi=250)
        for ref in self.ref_cases:
            ax.plot(self.listBU, ref.SERPENT_Keff, label=ref.case_name, marker = "D", linewidth=1, markersize=2, linestyle = "--")
        for test in self.test_cases:
            ax.plot(self.listBU, test.SERPENT_Keff, label=test.case_name, marker = "x", linewidth=1, markersize=2, linestyle = "--")
        ax.legend(loc="best")
        ax.set_xlabel("BU [MWj/t]")
        ax.set_ylabel("Keff Serpent2")
        ax.grid()
        ax.set_title(f"Comparison of Keffs for {self.print_name}", y=1.05)
        fig.savefig(self.savedir+f"Keffs_comparison_{self.save_name}")
        plt.close("all")
        return
    
    def compute_errors_Keff_to_ref(self, reference, test):
        tmp_err =[]
        for i in range(len(test.SERPENT_Keff)):
            tmp_err.append((test.SERPENT_Keff[i]-reference.SERPENT_Keff[i])*1e5)
        self.error_keffs_list.append(tmp_err)
        return

    def plot_errorK_to_ref(self):
        print('$$$ -------- POSTPROC.py : Serpent2 figures (Keff error) ')
        fig,ax = plt.subplots(dpi=250)
        for i in range(len(self.error_keffs_list)):
            ax.plot(self.listBU, self.error_keffs_list[i], label=f'{self.test_cases[i].case_name}', marker = "x", linewidth=1, markersize=4, linestyle = "--")
        ax.legend(loc="best")
        ax.set_xlabel("BU [MWj/t]")
        ax.set_ylabel(f"Difference on Keff [pcm]")
        ax.grid()
        ax.set_title(f"\u0394 Keff - {self.print_name} ", y=1.05)
        fig.savefig(self.savedir+f"Keff_error_{self.save_name}")
        plt.close("all")
        return
    
    
    def plot_comparisons_isotopeDens(self):
        print('$$$ -------- POSTPROC.py : Serpent2 figures (Iso dens) ')
        for iso in self.ref_cases[0].SERPENT_isotopes_dens.keys():
            fig,ax = plt.subplots(dpi=500)
            for ref in self.ref_cases:
                ax.plot(self.listBU, ref.SERPENT_isotopes_dens[iso], label=f'{ref.case_name} : {ref.library_used}', marker = "D", linewidth=1, markersize=2, linestyle = "--")
            for test in self.test_cases:
                ax.plot(self.listBU, test.SERPENT_isotopes_dens[iso], label=f'{test.case_name} : {test.library_used}', marker = "x", linewidth=1, markersize=2, linestyle = "--")
            ax.legend(loc="best")
            ax.set_xlabel("BU [MWj/t]")
            ax.set_ylabel(f"Isotope density of {iso} [a/b-cm]")
            ax.set_title(f"{iso} density for PyNjoy2016 tests, {self.print_name}")
            ax.grid()
            fig.savefig(self.savedir+f"{iso}_comparison_{self.save_name}")
            plt.tight_layout()
            plt.close("all")
        return

    def compute_error_isoDens(self, reference, test):
        # Temporarily suppress the warning for invalid division
        with np.errstate(divide='ignore', invalid='ignore'):
            error_per_isotope = {}
            for iso in reference.SERPENT_isotopes_dens.keys():
                ref_dens = reference.SERPENT_isotopes_dens[iso]
                test_dens = test.SERPENT_isotopes_dens[iso]
                relative_error = np.where(ref_dens != 0.0, (test_dens - ref_dens) * 100 / ref_dens, 0.0)
                error_per_isotope[iso] = relative_error
        return error_per_isotope
    
    def plot_error_isoDens(self):
        print('$$$ -------- POSTPROC.py : Serpent2 figures (Iso dens error) ')
        for iso in self.ref_cases[0].SERPENT_isotopes_dens.keys():
            fig,ax = plt.subplots(figsize=(8, 6))
            for i in range(len(self.error_iso_dens_for_tests)):
                ax.plot(self.listBU, self.error_iso_dens_for_tests[i][iso], label=f'{self.test_cases[i].case_name}', marker = "x", linewidth=1, markersize=4, linestyle = "--")
            ax.legend(loc="best")
            ax.set_xlabel("BU [MWj/t]")
            ax.set_ylabel(f"Relative difference on isotopic density [%]")
            ax.grid()
            ax.set_title(f"\u0394 {iso} - {self.print_name}")
            fig.savefig(self.savedir+f"Error_{iso}_dens_{self.save_name}")
            plt.tight_layout()
            plt.close("all")
        return
    def linear_model(self, BU, a, b):
        """
        Linear model for error in isotopic densities evolution in function of BU
        """
        return a*BU + b 

    def proportional_model(self, BU, a):
        """
        Proportional model for error in isotopic densities evolution in function of BU
        """
        return a*BU
    
    def constant_model(self, BU, a):
        """
        constant model for error in isotopic densities evolution in function of BU
        """
        return a*np.ones(len(BU))  

    def exponential_model(self, BU, a, b):
        """
        Exponential model for isotopic densities evolution in function of BU
        """
        return a*np.exp(b*BU)
    
    def exponenatial_linear_model(self, BU, a, b, c):
        """
        Expo-linear model for isotopic densities evolution in function of BU
        """
        return a*np.exp(b*BU)+c # test quality / applicability of this model
    
    def combination_of_exponentials_model(self, BU, a, b, c, d):
        """
        Linear combination of exponential model for isotopic densities evolution in function of BU
        """
        return a*np.exp(b*BU)+c*np.exp(d*BU) # expect a linear combination of rising/decresing exponentials for isotopic densities evolution in function of BU

    def combination_of_exponentials_model_with_constant(self, BU, a, b, c, d, e):
        """
        Linear combination of exponential model with added constant for isotopic densities evolution in function of BU
        """
        return a*np.exp(b*BU)+c*np.exp(d*BU)+e #check quality / applicability of this model

    def fit_error_isoDens(self, iso):
        """
        This function fits the error on isotopic densities to a polynomial function
        Purpose : support analysis and visualization of the error on isotopic densities
        Possibility to exclude the first n points of the BU list to fit the error on isotopic densities:
          - Allows for fitting points with residual=0 at the origin in the case of a proportional model,
          - Allows to exclude the initial transient of the isotopic densities evolution and get a better fit of the error after the first n points
        """
        print(f"$$$ --- Fitting error on isotopic densities for case {self.print_name} and isotope {iso}")
        fig,ax = plt.subplots(figsize=(8, 6))
        for i in range(len(self.error_iso_dens_for_tests)):
            BU_intervals = getBULists.get_SplitLists(self.listBU, self.BUsplit_indx_dict[iso])
            for j in range(len(BU_intervals)):
                BU_sublist = BU_intervals[j]
                print(f"error_model_to_fit dict = {self.error_models_to_fit}")
                print(f"error_model_to_fit dict [iso] = {self.error_models_to_fit[iso]}")
                model_to_fit = self.error_models_to_fit[iso][j] # associate a model to a sub-interval of the BU list
                print(f"$$$ --- Fitting error on isotopic densities for case {self.print_name} and isotope {iso} with model {model_to_fit}")
                print(f"$$$ --- BU sublist for fitting : {BU_sublist}")
                print(f"self.BUsplit_indx_dict[iso] = {self.BUsplit_indx_dict[iso]}")
                if self.BUsplit_indx_dict[iso]:
                    if len(self.BUsplit_indx_dict[iso]) > 1:
                        first_indx = self.BUsplit_indx_dict[iso][j]
                        last_indx = self.BUsplit_indx_dict[iso][j+1]
                    else:
                        if j == 0:
                            first_indx = 0
                            print(f"iso = {iso}, i = {i}, self.BUsplit_indx_dict[iso] = {self.BUsplit_indx_dict[iso]}")
                            last_indx = self.BUsplit_indx_dict[iso][0]
                        else:
                            first_indx = self.BUsplit_indx_dict[iso][0]
                            last_indx = len(self.listBU)
                else:
                    first_indx = 0
                    last_indx = len(self.listBU)
                if model_to_fit == "linear":
                    # Linear model
                    print(f"first index, last index = {first_indx}, {last_indx}")
                    print(f"length of BU sublist = {len(BU_sublist)}")
                    print(f"Length of isotopic densities = {len(self.error_iso_dens_for_tests[i][iso])}")
                    print(f"Length of subsets of errors = {len(self.error_iso_dens_for_tests[i][iso][first_indx:last_indx])}")
                    popt, pcov = curve_fit(self.linear_model, BU_sublist, self.error_iso_dens_for_tests[i][iso][first_indx:last_indx])
                    linear_fit_data = self.linear_model(BU_sublist, *popt)
                    chi_squared, reduced_chi_squared, r_squared = self.compute_metrics(linear_fit_data, self.error_iso_dens_for_tests[i][iso][first_indx:last_indx], popt)
                    print(f"$$ -- Linear fit : popt = {popt}, condition number = {np.linalg.cond(pcov)}")
                    print(f"$ - Quality of fit : chi_squared = {chi_squared}, reduced_chi_squared = {reduced_chi_squared}, r_squared = {r_squared}")
                    ax.plot(BU_sublist, linear_fit_data, label=f'Linear fit : error = {popt[0]:.3E}*BU+{popt[1]:.3E}', marker = ".", linewidth=0.5, linestyle = "solid")
                if model_to_fit == "proportional":
                    # Proportional model
                    popt, pcov = curve_fit(self.proportional_model, BU_sublist, self.error_iso_dens_for_tests[i][iso][first_indx:last_indx])
                    proportional_fit_data = self.proportional_model(BU_sublist, *popt)
                    chi_squared, reduced_chi_squared, r_squared = self.compute_metrics(proportional_fit_data, self.error_iso_dens_for_tests[i][iso][first_indx:last_indx], popt, print_info=True)
                    print(f"$$ -- Proportional fit : popt = {popt}, condition number = {np.linalg.cond(pcov)}")
                    print(f"$ - Quality of fit : chi_squared = {chi_squared}, reduced_chi_squared = {reduced_chi_squared}, r_squared = {r_squared}")
                    ax.plot(BU_sublist, proportional_fit_data, label=f'Proportional fit : error = {popt[0]:.3E}*BU', marker = ".", linewidth=0.5, linestyle = "solid")
                if model_to_fit == "constant":
                    # Constant model
                    popt, pcov = curve_fit(self.constant_model, BU_sublist, self.error_iso_dens_for_tests[i][iso][first_indx:last_indx])
                    constant_fit_data = self.constant_model(BU_sublist, *popt)
                    chi_squared, reduced_chi_squared, r_squared = self.compute_metrics(constant_fit_data, self.error_iso_dens_for_tests[i][iso][first_indx:last_indx], popt)
                    print(f"$$ -- Constant fit : popt = {popt}, condition number = {np.linalg.cond(pcov)}")
                    print(f"$ - Quality of fit : chi_squared = {chi_squared}, reduced_chi_squared = {reduced_chi_squared}, r_squared = {r_squared}")
                    ax.plot(BU_sublist, constant_fit_data, label=f'Constant fit : error = {popt[0]:.3E}', marker = ".", linewidth=0.5, linestyle = "solid")

            ax.plot(self.listBU, self.error_iso_dens_for_tests[i][iso], label=f'{self.test_cases[i].case_name} : {self.test_cases[i].library_used} - {self.ref_cases[0].library_used}', marker = "D", linewidth=1, markersize=4, linestyle = "--", color="red")
        ax.legend(loc="best")
        ax.set_xlabel("BU [MWj/t]")
        ax.set_ylabel(f"Relative difference on isotopic density [%]")
        ax.grid()
        ax.set_title(f"{self.print_name} - \u0394 {iso}")
        fig.savefig(self.savedir+f"Fit_Error_{iso}_dens_{self.save_name}")
        plt.tight_layout()
        plt.close("all")

        return
    
    def compute_metrics(self, fitted_data, ref_data, params, print_info = False):
        """
        Compute metrics to assess the quality of the fit
        """
        residuals = fitted_data - ref_data
        if print_info:
            print(f"residuals = {residuals}")
            print(f"fitted_data = {fitted_data}")
        chi_squared = np.sum((residuals ** 2) / fitted_data)
        reduced_chi_squared = chi_squared / (len(ref_data) - len(params))
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((ref_data - np.mean(ref_data))**2)
        r_squared = 1 - (ss_res / ss_tot)

        return chi_squared, reduced_chi_squared, r_squared
    
    def fit_isoDens(self):
        """
        This function fits the isotopic densities to an exponential function
        Purpose : support analysis and visualization of the isotopic densities
        Hopefully : analyse some of the physics behind the isotopic densities evolution and help identify differences between the reference and test cases
        """
        for iso in self.ref_cases[0].SERPENT_isotopes_dens.keys():
            fig,ax = plt.subplots(figsize=(8, 6))
            for i in range(len(self.test_cases)):
                # Exponential model
                popt, pcov = curve_fit(self.exponential_model, self.listBU, self.test_cases[i].SERPENT_isotopes_dens[iso])
                print(f"for isotopic density {iso} and test case {self.test_cases[i].case_name}, Exponential fit : popt = {popt}, condition number = {np.linalg.cond(pcov)}")
                ax.plot(self.listBU, self.exponential_model(self.listBU, *popt), label=f'Exponential fit : dens = {popt[0]}*exp({popt[1]}*BU)', marker = ".", linewidth=1, linestyle = "--")
                ax.plot(self.listBU, self.test_cases[i].SERPENT_isotopes_dens[iso], label=f'{self.test_cases[i].case_name} : {self.test_cases[i].library_used}', marker = "D", linewidth=1, markersize=4, linestyle = "--")
            ax.legend(loc="best")
            ax.set_xlabel("BU [MWj/t]")
            ax.set_ylabel(f"Isotope density of {iso} [a/b-cm]")
            ax.grid()
            ax.set_title(f"{iso} density for {self.print_name}")
            fig.savefig(self.savedir+f"Fit_{iso}_dens_{self.save_name}")
            plt.tight_layout()
            plt.close("all")
        return

# --- Main script

noGd_tests = []
Gd_tests = []

cell = "HOM_CELL"
all_cases = ["HOM_U5", "HOM_U5_U8", "HOM_UOX", "HOM_UOX_clad_noZr", "HOM_UOX_clad", "HOM_UOXGd_no155157", "HOM_UOX_Gd155", "HOM_UOX_Gd157", "HOM_UOXGd"]
cases_no_Gd = ["HOM_U5", "HOM_U5_U8", "HOM_UOX", "HOM_UOX_clad_noZr", "HOM_UOX_clad"]
cases_Gd = ["HOM_UOXGd_no155157", "HOM_UOX_Gd155", "HOM_UOX_Gd157", "HOM_UOXGd"]

for case in all_cases:

    # --- Get current working directory
    path=os.getcwd()

    # --- Create the directory to save the results 
    a=os.path.exists(f'Serpent2_PyNjoy2016_tests/{case}/')
    if a==False:
        os.makedirs(f'Serpent2_PyNjoy2016_tests/{case}/')

    SAVE_DIR=f'Serpent2_PyNjoy2016_tests/{case}/'

    a=os.path.exists(SAVE_DIR)
    if a==False:
        os.mkdir(SAVE_DIR)
    form = "png"

    path_to_data = f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/{cell}_study/{case}/"

    pynjoy_test_file_name=f'{case}_PyNjoy2016_mc' # name of the tests' outputs
    sss_test_file_name=f'{case}_oldlib_mc'
    path_to_tests = path_to_data


    isotopes_=['U235', 'U238', 'U234', 'Pu239', 'Pu240', 'Gd155','Gd157','Xe135','Sm149']

    HOM_PyNjoy = Serpent2_case(path_to_data, pynjoy_test_file_name, case, 750, "PyNjoy2016", isotopes_, SAVE_DIR)


    HOM_oldlib = Serpent2_case(path_to_tests, sss_test_file_name, case, 750, "oldlib", isotopes_, SAVE_DIR)
    HOM_oldlib.setisRef(True)
    if case in cases_no_Gd:
        noGd_tests.append(HOM_PyNjoy)
        noGd_tests.append(HOM_oldlib)
    elif case in cases_Gd:
        Gd_tests.append(HOM_PyNjoy)
        Gd_tests.append(HOM_oldlib)

    # --- Post-treatment of the Serpent2 test cases
    # --- Comparison of the two cases
    # -- Options for error fitting
    BU_list = getBULists.getLists("BOC_fine_autop5")[0]
    split_indxs_Gd = [len(BU_list)-3]
    split_indixs_U8_U4 = [10]
    split_lists_Gd = getBULists.get_SplitLists(BU_list, split_indxs_Gd)
    print(f"$$$ --- Split lists for Gd : {split_lists_Gd}")

    """
    Fit error on isotopic densities options: 
            BU_split_for_fit_dict = (dictionary) with (keys :) isotopes, (values :) lists of indices to split the BU list for fitting the error on sub-intervals
            error_models_to_fit_dict = (dictionary) with (keys :) isotopes, (values :) list of error models to fit as values, if empty, no fit is performed on the corresponding isotope/BU list
    BU_split_for_fit_dict, error_models_to_fit,
    """
    BU_split_for_fit_dict = {"U235": [10], "U238":split_indixs_U8_U4, "U234":split_indixs_U8_U4, "Gd155": split_indxs_Gd, "Gd157": split_indxs_Gd, "Pu239": [10], "Pu240": [10], "Xe135": [10], "Sm149": []}
    error_models_to_fit_dict = {"U235": ["none","constant"], "Gd155": ["linear","linear"], "Gd157": ["linear","linear"], "U238": ["linear","linear"], "U234": ["linear","linear"], "Pu239": ["none","constant"], "Pu240": ["none","constant"], "Xe135": ["none","constant"], "Sm149": ["proportional"]}
    
SAVE_DIR=f'Serpent2_PyNjoy2016_tests/PyNjoyvsOldlib/'

a=os.path.exists(SAVE_DIR)
if a==False:
    os.mkdir(SAVE_DIR)
form = "png"    
Serpent2_post_treatment(f"PyNjoy2016-oldlib_noGd", noGd_tests, SAVE_DIR)
Serpent2_post_treatment(f"PyNjoy2016-oldlib_Gd", Gd_tests, SAVE_DIR)