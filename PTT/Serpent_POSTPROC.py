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
import serpentTools
from serpentTools.settings import rc
import numpy as np
import os



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
        comparison_name = string to describe the comparison to plot
        comparision_list = list with a sequence of reference case followed by its associated test case, 
        save_dir = directory where to save figures
        """
        self.print_name = comparision_name
        self.ref_cases = [case for case in comparision_list if case.case_is_ref]
        print(f"ref_cases = {self.ref_cases}")
        self.test_cases = [case for case in comparision_list if not case.case_is_ref]           
        print(f"TEST CASES = {self.test_cases}")
        print(f"len(TEST CASES) = {len(self.test_cases)}")
        self.savedir = save_dir
        for ref in self.ref_cases:
            self.listBU = ref.SERPENT_BU
        self.save_name = self.print_name.replace(" ","_")
        print(f"SAVE NAME is = {self.save_name}")
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
                print(f"{self.test_cases[i]}")
                self.error_iso_dens_for_tests.append(self.compute_error_isoDens(self.ref_cases[i], self.test_cases[i])) # creation of the error_iso_dens_for_tests attribute
            print("ERROR ISO DENS FOR TESTS IS :")
            print(self.error_iso_dens_for_tests)
            print(f"length of error list attribute is {len(self.error_iso_dens_for_tests)}")
            self.plot_error_isoDens()

    def plot_comparisons_keff(self): 
        print('$$$ -------- POSTPROC.py : Serpent2 figures (Keff) ')
        fig,ax = plt.subplots(dpi=250)
        for ref in self.ref_cases:
            ax.plot(self.listBU, ref.SERPENT_Keff, label=ref.case_name, marker = "D", linewidth=1, markersize=2, linestyle = "--")
        for test in self.test_cases:
            print(test)
            ax.plot(self.listBU, test.SERPENT_Keff, label=test.case_name, marker = "x", linewidth=1, markersize=2, linestyle = "--")
        ax.legend(loc="best")
        ax.set_xlabel("BU (MWj/t)")
        ax.set_ylabel("Keff Serpent2")
        ax.grid()
        ax.set_title(f"Comparison of Keffs for {self.print_name}", y=1.05)
        fig.savefig(self.savedir+f"Keffs_comparison_{self.save_name}")
        plt.close("all")
        return
    
    def compute_errors_Keff_to_ref(self, reference, test):
        tmp_err =[]
        for i in range(len(test.SERPENT_Keff)):
            print("In compute errors")
            print(test.case_name)
            tmp_err.append((test.SERPENT_Keff[i]-reference.SERPENT_Keff[i])*1e5)
        self.error_keffs_list.append(tmp_err)
        return

    def plot_errorK_to_ref(self):
        print('$$$ -------- POSTPROC.py : Serpent2 figures (Keff error) ')
        fig,ax = plt.subplots(dpi=250)
        for i in range(len(self.error_keffs_list)):
            print(f"i = {i}")
            print(f"len(self.error_keffs_list) = {len(self.error_keffs_list)}")
            print(f"len(self.test_cases) = {len(self.test_cases)}")
            ax.plot(self.listBU, self.error_keffs_list[i], label=f'{self.test_cases[i].case_name} : {self.test_cases[i].library_used}-{self.ref_cases[0].library_used}', marker = "x", linewidth=1, markersize=4, linestyle = "--")
        ax.legend(loc="best")
        ax.set_xlabel("BU (MWj/t)")
        ax.set_ylabel(f"Difference on Keff (pcm)")
        ax.grid()
        ax.set_title(f"{self.print_name} Keff difference PyNjoy2016-oldlib.", y=1.05)
        fig.savefig(self.savedir+f"Keff_error_{self.save_name}")
        plt.close("all")
        return
    
    
    def plot_comparisons_isotopeDens(self):
        print('$$$ -------- POSTPROC.py : Serpent2 figures (Iso dens) ')
        for iso in self.ref_cases[0].SERPENT_isotopes_dens.keys():
            fig,ax = plt.subplots(dpi=500)
            for ref in self.ref_cases:
                ax.plot(self.listBU, ref.SERPENT_isotopes_dens[iso], label=f'{ref.case_name} : {ref.library_used}', marker = "D", linewidth=1, markersize=5, linestyle = "--")
            for test in self.test_cases:
                ax.plot(self.listBU, test.SERPENT_isotopes_dens[iso], label=f'{test.case_name} : {test.library_used}', marker = "x", linewidth=1, markersize=5, linestyle = "--")
            ax.legend(loc="best")
            ax.set_xlabel("BU (MWj/t)")
            ax.set_ylabel(f"Isotope density of {iso} (a/b*cm)")
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
            print(f"number of errors associated to tests = {len(self.error_iso_dens_for_tests)}")
            print(f"errors array is = {self.error_iso_dens_for_tests}")
            for i in range(len(self.error_iso_dens_for_tests)):
                ax.plot(self.listBU, self.error_iso_dens_for_tests[i][iso], label=f'{self.test_cases[i].case_name} : {self.test_cases[i].library_used} - {self.ref_cases[0].library_used}', marker = "x", linewidth=1, markersize=4, linestyle = "--")
            ax.legend(loc="best")
            ax.set_xlabel("BU (MWj/t)")
            ax.set_ylabel(f"Relative difference on isotopic density (%)")
            ax.grid()
            ax.set_title(f"{self.print_name} rel. diff. on {iso}")
            fig.savefig(self.savedir+f"Error_{iso}_dens_{self.save_name}")
            plt.tight_layout()
            plt.close("all")
        return

# --- Main script
# --- Get current working directory
path=os.getcwd()

# --- Create the directory to save the results 
a=os.path.exists('Serpent2_PyNjoy2016_tests')
if a==False:
    os.mkdir('Serpent2_PyNjoy2016_tests')

SAVE_DIR='Serpent2_PyNjoy2016_tests/'

a=os.path.exists(SAVE_DIR)
if a==False:
    os.mkdir(SAVE_DIR)
form = "png"


cell = "HOM_CELL"
case = "HOM_U5"
path_to_data = f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/{cell}_study/{case}/"

pynjoy_test_file_name=f'{case}_PyNjoy2016_mc' # name of the tests' outputs
sss_test_file_name=f'{case}_oldlib_mc'
path_to_tests = path_to_data


isotopes_=['U235', 'U238', 'U234', 'Pu239', 'Pu240', 'Gd155','Gd157','Xe135','Sm149']

HOM_U5_PyNjoy = Serpent2_case(path_to_data, pynjoy_test_file_name, "HOM_U5", 750, "PyNjoy2016", isotopes_, SAVE_DIR)

HOM_U5_oldlib = Serpent2_case(path_to_tests, sss_test_file_name, "HOM_U5", 750, "oldlib", isotopes_, SAVE_DIR)
HOM_U5_oldlib.setisRef(True)

Serpent2_post_treatment("HOM_U5 PyNjoy2016-oldlib", [HOM_U5_PyNjoy, HOM_U5_oldlib], SAVE_DIR)
