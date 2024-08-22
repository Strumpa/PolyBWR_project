# Python3 script used to retrive and analyse time + memory consumption of DRAGON5 modules.
# The aim is obtain the cost / precision optimum considering all different tracking modules and solution doors available in DRAGON5.
# The first topic of investigation is self-shelding methods and solution doors to the slowing down equation.

# Procedures called by SSH_tests.x2m are :
#       - SALSAL.c2m
#       - SALNXT.c2m
#       - NXTSAL.c2m
#       - SYBSYB.c2m
#       - SYBNXT.c2m
#       - SYBSAL.c2m
#       - NXTNXT.c2m
#  Every time, the first tracking module mentionned corresponds to the self-shielding tracking module.
#  When NXT or SALT are used as SSH tracking modules, the MERGMIX option is activated.
#  The latter give the opportunity of using an additional call to MCCGT to use a MOC solution door to the self-shielding.
#  Every time the procedure called treats UOX and UOX+Gd benchmark cells, comparing RSE and PT methods, 
#       - In the case of NXT or SALT ssh trackings, compare CP and MOC solution door to the SSH method,  (an2d_ss = 8, densur_ss = 15.0), (an2d = 18, densur=35.0)
#       - In the case of SYBILT use DP01 method on pincell, (iqua2_ss = 24, nseg_ss = 5) and (iqua2 = 30, nseg = 10).

import matplotlib.pyplot as plt
import numpy as np    

def parse_output(path_to_file):
    print("in Parse output")
    tests = []
    with open(path_to_file, 'r') as file:
            parsing_new_test = False
            for line in file:
                line = line.strip()  # Remove leading/trailing whitespace
                if "PROCEDURE FOR OPTIONS" in line:
                    parsing_new_test = True
                    procedure = line.split(" ")[1].split(".")[0]
                    benchmark_cell = line.split(" ")[5]
                    self_shielding_method = line.split(" ")[7]
                    if procedure not in ["SYBSAL", "SYBNXT", "SYBSYB"]:
                        self_shielding_solution_door = line.split(" ")[8]
                    else:
                        self_shielding_solution_door = "IC"
                    test = ssh_test(procedure, benchmark_cell, self_shielding_method, self_shielding_solution_door)
                if parsing_new_test and "TIME SPENT=" in line:
                    line = line.replace("        ", " ")
                    line = line.replace("       ", " ")
                    line = line.replace("      ", " ")
                    line = line.replace("     ", " ")
                    line = line.replace("    ", " ")
                    line = line.replace("  ", " ")
                    print(line)
                    print(line.split(" "))
                    module = line.split(" ")[1]
                    time_spent = float(line.split(" ")[5])
                    memory_used = float(line.split(" ")[8])
                    test.set_time_by_module(module, time_spent)
                    test.set_memory_by_module(module, memory_used)
                if f"MOC Kinf=" in line:
                    Kinf = float(line.split("MOC Kinf=  ")[-1].split("|")[0])
                    test.set_Kinf(Kinf)
                    parsing_new_test = False
                    tests.append(test)
                elif "PIJ Kinf=" in line:
                    Kinf = float(line.split("PIJ Kinf=  ")[-1].split("|")[0])
                    test.set_Kinf(Kinf)
                    parsing_new_test = False
                    tests.append(test)
    return tests

def analyse_tests(ssh_tests, Kinf_ref_UOX, Kinf_ref_UOX_Gd):
    precision_UOX = []
    time_UOX = []
    precision_UOX_Gd = []
    time_UOX_Gd = []

    SALSAL_UOX = {}
    precision_SALSAL_UOX = []
    SALNXT_UOX = {}
    precision_SALNXT_UOX = []
    NXTSAL_UOX = {}
    precision_NXTSAL_UOX = []
    SYBSYB_UOX = {}
    precision_SYBSYB_UOX = []
    SYBSAL_UOX = {}
    precision_SYBSAL_UOX = []
    SYBNXT_UOX = {}
    precision_SYBNXT_UOX = []
    NXTNXT_UOX = {}
    precision_NXTNXT_UOX = []

    
    SALSAL_Gd = {}
    precision_SALSAL_Gd = []
    SALNXT_Gd = {}
    precision_SALNXT_Gd = []
    NXTSAL_Gd = {}
    precision_NXTSAL_Gd = []
    SYBSYB_Gd = {}
    precision_SYBSYB_Gd = []
    SYBSAL_Gd = {}
    precision_SYBSAL_Gd = []
    SYBNXT_Gd = {}
    precision_SYBNXT_Gd = []
    NXTNXT_Gd = {}
    precision_NXTNXT_Gd = []

    for test in ssh_tests:
        test.total_time = sum(test.time_by_module.values())
        test.total_memory = sum(test.memory_by_module.values())
        test.filter_times()  
        if test.cell == "AT10_UOX":
            test.Error_on_Kinf = (test.Kinf - Kinf_ref_UOX)*1e5 # error on eigenvalue in pcm
            precision_UOX.append(test.Error_on_Kinf)
            time_UOX.append(test.total_time)
            if test.proc == "SALSAL":
                SALSAL_UOX[test.data] = test.time_by_module
                precision_SALSAL_UOX.append(test.Error_on_Kinf)
            elif test.proc == "SALNXT":
                SALNXT_UOX[test.data] = test.time_by_module
                precision_SALNXT_UOX.append(test.Error_on_Kinf)
            elif test.proc == "NXTSAL":
                NXTSAL_UOX[test.data] = test.time_by_module
                precision_NXTSAL_UOX.append(test.Error_on_Kinf)
            elif test.proc == "SYBSYB":
                SYBSYB_UOX[test.data] = test.time_by_module
                precision_SYBSYB_UOX.append(test.Error_on_Kinf)
            elif test.proc == "SYBNXT":
                SYBNXT_UOX[test.data] = test.time_by_module
                precision_SYBNXT_UOX.append(test.Error_on_Kinf)   
            elif test.proc == "SYBSAL":
                SYBSAL_UOX[test.data] = test.time_by_module
                precision_SYBSAL_UOX.append(test.Error_on_Kinf)
            elif test.proc == "NXTNXT":
                NXTNXT_UOX[test.data] = test.time_by_module  
                precision_NXTNXT_UOX.append(test.Error_on_Kinf)                            
        elif test.cell == "AT10_UOX_Gd2O3":
            test.Error_on_Kinf = (test.Kinf - Kinf_ref_UOX_Gd)*1e5 # error on eigenvalu in pc
            precision_UOX_Gd.append(test.Error_on_Kinf)
            time_UOX_Gd.append(test.total_time)
            if test.proc == "SALSAL":
                SALSAL_Gd[test.data] = test.time_by_module
                precision_SALSAL_Gd.append(test.Error_on_Kinf)
            elif test.proc == "SALNXT":
                SALNXT_Gd[test.data] = test.time_by_module
                precision_SALNXT_Gd.append(test.Error_on_Kinf)
            elif test.proc == "NXTSAL":
                NXTSAL_Gd[test.data] = test.time_by_module
                precision_NXTSAL_Gd.append(test.Error_on_Kinf)
            elif test.proc == "SYBSYB":
                SYBSYB_Gd[test.data] = test.time_by_module
                precision_SYBSYB_Gd.append(test.Error_on_Kinf)
            elif test.proc == "SYBNXT":
                SYBNXT_Gd[test.data] = test.time_by_module
                precision_SYBNXT_Gd.append(test.Error_on_Kinf)
            elif test.proc == "SYBSAL":
                SYBSAL_Gd[test.data] = test.time_by_module
                precision_SYBSAL_Gd.append(test.Error_on_Kinf)
            elif test.proc == "NXTNXT":
                NXTNXT_Gd[test.data] = test.time_by_module
                precision_NXTNXT_Gd.append(test.Error_on_Kinf)
    plot_precision_vs_time(time_UOX, precision_UOX, "UOX")
    plot_precision_vs_time(time_UOX_Gd, precision_UOX_Gd, "UOX_Gd")
    #for test in SALSAL:
    plot_module_times(SALSAL_UOX, "SALSAL_UOX_Cell", precision_SALSAL_UOX)
    plot_module_times(SALNXT_UOX, "SALNXT_UOX_Cell", precision_SALNXT_UOX)
    plot_module_times(NXTSAL_UOX, "NXTSAL_UOX_Cell", precision_NXTSAL_UOX)
    plot_module_times(SYBSYB_UOX, "SYBSYB_UOX_Cell", precision_SYBSYB_UOX)
    plot_module_times(SYBSAL_UOX, "SYBSAL_UOX_Cell", precision_SYBSAL_UOX)
    plot_module_times(SYBNXT_UOX, "SYBNXT_UOX_Cell", precision_SYBNXT_UOX)
    plot_module_times(NXTNXT_UOX, "NXTNXT_UOX_Cell", precision_NXTNXT_UOX)

    plot_module_times(SALSAL_Gd, "SALSAL_UOX_Gd_Cell", precision_SALSAL_Gd)
    plot_module_times(SALNXT_Gd, "SALNXT_UOX_Gd_Cell", precision_SALNXT_Gd)
    plot_module_times(NXTSAL_Gd, "NXTSAL_UOX_Gd_Cell", precision_NXTSAL_Gd)
    plot_module_times(SYBSYB_Gd, "SYBSYB_UOX_Gd_Cell", precision_SYBSYB_Gd)
    plot_module_times(SYBSAL_Gd, "SYBSAL_UOX_Gd_Cell", precision_SYBSAL_Gd)
    plot_module_times(SYBNXT_Gd, "SYBNXT_UOX_Gd_Cell", precision_SYBNXT_Gd)
    plot_module_times(NXTNXT_Gd, "NXTNXT_UOX_Gd_Cell", precision_NXTNXT_Gd)
    

    return




def plot_module_times(test_cases, test_name, precision_):
    """
    Plots histograms representing the time spent in each module for different test cases.
    
    Parameters:
    test_cases (dict): A dictionary where keys are test case names and values are dictionaries
                       with module names as keys and time spent as values.
    """
    # Set the width of each bar
    bar_width = 0.2

    # Extract module names
    modules = sorted(set(module for times in test_cases.values() for module in times))
    
    # Number of test cases
    num_test_cases = len(test_cases)

    # Generate a color for each test case
    colors = plt.cm.viridis(np.linspace(0, 1, num_test_cases))

    # Create a figure and a set of subplots
    fig, ax = plt.subplots()

    # Bar positions
    indices = np.arange(len(modules))

    for i, (test_case, times) in enumerate(test_cases.items()):
        module_times = [times.get(module, 0) for module in modules]
        ax.bar(indices + i * bar_width, module_times, bar_width, label=test_case+f" error = {precision_[i]: .0f} pcm", color=colors[i])

    # Adding labels
    ax.set_xlabel('Modules')
    ax.set_ylabel('Time Spent (s)')
    ax.set_title(f'Time Spent in Each Module for {test_name}')
    ax.set_xticks(indices + bar_width * (num_test_cases - 1) / 2)
    ax.set_xticklabels(modules, rotation=45, ha='right')
    ax.legend()

    plt.tight_layout()
    plt.savefig(f'Module_Times_{test_name}.png')
    plt.show()

def plot_histogram_time_by_module(ssh_tests):
    # plot histogram of time spent by module
    # each module is a bar, each test is a color
    # x axis : modules
    # y axis : time spent
    times = [test.times for test in ssh_tests]
    bins = len(ssh_tests[0].modules)
    print(times)
    print(bins)
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'crimson']
    fig, ax = plt.subplots()
    ax.hist(times, bins, color = colors, label = [test.data for test in ssh_tests], density=True, histtype='bar')
    ax.legend()
    
    plt.xlabel('Module')
    plt.ylabel(f'Time spent (s)')
    plt.title(f"Time spent by module for {ssh_tests[0].proc}")
    plt.savefig(f"test_histogram_{ssh_tests[0].proc}.png")
    plt.show()
    return




def plot_precision_vs_time(time, precision, bechmark_cell):
    plt.scatter(time, precision)
    plt.ylabel('Error on Kinf (pcm)')
    plt.xlabel('Total time spent (s)')
    plt.title(f'Error on Kinf vs Total time spent for benchmark cell {bechmark_cell}')
    plt.savefig(f"precision_time_{bechmark_cell}.png")
    plt.show()
    plt.close()
    

class ssh_test:
    def __init__(self, procedure, cell, ssh_opt, ssh_sol):
        self.proc = procedure
        self.cell = cell
        self.ssh_opt = ssh_opt
        self.ssh_sol = ssh_sol

        if self.proc == "SYBSYB":
            print(f"proc is {self.proc}")
            self.ssh_trk = "SYBILT"
            self.flx_trk = "SYBILT"
        elif self.proc == "SYBSAL":
            self.ssh_trk = "SYBILT"
            self.flx_trk = "SALT"
        elif self.proc == "SYBNXT":
            self.ssh_trk = "SYBILT"
            self.flx_trk = "NXT"
        elif self.proc == "NXTNXT":
            self.ssh_trk = "NXT"
            self.flx_trk = "NXT" 
        elif self.proc == "SALNXT":
            self.ssh_trk = "SALT"
            self.flx_trk = "NXT"       
        elif self.proc == "NXTSAL":
            self.ssh_trk = "NXT"
            self.flx_trk = "SALT"     

        self.data = f"{self.ssh_opt}, {self.ssh_sol}"
        if self.proc == "SYBSYB":
            print(f"data is = {self.data}")
        self.time_by_module = {}   
        self.memory_by_module = {}
        return
    
    def set_time_by_module(self, module, time_spent):
        if module in self.time_by_module.keys():
            self.time_by_module[module] = self.time_by_module[module]+time_spent
            print(f"module {module} was already in time_by_module dictionnary !")
        else:
            self.time_by_module[module] = time_spent
            print(f"adding a pair to time_by_module : module = {module}, time = {time_spent}")
        return

    def filter_times(self, threshold=1):
        tmp_time_by_module = self.time_by_module.copy()
        for module in tmp_time_by_module.keys():
            if tmp_time_by_module[module] < threshold:
                del self.time_by_module[module]    

    def set_memory_by_module(self, module, memory_used):
        if module in self.memory_by_module.keys():
            self.memory_by_module[module] = self.memory_by_module[module]+memory_used
            print(f"module {module} was already in memory_by_module dictionnary !")
        else:
            self.memory_by_module[module] = memory_used
            print(f"adding a pair to memory_by_module : module = {module}, memory = {memory_used}")
        return
    
    def set_Kinf(self, eigenvalue):
        self.Kinf = eigenvalue
        print(f"Adding kinf= {self.Kinf} to test {self.proc}, {self.cell}, {self.ssh_opt}, {self.ssh_sol}")
        return
    
    def plot_time_consumption(self):
        print(f"Test with {self.Kinf} to test {self.proc}, {self.cell}, {self.ssh_opt}, {self.ssh_sol}")
        #print(self.time_by_module)
        # plot histogram of time spent by module

        plt.bar(self.time_by_module.keys(), self.time_by_module.values())
        plt.xlabel('Module')
        plt.ylabel(f'Time spent (s)')
        plt.title(f"Time spent by module for {self.proc}, {self.cell}, {self.ssh_opt}, {self.ssh_sol}")
        plt.show()
        plt.savefig(f"test_histogram_{self.proc}_{self.ssh_opt}_{self.ssh_sol}.png")
        plt.close()

    def plot_memory_consumption(self):
        print(f"Test with {self.Kinf} to test {self.proc}, {self.cell}, {self.ssh_opt}, {self.ssh_sol}")
        print(self.memory_by_module)




path = "../Version5_wc/Dragon/Linux_aarch64/SSH_tests.result"
#path = "/tmp/rundir36/SSH_tests.result"
ssh_tests = parse_output(path)
for test in ssh_tests:
    print(f"Test with {test.Kinf} to test {test.proc}, {test.cell}, {test.ssh_opt}, {test.ssh_sol}")

analyse_tests(ssh_tests, 1.24378E+00, 4.29394E-01)