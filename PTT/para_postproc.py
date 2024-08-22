# Python3 script to post-process the results obtained for Dragon5 parallelization
# Author: R. Guasch
# Date : 11 July 2024

import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Function to read the output file of the parallelization test


###### Parallelization test

def parse_output(path_to_file, omp_threads):
    print("in Parse output")
    tests = []
    procedure = path_to_file.split("/")[-1].split(".")[0]
    if "3x3" in procedure:
        benchmark_case = "AT10_3x3"
    elif "2x2" in procedure:
        benchmark_case = "AT10_2x2"
    elif "ASSBLY" in procedure:
        benchmark_case = "AT10_ASSBLY"
    else:
        benchmark_case = "AT10_cell"
    with open(path_to_file, 'r') as file:
            parsing_new_test = False
            for line in file:
                line = line.strip()  # Remove leading/trailing whitespace
                if f"Perform {benchmark_case} calculation" in line and "ECHO" not in line:
                    parsing_new_test = True
                    print(line)
                    print(line.split(" "))
                    #procedure = line.split(" ")[1].split(".")[0]
                    #benchmark_cell = line.split(" ")[5]
                    #self_shielding_method = line.split(" ")[7]
                    #if procedure not in ["SYBSAL", "SYBNXT", "SYBSYB"]:
                    #    self_shielding_solution_door = line.split(" ")[8]
                    #else:
                    #    self_shielding_solution_door = "IC"
                    ssh_opt = line.split(" ")[12].split("/")[0]
                    ssh_sol = line.split(" ")[12].split("/")[1].split(",")[0]
                    flx_sol = line.split(" ")[13]
                    batch_ssh = int(line.split(" ")[6])
                    batch_flx = int(line.split(" ")[10])
                    test = para_test(procedure, omp_threads, benchmark_case, ssh_opt, ssh_sol, flx_sol, batch_ssh, batch_flx)
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
                if f"MOC on 3x3_UOX_Gd: Kinf=" in line:
                    Kinf = float(line.split("Kinf=  ")[-1].split("|")[0])
                    test.set_Kinf(Kinf)
                    parsing_new_test = False
                    tests.append(test)
                elif "CP on 3x3_UOX_Gd: Kinf=" in line:
                    Kinf = float(line.split("Kinf=  ")[-1].split("|")[0])
                    test.set_Kinf(Kinf)
                    parsing_new_test = False
                    tests.append(test)

    return tests

def analyse_tests(ssh_tests, Kinf_ref):
    precision = []
    time = []

    RSEMOCMOC = {}
    RSEMOCMOC_list = []
    precision_RSEMOCMOC = []
    RSECPMOC = {}
    RSECPMOC_list = []
    precision_RSECPMOC = []
    RSECPCP = {}
    RSECPCP_list = []
    precision_RSECPCP = []

    min_MOCMOC_ssh, min_MOCMOC_flx = 1e6, 1e6
    min_CPCP_ssh, min_CPCP_flx = 1e6, 1e6
    min_CPMOC_ssh, min_CPMOC_flx = 1e6, 1e6
    for test in ssh_tests:
        test.total_time = sum(test.time_by_module.values())
        #test.set_time_by_module("total", test.total_time)
        test.total_memory = sum(test.memory_by_module.values())
        if "SALT:" in test.time_by_module.keys():
            test.tracking_module = "SALT"
        elif "NXT:" in test.time_by_module.keys():
            test.tracking_module = "NXT"
        elif "SYBILT:" in test.time_by_module.keys():
            test.tracking_module = "SYBILT"
        test.filter_times()  
        
        test.Error_on_Kinf = (test.Kinf - Kinf_ref)*1e5 # error on eigenvalue in pcm
        precision.append(test.Error_on_Kinf)
        time.append(test.total_time)

        if test.ssh_opt == "RSE" and test.ssh_sol == "MOC" and test.flx_sol == "MOC":
            RSEMOCMOC[test.data] = test.time_by_module
            RSEMOCMOC_list.append(test)
            if test.time_by_module["USS:"] < min_MOCMOC_ssh:
                min_MOCMOC_ssh = test.time_by_module["USS:"]
                best_MOCMOC_ssh = test
            if test.time_by_module["FLU:"]+test.time_by_module["ASM:"] < min_MOCMOC_flx:
                min_MOCMOC_flx = test.time_by_module["FLU:"]+test.time_by_module["ASM:"]
                best_MOCMOC_flx = test
            precision_RSEMOCMOC.append(test.Error_on_Kinf)
        elif test.ssh_opt == "RSE" and test.ssh_sol == "CP" and test.flx_sol == "MOC":
            RSECPMOC[test.data] = test.time_by_module
            RSECPMOC_list.append(test)
            if test.time_by_module["USS:"] < min_CPMOC_ssh:
                min_CPMOC_ssh = test.time_by_module["USS:"]
                best_CPMOC_ssh = test
            if test.time_by_module["FLU:"]+test.time_by_module["ASM:"] < min_CPMOC_flx:
                min_CPMOC_flx = test.time_by_module["FLU:"]+test.time_by_module["ASM:"]
                best_CPMOC_flx = test
            precision_RSECPMOC.append(test.Error_on_Kinf)
        elif test.ssh_opt == "RSE" and test.ssh_sol == "CP" and test.flx_sol == "CP":
            RSECPCP[test.data] = test.time_by_module
            RSECPCP_list.append(test)
            if test.time_by_module["USS:"] < min_CPCP_ssh:
                min_CPCP_ssh = test.time_by_module["USS:"]
                best_CPCP_ssh = test
            if test.time_by_module["ASM:"] < min_CPCP_flx:
                min_CPCP_flx = test.time_by_module["ASM:"]
                best_CPCP_flx = test
            precision_RSECPCP.append(test.Error_on_Kinf)
                            
    plot_module_times(RSEMOCMOC, tests[0].tracking_module, f"{tests[0].benchmark} MOC/MOC, omp = {tests[0].omp_threads}", f"{tests[0].benchmark}_{tests[0].tracking_module}_Module_Times_RSEMOCMOC_omp{tests[0].omp_threads}")
    plot_module_times(RSECPMOC, tests[0].tracking_module, f"{tests[0].benchmark} CP/MOC, omp = {tests[0].omp_threads}",  f"{tests[0].benchmark}_{tests[0].tracking_module}_Module_Times_RSECPMOC_omp{tests[0].omp_threads}")
    plot_module_times(RSECPCP, tests[0].tracking_module, f"{tests[0].benchmark} CP/CP, omp = {tests[0].omp_threads}", f"{tests[0].benchmark}_{tests[0].tracking_module}_Module_Times_RSECPCP_omp{tests[0].omp_threads}")

    print(f"$$ Outputs for {tests[0].benchmark} benchmark case, with {tests[0].omp_threads} omp threads $$")
    print("Best tests are :")
    print(f"Best MOC/MOC test for ssh : BATCH ssh = {best_MOCMOC_ssh.batch_ssh}, with time spent in USS = {min_MOCMOC_ssh}, and omp = {best_MOCMOC_ssh.omp_threads}")
    print(f"Best MOC/MOC test for flx : BATCH flx = {best_MOCMOC_flx.batch_flx}, with time spent in FLU+ASM = {min_MOCMOC_flx} and omp = {best_MOCMOC_flx.omp_threads}")
    print(f"Best CP/MOC test for ssh : BATCH ssh = {best_CPMOC_ssh.batch_ssh}, with time spent in USS = {min_CPMOC_ssh} and omp = {best_CPMOC_ssh.omp_threads}")
    print(f"Best CP/MOC test for flx : BATCH flx = {best_CPMOC_flx.batch_flx}, with time spent in FLU+ASM = {min_CPMOC_flx} and omp = {best_CPMOC_flx.omp_threads}")
    print(f"Best CP/CP test for ssh : BATCH ssh = {best_CPCP_ssh.batch_ssh}, with time spent in USS = {min_CPCP_ssh} and omp = {best_CPCP_ssh.omp_threads}")
    print(f"Best CP/CP test for flx : BATCH flx = {best_CPCP_flx.batch_flx}, with time spent in ASM = {min_CPCP_flx} and omp = {best_CPCP_flx.omp_threads}")


    return RSEMOCMOC_list, RSECPMOC_list, RSECPCP_list

def time_vs_omp_for_batch(benchmark, MOCMOC_tests, CPMOC_tests, CPCP_tests, batch_ssh, batch_flx):
    print("In time_vs_omp_for_batch")
    # filter tests
    MOCMOC_for_batch_SALT = [test for test in MOCMOC_tests if test.batch_ssh == batch_ssh and test.batch_flx == batch_flx and test.tracking_module == "SALT"]
    CPMOC_for_batch_SALT = [test for test in CPMOC_tests if test.batch_ssh == batch_ssh and test.batch_flx == batch_flx and test.tracking_module == "SALT"]
    CPCP_for_batch_SALT = [test for test in CPCP_tests if test.batch_ssh == batch_ssh and test.batch_flx == batch_flx and test.tracking_module == "SALT"]
    omp_threads = [test.omp_threads for test in MOCMOC_for_batch_SALT]
    print(omp_threads)
    time_spent_MOCMOC_SALT = [test.total_time for test in MOCMOC_for_batch_SALT]
    time_spent_CPMOC_SALT = [test.total_time for test in CPMOC_for_batch_SALT]
    time_spent_CPCP_SALT = [test.total_time for test in CPCP_for_batch_SALT]
    plt.plot(omp_threads, time_spent_MOCMOC_SALT, label="SALT: MOC/MOC", color='r', marker='o')
    plt.plot(omp_threads, time_spent_CPMOC_SALT, label="SALT: CP/MOC", color='b', marker='x')
    plt.plot(omp_threads, time_spent_CPCP_SALT, label="SALT: CP/CP", color='g', marker='s')
    plt.legend()
    plt.xlabel('OMP threads')
    plt.ylabel('Total time spent (s)')
    plt.title(f"Total time vs OMP threads for batch = {batch_flx}")
    plt.savefig(f"figures_time_opti/{benchmark}_time_vs_omp_{batch_ssh}_{batch_flx}.png")
    plt.show()
    plt.close()

    # Now plot time spent in USS and FLU+ASM for MOC/MOC and CP/MOC
    time_spent_MOCMOC_USS_SALT = [test.time_by_module["USS:"] for test in MOCMOC_for_batch_SALT]
    time_spent_MOCMOC_FLU_SALT = [test.time_by_module["FLU:"] for test in MOCMOC_for_batch_SALT]
    time_spent_CPMOC_USS_SALT = [test.time_by_module["USS:"] for test in CPMOC_for_batch_SALT]
    time_spent_CPMOC_FLU_SALT = [test.time_by_module["FLU:"] for test in CPMOC_for_batch_SALT]
    plt.plot(omp_threads, time_spent_MOCMOC_USS_SALT, label="SALT: MOC USS: ", color='r', marker='o')
    plt.plot(omp_threads, time_spent_MOCMOC_FLU_SALT, label="SALT: MOC FLU:", color='b', marker='x')
    plt.plot(omp_threads, time_spent_CPMOC_USS_SALT, label="SALT: CP USS:", color='g', marker='s')
    #plt.plot(omp_threads, time_spent_CPMOC_FLU_SALT, label="CP->MOC FLU:", color='y', marker='d')
    plt.legend()
    plt.xlabel('OMP threads')
    plt.ylabel('Time spent / Module(s)')
    plt.title(f"Time / Module vs OMP threads for {benchmark} batch = {batch_ssh}, {batch_flx}")
    plt.savefig(f"figures_time_opti/{benchmark}_time_vs_omp_{batch_ssh}_{batch_flx}_USS_FLU_ASM.png")
    plt.show()
    plt.close()


    # now plot time spent in USS and ASM for CP/CP
    time_spent_CPCP_USS_SALT = [test.time_by_module["USS:"] for test in CPCP_for_batch_SALT]
    time_spent_CPCP_ASM_SALT = [test.time_by_module["ASM:"] for test in CPCP_for_batch_SALT]
    plt.plot(omp_threads, time_spent_CPCP_USS_SALT, label="SALT: CP USS:", color='r', marker='o')
    plt.plot(omp_threads, time_spent_CPCP_ASM_SALT, label="SALT: CP ASM:", color='b', marker='x')
    plt.legend()
    plt.xlabel('OMP threads')
    plt.ylabel('Time spent / Module (s)')
    plt.title(f"Time / Module vs OMP threads for for {benchmark} batch = {batch_ssh}, {batch_flx}")
    plt.savefig(f"figures_time_opti/{benchmark}_time_vs_omp_{batch_ssh}_{batch_flx}_USS_ASM_CPCP.png")
    plt.show()
    plt.close()




def plot_module_times(test_cases, tracking_module, test_name, save_name):
    """
    Plots histograms representing the time spent in each module for different test cases.
    
    Parameters:
    test_cases (dict): A dictionary where keys are test case names and values are dictionaries
                       with module names as keys and time spent as values.
    """
    # Set the width of each bar
    bar_width = 0.1

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
        ax.bar(indices + i * bar_width, module_times, bar_width, label=test_case, color=colors[i])

    # Adding labels
    ax.set_xlabel('Modules')
    ax.set_ylabel('Time Spent (s)')
    ax.set_title(f'Time/Module for {test_name}, {tracking_module}')
    ax.set_xticks(indices + bar_width * (num_test_cases - 1) / 2)
    ax.set_xticklabels(modules, rotation=45, ha='right')
    ax.legend(loc="best", fontsize=9)

    plt.tight_layout()
    plt.savefig(f'figures_time_opti/{save_name}.png')
    plt.show()
    plt.close()






def create_dictionary_for_plot(tests_list):
    """
    In order to use histogram plotting function, need to create a dictionary where keys are test case names and values are dictionaries
    """
    test_dict = {}
    for test in tests_list:
        test_dict[test.data] = test.time_by_module
    return test_dict


class para_test:
    def __init__(self, procedure, omp_threads, benchmark, ssh_opt, ssh_sol, flx_sol, batch_ssh, batch_flx):
        self.proc = procedure
        self.omp_threads = omp_threads
        self.benchmark = benchmark
        self.ssh_opt = ssh_opt
        self.ssh_sol = ssh_sol
        self.flx_sol = flx_sol
        self.batch_ssh = batch_ssh
        self.batch_flx = batch_flx

        self.data = f"BATCH ssh: {self.batch_ssh}, flx: {self.batch_flx}"

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
        print(f"Adding kinf= {self.Kinf} to test {self.proc}, {self.benchmark}, {self.ssh_opt}, {self.ssh_sol}")
        return
    
    def plot_time_consumption(self):
        print(f"Test with {self.Kinf} to test {self.proc}, {self.benchmark}, {self.ssh_opt}, {self.ssh_sol}")
        #print(self.time_by_module)
        # plot histogram of time spent by module

        plt.bar(self.time_by_module.keys(), self.time_by_module.values())
        plt.xlabel('Module')
        plt.ylabel(f'Time spent (s)')
        plt.title(f"Time spent by module for {self.proc}, {self.benchmark}, {self.ssh_opt}, {self.ssh_sol}")
        plt.show()
        plt.savefig(f"figures_time_opti/test_histogram_{self.proc}_{self.ssh_sol}_{self.batch_ssh}_{self.flx_sol}_{self.batch_flx}.png")
        plt.close()

    def plot_memory_consumption(self):
        print(f"Test with {self.Kinf} to test {self.proc}, {self.benchmark}, {self.ssh_opt}, {self.ssh_sol}")
        print(self.memory_by_module)


ompthread = [10, 20, 30, 40, 50]
MOCMOC = []
CPMOC = []
CPCP = []
for omp in ompthread:
    path_to_results = f"/home/p117902/working_dir/PolyBWR_project/Version5_wc/Dragon/Linux_aarch64/AT10_3x3_para_{omp}.result"
    tests = parse_output(path_to_results, omp)
    results = analyse_tests(tests, 0.89203)
    MOCMOC.extend(results[0])
    CPMOC.extend(results[1])
    CPCP.extend(results[2])

common_batch_lengths = [20,50,100,200,350,500]
for batch in common_batch_lengths:
    time_vs_omp_for_batch("AT10_3x3", MOCMOC, CPMOC, CPCP, batch, batch)
#time_vs_omp_for_batch(MOCMOC, CPMOC, CPCP, 100, 100)

#compare_para_settings(all_tests, ["AT10_3x3"], "MOC", "MOC")
#compare_para_settings(all_tests, ["AT10_3x3"], "CP", "CP")
#compare_para_settings(all_tests, ["AT10_3x3"], "CP", "MOC")