# Class intended to replace the compareD5_S2.py classes and make comparisons more general


import numpy as np
import matplotlib.pyplot as plt
import os
import sys

class compare_D5_S2_rates_XS:
    def __init__(self, case_name, study_type, D5_case, S2_case, S2_libs, compo_keywords, isotopes, self_shielding_methods, save_path):
        self.E0 = 1.0E+07 # eV --> 10 MeV
        self.case_name = case_name
        self.study_type = study_type
        self.D5_case = D5_case
        self.S2_case = S2_case

        # Data used to initialize the nested dictionaries for the D5-S2 differences
        self.S2_libs = S2_libs
        self.compo_keywords = compo_keywords
        self.isotopes = isotopes
        self.self_shielding_methods = self_shielding_methods

        # Used for plot exporation
        self.save_path = save_path
        reactions = ["ngamma"]
        # Initialize the nested dictionaries for the D5-S2 differences
        self.D5_S2_rates_diff = self.initialize_nested_dict(S2_libs, compo_keywords, isotopes, reactions, self_shielding_methods)
        self.D5_S2_XS_diff =  self.initialize_nested_dict(S2_libs, compo_keywords, isotopes, reactions, self_shielding_methods)

        return  

    def initialize_nested_dict(self, S2_libs, keywords, isotopes, reactions, ssh_methods):
        """
        Initializes a nested dictionary with the structure:
        S2_libs -> keywords -> isotopes -> reactions -> ssh_methods -> empty list.
        
        Parameters:
            S2_libs (list): List of keys for the first level.
            keywords (list): List of keys for the second level.
            isotopes (list): List of keys for the third level.
            reactions (list): List of keys for the fourth level.
            ssh_methods (list): List of keys for the innermost level.
            
        Returns:
            dict: Nested dictionary with the specified structure.
        """
        nested_dict = {}
        
        for lib in S2_libs:
            nested_dict[lib] = {}
            for keyword in keywords:
                nested_dict[lib][keyword] = {}
                for isotope in isotopes:
                    nested_dict[lib][keyword][isotope] = {}
                    for reaction in reactions:
                        nested_dict[lib][keyword][isotope][reaction] = {}
                        for ssh_method in ssh_methods:
                            nested_dict[lib][keyword][isotope][reaction][ssh_method] = []
        
        return nested_dict
       
    
    def compare_reaction_rates_and_XS(self, reaction_id, keyword, SSH_methods):
        iso = reaction_id.split("_")[0]
        reaction = reaction_id.split("_")[1]
        for library in self.S2_libs:
            # for now only consider BU step = 0
            for SSH in SSH_methods:
                if self.study_type == "EMESH":
                    self.D5_S2_rates_diff[library][keyword][iso][reaction][SSH] = (self.D5_case.reaction_data_pair[reaction_id][keyword][SSH]["rates"] - self.S2_case.rates[keyword][f"{iso}_{reaction}_{library}_0"])*100/self.S2_case.rates[keyword][f"{iso}_{reaction}_{library}_0"]
                    self.D5_S2_XS_diff[library][keyword][iso][reaction][SSH] = (self.D5_case.reaction_data_pair[reaction_id][keyword][SSH]["XS"] - self.S2_case.XS[keyword][f"{iso}_{reaction}_{library}_0"])*100/self.S2_case.XS[keyword][f"{iso}_{reaction}_{library}_0"]
                else:
                    self.D5_S2_rates_diff[library][keyword][iso][reaction][SSH] = (self.D5_case.reaction_data_pair[reaction_id][keyword][SSH]["rates"] - self.S2_case.rates["SHEM295"][f"{iso}_{reaction}_{library}_0"])*100/self.S2_case.rates["SHEM295"][f"{iso}_{reaction}_{library}_0"]
                    self.D5_S2_XS_diff[library][keyword][iso][reaction][SSH] = (self.D5_case.reaction_data_pair[reaction_id][keyword][SSH]["XS"] - self.S2_case.XS["SHEM295"][f"{iso}_{reaction}_{library}_0"])*100/self.S2_case.XS["SHEM295"][f"{iso}_{reaction}_{library}_0"]
        return
    
    def plot_XS_D5_S2(self, reaction_id, keyword, mesh_name, SSH_methods, bu_step = 0):
        iso = reaction_id.split("_")[0]
        reaction = reaction_id.split("_")[1]
        if reaction == "ngamma":
            reaction_print = "$(n,\\gamma)$"
        u_mesh = self.D5_case.mesh_objects[mesh_name].lethargyMesh
        plt.figure(figsize=(10, 6))
        lib_id = ""
        for library in self.S2_libs:
            lib_id += library
            u = []
            S2_XS_to_plot = []
            S2_XS = self.S2_case.XS[mesh_name][f"{iso}_{reaction}_{library}_{bu_step}"]
            for i in range(len(S2_XS)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                S2_XS_to_plot.extend([S2_XS[i],S2_XS[i]])
            plt.step(u, S2_XS_to_plot, where='post', label=f"S2: $\\sigma$ {reaction_print} {library} for {iso}")
        ssh_id = "" 
        for SSH in SSH_methods:    
            ssh_id += SSH
            u = []   
            D5_XS_to_plot = []     
            D5_XS = self.D5_case.reaction_data_pair[reaction_id][keyword][SSH]["XS"]
            for i in range(len(D5_XS)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                D5_XS_to_plot.extend([D5_XS[i],D5_XS[i]])
            plt.step(u, D5_XS_to_plot, where='post', label=f"{SSH}: $\\sigma$ {reaction_print} for {iso}")
        plt.xlabel("Lethargy")
        plt.ylabel("Cross section [barn]")
        plt.yscale("log")
        plt.title(f"Comparison of {iso} {reaction_print} cross sections between DRAGON5 and SERPENT2")
        plt.grid()
        plt.legend()
        plt.savefig(f"{self.save_path}/{iso}_{reaction}_XS_D5_{keyword}_{ssh_id}_S2_{lib_id}_{mesh_name}.png")
        plt.close()
        return
    
    def plot_rates_D5_S2(self, reaction_id, keyword, mesh_name, SSH_methods, bu_step = 0):
        iso = reaction_id.split("_")[0]
        reaction = reaction_id.split("_")[1]
        if reaction == "ngamma":
            reaction_print = "$(n,\\gamma)$"
        u_mesh = self.D5_case.mesh_objects[mesh_name].lethargyMesh
        plt.figure(figsize=(10, 6))
        lib_id = ""
        for library in self.S2_libs:
            lib_id += library
            u = []
            S2_rates_to_plot = []
            S2_rates = self.S2_case.rates[mesh_name][f"{iso}_{reaction}_{library}_{bu_step}"]
            for i in range(len(S2_rates)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                S2_rates_to_plot.extend([S2_rates[i],S2_rates[i]])
            plt.step(u, S2_rates_to_plot, where='post', label=f"S2: $\\tau$ {reaction_print} {library} for {iso}")
        ssh_id = ""
        for SSH in SSH_methods:   
            ssh_id += SSH 
            u = []   
            D5_rates_to_plot = []     
            D5_rates = self.D5_case.reaction_data_pair[reaction_id][keyword][SSH]["rates"]
            for i in range(len(D5_rates)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                D5_rates_to_plot.extend([D5_rates[i],D5_rates[i]])
            plt.step(u, D5_rates_to_plot, where='post', label=f"{SSH}: $\\tau$ {reaction_print} for {iso}")
        plt.xlabel("Lethargy")
        plt.ylabel(f"$\\tau$ {reaction_print} {iso}")
        plt.yscale("log")
        plt.title(f"Comparison of {iso} {reaction_print} reaction rates between DRAGON5 and SERPENT2")
        plt.grid()
        plt.legend()
        plt.savefig(f"{self.save_path}/{iso}_{reaction}_rates_D5_{keyword}_{ssh_id}_S2_{lib_id}_{mesh_name}.png")
        plt.close()
        return
    
    def plot_diff_rates_D5_S2(self, reaction_id, keyword, mesh_name, library, SSH_methods):
        iso = reaction_id.split("_")[0]
        reaction = reaction_id.split("_")[1]
        u_mesh = self.D5_case.mesh_objects[mesh_name].lethargyMesh
        if reaction == "ngamma":
            reaction_print = "$(n,\\gamma)$"
        ssh_id = ""
        plt.figure(figsize=(10, 6))
        for SSH in SSH_methods:
            ssh_id += SSH 
            u = []   
            D5_S2_rates_diff_to_plot = []     
            D5_S2_rates_diff = self.D5_S2_rates_diff[library][keyword][iso][reaction][SSH]
            for i in range(len(D5_S2_rates_diff)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                D5_S2_rates_diff_to_plot.extend([D5_S2_rates_diff[i],D5_S2_rates_diff[i]])
            plt.step(u, D5_S2_rates_diff_to_plot, label=f"{SSH} - {library} : {iso} {reaction_print}")
        plt.xlabel("Lethargy")
        plt.ylabel(f"$\\Delta \\tau$ {iso} {reaction_print} [%]")
        plt.title(f"Relative difference on $\\tau$ {reaction_print} {iso} between DRAGON5 {keyword} and SERPENT2 ({library}) on {mesh_name}")
        plt.grid()
        plt.legend()
        plt.savefig(f"{self.save_path}/{iso}_{reaction}_rates_diff_{keyword}_{mesh_name}_{library}_{ssh_id}.png")
        plt.close()
        return
    
    def plot_diff_XS_D5_S2(self, reaction_id, keyword, mesh_name, library, SSH_methods):
        iso = reaction_id.split("_")[0]
        reaction = reaction_id.split("_")[1]
        u_mesh = self.D5_case.mesh_objects[mesh_name].lethargyMesh
        if reaction == "ngamma":
            reaction_print = "$(n,\\gamma)$"
        ssh_id = ""
        plt.figure(figsize=(10, 6))
        for SSH in SSH_methods:
            ssh_id += SSH 
            u = []   
            D5_S2_XS_diff_to_plot = []     
            D5_S2_XS_diff = self.D5_S2_XS_diff[library][keyword][iso][reaction][SSH]
            for i in range(len(D5_S2_XS_diff)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                D5_S2_XS_diff_to_plot.extend([D5_S2_XS_diff[i],D5_S2_XS_diff[i]])
            plt.step(u, D5_S2_XS_diff_to_plot, label=f"{SSH} - {library} : {iso} {reaction_print}")
        plt.xlabel("Lethargy")
        plt.ylabel(f"$\\Delta \\sigma$ {reaction_print} {iso}  [%]")
        plt.title(f"Relative difference on $\\sigma$ {reaction_print} {iso} between DRAGON5 and SERPENT2 ({library}) on {mesh_name}")
        plt.grid()
        plt.legend()
        plt.savefig(f"{self.save_path}/{iso}_{reaction}_XS_diff_{keyword}_{mesh_name}_{library}_{ssh_id}.png")
        plt.close()
        return
    
    def renorm_rates(self, reaction_id, keyword, mesh_name):
        """
        Since DRAGON5 and SERPENT2 do not use the same normalization for fluxes, 
        this function renormalizes the rates in order to allow for a direct D5-S2 comparison.
        """
        # Renormalize DRAGON5 and Serpent2 rates to 1
        # For all Serpent2 libraries and DRAGON5 SSH methods
        #D5_renorm = self.D5_case.reaction_data_pair[reaction_id]["SHEM281"]["Autosecol"]["rates"][0]/self.S2_case.rates["SHEM281"][f"{reaction_id}_oldlib_0"][0]
        for SSH in self.D5_case.reaction_data_pair[reaction_id][keyword].keys():
            self.D5_case.reaction_data_pair[reaction_id][keyword][SSH]["rates"] = self.D5_case.reaction_data_pair[reaction_id][keyword][SSH]["rates"]/np.sum(self.D5_case.reaction_data_pair[reaction_id][keyword][SSH]["rates"])
        for library in self.S2_libs:
            self.S2_case.rates[mesh_name][f"{reaction_id}_{library}_0"] = self.S2_case.rates[mesh_name][f"{reaction_id}_{library}_0"]/np.sum(self.S2_case.rates[mesh_name][f"{reaction_id}_{library}_0"])
        return

    def find_top_differences_rates(self, library, reaction_id, keyword, mesh_name, top_n=5):
        """
        Finds the indices of the top N maximal absolute differences in a nested dictionary.
        
        Parameters:
            diff (dict): Nested dictionary containing differences as 
                        diff[library][mesh][reaction_id][SSH].
            library (str): The library to extract data from.
            reaction_id (str): The reaction ID to extract data from.
            keyword (str): The keyword/specific COMPO calculation to retreive data from.
            mesh_name (str): The identifier for energy mesh.
            ssh (str): The SSH key to extract data from.
            top_n (int): Number of top differences to return (default: 5).
        
        Returns:
            list: Indices of the top N absolute differences.
        """
        print("$$$--- Begin find_top_differences_rates ---$$$")
        diff_data = []
        iso = reaction_id.split("_")[0]
        reaction = reaction_id.split("_")[1]
        # Retrieve the differences for the given parameters
    
        for SSH in self.D5_S2_rates_diff[library][keyword][iso][reaction].keys():
            differences = self.D5_S2_rates_diff[library][keyword][iso][reaction][SSH]
    
            # Compute the absolute values of the differences
            absolute_differences = np.abs(differences)
    
            # Get the indices of the top N maximal differences
            top_indices = sorted(range(len(absolute_differences)), 
                                key=lambda i: absolute_differences[i], 
                                reverse=True)[:top_n]
            
            # Print the top N maximal differences
            print("\n")
            print(f"Top {top_n} maximal differences on reaction rates for {reaction_id} in {library} library, {keyword} keyword {mesh_name} mesh, {SSH} SSH:")
            for i, index in enumerate(top_indices):
                print(f"Error is {differences[index]} % in group {index+1}")
                print(f"With lethargy bounds u_min = {self.D5_case.mesh_objects[mesh_name].lethargyMesh[index]} and u_max = {self.D5_case.mesh_objects[mesh_name].lethargyMesh[index+1]}")
                print(f"and energy bounds E_min = {self.E0*np.exp(-self.D5_case.mesh_objects[mesh_name].lethargyMesh[index+1])} eV and E_max = {self.E0*np.exp(-self.D5_case.mesh_objects[mesh_name].lethargyMesh[index])} eV")
                print("\n")
                u_min, u_max = self.D5_case.mesh_objects[mesh_name].lethargyMesh[index], self.D5_case.mesh_objects[mesh_name].lethargyMesh[index+1]
                E_min, E_max = self.E0*np.exp(-u_max), self.E0*np.exp(-u_min)
                diff_data.append({"SSH Method": SSH, "Error (%)": differences[index], "Group": index+1, "u_min": u_min, "u_max": u_max, "E_min": E_min, "E_max": E_max})
            
        print("$$$--- End find_top_differences_rates ---$$$")
        return diff_data

    
    def find_top_differences_XS(self, library, reaction_id, keyword, mesh_name, top_n=5):
        """
        Finds the indices of the top N maximal absolute differences in a nested dictionary.
        
        Parameters:
            diff (dict): Nested dictionary containing differences as 
                        diff[library][mesh][reaction_id][SSH].
            library (str): The library to extract data from.
            reaction_id (str): The reaction ID to extract data from.
            keyword (str): The keyword/specific COMPO calculation to retreive data from.
            mesh_name (str): The identifier for energy mesh.
            ssh (str): The SSH key to extract data from.
            top_n (int): Number of top differences to return (default: 5).
        
        Returns:
            list: Indices of the top N absolute differences.
        """
        print("$$$--- Begin find_top_differences_XS ---$$$")
        diff_data = []
        iso = reaction_id.split("_")[0]
        reaction = reaction_id.split("_")[1]
        # Retrieve the differences for the given parameters
    
        for SSH in self.D5_S2_XS_diff[library][keyword][iso][reaction].keys():
            differences = self.D5_S2_XS_diff[library][keyword][iso][reaction][SSH]
    
            # Compute the absolute values of the differences
            absolute_differences = np.abs(differences)
    
            # Get the indices of the top N maximal differences
            top_indices = sorted(range(len(absolute_differences)), 
                                key=lambda i: absolute_differences[i], 
                                reverse=True)[:top_n]
            
            # Print the top N maximal differences
            print("\n")
            print(f"Top {top_n} maximal differences on reaction XS for {reaction_id} in {library} library, {keyword} keyword {mesh_name} mesh, {SSH} SSH:")
            for i, index in enumerate(top_indices):
                print(f"Error is {differences[index]} % in group {index+1}")
                print(f"With lethargy bounds u_min = {self.D5_case.mesh_objects[mesh_name].lethargyMesh[index]} and u_max = {self.D5_case.mesh_objects[mesh_name].lethargyMesh[index+1]}")
                print(f"and energy bounds E_min = {self.E0*np.exp(-self.D5_case.mesh_objects[mesh_name].lethargyMesh[index+1])} eV and E_max = {self.E0*np.exp(-self.D5_case.mesh_objects[mesh_name].lethargyMesh[index])} eV")
                print("\n")
                u_min, u_max = self.D5_case.mesh_objects[mesh_name].lethargyMesh[index], self.D5_case.mesh_objects[mesh_name].lethargyMesh[index+1]
                E_min, E_max = self.E0*np.exp(-u_max), self.E0*np.exp(-u_min)
                diff_data.append({"SSH Method": SSH, "Error (%)": differences[index], "Group": index+1, "u_min": u_min, "u_max": u_max, "E_min": E_min, "E_max": E_max})
            
        print("$$$--- End find_top_differences_XS ---$$$")
        return diff_data
    
    def plot_zoom_XS_and_errors(self, library, rection_id, keyword, mesh_name, SSH_methods, grmin, grmax):
        """
        use subplot to plot 1) XS (D5 and S2) and 2) relative errors (D5-S2)/S2 for each SSH method
        grmin and grmax are the group indices to zoom in
        """
        ssh_id = ""
        iso = rection_id.split("_")[0]
        reaction = rection_id.split("_")[1]
        if reaction == "ngamma":
            reaction_print = "$(n,\\gamma)$"
        # create lethargy mesh to be used for step plot:
        u_mesh = self.D5_case.mesh_objects[mesh_name].lethargyMesh[grmin-1:grmax+1]
        u = []
        for i in range(len(u_mesh)-1): # -1 since we have one less interval than mesh points
            u.extend([u_mesh[i],u_mesh[i+1]])
        #plt.figure(figsize=(10, 6))
        number_of_subplot = len(SSH_methods)+1
        fig,ax = plt.subplots(number_of_subplot,1,figsize=(20, 3*number_of_subplot))
        lib_id = library
        S2_XS = self.S2_case.XS[mesh_name][f"{iso}_{reaction}_{library}_0"][grmin-1:grmax] # correct for python indexing
        S2_XS_to_plot = []
        for i in range(len(S2_XS)):
            S2_XS_to_plot.extend([S2_XS[i],S2_XS[i]])
        ax[0].step(u, S2_XS_to_plot, where='post', label=f"S2: $\\sigma$ {reaction_print} {library} for {iso}")
        for SSH in SSH_methods:
            D5_XS_to_plot = []     
            D5_XS = self.D5_case.reaction_data_pair[rection_id][keyword][SSH]["XS"][grmin-1:grmax]
            for i in range(len(D5_XS)):
                D5_XS_to_plot.extend([D5_XS[i],D5_XS[i]])
            ax[0].step(u, D5_XS_to_plot, where='post', label=f"{SSH}: $\\sigma$ {reaction_print} for {iso}")
        ax[0].set_xlabel("Lethargy")
        ax[0].set_ylabel("Cross section [barn]")
        ax[0].set_yscale("log")
        ax[0].grid()
        ax[0].legend()
        ax[0].set_title(f"Comparison of {iso} {reaction_print} cross sections between DRAGON5 {keyword} and SERPENT2")

        for i, SSH in enumerate(SSH_methods):
            ssh_id += SSH
            D5_S2_XS_diff_to_plot = []     
            D5_S2_XS_diff = self.D5_S2_XS_diff[library][keyword][iso][reaction][SSH][grmin-1:grmax]
            for j in range(len(D5_S2_XS_diff)):
                D5_S2_XS_diff_to_plot.extend([D5_S2_XS_diff[j],D5_S2_XS_diff[j]])
            ax[i+1].step(u, D5_S2_XS_diff_to_plot, label=f"{SSH} - {library} : {iso} {reaction_print}")
            ax[i+1].set_xlabel("Lethargy")
            ax[i+1].set_ylabel(f"$\\Delta \\sigma$ {reaction_print} {iso}  [%]")
            ax[i+1].set_title(f"Relative difference on $\\sigma$ {reaction_print} {iso} between DRAGON5 {SSH} and SERPENT2 ({library}) on {mesh_name}")
            ax[i+1].grid()
            ax[i+1].legend()
        
        plt.tight_layout()
        #plt.title(f"Comparison of {iso} {reaction_print} cross sections between DRAGON5 and SERPENT2")
        plt.savefig(f"{self.save_path}/{iso}_{reaction}_XS_D5_{keyword}_{ssh_id}_S2_{lib_id}_{mesh_name}_zoomed_{grmin}_{grmax}.png")
        plt.close()
        return
    
    def plot_zoom_rates_and_errors(self, library, rection_id, keyword, mesh_name, SSH_methods, grmin, grmax):
        """
        use subplot to plot 1) rates (D5 and S2) and 2) relative errors (D5-S2)/S2 for each SSH method
        grmin and grmax are the group indices to zoom in
        """
        ssh_id = ""
        iso = rection_id.split("_")[0]
        reaction = rection_id.split("_")[1]
        if reaction == "ngamma":
            reaction_print = "$(n,\\gamma)$"
        # create lethargy mesh to be used for step plot:
        u_mesh = self.D5_case.mesh_objects[mesh_name].lethargyMesh[grmin-1:grmax+1]
        u = []
        for i in range(len(u_mesh)-1):
            u.extend([u_mesh[i],u_mesh[i+1]])
        #plt.figure(figsize=(10, 6))
        number_of_subplot = len(SSH_methods)+1
        fig,ax = plt.subplots(number_of_subplot,1,figsize=(20, 3*number_of_subplot))
        lib_id = library
        S2_rates = self.S2_case.rates[mesh_name][f"{iso}_{reaction}_{library}_0"][grmin-1:grmax]
        S2_rates_to_plot = []
        for i in range(len(S2_rates)):
            S2_rates_to_plot.extend([S2_rates[i],S2_rates[i]])
        ax[0].step(u, S2_rates_to_plot, where='post', label=f"S2: $\\tau$ {reaction_print} {library} for {iso}")
        for SSH in SSH_methods:
            D5_rates_to_plot = []     
            D5_rates = self.D5_case.reaction_data_pair[rection_id][keyword][SSH]["rates"][grmin-1:grmax]
            for i in range(len(D5_rates)):
                D5_rates_to_plot.extend([D5_rates[i],D5_rates[i]])
            ax[0].step(u, D5_rates_to_plot, where='post', label=f"{SSH}: $\\tau$ {reaction_print} for {iso}")
        ax[0].set_xlabel("Lethargy")
        ax[0].set_ylabel("Reaction rates")
        ax[0].set_yscale("log")
        ax[0].grid()
        ax[0].legend()
        ax[0].set_title(f"Comparison of {iso} {reaction_print} reaction rates between DRAGON5 {keyword} and SERPENT2")

        for i, SSH in enumerate(SSH_methods):
            ssh_id += SSH
            D5_S2_rates_diff_to_plot = []     
            D5_S2_rates_diff = self.D5_S2_rates_diff[library][keyword][iso][reaction][SSH][grmin-1:grmax]
            for j in range(len(D5_S2_rates_diff)):
                D5_S2_rates_diff_to_plot.extend([D5_S2_rates_diff[j],D5_S2_rates_diff[j]])
            ax[i+1].step(u, D5_S2_rates_diff_to_plot, label=f"{SSH} - {library} : {iso} {reaction_print}")
            ax[i+1].set_xlabel("Lethargy")
            ax[i+1].set_ylabel(f"$\\Delta \\tau$ {reaction_print} {iso}  [%]")
            ax[i+1].set_title(f"Relative difference on $\\tau$ {reaction_print} {iso} between DRAGON5 {SSH} and SERPENT2 ({library}) on {mesh_name}")
            ax[i+1].grid()
            ax[i+1].legend()
        plt.tight_layout()
        #plt.title(f"Comparison of {iso} {reaction_print} reaction rates between DRAGON5 and SERPENT2")
        plt.savefig(f"{self.save_path}/{iso}_{reaction}_rates_D5_{keyword}_{ssh_id}_S2_{lib_id}_{mesh_name}_zoomed_{grmin}_{grmax}.png")
        plt.close()
        return