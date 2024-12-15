# Python3 class to unify post treatment of results from HOM_UOX_Gd157 calculations
# Author: R. Guasch
# Date: 2024-11-14
# Purpose: Unify post treatment of results from HOM_UOX_Gd157 and facilitate D5-AUTO and D5-S2 comparisons
#

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import serpentTools as st

def initialize_nested_dict(level1_keys, level2_keys, subkeys, default_value=[]):
        """
        Initializes a nested dictionary with a specified structure.

        Parameters:
            level1_keys (list): Keys for the first level of the dictionary.
            level2_keys (list): Keys for the second level of the dictionary.
            subkeys (list): Subkeys for the innermost level.
            default_value: The default value for innermost subkeys. Defaults to an empty list.

        Returns:
            dict: A nested dictionary with the specified structure.
        """
        # Create the innermost dictionary structure using subkeys
        innermost_dict = {key: default_value.copy() if isinstance(default_value, list) else default_value for key in subkeys}

        # Create the second level dictionary structure using level2_keys
        second_level_dict = {key: innermost_dict.copy() for key in level2_keys}

        # Create the top-level dictionary structure using level1_keys
        nested_dict = {key: second_level_dict.copy() for key in level1_keys}

        return nested_dict

class postTreatment_rates_XS_D5_EMESHES:
    def __init__(self, case_name, mesh_objects, self_shielding_methods, save_path):
        self.case_name = case_name
        self.mesh_objects = mesh_objects
        self.self_shielding_methods = self_shielding_methods
        self.save_path = save_path
        self.reaction_data_pair = {}
        
        self.study_type = "MESHES"
        # dictionaries to store the relative differences between self shielded cross sections, reaction rates and fluxes, identified by reaction_id, mesh_name and SSH method
        mesh_names = self.mesh_objects.keys()
        self.delta_XS = initialize_nested_dict(["Gd157_ngamma", "U238_ngamma"], mesh_names, self.self_shielding_methods, default_value=[])
        self.delta_Rates = initialize_nested_dict(["Gd157_ngamma", "U238_ngamma"], mesh_names, self.self_shielding_methods, default_value=[])
       
        self.delta_Fluxes = {}

        return

    def set_reaction_data(self, reaction_id, reaction_data):
        self.reaction_data_pair[reaction_id] = reaction_data
        return
    
    def set_fluxes(self, fluxes):
        self.fluxes = fluxes
        return
    
    def plot_fluxes(self, mesh_name):
        plt.figure(figsize=(10, 6))
        u_mesh = self.mesh_objects[mesh_name].lethargyMesh
        for SSH in self.fluxes[mesh_name].keys():
            flux = self.fluxes[mesh_name][SSH] 
            u = []
            flux_to_plot = [] 
            for i in range(len(flux)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                flux_to_plot.extend([flux[i],flux[i]])
            plt.step(u, flux_to_plot, where='post', label=f"{SSH} multigroup flux")
        plt.legend()
        plt.yscale('log')
        plt.ylabel('Flux')
        plt.xlabel('Lethargy')
        plt.title(f"Fluxes for {self.case_name}")
        plt.savefig(f"{self.save_path}/{self.case_name}_fluxes_{mesh_name}.png")
        plt.close()
        return
    
    def plot_XS(self, mesh_name, reaction_id):
        plt.figure(figsize=(10, 6))
        iso = reaction_id.split("_")[0]
        reaction = reaction_id.split("_")[1]
        u_mesh = self.mesh_objects[mesh_name].lethargyMesh
        if reaction == "ngamma":
            reaction = "$(n,\gamma)$"
        for SSH in self.reaction_data_pair[reaction_id][mesh_name].keys():
            print(f"plotting XS for {reaction_id} {mesh_name} {SSH}")
            XS = self.reaction_data_pair[reaction_id][mesh_name][SSH]["XS"]
            print(f"XS = {XS}")
            print(len(XS))
            u = []
            XS_to_plot = [] 
            for i in range(len(XS)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                XS_to_plot.extend([XS[i],XS[i]])
            plt.step(u, XS_to_plot, where='post', label=f"{SSH} {reaction} cross section")
        plt.legend()
        plt.xlabel('Lethargy')
        plt.ylabel(f"$\sigma$ {reaction} (barns)")
        plt.yscale('log')
        plt.title(f"Cross sections for {iso} {reaction}")
        plt.savefig(f"{self.save_path}/{self.case_name}_{reaction_id}_XS_{mesh_name}.png")
        plt.close()
        return
    
    def plot_reaction_rates(self, mesh_name, reaction_id):
        plt.figure(figsize=(10, 6))
        iso = reaction_id.split("_")[0]
        reaction = reaction_id.split("_")[1]
        u_mesh = self.mesh_objects[mesh_name].lethargyMesh
        if reaction == "ngamma":
            reaction = "$(n,\gamma)$"
        for SSH in self.reaction_data_pair[reaction_id][mesh_name].keys():
            rates = self.reaction_data_pair[reaction_id][mesh_name][SSH]["rates"]
            u = []
            rates_to_plot = [] 
            for i in range(len(rates)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                rates_to_plot.extend([rates[i],rates[i]])
            plt.step(u, rates_to_plot, where='post', label=f"{SSH} {reaction} reaction rates")
        plt.legend()
        plt.xlabel('Lethargy')
        plt.ylabel(f"{reaction} rates")
        plt.yscale('log')
        plt.grid()
        plt.title(f"Reaction rates for {iso} {reaction}")
        plt.savefig(f"{self.save_path}/{self.case_name}_{reaction_id}_rates_{mesh_name}.png")
        plt.close()
        return

    def compute_relative_differences_XS(self, reaction_id, mesh_name, reference_SSH_method, SSH_methods_to_compare):
        """
        For a given reaction, energy mesh, compute the relative differences on self shielded cross sections between a reference self shielding method and the others.
        """
        XS_ref = self.reaction_data_pair[reaction_id][mesh_name][reference_SSH_method]["XS"]
        self.delta_XS[reaction_id] = {mesh_name: {}}
        for SSH in SSH_methods_to_compare:
            print(f"computing delta XS for {reaction_id} {mesh_name} {SSH} - {reference_SSH_method}")
            XS = self.reaction_data_pair[reaction_id][mesh_name][SSH]["XS"]
            delta_XS = (np.array(XS) - np.array(XS_ref))*100/np.array(XS_ref)
            self.delta_XS[reaction_id][mesh_name][SSH] = delta_XS

        print(f"self.delta_XS = {self.delta_XS}")
        return
    def compute_relative_differences_Rates(self, reaction_id, mesh_name, reference_SSH_method, SSH_methods_to_compare):
        """
        For a given reaction, energy mesh, compute the relative differences on reaction rates between a reference self shielding method and the others.
        """
        rates_ref = self.reaction_data_pair[reaction_id][mesh_name][reference_SSH_method]["rates"]
        self.delta_Rates[reaction_id] = {mesh_name: {}}
        for SSH in SSH_methods_to_compare:
            rates = self.reaction_data_pair[reaction_id][mesh_name][SSH]["rates"]
            
            delta_Rates = (np.array(rates) - np.array(rates_ref))*100/np.array(rates_ref)
            self.delta_Rates[reaction_id][mesh_name][SSH] = delta_Rates
        return
    
    def compute_relative_differences_Fluxes(self, mesh_name, reference_SSH_method, SSH_methods_to_compare):
        """
        For a given energy mesh, compute the relative differences on fluxes between a reference self shielding method and the others.
        """
        flux_ref = self.fluxes[mesh_name][reference_SSH_method]
        self.delta_Fluxes[mesh_name] = {}
        for SSH in SSH_methods_to_compare:
            flux = self.fluxes[mesh_name][SSH]
            delta_Fluxes = (np.array(flux) - np.array(flux_ref))*100/np.array(flux_ref)
            self.delta_Fluxes[mesh_name][SSH] = delta_Fluxes
        return
    
    def plot_histogram_relative_differences_Fluxes(self, mesh_name, SSH_methods_to_compare, grmin, grmax):
        colors = {"RSE":'skyblue', "PT": 'red', "SUBG":"green"}
        save_ssh_id =""
        for SSH in SSH_methods_to_compare:
            save_ssh_id += SSH
        energy_groups = list(range(grmin, grmax + 1))
        plt.figure(figsize=(10, 6))
        for SSH in SSH_methods_to_compare:
            delta_Fluxes = self.delta_Fluxes[mesh_name][SSH]
            plt.bar(energy_groups, delta_Fluxes, edgecolor='black', alpha=0.5, color=colors[SSH], label=f"{SSH} - AUTO")
        plt.ylabel('Relative difference on fluxes (%)')
        plt.xticks(range(grmin, grmax+1, 20))  # Set x-ticks with an interval of 20 for readability
        plt.xlabel('Energy group')
        plt.title(f"Relative differences on fluxes for {mesh_name}")
        plt.legend()
        plt.savefig(f"{self.save_path}/{self.case_name}_delta_Fluxes_{mesh_name}_{save_ssh_id}.png")
        plt.close()
        return
        

    def plot_histogram_relative_differences_Rates(self, reaction_id, mesh_name, SSH_methods_to_compare, grmin, grmax):
        colors = {"RSE":'skyblue', "PT": 'red', "SUBG":"green"}
        iso = reaction_id.split("_")[0]
        reaction = reaction_id.split("_")[1]
        if reaction == "ngamma":
            reaction = "$(n,\gamma)$"
        save_ssh_id =""
        energy_groups = list(range(grmin, grmax + 1))
        plt.figure(figsize=(10, 6))
        for SSH in SSH_methods_to_compare:
            save_ssh_id += SSH
            delta_Rates = self.delta_Rates[reaction_id][mesh_name][SSH]
            plt.bar(energy_groups, delta_Rates, edgecolor='black', alpha=0.5, color=colors[SSH], label=f"{iso} {reaction} {SSH} - AUTO")
        plt.ylabel(f'Relative difference on {reaction} rates (%)')
        plt.xticks(range(grmin, grmax+1, 20))  # Set x-ticks with an interval of 20 for readability
        plt.xlabel('Energy group')
        plt.title(f"Relative differences on {iso} {reaction} rates for {mesh_name}")
        plt.legend()
        plt.savefig(f"{self.save_path}/{self.case_name}_{reaction_id}_delta_Rates_hist_{mesh_name}_{save_ssh_id}.png")
        plt.close()
        return

    
    def plot_histogram_relative_differences_XS(self, reaction_id, mesh_name, SSH_methods_to_compare, grmin, grmax):
        colors = {"RSE":'skyblue', "PT": 'red', "SUBG":"green"}
        iso = reaction_id.split("_")[0]
        reaction = reaction_id.split("_")[1]
        if reaction == "ngamma":
            reaction = "$(n,\gamma)$"
        elif reaction == "abs":
            reaction = "a"
        
        save_ssh_id =""
        energy_groups = list(range(grmin, grmax + 1))
        plt.figure(figsize=(10, 6))
        for SSH in SSH_methods_to_compare:
            save_ssh_id += SSH
            delta_XS = self.delta_XS[reaction_id][mesh_name][SSH]
            plt.bar(energy_groups, delta_XS, edgecolor='black', alpha=0.5, color=colors[SSH], label=f"{iso} {reaction} {SSH} - AUTO")
        plt.ylabel('Relative difference on XS (%)')
        plt.xticks(range(grmin, grmax+1, 20))  # Set x-ticks with an interval of 20 for readability
        plt.xlabel('Energy group')
        plt.title(f"Relative differences on {iso} $\sigma${reaction} for {mesh_name}")
        plt.legend()
        plt.savefig(f"{self.save_path}/{self.case_name}_{reaction_id}_delta_XS_hist_{mesh_name}_{save_ssh_id}.png")
        plt.close()
        return
    
    def plot_relative_differences_Rates(self, reaction_id, mesh_name, SSH_methods_to_compare):
        """
        Plot relative differences on reaction rates for a given reaction and energy mesh
        """
        plt.figure(figsize=(10, 6))
        iso = reaction_id.split("_")[0]
        reaction = reaction_id.split("_")[1]
        if reaction == "ngamma":
            reaction = "$(n,\gamma)$"
        for SSH in SSH_methods_to_compare:
            delta_Rates = self.delta_Rates[reaction_id][mesh_name][SSH]
            u_mesh = self.mesh_objects[mesh_name].lethargyMesh
            u = []
            delta_rates_to_plot = [] 
            for i in range(len(delta_Rates)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                delta_rates_to_plot.extend([delta_Rates[i],delta_Rates[i]])
            plt.step(u, delta_rates_to_plot, where='post', label=f"{iso} : $\\Delta \\tau$ {reaction} {SSH}")
        plt.xlabel('Lethargy')
        plt.ylabel(f"$\\Delta \\tau$ {reaction} (%)")
        plt.legend()
        plt.grid()
        plt.title(f"{iso} : $\\Delta \\tau$ {reaction} (D5 - AUTO) {mesh_name}")
        plt.savefig(f"{self.save_path}/{self.case_name}_{reaction_id}_delta_Rates_{mesh_name}.png")
        plt.close()

        return
    
    def plot_relative_differences_XS(self, reaction_id, mesh_name, SSH_methods_to_compare):
        """
        Plot relative differences in cross sections for a given reaction and energy mesh
        This is not a regular histogram as the energy groups are not equally spaced : each energy group is represented by a bar of lethargy width
        """
        plt.figure(figsize=(10, 6))
        iso = reaction_id.split("_")[0]
        reaction = reaction_id.split("_")[1]
        if reaction == "ngamma":
            reaction = "$(n,\gamma)$"
        elif reaction == "abs":
            reaction = "a"
        for SSH in SSH_methods_to_compare:
            delta_XS = self.delta_XS[reaction_id][mesh_name][SSH]
            u_mesh = self.mesh_objects[mesh_name].lethargyMesh
            u = []
            delta_XS_to_plot = [] 
            for i in range(len(delta_XS)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                delta_XS_to_plot.extend([delta_XS[i],delta_XS[i]])
            plt.step(u, delta_XS_to_plot, where='post', label=f"{iso} : $\\Delta \\sigma$ {reaction} {SSH}")
        plt.xlabel('Lethargy')
        plt.ylabel(f"$\\Delta \\sigma$ {reaction} (%)")
        plt.legend()
        plt.grid()
        plt.title(f"{iso} : $\\Delta \\sigma$ {reaction} (D5 - AUTO) {mesh_name}")
        plt.savefig(f"{self.save_path}/{self.case_name}_{reaction_id}_delta_XS_{mesh_name}.png")
        plt.close()
        return
    
class postTreatment_rates_XS_D5_AUTOlib:
    def __init__(self, case_name, mesh_objects, self_shielding_methods, save_path):
        self.case_name = case_name
        self.mesh_objects = mesh_objects
        self.self_shielding_methods = self_shielding_methods
        self.save_path = save_path
        self.reaction_data_pair = {}
        
        self.study_type = "AUTOLIB"
        # dictionaries to store the relative differences between self shielded cross sections, reaction rates and fluxes, identified by reaction_id, mesh_name and SSH method
        self.delta_XS = initialize_nested_dict(["Gd157_ngamma", "U238_ngamma"], ["J311_295", "J311_295F", "J311_295FF"], self.self_shielding_methods, default_value=[])
        self.delta_rates = initialize_nested_dict(["Gd157_ngamma", "U238_ngamma"], ["J311_295", "J311_295F", "J311_295FF"], self.self_shielding_methods, default_value=[])
        self.delta_Fluxes = {}

        return

    def set_reaction_data(self, reaction_id, reaction_data):
        self.reaction_data_pair[reaction_id] = reaction_data
        return
    
    def set_fluxes(self, fluxes):
        self.fluxes = fluxes
        return
    
    def plot_fluxes(self, autolib_name):
        plt.figure(figsize=(10, 6))
        u_mesh = self.mesh_objects["SHEM295"].lethargyMesh
        for SSH in self.fluxes[autolib_name].keys():
            flux = self.fluxes[autolib_name][SSH] 
            u = []
            flux_to_plot = [] 
            for i in range(len(flux)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                flux_to_plot.extend([flux[i],flux[i]])
            plt.step(u, flux_to_plot, where='post', label=f"{SSH} multigroup flux")
        plt.legend()
        plt.yscale('log')
        plt.ylabel('Flux')
        plt.xlabel('Lethargy')
        plt.title(f"Fluxes for {self.case_name}")
        plt.savefig(f"{self.save_path}/{self.case_name}_fluxes_{autolib_name}.png")
        plt.close()
        return
    
    def plot_XS(self, autolib_name, reaction_id):
        plt.figure(figsize=(10, 6))
        iso = reaction_id.split("_")[0]
        reaction = reaction_id.split("_")[1]
        u_mesh = self.mesh_objects["SHEM295"].lethargyMesh
        if reaction == "ngamma":
            reaction = "$(n,\gamma)$"
        for SSH in self.reaction_data_pair[reaction_id][autolib_name].keys():
            XS = self.reaction_data_pair[reaction_id][autolib_name][SSH]["XS"]
            u = []
            XS_to_plot = [] 
            for i in range(len(XS)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                XS_to_plot.extend([XS[i],XS[i]])
            plt.step(u, XS_to_plot, where='post', label=f"{SSH} {reaction} cross section")
        plt.legend()
        plt.xlabel('Lethargy')
        plt.ylabel(f"$\sigma$ {reaction} (barns)")
        plt.yscale('log')
        plt.title(f"Cross sections for {iso} {reaction}")
        plt.savefig(f"{self.save_path}/{self.case_name}_{reaction_id}_XS_{autolib_name}.png")
        plt.close()
        return
    
    def plot_reaction_rates(self, autolib_name, reaction_id):
        plt.figure(figsize=(10, 6))
        iso = reaction_id.split("_")[0]
        reaction = reaction_id.split("_")[1]
        u_mesh = self.mesh_objects["SHEM295"].lethargyMesh
        if reaction == "ngamma":
            reaction = "$(n,\gamma)$"
        for SSH in self.reaction_data_pair[reaction_id][autolib_name].keys():
            rates = self.reaction_data_pair[reaction_id][autolib_name][SSH]["rates"]
            u = []
            rates_to_plot = [] 
            for i in range(len(rates)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                rates_to_plot.extend([rates[i],rates[i]])
            plt.step(u, rates_to_plot, where='post', label=f"{SSH} {reaction} reaction rates")
        plt.legend()
        plt.xlabel('Lethargy')
        plt.ylabel(f"{reaction} rates")
        plt.yscale('log')
        plt.grid()
        plt.title(f"Reaction rates for {iso} {reaction}")
        plt.savefig(f"{self.save_path}/{self.case_name}_{reaction_id}_rates_{autolib_name}.png")
        plt.close()
        return

class postTreatment_rates_XS_D5_IRSET:
    def __init__(self, case_name, mesh_objects, self_shielding_methods, save_path):
        self.case_name = case_name
        self.mesh_objects = mesh_objects
        self.self_shielding_methods = self_shielding_methods
        self.save_path = save_path
        self.reaction_data_pair = {}
        
        self.study_type = "IRSET"
        # dictionaries to store the relative differences between self shielded cross sections, reaction rates and fluxes, identified by reaction_id, mesh_name and SSH method
        self.delta_XS = initialize_nested_dict(["Gd157_ngamma", "U238_ngamma"], ["noIRSET", "noIRSETF", "noIRSETFF", "IRSET_101", "IRSET_101F", "IRSET_101FF", "IRSETNone", "IRSETNoneF", "IRSETNoneFF"], self.self_shielding_methods, default_value=[])
        self.delta_rates = initialize_nested_dict(["Gd157_ngamma", "U238_ngamma"], ["noIRSET", "noIRSETF", "noIRSETFF", "IRSET_101", "IRSET_101F", "IRSET_101FF", "IRSETNone", "IRSETNoneF", "IRSETNoneFF"], self.self_shielding_methods, default_value=[])
        self.delta_Fluxes = {}

        return


    def set_reaction_data(self, reaction_id, reaction_data):
        self.reaction_data_pair[reaction_id] = reaction_data
        return
    
    def set_fluxes(self, fluxes):
        self.fluxes = fluxes
        return
    
    def plot_fluxes(self, IRSET_name):
        plt.figure(figsize=(10, 6))
        u_mesh = self.mesh_objects["SHEM295"].lethargyMesh
        for SSH in self.fluxes[IRSET_name].keys():
            flux = self.fluxes[IRSET_name][SSH] 
            u = []
            flux_to_plot = [] 
            for i in range(len(flux)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                flux_to_plot.extend([flux[i],flux[i]])
            plt.step(u, flux_to_plot, where='post', label=f"{SSH} multigroup flux")
        plt.legend()
        plt.yscale('log')
        plt.ylabel('Flux')
        plt.xlabel('Lethargy')
        plt.title(f"Fluxes for {self.case_name}")
        plt.savefig(f"{self.save_path}/{self.case_name}_fluxes_{IRSET_name}.png")
        plt.close()
        return
    
    def plot_XS(self, IRSET_name, reaction_id):
        plt.figure(figsize=(10, 6))
        iso = reaction_id.split("_")[0]
        reaction = reaction_id.split("_")[1]
        u_mesh = self.mesh_objects["SHEM295"].lethargyMesh
        if reaction == "ngamma":
            reaction = "$(n,\gamma)$"
        for SSH in self.reaction_data_pair[reaction_id][IRSET_name].keys():
            XS = self.reaction_data_pair[reaction_id][IRSET_name][SSH]["XS"]
            u = []
            XS_to_plot = [] 
            for i in range(len(XS)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                XS_to_plot.extend([XS[i],XS[i]])
            plt.step(u, XS_to_plot, where='post', label=f"{SSH} {reaction} cross section")
        plt.legend()
        plt.xlabel('Lethargy')
        plt.ylabel(f"$\sigma$ {reaction} (barns)")
        plt.yscale('log')
        plt.title(f"Cross sections for {iso} {reaction}")
        plt.savefig(f"{self.save_path}/{self.case_name}_{reaction_id}_XS_{IRSET_name}.png")
        plt.close()
        return
    
    def plot_reaction_rates(self, IRSET_name, reaction_id):
        plt.figure(figsize=(10, 6))
        iso = reaction_id.split("_")[0]
        reaction = reaction_id.split("_")[1]
        u_mesh = self.mesh_objects["SHEM295"].lethargyMesh
        if reaction == "ngamma":
            reaction = "$(n,\gamma)$"
        for SSH in self.reaction_data_pair[reaction_id][IRSET_name].keys():
            rates = self.reaction_data_pair[reaction_id][IRSET_name][SSH]["rates"]
            u = []
            rates_to_plot = [] 
            for i in range(len(rates)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                rates_to_plot.extend([rates[i],rates[i]])
            plt.step(u, rates_to_plot, where='post', label=f"{SSH} {reaction} reaction rates")
        plt.legend()
        plt.xlabel('Lethargy')
        plt.ylabel(f"{reaction} rates")
        plt.yscale('log')
        plt.grid()
        plt.title(f"Reaction rates for {iso} {reaction}")
        plt.savefig(f"{self.save_path}/{self.case_name}_{reaction_id}_rates_{IRSET_name}.png")
        plt.close()
        return
    
class postTreatment_rates_XS_IRtest:
    def __init__(self, case_name, mesh_objects, self_shielding_methods, save_path):
        self.case_name = case_name
        self.mesh_objects = mesh_objects
        self.self_shielding_methods = self_shielding_methods
        self.save_path = save_path
        self.reaction_data_pair = {}
        
        self.study_type = "IRTEST"
        # dictionaries to store the relative differences between self shielded cross sections, reaction rates and fluxes, identified by reaction_id, mesh_name and SSH method
        self.delta_XS = initialize_nested_dict(["Gd157_ngamma", "U238_ngamma"], ["noIRTEST", "noIRTESTF", "noIRTESTFF", "IRTEST_101", "IRTEST_101F", "IRTEST_101FF", "IRTESTNone", "IRTESTNoneF", "IRTESTNoneFF"], self.self_shielding_methods, default_value=[])
        self.delta_rates = initialize_nested_dict(["Gd157_ngamma", "U238_ngamma"], ["noIRTEST", "noIRTESTF", "noIRTESTFF", "IRTEST_101", "IRTEST_101F", "IRTEST_101FF", "IRTESTNone", "IRTESTNoneF", "IRTESTNoneFF"], self.self_shielding_methods, default_value=[])
        self.delta_Fluxes = {}

        return


class postTreatment_rates_XS_S2:
    def __init__(self, case_name, mesh_objects, S2_libraries, BU_steps_to_treat, save_path):
        self.case_name = case_name
        self.mesh_objects = mesh_objects
        self.libraries = S2_libraries
        self.BU_steps_to_treat = BU_steps_to_treat
        self.rates = {"SHEM281":{}, "SHEM295":{}, "SHEM315":{}}
        self.XS = {"SHEM281":{}, "SHEM295":{}, "SHEM315":{}}
        self.delta_XS = {"Gd157_ngamma":{"SHEM281":[], "SHEM295":[], "SHEM315":[]},
                         "Gd157_abs":{"SHEM281":[], "SHEM295":[], "SHEM315":[]},
                         "U238_ngamma":{"SHEM281":[], "SHEM295":[], "SHEM315":[]}}  # dictionary to store the relative differences between cross sections, identified by reaction_id, mesh_name and library
        self.delta_Rates = {"Gd157_ngamma":{"SHEM281":[], "SHEM295":[], "SHEM315":[]},
                            "Gd157_abs":{"SHEM281":[], "SHEM295":[], "SHEM315":[]},
                            "U238_ngamma":{"SHEM281":[], "SHEM295":[], "SHEM315":[]}} # dictionary to store the relative differences between reaction rates, identified by reaction_id, mesh_name and library
        self.save_path = save_path

        return
    def parse_S2_outputs(self, path_S2):
        """
        Parse the output file from S2 calculation
        """
        #path_S2 = f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/HOM_CELL_study/{self.case_name}/XS_rates_study"
        for bu_step in self.BU_steps_to_treat:
            for library in self.libraries:
                for mesh_name in self.mesh_objects.keys(): 
                    print(f"parsing {mesh_name} {library} {bu_step}")
                    Gd_nGamma_rates, Gd_abs_rates, U238_nGamma_rates, U238_abs_rates = self.parse_Serpent_detector(path_S2, mesh_name, library, bu_step)
                    self.rates[mesh_name][f"Gd157_ngamma_{library}_{bu_step}"] = Gd_nGamma_rates
                    self.rates[mesh_name][f"Gd157_abs_{library}_{bu_step}"] = Gd_abs_rates
                    self.rates[mesh_name][f"U238_ngamma_{library}_{bu_step}"] = U238_nGamma_rates
                    self.rates[mesh_name][f"U238_abs_{library}_{bu_step}"] = U238_abs_rates
                    if bu_step == 0:
                        Gd157_XS_102, Gd157_XS_101, U238_XS_102, U238_XS_101 = self.parse_Serpent_microdepletion(path_S2, mesh_name, library, bu_step)
                        self.N_Gd157 = Gd157_XS_102[0]
                        self.N_U238 = U238_XS_102[0]
                        self.XS[mesh_name][f"Gd157_ngamma_{library}_{bu_step}"] = Gd157_XS_102[1:]
                        self.XS[mesh_name][f"Gd157_abs_{library}_{bu_step}"] = Gd157_XS_101[1:]
                        self.XS[mesh_name][f"U238_ngamma_{library}_{bu_step}"] = U238_XS_102[1:]
                        self.XS[mesh_name][f"U238_abs_{library}_{bu_step}"] = U238_XS_101[1:]
        return
    

    def parse_Serpent_detector(self, path_to_serpent_results, mesh_name, library, bu_step):
        # Read detector file
        if mesh_name == "SHEM295":
            # 26/11/2024 : remember to change so that SHEM295 is not a different case, and detectors are read from _rates_ files
            det = st.read(f"{path_to_serpent_results}/HOM_UOX_Gd157_XS_{library}_mc_det{bu_step}.m")
        else:
            det = st.read(f"{path_to_serpent_results}/HOM_UOX_Gd157_XS_{mesh_name}_{library}_mc_det{bu_step}.m")
        # Get detector names
        Gd_det = det.detectors["Gd_det"].tallies
        n_groups = Gd_det.shape[0]
        n_reactions = Gd_det.shape[1]
        print(f"Number of energy groups: {n_groups}")
        print(f"Number of reactions: {n_reactions}")
        Gd_nGamma_rates = []
        Gd_abs_rates = []
        U238_nGamma_rates = []
        U238_abs_rates = []
        for group in range(n_groups):
            Gd_nGamma_rates.append(Gd_det[group, 0])
            Gd_abs_rates.append(Gd_det[group, 1])
            U238_nGamma_rates.append(Gd_det[group, 2])
            U238_abs_rates.append(Gd_det[group, 3])
        Gd_nGamma_rates_reversed = np.array(Gd_nGamma_rates[::-1])
        Gd_abs_rates_reversed = np.array(Gd_abs_rates[::-1])
        U238_nGamma_rates_reversed = np.array(U238_nGamma_rates[::-1])
        U238_abs_rates_reversed = np.array(U238_abs_rates[::-1])
        return Gd_nGamma_rates_reversed, Gd_abs_rates_reversed, U238_nGamma_rates_reversed, U238_abs_rates_reversed
    
    def parse_Serpent_microdepletion(self, path_to_serpent_results, mesh_name, library, bu_step=0):
        # Read mdep file
        mdep = st.read(f"{path_to_serpent_results}/HOM_UOX_Gd157_XS_{mesh_name}_{library}_mc_mdx{bu_step}.m")
        mdep_scores = mdep.xsVal

        for key in mdep_scores['1'].keys():
            Gd157_XS_102 = mdep_scores['1'][key]
        for key in mdep_scores['2'].keys():
            Gd157_XS_101 = mdep_scores['2'][key]
        for key in mdep_scores['3'].keys():
            U238_XS_102 = mdep_scores['3'][key]
        for key in mdep_scores['4'].keys():
            U238_XS_101 = mdep_scores['4'][key]
        return  Gd157_XS_102, Gd157_XS_101, U238_XS_102, U238_XS_101
    
    def plot_XS(self, mesh_name, reaction_id, bu_step):
        plt.figure(figsize=(10, 6))
        iso = reaction_id.split("_")[0]
        reaction = reaction_id.split("_")[1]
        u_mesh = self.mesh_objects[mesh_name].lethargyMesh
        if reaction == "ngamma":
            reaction = "$(n,\gamma)$"
        for library in self.libraries:
            XS = self.XS[mesh_name][f"{reaction_id}_{library}_{bu_step}"]
            u = []
            XS_to_plot = [] 
            for i in range(len(XS)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                XS_to_plot.extend([XS[i],XS[i]])
            plt.step(u, XS_to_plot, where='post', label=f"{library} {reaction} cross section")
        plt.legend()
        plt.xlabel('Lethargy')
        plt.ylabel(f"$\sigma$ {reaction} (barns)")
        plt.yscale('log')
        plt.grid()
        plt.title(f"Cross sections for {iso} {reaction}")
        plt.savefig(f"{self.save_path}/{self.case_name}_{reaction_id}_XS_{mesh_name}_SERPENT2_{bu_step}.png")
        plt.close()
        return
    
    def plot_reaction_rates(self, mesh_name, reaction_id, bu_step):
        plt.figure(figsize=(10, 6))
        iso = reaction_id.split("_")[0]
        reaction = reaction_id.split("_")[1]
        u_mesh = self.mesh_objects[mesh_name].lethargyMesh
        if reaction == "ngamma":
            reaction = "$(n,\gamma)$"
        for library in self.libraries:
            rates = self.rates[mesh_name][f"{reaction_id}_{library}_{bu_step}"]
            u = []
            rates_to_plot = [] 
            for i in range(len(rates)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                rates_to_plot.extend([rates[i],rates[i]])
            plt.step(u, rates_to_plot, where='post', label=f"{library} {reaction} reaction rates")
        plt.legend()
        plt.xlabel('Lethargy')
        plt.ylabel(f"{reaction} rates")
        plt.yscale('log')
        plt.grid()
        plt.title(f"Reaction rates for {iso} {reaction}")
        plt.savefig(f"{self.save_path}/{self.case_name}_{reaction_id}_rates_{mesh_name}_SERPENT2_{bu_step}.png")
        plt.close()
        return
    
    def compute_relative_differences_XS(self, reaction_id, mesh_name, library, reference_library):
        """
        For a given reaction, energy mesh, compute the relative differences on cross sections between a reference library and the others.
        """
        for bu_step in self.BU_steps_to_treat:
            XS_ref = self.XS[mesh_name][f"{reaction_id}_{reference_library}_{bu_step}"]
            XS = self.XS[mesh_name][f"{reaction_id}_{library}_{bu_step}"]
            delta_XS = (np.array(XS) - np.array(XS_ref))*100/np.array(XS_ref)
            self.delta_XS[reaction_id][mesh_name].append(np.array(delta_XS))
        return
    def compute_relative_differences_Rates(self, reaction_id, mesh_name, library, reference_library):
        """
        For a given reaction, energy mesh, compute the relative differences on reaction rates between a reference library and the others.
        """
        for bu_step in self.BU_steps_to_treat:
            rates_ref = self.rates[mesh_name][f"{reaction_id}_{reference_library}_{bu_step}"]
            rates = self.rates[mesh_name][f"{reaction_id}_{library}_{bu_step}"]
            delta_Rates = (np.array(rates) - np.array(rates_ref))*100/np.array(rates_ref)
            self.delta_Rates[reaction_id][mesh_name].append(np.array(delta_Rates))
        return 
    
    def plot_relative_differences_Rates(self, reaction_id, mesh_name, bu_step):
        """
        Plot relative differences on reaction rates for a given reaction and energy mesh at a given burnup step
        This method is intended to be used after calling compute_relative_differences_rates
        """
        plt.figure(figsize=(10, 6))
        iso = reaction_id.split("_")[0]
        reaction = reaction_id.split("_")[1]
        if reaction == "ngamma":
            reaction = "$(n,\gamma)$"
        # return index of differences list corresponding to the bu step to plot
        idx = self.BU_steps_to_treat.index(bu_step)
        # retrieve corresponding differences
        delta_Rates = self.delta_Rates[reaction_id][mesh_name][idx]
        u_mesh = self.mesh_objects[mesh_name].lethargyMesh
        u = []
        delta_rates_to_plot = [] 
        for i in range(len(delta_Rates)):
            u.extend([u_mesh[i],u_mesh[i+1]])
            delta_rates_to_plot.extend([delta_Rates[i],delta_Rates[i]])
        plt.step(u, delta_rates_to_plot, where='post', label=f"{iso} : $\\Delta \\tau$ {reaction}")
        plt.xlabel('Lethargy')
        plt.ylabel(f"$\\Delta \\tau$ {reaction} (%)")
        plt.legend()
        plt.grid()
        #plt.yscale('log')
        plt.title(f"{iso} : $\\Delta \\tau$ {reaction} (PyNjoy-oldlib) {mesh_name} at BU step {bu_step}")
        plt.savefig(f"{self.save_path}/{self.case_name}_{reaction_id}_delta_Rates_{mesh_name}_SERPENT2_{bu_step}.png")
        plt.close()
        return
    
    def plot_relative_differences_XS(self, reaction_id, mesh_name, bu_step):
        """
        Plot relative differences on cross sections for a given reaction and energy mesh at a given burnup step
        This method is intended to be used after calling compute_relative_differences_XS
        """
        plt.figure(figsize=(10, 6))
        iso = reaction_id.split("_")[0]
        reaction = reaction_id.split("_")[1]
        if reaction == "ngamma":
            reaction = "$(n,\gamma)$"
        # return index of differences list corresponding to the bu step to plot
        idx = self.BU_steps_to_treat.index(bu_step)
        # retrieve corresponding differences
        delta_XS = self.delta_XS[reaction_id][mesh_name][idx]
        u_mesh = self.mesh_objects[mesh_name].lethargyMesh
        u = []
        delta_XS_to_plot = [] 
        for i in range(len(delta_XS)):
            u.extend([u_mesh[i],u_mesh[i+1]])
            delta_XS_to_plot.extend([delta_XS[i],delta_XS[i]])
        plt.step(u, delta_XS_to_plot, where='post', label=f"{iso} : $\\Delta \\sigma$ {reaction}")
        plt.xlabel('Lethargy')
        plt.ylabel(f"$\\Delta \\sigma$ {reaction} (%)")
        plt.legend()
        plt.grid()
        plt.title(f"{iso} : $\\Delta \\sigma$ {reaction} (PyNjoy-oldlib) {mesh_name} at BU step {bu_step}")
        plt.savefig(f"{self.save_path}/{self.case_name}_{reaction_id}_delta_XS_{mesh_name}_SERPENT2_{bu_step}.png")
        plt.close()
        return
    



        