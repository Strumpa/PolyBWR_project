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

class postTreatment_rates_XS_D5:
    def __init__(self, case_name, mesh_objects, self_shielding_methods, save_path):
        self.case_name = case_name
        self.mesh_objects = mesh_objects
        self.self_shielding_methods = self_shielding_methods
        self.save_path = save_path
        self.reaction_data_pair = {}
        # dictionaries to store the relative differences between self shielded cross sections, reaction rates and fluxes, identified by reaction_id, mesh_name and SSH method
        self.delta_XS = {"Gd157_ngamma":{"SHEM281":{"RSE":None, "PT":None, "SUBG":None}, "SHEM295":{"RSE":None, "PT":None, "SUBG":None}, "SHEM315":{"RSE":None, "PT":None, "SUBG":None}},
                         "Gd157_abs":{"SHEM281":{"RSE":None, "PT":None, "SUBG":None}, "SHEM295":{"RSE":None, "PT":None, "SUBG":None}, "SHEM315":{"RSE":None, "PT":None, "SUBG":None}},
                         "U8_ngamma":{"SHEM281":{"RSE":None, "PT":None, "SUBG":None}, "SHEM295":{"RSE":None, "PT":None, "SUBG":None}, "SHEM315":{"RSE":None, "PT":None, "SUBG":None}}}
        self.delta_Rates = {"Gd157_ngamma":{"SHEM281":{"RSE":None, "PT":None, "SUBG":None}, "SHEM295":{"RSE":None, "PT":None, "SUBG":None}, "SHEM315":{"RSE":None, "PT":None, "SUBG":None}},
                         "Gd157_abs":{"SHEM281":{"RSE":None, "PT":None, "SUBG":None}, "SHEM295":{"RSE":None, "PT":None, "SUBG":None}, "SHEM315":{"RSE":None, "PT":None, "SUBG":None}},
                         "U8_ngamma":{"SHEM281":{"RSE":None, "PT":None, "SUBG":None}, "SHEM295":{"RSE":None, "PT":None, "SUBG":None}, "SHEM315":{"RSE":None, "PT":None, "SUBG":None}}}
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
        for SSH in self.self_shielding_methods:
            flux = self.fluxes[mesh_name][SSH] 
            u = []
            flux_to_plot = [] 
            for i in range(len(flux)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                flux_to_plot.extend([flux[i],flux[i]])
            plt.step(u, flux_to_plot, where='post', label=f"{SSH} multigroup flux")
        plt.legend()
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
        for SSH in self.self_shielding_methods:
            XS = self.reaction_data_pair[reaction_id][mesh_name][SSH]["XS"]
            u = []
            XS_to_plot = [] 
            for i in range(len(XS)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                XS_to_plot.extend([XS[i],XS[i]])
            plt.step(u, XS_to_plot, where='post', label=f"{SSH} {reaction} cross section")
        plt.legend()
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
        for SSH in self.self_shielding_methods:
            rates = self.reaction_data_pair[reaction_id][mesh_name][SSH]["rates"]
            u = []
            rates_to_plot = [] 
            for i in range(len(rates)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                rates_to_plot.extend([rates[i],rates[i]])
            plt.step(u, rates_to_plot, where='post', label=f"{SSH} {reaction} reaction rates")
        plt.legend()
        plt.yscale('log')
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
            print(f"delta XS = {delta_XS}")
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
        for SSH in SSH_methods_to_compare:
            save_ssh_id += SSH
        energy_groups = list(range(grmin, grmax + 1))
        plt.figure(figsize=(10, 6))
        for SSH in SSH_methods_to_compare:
            delta_Rates = self.delta_Rates[reaction_id][mesh_name][SSH]
            plt.bar(energy_groups, delta_Rates, edgecolor='black', alpha=0.5, color=colors[SSH], label=f"{iso} {reaction} {SSH} - AUTO")
        plt.ylabel(f'Relative difference on {reaction} rates (%)')
        plt.xticks(range(grmin, grmax+1, 20))  # Set x-ticks with an interval of 20 for readability
        plt.xlabel('Energy group')
        plt.title(f"Relative differences on {iso} {reaction} rates for {mesh_name}")
        plt.legend()
        plt.savefig(f"{self.save_path}/{self.case_name}_{reaction_id}_delta_Rates_{mesh_name}_{save_ssh_id}.png")
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
        for SSH in SSH_methods_to_compare:
            save_ssh_id += SSH
        energy_groups = list(range(grmin, grmax + 1))
        plt.figure(figsize=(10, 6))
        for SSH in SSH_methods_to_compare:
            delta_XS = self.delta_XS[reaction_id][mesh_name][SSH]
            plt.bar(energy_groups, delta_XS, edgecolor='black', alpha=0.5, color=colors[SSH], label=f"{iso} {reaction} {SSH} - AUTO")
        plt.ylabel('Relative difference on XS (%)')
        plt.xticks(range(grmin, grmax+1, 20))  # Set x-ticks with an interval of 20 for readability
        plt.xlabel('Energy group')
        plt.title(f"Relative differences on {iso} $\sigma${reaction} for {mesh_name}")
        plt.legend()
        plt.savefig(f"{self.save_path}/{self.case_name}_{reaction_id}_delta_XS_{mesh_name}_{save_ssh_id}.png")
        plt.close()
        return
    
class postTreatment_rates_XS_S2:
    def __init__(self, case_name, mesh_objects, S2_libraries, BU_steps_to_treat):
        self.case_name = case_name
        self.mesh_objects = mesh_objects
        self.libraries = S2_libraries
        self.BU_steps_to_treat = BU_steps_to_treat
        self.rates = {"SHEM281":{}, "SHEM295":{}, "SHEM315":{}}
        self.XS = {"SHEM281":{}, "SHEM295":{}, "SHEM315":{}}

        return
    def parse_S2_outputs(self):
        """
        Parse the output file from S2 calculation
        """
        path_S2 = f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/HOM_CELL_study/{self.case_name}/XS_study"
        for bu_step in self.BU_steps_to_treat:
            for library in self.libraries:
                for mesh_name in self.mesh_objects.keys():
                    Gd_nGamma_rates, Gd_abs_rates, U8_nGamma_rates, U8_abs_rates = self.parse_Serpent_detector(path_S2, library, bu_step)
                    self.rates[mesh_name][f"Gd157_ngamma_{library}_{bu_step}"] = Gd_nGamma_rates
                    self.rates[mesh_name][f"Gd157_abs_{library}_{bu_step}"] = Gd_abs_rates
                    self.rates[mesh_name][f"U8_ngamma_{library}_{bu_step}"] = U8_nGamma_rates
                    self.rates[mesh_name][f"U8_abs_{library}_{bu_step}"] = U8_abs_rates
                    Gd157_XS_102, Gd157_XS_101, U8_XS_102, U8_XS_101 = self.parse_Serpent_microdepletion(path_S2, library, bu_step)
                    self.XS[mesh_name][f"Gd157_ngamma_{library}_{bu_step}"] = Gd157_XS_102
                    self.XS[mesh_name][f"Gd157_abs_{library}_{bu_step}"] = Gd157_XS_101
                    self.XS[mesh_name][f"U8_ngamma_{library}_{bu_step}"] = U8_XS_102
                    self.XS[mesh_name][f"U8_abs_{library}_{bu_step}"] = U8_XS_101
        return
    

    def parse_Serpent_detector(path_to_serpent_results, library, bu_step):
        # Read detector file
        det = st.read(f"{path_to_serpent_results}/HOM_UOX_Gd157_XS_{library}_mc_det{bu_step}.m")
        # Get detector names
        Gd_det = det.detectors["Gd_det"].tallies
        print(Gd_det.shape)
        n_groups = Gd_det.shape[0]
        n_reactions = Gd_det.shape[1]
        print(f"Number of energy groups: {n_groups}")
        print(f"Number of reactions: {n_reactions}")
        Gd_nGamma_rates = []
        Gd_abs_rates = []
        U8_nGamma_rates = []
        U8_abs_rates = []
        for group in range(n_groups):
            Gd_nGamma_rates.append(Gd_det[group, 0])
            Gd_abs_rates.append(Gd_det[group, 1])
            U8_nGamma_rates.append(Gd_det[group, 2])
            U8_abs_rates.append(Gd_det[group, 3])
        Gd_nGamma_rates_reversed = np.array(Gd_nGamma_rates[::-1])
        Gd_abs_rates_reversed = np.array(Gd_abs_rates[::-1])
        U8_nGamma_rates_reversed = np.array(U8_nGamma_rates[::-1])
        U8_abs_rates_reversed = np.array(U8_abs_rates[::-1])
        return Gd_nGamma_rates_reversed, Gd_abs_rates_reversed, U8_nGamma_rates_reversed, U8_abs_rates_reversed
    
    def parse_Serpent_microdepletion(path_to_serpent_results, library, bu_step):
        # Read mdep file
        mdep = st.read(f"{path_to_serpent_results}/HOM_UOX_Gd157_XS_{library}_mc_mdx{bu_step}.m")

        mdep_scores = mdep.xsVal

        for key in mdep_scores['1'].keys():
            Gd157_XS_102 = mdep_scores['1'][key]
        for key in mdep_scores['2'].keys():
            Gd157_XS_101 = mdep_scores['2'][key]
        for key in mdep_scores['3'].keys():
            U8_XS_102 = mdep_scores['3'][key]
        for key in mdep_scores['4'].keys():
            U8_XS_101 = mdep_scores['4'][key]
        return  Gd157_XS_102, Gd157_XS_101, U8_XS_102, U8_XS_101


        