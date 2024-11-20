# Python3 class to compare DRAGON5 to SERPENT2 results 
# Author: R. Guasch
# Date: 2024-11-18
# Purpose: Post-treat results for HOM_UOX_Gd155/157 calculations : compare DRAGON5 to SERPENT2 results
#        - compare XS and rates from DRAGON5 to SERPENT2 results
#        - compare iso dens and keffs from DRAGON5 to SERPENT2 results

import numpy as np
import matplotlib.pyplot as plt
import os
import sys


class compare_D5_S2_rates_XS:
    def __init__(self, case_name, D5_case, S2_case, save_path):
        self.case_name = case_name
        self.D5_case = D5_case
        self.S2_case = S2_case
        self.S2_lib = ["oldlib","PyNjoy2016"]
        
        self.save_path = save_path
        self.D5_S2_rates_diff = {"oldlib":{"SHEM281":{"U238":{"ngamma":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None},"abs":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None}},
                                                 "Gd157":{"ngamma":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None},"abs":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None}}},
                                    "SHEM295":{"U238":{"ngamma":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None},"abs":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None}},
                                            "Gd157":{"ngamma":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None},"abs":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None}}},
                                    "SHEM315":{"U238":{"ngamma":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None},"abs":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None}},
                                            "Gd157":{"ngamma":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None},"abs":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None}}}},
                            "PyNjoy2016":{"SHEM281":{"U238":{"ngamma":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None},"abs":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None}},
                                                 "Gd157":{"ngamma":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None},"abs":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None}}},
                                    "SHEM295":{"U238":{"ngamma":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None},"abs":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None}},
                                            "Gd157":{"ngamma":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None},"abs":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None}}},
                                    "SHEM315":{"U238":{"ngamma":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None},"abs":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None}},
                                            "Gd157":{"ngamma":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None},"abs":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None}}}}}
        self.D5_S2_XS_diff =  {"oldlib":{"SHEM281":{"U238":{"ngamma":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None},"abs":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None}},
                                                 "Gd157":{"ngamma":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None},"abs":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None}}},
                                    "SHEM295":{"U238":{"ngamma":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None},"abs":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None}},
                                            "Gd157":{"ngamma":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None},"abs":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None}}},
                                    "SHEM315":{"U238":{"ngamma":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None},"abs":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None}},
                                            "Gd157":{"ngamma":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None},"abs":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None}}}},
                            "PyNjoy2016":{"SHEM281":{"U238":{"ngamma":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None},"abs":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None}},
                                                 "Gd157":{"ngamma":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None},"abs":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None}}},
                                    "SHEM295":{"U238":{"ngamma":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None},"abs":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None}},
                                            "Gd157":{"ngamma":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None},"abs":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None}}},
                                    "SHEM315":{"U238":{"ngamma":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None},"abs":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None}},
                                            "Gd157":{"ngamma":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None},"abs":{"Autosecol":None, "RSE": None, "PT":None, "SUBG":None}}}}}
        return            
    
    def compare_reaction_rates_and_XS(self, reaction_id, mesh, SSH_methods):
        iso = reaction_id.split("_")[0]
        reaction = reaction_id.split("_")[1]
        for library in self.S2_lib:
            # for now only consider BU step = 0
            for SSH in SSH_methods:
                self.D5_S2_rates_diff[library][mesh][iso][reaction][SSH] = (self.D5_case.reaction_data_pair[reaction_id][mesh][SSH]["rates"] - self.S2_case.rates[mesh][f"{iso}_{reaction}_{library}_0"])*100/self.S2_case.rates[mesh][f"{iso}_{reaction}_{library}_0"]
                self.D5_S2_XS_diff[library][mesh][iso][reaction][SSH] = (self.D5_case.reaction_data_pair[reaction_id][mesh][SSH]["XS"] - self.S2_case.XS[mesh][f"{iso}_{reaction}_{library}_0"])*100/self.S2_case.XS[mesh][f"{iso}_{reaction}_{library}_0"]
        return
    
    def plot_XS_D5_S2(self, reaction_id, mesh, SSH_methods, bu_step = 0):
        iso = reaction_id.split("_")[0]
        reaction = reaction_id.split("_")[1]
        if reaction == "ngamma":
            reaction_print = "$(n,\\gamma)$"
        u_mesh = self.D5_case.mesh_objects[mesh].lethargyMesh
        plt.figure(figsize=(10, 6))
        lib_id = ""
        for library in self.S2_lib:
            lib_id += library
            u = []
            S2_XS_to_plot = []
            S2_XS = self.S2_case.XS[mesh][f"{iso}_{reaction}_{library}_{bu_step}"]
            for i in range(len(S2_XS)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                S2_XS_to_plot.extend([S2_XS[i],S2_XS[i]])
            plt.step(u, S2_XS_to_plot, where='post', label=f"S2: $\\sigma$ {reaction_print} {library} for {iso}")
        ssh_id = "" 
        for SSH in SSH_methods:    
            ssh_id += SSH
            u = []   
            D5_XS_to_plot = []     
            D5_XS = self.D5_case.reaction_data_pair[reaction_id][mesh][SSH]["XS"]
            for i in range(len(D5_XS)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                D5_XS_to_plot.extend([D5_XS[i],D5_XS[i]])
            plt.step(u, D5_XS_to_plot, where='post', label=f"{SSH}: $\\sigma$ {reaction_print} for {iso}")
        plt.xlabel("Lethargy")
        plt.ylabel("Cross section [barn]")
        plt.yscale("log")
        plt.title(f"Comparison of {iso}_{reaction} cross sections between DRAGON5 and SERPENT2")
        plt.grid()
        plt.legend()
        plt.savefig(f"{self.save_path}/{iso}_{reaction}_XS_D5_{ssh_id}_S2_{lib_id}.png")
        plt.close()
        return
    
    def plot_rates_D5_S2(self, reaction_id, mesh, SSH_methods, bu_step = 0):
        iso = reaction_id.split("_")[0]
        reaction = reaction_id.split("_")[1]
        if reaction == "ngamma":
            reaction_print = "$(n,\\gamma)$"
        u_mesh = self.D5_case.mesh_objects[mesh].lethargyMesh
        plt.figure(figsize=(10, 6))
        lib_id = ""
        for library in self.S2_lib:
            lib_id += library
            u = []
            S2_rates_to_plot = []
            S2_rates = self.S2_case.rates[mesh][f"{iso}_{reaction}_{library}_{bu_step}"]
            for i in range(len(S2_rates)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                S2_rates_to_plot.extend([S2_rates[i],S2_rates[i]])
            plt.step(u, S2_rates_to_plot, where='post', label=f"S2: $\\tau$ {reaction_print} {library} for {iso}")
        ssh_id = ""
        for SSH in SSH_methods:   
            ssh_id += SSH 
            u = []   
            D5_rates_to_plot = []     
            D5_rates = self.D5_case.reaction_data_pair[reaction_id][mesh][SSH]["rates"]
            for i in range(len(D5_rates)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                D5_rates_to_plot.extend([D5_rates[i],D5_rates[i]])
            plt.step(u, D5_rates_to_plot, where='post', label=f"{SSH}: $\\tau$ {reaction_print} for {iso}")
        plt.xlabel("Lethargy")
        plt.ylabel(f"$\\tau$ {reaction_print} {iso}")
        plt.yscale("log")
        plt.title(f"Comparison of {iso} {reaction_print} cross sections between DRAGON5 and SERPENT2")
        plt.grid()
        plt.legend()
        plt.savefig(f"{self.save_path}/{iso}_{reaction}_rates_D5_{ssh_id}_S2_{lib_id}.png")
        plt.close()
        return
    
    def plot_diff_rates_D5_S2(self, reaction_id, mesh_name, library, SSH_methods):
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
            D5_S2_rates_diff = self.D5_S2_rates_diff[library][mesh_name][iso][reaction][SSH]
            for i in range(len(D5_S2_rates_diff)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                D5_S2_rates_diff_to_plot.extend([D5_S2_rates_diff[i],D5_S2_rates_diff[i]])
            plt.step(u, D5_S2_rates_diff_to_plot, label=f"{SSH} - {library} : {iso} {reaction_print}")
        plt.xlabel("Lethargy")
        plt.ylabel(f"$\\Delta \\tau$ {iso} {reaction_print} [%]")
        plt.title(f"Relative difference on $\\tau$ {reaction_print} {iso} between DRAGON5 and SERPENT2 ({library}) on {mesh_name}")
        plt.grid()
        plt.legend()
        plt.savefig(f"{self.save_path}/{iso}_{reaction}_rates_diff_{mesh_name}_{library}_{ssh_id}.png")
        plt.close()
        return
    
    def plot_diff_XS_D5_S2(self, reaction_id, mesh_name, library, SSH_methods):
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
            D5_S2_XS_diff = self.D5_S2_XS_diff[library][mesh_name][iso][reaction][SSH]
            for i in range(len(D5_S2_XS_diff)):
                u.extend([u_mesh[i],u_mesh[i+1]])
                D5_S2_XS_diff_to_plot.extend([D5_S2_XS_diff[i],D5_S2_XS_diff[i]])
            plt.step(u, D5_S2_XS_diff_to_plot, label=f"{SSH} - {library} : {iso} {reaction_print}")
        plt.xlabel("Lethargy")
        plt.ylabel(f"$\\Delta \\sigma$ {reaction_print} {iso}  [%]")
        plt.title(f"Relative difference on $\\sigma$ {reaction_print} {iso} between DRAGON5 and SERPENT2 ({library}) on {mesh_name}")
        plt.grid()
        plt.legend()
        plt.savefig(f"{self.save_path}/{iso}_{reaction}_XS_diff_{mesh_name}_{library}_{ssh_id}.png")
        plt.close()
        return
    
    def renorm_rates(self, reaction_id, mesh_name):
        """
        Since DRAGON5 and SERPENT2 do not use the same normalization for fluxes, 
        this function renormalizes the rates in order to allow for a direct D5-S2 comparison.
        """
        # Renormalize DRAGON5 and Serpent2 rates to 1
        # For all Serpent2 libraries and DRAGON5 SSH methods
        #D5_renorm = self.D5_case.reaction_data_pair[reaction_id]["SHEM281"]["Autosecol"]["rates"][0]/self.S2_case.rates["SHEM281"][f"{reaction_id}_oldlib_0"][0]
        for SSH in self.D5_case.reaction_data_pair[reaction_id][mesh_name].keys():
            self.D5_case.reaction_data_pair[reaction_id][mesh_name][SSH]["rates"] = self.D5_case.reaction_data_pair[reaction_id][mesh_name][SSH]["rates"]/np.sum(self.D5_case.reaction_data_pair[reaction_id][mesh_name][SSH]["rates"])
        for library in self.S2_lib:
            self.S2_case.rates[mesh_name][f"{reaction_id}_{library}_0"] = self.S2_case.rates[mesh_name][f"{reaction_id}_{library}_0"]/np.sum(self.S2_case.rates[mesh_name][f"{reaction_id}_{library}_0"])
        return