# Python3 post-treatment script for BWR 2x2 case
# Author: R. Guasch
# Date: 1st August 2024
# Purpose: Post-treatment for global (Keff, Ni) and local (Pin-wise fission rates) results.

import os
import shutil # move files or folders
import lifo
import lcm
import cle2000
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import serpentTools as st
from serpentTools.settings import rc
from getLists import *

def POSTPROC_2x2(pyCOMPO, ListCOMPO, listSSH, BU_points, case_name, ssh_method, ssh_sol, flx_sol, CALC_opt):
    BWR_2x2_case_instance = BWR_2x2_case(pyCOMPO, ListCOMPO, listSSH, BU_points, case_name, ssh_method, ssh_sol, flx_sol, CALC_opt)
    return BWR_2x2_case_instance

class BWR_2x2_case:
    def __init__(self, pyCOMPO, ListCOMPO, listSSH, BU_points, case_name, ssh_method, ssh_sol, flx_sol, CALC_opt):
        self.pyCOMPO = pyCOMPO
        self.ListCOMPO = ListCOMPO
        self.listSSH = listSSH
        self.BU_points = BU_points
        self.case_name = case_name
        self.ssh_method = ssh_method
        self.ssh_sol = ssh_sol
        self.flx_sol = flx_sol
        self.CALC_opt = CALC_opt

        self.path=os.getcwd()

        self.SAVE_DIR=f'BWRresults_PyGan_{case_name}/{BU_points}_postprocess/'

        a=os.path.exists(self.SAVE_DIR)
        if a==False:
            os.makedirs(self.SAVE_DIR)

        self.isotopes_to_compare = ['U235','U236','U238','Pu239','Pu240','Pu241','Pu242','Gd155','Gd157','Xe135','Sm149']
        self.fissionable_isotopes = ['U235 ', 'U238 ', 'Pu239', 'Pu241']
        self.nCELL = 4
        # Load Serpent2 reference case
        self.load_serpent2_reference_case()

        # Load DRAGON5 results
        self.load_dragon5_results()

        self.plot_COMPARISON()
        self.compute_errors()
        self.plot_ERRORS()

        return

    def load_serpent2_reference_case(self):
        print("Loading Serpent2 reference case")
        if self.case_name == "2x2_UOX":
            # Load Serpent2 reference case : Keff results
            serpent2_case = st.read("/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/PyNjoy2016_results/AT10_2x2/AT10_2x2_UOX_mc_res.m")
            # Load Serpent2 reference case : depletion data
            self.serpent_dep = st.read("/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/PyNjoy2016_results/AT10_2x2/AT10_2x2_UOX_mc_dep.m")

        elif self.case_name == "2x2_UOX_Gd":
            # Load Serpent2 reference case : Keff results
            serpent2_case = st.read("/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/PyNjoy2016_results/AT10_2x2/AT10_2x2_UOX_Gd_mc_res.m")
            # Load Serpent2 reference case : depletion data
            self.serpent_dep = st.read("/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/PyNjoy2016_results/AT10_2x2/AT10_2x2_UOX_Gd_mc_dep.m")


        serpent_keff=serpent2_case.resdata["absKeff"]
            
        self.fuel=self.serpent_dep.materials['total']
        serpent_BU=self.fuel.burnup
        np.savetxt(f'serpent_BU_{self.case_name}.txt',serpent_BU)
        self.serpent_BU=np.loadtxt(f'serpent_BU_{self.case_name}.txt',dtype=float)

        serpent_iso_dens = self.fuel.toDataFrame("adens",names=self.isotopes_to_compare)
        np.savetxt(f'serpent_ISOTOPESDENS_{self.case_name}.txt',serpent_iso_dens)
        SERPENT_ISOTOPESDENS=np.loadtxt(f'serpent_ISOTOPESDENS_{self.case_name}.txt',dtype=float)
        serpent_iso_dens = np.transpose(SERPENT_ISOTOPESDENS)

        self.Serpent_Ni = [ serpent_iso_dens[0,:],
                            serpent_iso_dens[1,:],
                            serpent_iso_dens[2,:],
                            serpent_iso_dens[3,:],
                            serpent_iso_dens[4,:],
                            serpent_iso_dens[5,:],
                            serpent_iso_dens[6,:],
                            serpent_iso_dens[7,:],
                            serpent_iso_dens[8,:],
                            serpent_iso_dens[9,:],
                            serpent_iso_dens[10,:]
                        ]

        self.renomalize_serpent_BU(serpent_keff)
        

        return

    def renomalize_serpent_BU(self, serpent_keff):
        # Modification unite BU pour match avec DRAGON
        self.serpent_keff=np.zeros(len(self.serpent_BU))    
        for k in range(len(self.serpent_BU)):
            self.serpent_BU[k]=1000*self.serpent_BU[k]
            self.serpent_keff[k]=serpent_keff[k][0]



    def load_dragon5_results(self):
        print("Loading DRAGON5 results")
        # Load DRAGON5 results
        self.lenBU_DRAGON=np.shape(self.ListCOMPO)[0]
         
        ISOTOPES=self.pyCOMPO['EDIBU_HOM']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESDENS']
        print(f"Dragon isotopes = {ISOTOPES}")
        self.lenISOT_DRAGON=np.shape(ISOTOPES)[0]

        self.DRAGON_BU=self.ListCOMPO
        DRAGON_ISOTOPESDENS=np.zeros((self.lenISOT_DRAGON,self.lenBU_DRAGON))
        self.DRAGON_Keff=np.zeros(self.lenBU_DRAGON)

        #print("$$$ ---------------- DRAGON_ISOTOPESDENS shape = ",lenISOT_DRAGON,lenBU_DRAGON)

        for k in range(self.lenBU_DRAGON):
            self.DRAGON_Keff[k]=self.pyCOMPO['EDIBU_HOM']['MIXTURES'][0]['CALCULATIONS'][k]['K-EFFECTIVE']
            for j in range(self.lenISOT_DRAGON):
                print("$$$ ---------------- ISOTOPES ",j,"/",self.lenISOT_DRAGON," = ",self.pyCOMPO['EDIBU']['MIXTURES'][0]['CALCULATIONS'][k]['ISOTOPESDENS'][j])
                DRAGON_ISOTOPESDENS[j][k]=self.pyCOMPO['EDIBU_HOM']['MIXTURES'][0]['CALCULATIONS'][k]['ISOTOPESDENS'][j]

        print('$$$ ---------------- DRAGON_BU =',self.DRAGON_BU)
        print('$$$ ---------------- LEN_DRAGON_BU =',len(self.DRAGON_BU))
        print("$$$ ---------------- DRAGON_Keff = ",self.DRAGON_Keff)    
        print("$$$ ---------------- DRAGON_ISOTOPESDENS = ",DRAGON_ISOTOPESDENS)

        # --------- Liste des isotopes recuperes dans la Multicompo
        isotopes2=[]
        isotopes=[]
        for k in range(self.lenISOT_DRAGON):
            isotopes2=isotopes2+[self.pyCOMPO["EDIBU_HOM"]['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][k]['ALIAS']]
        for k in range(self.lenISOT_DRAGON):
            if isotopes2[k][0]=='U':
                isotopes=isotopes+[isotopes2[k][0:4]]
            else:
                isotopes=isotopes+[isotopes2[k][0:5]]

        indices=np.zeros(len(self.isotopes_to_compare))
        for n in range(len(self.isotopes_to_compare)):
            for m in range(len(isotopes)):
                if self.isotopes_to_compare[n]==isotopes[m]:
                    indices[n]=m

        print("$$$ ---------------- DRAGON isotopes = ",isotopes)
        print("$$$ ---------------- DRAGON isotopes souhaites = ",self.isotopes_to_compare)
        print("$$$ ---------------- indices correspondances = ",indices)
            
        # Store in DRAGON_ALL
        self.DRAGON_Ni=[
               DRAGON_ISOTOPESDENS[int(indices[0]),:],
               DRAGON_ISOTOPESDENS[int(indices[1]),:],
               DRAGON_ISOTOPESDENS[int(indices[2]),:],
               DRAGON_ISOTOPESDENS[int(indices[3]),:],
               DRAGON_ISOTOPESDENS[int(indices[4]),:],
               DRAGON_ISOTOPESDENS[int(indices[5]),:],
               DRAGON_ISOTOPESDENS[int(indices[6]),:],
               DRAGON_ISOTOPESDENS[int(indices[7]),:],
               DRAGON_ISOTOPESDENS[int(indices[8]),:],
               DRAGON_ISOTOPESDENS[int(indices[9]),:],
               DRAGON_ISOTOPESDENS[int(indices[10]),:]
           ]
        return

    
    def plot_COMPARISON(self):
        # Plot Keff
        plt.figure()
        print(f"$$$ ---------------- serpent_BU = {self.serpent_BU}")
        print(f"$$$ ---------------- serpent_keff = {self.serpent_keff}")
        print(f"$$$ ---------------- DRAGON_BU = {self.DRAGON_BU}")
        print(f"$$$ ---------------- DRAGON_Keff = {self.DRAGON_Keff}")
        plt.plot(self.serpent_BU, self.serpent_keff, '2-',label='Serpent2', linewidth=1)
        plt.plot(self.DRAGON_BU, self.DRAGON_Keff, '2-',label='DRAGON5', linewidth=1)
        plt.xlabel('Burnup [MWd/t]')
        plt.ylabel('Keff')
        plt.legend()
        plt.grid()

        save_name=f'{self.case_name}_{self.BU_points}_{self.CALC_opt}_{self.ssh_method}_COMPARE_Keff'
        fig_name=f'{self.case_name} {self.BU_points} {self.CALC_opt} - Keff'

        plt.title(fig_name)
        os.chdir(self.path+'/'+self.SAVE_DIR)
        plt.savefig(save_name+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
        os.chdir(self.path)
        plt.close('all')

        # Plot isotopes
        for i in range(len(self.isotopes_to_compare)):
            plt.figure()
            plt.plot(self.serpent_BU, self.Serpent_Ni[i], label='Serpent2', )
            plt.plot(self.DRAGON_BU, self.DRAGON_Ni[i], label='DRAGON5')
            plt.xlabel('Burnup [MWd/t]')
            plt.ylabel('Density [atom/b-cm]')
            plt.legend()
            plt.grid()
            save_name = f'{self.isotopes_to_compare[i]}_{self.case_name}_{self.BU_points}_{self.CALC_opt}_{self.ssh_method}_COMPARE'
            fig_name = f'{self.case_name} {self.BU_points} {self.CALC_opt} {self.isotopes_to_compare[i]} - Density'
            plt.title(fig_name)
            
            os.chdir(self.path+'/'+self.SAVE_DIR)
            plt.savefig(save_name+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
            os.chdir(self.path)
            plt.close('all')

        return
    
    def compute_errors(self):
        # Compute errors
        self.err_Keff = np.zeros(self.lenBU_DRAGON)
        self.err_Ni = np.zeros((len(self.isotopes_to_compare), self.lenBU_DRAGON))

        for k in range(self.lenBU_DRAGON):
            self.err_Keff[k] = (self.DRAGON_Keff[k] - self.serpent_keff[k])*1e5
            for j in range(len(self.isotopes_to_compare)):
                if self.Serpent_Ni[j][k] == 0:
                    self.err_Ni[j][k] = 0
                else:
                    self.err_Ni[j][k] = (self.DRAGON_Ni[j][k] - self.Serpent_Ni[j][k])*100 / self.Serpent_Ni[j][k]

        return
    
    def plot_ERRORS(self):
        # Plot Keff error
        plt.figure()
        plt.plot(self.DRAGON_BU, self.err_Keff, '2-',linewidth=1)
        plt.plot([0,self.ListCOMPO[-1]],[300,300],'r-.') # limite +300pcm
        plt.plot([0,self.ListCOMPO[-1]],[-300,-300],'r-.') # limite -300pcm
        plt.xlabel('Burnup [MWd/t]')
        plt.ylabel('Keff error [pcm]')
        plt.grid()
        save_name = f'{self.case_name}_{self.BU_points}_{self.CALC_opt}_{self.ssh_method}_ERROR_Keff'
        fig_name = f'{self.case_name} {self.BU_points} {self.CALC_opt} - Keff error'
        plt.title(fig_name)
        os.chdir(self.path+'/'+self.SAVE_DIR)
        plt.savefig(save_name+'.png',bbox_inches = 'tight')
        os.chdir(self.path)
        plt.close('all')

        # Plot isotopes errors
        for i in range(len(self.isotopes_to_compare)):
            plt.figure()
            plt.plot(self.DRAGON_BU, self.err_Ni[i])
            plt.plot([0,self.ListCOMPO[-1]],[2,2],'r-.') # limite +2%
            plt.plot([0,self.ListCOMPO[-1]],[-2,-2],'r-.') # limite -2%
            plt.xlabel('Burnup [MWd/t]')
            plt.ylabel(f'{self.isotopes_to_compare[i]} relative error [%]')
            plt.grid()
            save_name = f'{self.isotopes_to_compare[i]}_{self.case_name}_{self.BU_points}_{self.CALC_opt}_{self.ssh_method}_ERROR'
            fig_name = f'{self.case_name} {self.BU_points} {self.CALC_opt} {self.isotopes_to_compare[i]} - Density error'
            plt.title(fig_name)
            os.chdir(self.path+'/'+self.SAVE_DIR)
            plt.savefig(save_name+'.png',bbox_inches = 'tight')
            os.chdir(self.path)
            plt.close('all')

        return
    
    def load_Serpent2_fission_rates(self):
        if self.case_name == "2x2_UOX":
            detFile_0 = '/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/PyNjoy2016_results/A10_2x2/AT10_2x2_UOX_mc_res_det0.m'
            detFile_48 = '/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/PyNjoy2016_results/A10_2x2/AT10_2x2_UOX_mc_res_det48.m'
        elif self.case_name == "2x2_UOX_Gd":
            detFile_0 = '/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/PyNjoy2016_results/A10_2x2/AT10_2x2_UOX_Gd_mc_res_det0.m'
            detFile_48 = '/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/PyNjoy2016_results/A10_2x2/AT10_2x2_UOX_Gd_mc_res_det48.m'
        self.det_BU0 = st.read(detFile_0, reader = 'det').detectors
        self.det_BU_last = st.read(detFile_48, reader = 'det').detectors

    def compute_D5_production(self, BU):
        self.prodD5 = 0.0
        for isotope in self.fissionable_isotopes:
            for mix in range(self.nCELL):
                NWT0 = self.pyCOMPO['EDIBU']['MIXTURES'][mix]['CALCULATIONS'][BU]['ISOTOPESLIST'][isotope]['NWT0']
                N = self.pyCOMPO['EDIBU']['MIXTURES'][mix]['CALCULATIONS'][BU]['ISOTOPESDENS'][isotope]
                vol = self.pyCOMPO['EDIBU']['MIXTURES'][mix]['CALCULATIONS'][BU]['ISOTOPESVOL'][isotope]
                NFTOT = self.pyCOMPO['EDIBU']['MIXTURES'][mix]['CALCULATIONS'][BU]['ISOTOPESLIST'][isotope]['NFTOT']
#                NGAMMA = pyCOMPO['EDIBU_2gr']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NG']
                for gr in range(2):
#                    prodD5 += 2*NWT0[gr]*(NFTOT[gr]+NGAMMA[gr])*N*vol
                    self.prodD5 += 2*NWT0[gr]*NFTOT[gr]*N*vol

    def compute_S2_production(self,BU):
        if BU == 0:
            detFile = self.det_BU0
        elif BU == -1:
            detFile = self.det_BU_last
        self.prodS2 = 0.0
        for groupS2 in range(2):
            for n_iso, isotope in enumerate(self.fissionable_isotopes[0:4]):
                for n,name_cell in enumerate(detFile):
                    print(f"$$$ ---------------- n = {n}")



                    




