#####################################################################
#                                                                   #
# Description : PyGan script for evo on homogeneous cell            #
# Author      : R. Guasch adapted for BWR_project                   #
# Date        : 2024                                                #
# Purpose     : Post-processing for global values (Keff,Ni)         #
#               To compare old serpent2 lib and PyNjoy2016          #
#                                                                   #
#####################################################################
#
import os
import shutil # move files or folders
import lifo
import lcm
import cle2000
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import serpentTools
from serpentTools.settings import rc

def plot_DRAGON(DRAGON_results, isotopes_SOUHAITES, SIZE, path, SAVE_DIR, ssh_method, correlation_name, energy_mesh):

    if correlation_name == "_CORR":
        CORR = " CORR " 
    elif correlation_name == "_noCORR":
        CORR = " noCORR "
    else:
        CORR = " "
        correlation_name = ""
    print('$$$ -------- POSTPROC.py : DRAGON5 figures ')
    for k in range(len(isotopes_SOUHAITES)+1):
        # UOX tests :
        legends=[]
        plt.figure(figsize=SIZE)
        for test_name in DRAGON_results.keys():
            plt.plot(DRAGON_results[test_name][0], DRAGON_results[test_name][k+1],'2-',linewidth=1)
            legends.append(f"{test_name} : D5")
        plt.xlabel('BU [MWd/t]')
        plt.grid()
        plt.legend(legends)

        if k == 0: # Comparaison des Keff
            plt.ylabel('Keff')
            save_name=f'HOM_UOX_Gd157_{energy_mesh}_{ssh_method}{correlation_name}_DRAGON5_Keff_PCC'
            fig_name=f'HOM UOX + Gd157 {CORR}{ssh_method} - Keff on {energy_mesh}'
        else : # Erreur sur isotopes
            plt.ylabel('Isotopic density [atom/b-cm]')
            save_name=f'HOM_UOX_Gd157_{energy_mesh}_{ssh_method}{correlation_name}_{isotopes_SOUHAITES[k-1]}_DRAGON5_PCC'
            fig_name=f'HOM UOX + Gd157 {CORR} {ssh_method} - {isotopes_SOUHAITES[k-1]} on {energy_mesh}'

        plt.title(fig_name)
        os.chdir(path+'/'+SAVE_DIR)
        plt.savefig(save_name+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
        os.chdir(path)
        plt.close('all')

    return

def plot_Serpent2(SERPENT_results, isotopes_SOUHAITES, SIZE, path, SAVE_DIR):
    print('$$$ -------- POSTPROC.py : Serpent2 figures ')

    for k in range(len(isotopes_SOUHAITES)+1):
        legends = []
        plt.figure(figsize=SIZE)
        for test_name in SERPENT_results.keys():
            lib_name = test_name.split('_')[0]
            pcc = test_name.split('_')[1]
            legends.append(f"{lib_name} : {pcc}")
            plt.plot(SERPENT_results[test_name][0],SERPENT_results[test_name][k+1],'2-',linewidth=1)
        plt.xlabel('BU [MWd/t]')
        plt.grid()
            
        plt.legend(legends)
        if k == 0: # Keff
            plt.ylabel('Keff')
            save_name=f'HOM_UOX_Gd157_Serpent2_Keff_PCC'
            fig_name=f'HOM CELL : Gd157 test - Keff Serpent2'
        else : # isotopes
            plt.ylabel('Isotopic density [atom/b-cm]')
            save_name=f'HOM_UOX_Gd157_Serpent2_{isotopes_SOUHAITES[k-1]}_PCC'
            fig_name=f'HOM_CELL : Gd157 test  - {isotopes_SOUHAITES[k-1]}'

        plt.title(fig_name)
        os.chdir(path+'/'+SAVE_DIR)
        plt.savefig(save_name+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
        os.chdir(path)
        plt.close('all')

    return

def plot_Comparison(DRAGON_results, SERPENT_results, isotopes_SOUHAITES, SIZE, path, SAVE_DIR, ssh_method, correlation_name, energy_mesh):

    if correlation_name == "_CORR":
        CORR = " CORR " 
    elif correlation_name == "_noCORR":
        CORR = " noCORR "
    else:
        CORR = " "
        correlation_name = ""
    print('$$$ -------- POSTPROC.py : Comparison DRAGON5 / Serpent2 figures ')

    for k in range(len(isotopes_SOUHAITES)+1):
        fig, ax = plt.subplots(figsize=SIZE)
        for test_name in DRAGON_results.keys():
            ax.plot(DRAGON_results[test_name][0],DRAGON_results[test_name][k+1],'2-',linewidth=1, label=f"{test_name} : D5")
        for test_name in SERPENT_results.keys():
            lib_name = test_name.split('_')[0]
            pcc = test_name.split('_')[1]        
            ax.plot(SERPENT_results[test_name][0],SERPENT_results[test_name][k+1],'2-',linewidth=1, label=f"{lib_name} {pcc} : S2")


        if k == 0: # Comparaison des Keff
            ax.set_ylabel('Keff')
            save_name=f'HOM_UOX_Gd157_{energy_mesh}_{ssh_method}{correlation_name}_COMP_Keff_PCC'
            fig_name=f'HOM_UOX_Gd157 edep/PCC {ssh_method}{CORR}- Keff on {energy_mesh}'
        else : # comparaison isotopes
            ax.set_ylabel('Isotopic density [atom/b-cm]')
            save_name=f'HOM_UOX_Gd157_{energy_mesh}_{ssh_method}{correlation_name}_COMP_{isotopes_SOUHAITES[k-1]}'
            fig_name=f'HOM_UOX_Gd157 edep/PCC {ssh_method}{CORR}- {isotopes_SOUHAITES[k-1]} on {energy_mesh}'

        ax.set_xlabel('BU [MWd/t]')
        ax.grid()
        ax.legend(loc="best")
        ax.set_title(fig_name)
        os.chdir(path+'/'+SAVE_DIR)
        fig.savefig(save_name+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
        os.chdir(path)
        plt.close('all')

    return

def plot_errors(D5_test_name, ERRORS, isotopes_SOUHAITES, SIZE, path, SAVE_DIR, ssh_method, correlation_name, ListAUTOP, ListeCOMPO, energy_mesh):

    if correlation_name == "_CORR":
        CORR = " CORR " 
    elif correlation_name == "_noCORR":
        CORR = " noCORR "
    else:
        CORR = " "
        correlation_name = ""

    print('$$$ -------- POSTPROC.py : ERROR DRAGON5-Serpent2 figures ')

    for k in range(len(isotopes_SOUHAITES)+1):
        fig,ax = plt.subplots(figsize=SIZE)
        for test_name in ERRORS.keys():
            #D5_test_name = test_name.split('_')[0]+"_"+test_name.split('_')[1]+"_"+test_name.split('_')[2]
            ax.plot(ERRORS[test_name][0],ERRORS[test_name][k+1],'2-',linewidth=1, label=f"{test_name}")
            ax.set_xlabel('BU [MWd/t]')
            ax.legend(loc="best")
            
        if k == 0: # Erreur sur Keff
            for step in ListAUTOP:
                if step <= ListeCOMPO[-1]:
                    ax.axvline(x=step,marker='o',color='red')
            ax.plot([0,ListeCOMPO[-1]],[300,300],'r-.') # limite +300pcm
            ax.plot([0,ListeCOMPO[-1]],[-300,-300],'r-.') # limite -300pcm
            ax.set_ylabel('\u0394 Keff (pcm)')
            save_name=f'{D5_test_name}_{energy_mesh}_{ssh_method}{correlation_name}_ERROR_Keff_PCC'
            fig_name=f'HOM CELL Gd157 {energy_mesh} {ssh_method}{CORR}- \u0394 Keff'
        else : # Erreur sur isotopes
            ax.plot([0,ListeCOMPO[-1]],[2,2],'r-.') # limite +2%
            ax.plot([0,ListeCOMPO[-1]],[-2,-2],'r-.') # limite -2%
            ax.set_ylabel('Relative error (%)')

            save_name=f'{D5_test_name}_{energy_mesh}_{ssh_method}{correlation_name}_ERROR_{isotopes_SOUHAITES[k-1]}_PCC'
            fig_name=f'HOM CELL Gd157 {energy_mesh} {ssh_method}{CORR}- \u0394 {isotopes_SOUHAITES[k-1]}'

        ax.grid()
        ax.set_title(fig_name)

        os.chdir(path+'/'+SAVE_DIR)
        fig.savefig(save_name+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
        os.chdir(path)
        plt.close('all')

    return

def POSTPROC_Gd157(pyCOMPOs,ListeCOMPO,ListeAUTOP,name_geom,bu_autop_list,VISU_param,Nmin,S2_libs,S2_PCC,S2_EDEP,ssh_module,ssh_method,CORR,energy_mesh):

    ######################################################################
    #  POST-PROCESSING OF DRAGON5 AND SERPENT2 RESULTS - ERROR PLOTTING  #
    ######################################################################

    # Retrive visualisation parameters
    visu_DRAGON=VISU_param[0]
    visu_SERPENT=VISU_param[1]
    visu_COMP=VISU_param[2]
    visu_DELTA=VISU_param[3]

    # --- Figure size
    SIZE=(8,6)
    # --- Creating SAVEDIR if not existing
    path=os.getcwd()
    a=os.path.exists('BWRresults_PyGan_'+name_geom)
    if a==False:
        os.mkdir('BWRresults_PyGan_'+name_geom)

    SAVE_DIR=f'BWRresults_PyGan_{name_geom}/{bu_autop_list}_postprocess/{ssh_module}_{ssh_method}/'
    print(SAVE_DIR)
    a=os.path.exists(SAVE_DIR)
    print(a)
    if a==False:
        os.makedirs(SAVE_DIR)


    # --- List of isotopes to be post-treated

    isotopes_SOUHAITES=['U235','U238','Pu239','Pu240','Gd155','Gd156','Gd157','Gd158','Xe135','Sm149', 'Sm154','Eu155', 'Eu156', 'Eu157']
    
    if ssh_module == "AUTO":
        ssh_method = "Autosecol"

    if CORR == "CORR":
        CORR = " CORR " 
        correlation_name = "_CORR"
    elif CORR == "noCORR":
        CORR = " noCORR "
        correlation_name = "_noCORR"
    else:
        CORR = " "
        correlation_name = ""
    ################################################################
    #   RETRIVE ALL RESULTS FROM PYCOMPO and SERPENT2 outputs
    ################################################################

    # ------------------------------
    #   RETRIEVE DRAGON RESULTS 
    # ------------------------------
    if visu_DRAGON==1 or visu_COMP==1 or visu_DELTA==1 :

        if 'ASS' in name_geom :
            DIR='EDIBU_HOM'
        else :
            DIR='EDIBU'   

        lenBU_DRAGON=np.shape(ListeCOMPO)[0]
        DRAGON_results_ALL = {}
        print(pyCOMPOs.keys())
        for test_name in pyCOMPOs.keys(): 
            ISOTOPES=pyCOMPOs[test_name][DIR]['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESDENS']
            print(f"Dragon isotopes = {ISOTOPES}")
            #lenISOT_DRAGON=np.shape(ISOTOPES)[0]-1
            lenISOT_DRAGON=np.shape(ISOTOPES)[0]

            DRAGON_BU=ListeCOMPO
            DRAGON_ISOTOPESDENS=np.zeros((lenISOT_DRAGON,lenBU_DRAGON))
            DRAGON_Keff=np.zeros(lenBU_DRAGON)

            #print("$$$ ---------------- DRAGON_ISOTOPESDENS shape = ",lenISOT_DRAGON,lenBU_DRAGON)

            for k in range(lenBU_DRAGON):
                DRAGON_Keff[k]=pyCOMPOs[test_name][DIR]['MIXTURES'][0]['CALCULATIONS'][k]['K-EFFECTIVE']
                #print("$$$ ---------------- ISOTOPES BU step ",k,"/",lenBU_DRAGON," = ",COMPO_py['EDIBU_HOM']['MIXTURES'][0]['CALCULATIONS'][k]['ISOTOPESDENS'])    
                for j in range(lenISOT_DRAGON):
                    #print("$$$ ---------------- ISOTOPES ",j,"/",lenISOT_DRAGON," = ",COMPO_py['EDIBU']['MIXTURES'][0]['CALCULATIONS'][k]['ISOTOPESDENS'][j])
                    DRAGON_ISOTOPESDENS[j][k]=pyCOMPOs[test_name][DIR]['MIXTURES'][0]['CALCULATIONS'][k]['ISOTOPESDENS'][j]


            # --------- List of isotopes from the DRAGON Multicompo results
            isotopes2=[]
            isotopes=[]
            for k in range(lenISOT_DRAGON):
                isotopes2=isotopes2+[pyCOMPOs[test_name][DIR]['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][k]['ALIAS']]
            for k in range(lenISOT_DRAGON):
                if isotopes2[k][0]=='U':
                    isotopes=isotopes+[isotopes2[k][0:4]]
                else:
                    isotopes=isotopes+[isotopes2[k][0:5]]

            indices=np.zeros(len(isotopes_SOUHAITES))
            for n in range(len(isotopes_SOUHAITES)):
                for m in range(len(isotopes)):
                    if isotopes_SOUHAITES[n]==isotopes[m]:
                        indices[n]=m

            # Store in DRAGON_ALL
            DRAGON_test=[
                DRAGON_BU,
                DRAGON_Keff,
            ]
            for k in range(len(isotopes_SOUHAITES)):
                DRAGON_test.append(DRAGON_ISOTOPESDENS[int(indices[k])])
            #print("$$$ ---------------- DRAGON_ALL",DRAGON_AL)
            DRAGON_results_ALL[test_name] = DRAGON_test
        print(f"DRAGON_RESULTS.keys = {DRAGON_results_ALL.keys()}")
    # -------------------------------
    #   RETRIVE SERPENT2 RESULTS 
    # -------------------------------
    if visu_SERPENT==1 or visu_COMP==1 or visu_DELTA==1 :
        SERPENT_results_ALL = {}
        # Importing Keff and isotopes densities from Serpent2
        #for test_name in DRAGON_results_ALL.keys():
        test_name_file = "HOM_UOX_Gd157"
        
        for lib_name in S2_libs:
            for PCC in S2_PCC:
                for edep in S2_EDEP:
                    # --- Keff
                    # 
                    res=serpentTools.read(f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/HOM_CELL_study/{test_name_file}/ENERGY_DEPOSITION_study/{test_name_file}_{lib_name}_edep{edep}_pcc{PCC}_mc_res.m")
                    serpent_keff=res.resdata["absKeff"]
                    np.savetxt(f'serpent_keff_{lib_name}.txt',serpent_keff)
                    SERPENT_keff=np.loadtxt(f'serpent_keff_{lib_name}.txt',dtype=float)
                        
                    # --- BU
                    depl = serpentTools.read(f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/HOM_CELL_study/{test_name_file}/ENERGY_DEPOSITION_study/{test_name_file}_{lib_name}_edep{edep}_pcc{PCC}_mc_dep.m")
                    fuel=depl.materials['total']
                    serpent_BU=fuel.burnup
                    np.savetxt(f'serpent_BU_{lib_name}.txt',serpent_BU)
                    SERPENT_BU=np.loadtxt(f'serpent_BU_{lib_name}.txt',dtype=float)
                    
                    # --- ISOTOPES DENSITIES
                    serpent_ISOTOPESDENS=fuel.toDataFrame("adens",names=isotopes_SOUHAITES)
                    np.savetxt(f'serpent_ISOTOPESDENS_{lib_name}.txt',serpent_ISOTOPESDENS)
                    SERPENT_ISOTOPESDENS=np.loadtxt(f'serpent_ISOTOPESDENS_{lib_name}.txt',dtype=float)
                    SERPENT_ISOTOPESDENS=np.transpose(SERPENT_ISOTOPESDENS)

                    Ls=np.shape(SERPENT_ISOTOPESDENS)
                    lenISOT_SERPENT=Ls[0]
                    lenBU_SERPENT=Ls[1]
                    print('$$$ ---------------- SERPENT_ISOTOPESDENS shape =',Ls)
                    
                    # Rescaling Serpent2 BU to match DRAGON Burnup units.
                    serpent_keff=np.zeros(lenBU_SERPENT)    
                    for k in range(lenBU_SERPENT):
                        SERPENT_BU[k]=1000*SERPENT_BU[k]
                        serpent_keff[k]=SERPENT_keff[k][0]


                    SERPENT_lib_edep_pcc_res=[
                        SERPENT_BU,
                        serpent_keff,
                        ]
                    for k in range(len(isotopes_SOUHAITES)):
                        SERPENT_lib_edep_pcc_res.append(SERPENT_ISOTOPESDENS[k,:])
                    #print("$$$ ---------------- SERPENT_ALL",SERPENT_ALL)
                    SERPENT_results_ALL[f"{lib_name}_{edep}_{PCC}"] = SERPENT_lib_edep_pcc_res
                    print(f"SERPENT_RESULTS = {SERPENT_results_ALL}") 
            print(f"SERPENT_RESULTS.keys = {SERPENT_results_ALL.keys()}")
    # Structure of S2 results = [S2 library used _ S2 energy deposition mode _ S2 PCC][BU, Keff, Isotopes densities]
    # -------------------------------
    #   ERROR MATRICES COMPUTATION 
    # -------------------------------
    
    if visu_DELTA==1 :
        ERRORS_ALL = {}
        print(DRAGON_results_ALL.keys())
        for D5_test_name in DRAGON_results_ALL.keys():
            print(f"D5_test_name = {D5_test_name}")
            for lib_edep_pcc in SERPENT_results_ALL.keys():
                print(f"lib_edep_pcc = {lib_edep_pcc}")
                ERROR=np.zeros((len(isotopes_SOUHAITES)+2,lenBU_DRAGON))
                LE=np.shape(ERROR)
                print('$$$ ------------------------ ERROR shape=',LE)
                for k in range(len(isotopes_SOUHAITES)+2):
                    for j in range(Nmin,lenBU_DRAGON):
                        #print('$$$ ----------------------- k=',k,'    j=',j)
                        #print('$$$ ----------------------- SERPENT_ALL[k][j]=',SERPENT_ALL[k][j])
                        if k==0: # Burnup points
                            ERROR[k][j-Nmin]=SERPENT_results_ALL[lib_edep_pcc][k][j]
                        elif k==1:  # Keff values ---> compute errors in pcm 
                                ERROR[k][j-Nmin]=1.0E+5*(DRAGON_results_ALL[D5_test_name][k][j]-SERPENT_results_ALL[lib_edep_pcc][k][j])
                        else: # Isotopic compositions --> compute relative errors in %   
                            if SERPENT_results_ALL[lib_edep_pcc][k][j]==0 :
                                ERROR[k][j-Nmin]=0
                            else:
                                ERROR[k][j-Nmin]=100*(DRAGON_results_ALL[D5_test_name][k][j]-SERPENT_results_ALL[lib_edep_pcc][k][j])/SERPENT_results_ALL[lib_edep_pcc][k][j]
                ERRORS_ALL[f"{D5_test_name}_{lib_edep_pcc}"] = ERROR            
        print(f"ERRORS = {ERRORS_ALL}")
        print(f"ERRORS_ALL.keys = {ERRORS_ALL.keys()}")
        print(f"Len of BU DRAGON = {lenBU_DRAGON}")
        print(f"SERPENT_RESULTS = {SERPENT_results_ALL}")
        print(f"DRAGON_RESULTS = {DRAGON_results_ALL}")


    ################################################################
    #             PLOT AND SAVE FIGURES
    ################################################################
            

    if visu_DRAGON==1:
        plot_DRAGON(DRAGON_results_ALL, isotopes_SOUHAITES, SIZE, path, SAVE_DIR, DEPL_SOL, ssh_method, correlation_name, sat_name)

    if visu_SERPENT==1:
        plot_Serpent2(SERPENT_results_ALL, isotopes_SOUHAITES, SIZE, path, SAVE_DIR)

    if visu_COMP==1:
        print('$$$ -------- POSTPROC.py : COMPARISON FIGURES ')
        plot_Comparison(DRAGON_results_ALL, SERPENT_results_ALL, isotopes_SOUHAITES, SIZE, path, SAVE_DIR, DEPL_SOL, ssh_method, correlation_name, sat_name)
        
    if visu_DELTA==1 :
        print('$$$ -------- POSTPROC.py : ERROR FIGURES ')
        ERRORS_D5_PCC0 = {}
        ERRORS_D5_PCC1 = {}
        ERRORS_D5_PCC2 = {}
        ERRORS_D5_PCC3 = {}
        ERRORS_D5_PCC3b = {}

        for test_name in ERRORS_ALL.keys():
            D5_pcc = test_name.split('_')[2]
            if D5_pcc == "PCC0":
                ERRORS_D5_PCC0[test_name] = ERRORS_ALL[test_name]
            elif D5_pcc == "PCC1":
                ERRORS_D5_PCC1[test_name] = ERRORS_ALL[test_name]
            elif D5_pcc == "PCC2":
                ERRORS_D5_PCC2[test_name] = ERRORS_ALL[test_name]
            elif D5_pcc == "PCC3":
                ERRORS_D5_PCC3[test_name] = ERRORS_ALL[test_name]
            elif D5_pcc == "PCC3b":
                ERRORS_D5_PCC3b[test_name] = ERRORS_ALL[test_name]

        plot_errors("HOM_UOX_Gd157_PCC0", ERRORS_D5_PCC0, isotopes_SOUHAITES, SIZE, path, SAVE_DIR, DEPL_SOL, ssh_method, correlation_name, sat_name, ListeAUTOP, ListeCOMPO)
        plot_errors("HOM_UOX_Gd157_PCC1", ERRORS_D5_PCC1, isotopes_SOUHAITES, SIZE, path, SAVE_DIR, DEPL_SOL, ssh_method, correlation_name, sat_name, ListeAUTOP, ListeCOMPO)
        plot_errors("HOM_UOX_Gd157_PCC2", ERRORS_D5_PCC2, isotopes_SOUHAITES, SIZE, path, SAVE_DIR, DEPL_SOL, ssh_method, correlation_name, sat_name, ListeAUTOP, ListeCOMPO)
        plot_errors("HOM_UOX_Gd157_PCC3", ERRORS_D5_PCC3, isotopes_SOUHAITES, SIZE, path, SAVE_DIR, DEPL_SOL, ssh_method, correlation_name, sat_name, ListeAUTOP, ListeCOMPO)
        plot_errors("HOM_UOX_Gd157_PCC3b", ERRORS_D5_PCC3b, isotopes_SOUHAITES, SIZE, path, SAVE_DIR, DEPL_SOL, ssh_method, correlation_name, sat_name, ListeAUTOP, ListeCOMPO)
























