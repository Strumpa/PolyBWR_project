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


def MULTI_SERP_POSTPROC(pyCOMPOs,ListeCOMPO,ListeAUTOP,name_geom,bu_autop_list,VISU_param,Nmin,S2_libs,ssh_module,ssh_method,CORR,DEPL_SOL,SAT, set_edep_to_treat):

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
    print(SAT)
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


    if SAT == "SAT":
        SAT=" SAT "
        sat_name = "_SAT"
    elif SAT == "SATOFF":
        SAT = " SATOFF "
        sat_name = "_SATOFF"
    elif SAT == "NSAT":
        SAT = " NSAT "
        sat_name = "_NSAT"
    else:
        SAT = " "
        sat_name = ""

    if CORR == "CORR":
        CORR = " CORR " 
        correlation_name = "_CORR"
    elif CORR == "noCORR":
        CORR = " noCORR "
        correlation_name = "_noCORR"
    else:
        CORR = " "
        correlation_name = ""

    list_tests = pyCOMPOs.keys()
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

    # -------------------------------
    #   RETRIVE SERPENT2 RESULTS 
    # -------------------------------
    if visu_SERPENT==1 or visu_COMP==1 or visu_DELTA==1 :
        SERPENT_results_ALL = {}
        # Importing Keff and isotopes densities from Serpent2
        for test_name in DRAGON_results_ALL.keys():
            if "HOM_UOX_Gd157" in test_name: 
                test_name_file = "HOM_UOX_Gd157"
            elif "HOM_UOX_Gd155" in test_name:
                test_name_file = "HOM_UOX_Gd155"
            elif "HOM_UOXGd" in test_name:
                test_name_file = "HOM_UOXGd"
            SERPENT_test_RESULTS = {}
            for lib_name in S2_libs:
                for edep_id in set_edep_to_treat[lib_name]:
                    if edep_id == 0:
                        edep = "_edep0"
                    elif edep_id == 1:
                        edep = "_edep1"
                    elif edep_id == 2:
                        edep = "_edep2"
                    # --- Keff
                    # 
                    if "NO_NG_toGd158" in test_name:
                        res=serpentTools.read(f"/home/p117902/working_dir/Serpent2_tests/Linux_aarch64/{test_name_file}_{lib_name}_noNG_toGd158_mc_res.m")
                    else:
                        res=serpentTools.read(f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/HOM_CELL_study/{test_name_file}/{test_name_file}_{lib_name}{edep}_mc_res.m")
                    serpent_keff=res.resdata["absKeff"]
                    np.savetxt(f'serpent_keff_{lib_name}.txt',serpent_keff)
                    SERPENT_keff=np.loadtxt(f'serpent_keff_{lib_name}.txt',dtype=float)
                        
                    # --- BU
                    if "NO_NG_toGd158" in test_name:
                        depl = serpentTools.read(f"/home/p117902/working_dir/Serpent2_tests/Linux_aarch64/{test_name_file}_{lib_name}_noNG_toGd158_mc_dep.m")
                    else:
                        depl = serpentTools.read(f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/HOM_CELL_study/{test_name_file}/{test_name_file}_{lib_name}{edep}_mc_dep.m")
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


                    SERPENT_lib_res=[
                        SERPENT_BU,
                        serpent_keff,
                        ]
                    for k in range(len(isotopes_SOUHAITES)):
                        SERPENT_lib_res.append(SERPENT_ISOTOPESDENS[k,:])
                    #print("$$$ ---------------- SERPENT_ALL",SERPENT_ALL)
                    SERPENT_test_RESULTS[f"{lib_name}_edep{edep_id}"] = SERPENT_lib_res
                    print(f"SERPENT_RESULTS = {SERPENT_test_RESULTS}")
                SERPENT_results_ALL[test_name] = SERPENT_test_RESULTS   
            print(f"SERPENT_RESULTS = {SERPENT_results_ALL}") 

    # -------------------------------
    #   ERROR MATRICES COMPUTATION 
    # -------------------------------
    
    if visu_DELTA==1 :
        ERRORS_ALL = {}
        print(DRAGON_results_ALL.keys())
        for test_name in DRAGON_results_ALL.keys():
            ERRORS = {}
            print(f"Test name = {test_name}")
            for lib_name in S2_libs:
                for edep_id in set_edep_to_treat[lib_name]:
                    if edep_id == 0:
                        edep = "_edep0"
                    elif edep_id == 1:
                        edep = "_edep1"
                    elif edep_id == 2:
                        edep = "_edep2"
                    ERROR=np.zeros((len(isotopes_SOUHAITES)+2,lenBU_DRAGON))
                    LE=np.shape(ERROR)
                    print('$$$ ------------------------ ERROR shape=',LE)
                    for k in range(len(isotopes_SOUHAITES)+2):
                        for j in range(Nmin,lenBU_DRAGON):
                            #print('$$$ ----------------------- k=',k,'    j=',j)
                            #print('$$$ ----------------------- SERPENT_ALL[k][j]=',SERPENT_ALL[k][j])
                            if k==0: # Burnup points
                                ERROR[k][j-Nmin]=SERPENT_results_ALL[test_name][f"{lib_name}_edep{edep_id}"][k][j]
                            elif k==1:  # Keff values ---> compute errors in pcm 
                                    ERROR[k][j-Nmin]=1.0E+5*(DRAGON_results_ALL[test_name][k][j]-SERPENT_results_ALL[test_name][f"{lib_name}_edep{edep_id}"][k][j])
                            else: # Isotopic compositions --> compute relative errors in %   
                                if SERPENT_results_ALL[test_name][f"{lib_name}_edep{edep_id}"][k][j]==0 :
                                    ERROR[k][j-Nmin]=0
                                else:
                                    ERROR[k][j-Nmin]=100*(DRAGON_results_ALL[test_name][k][j]-SERPENT_results_ALL[test_name][f"{lib_name}_edep{edep_id}"][k][j])/SERPENT_results_ALL[test_name][f"{lib_name}_edep{edep_id}"][k][j]
                    ERRORS[f"{lib_name}_edep{edep_id}"] = ERROR
            ERRORS_ALL[test_name] = ERRORS               
        print(f"ERRORS = {ERRORS}")
        print(f"ERRORS_ALL = {ERRORS_ALL}")
        print(f"SERPENT_RESULTS = {SERPENT_results_ALL}")
        print(f"DRAGON_RESULTS = {DRAGON_results_ALL}")


    ################################################################
    #             PLOT AND SAVE FIGURES
    ################################################################
            
    colors = {"DRAGON5": "black", "oldlib": "red", "PyNjoy2016": "blue"}
    UOX_tests = []
    Gd_tests = []
    for test in list_tests:
        if "Gd" in test:
            Gd_tests.append(test)
        else:
            UOX_tests.append(test)
    if visu_DRAGON==1:
        print('$$$ -------- POSTPROC.py : DRAGON5 figures ')
        for k in range(len(isotopes_SOUHAITES)+1):
            if UOX_tests:
                # UOX tests :
                legends=[]
                plt.figure(figsize=SIZE)
                for test_name in UOX_tests:
                    plt.plot(DRAGON_results_ALL[test_name][0], DRAGON_results_ALL[test_name][k+1],'2-',linewidth=1)
                    legends.append(f"{test_name} : D5")
                plt.xlabel('BU [MWd/t]')
                plt.grid()
                plt.legend(legends)

                if k == 0: # Comparaison des Keff
                    plt.ylabel('Keff')
                    save_name=f'UOX_tests_{DEPL_SOL}{sat_name}_{ssh_method}{correlation_name}_DRAGON5_Keff'
                    fig_name=f'HOM_CELL : UOX tests {DEPL_SOL}{SAT}{CORR}{ssh_method} - Keff'
                else : # Erreur sur isotopes
                    plt.ylabel('Isotopic density [atom/b-cm]')
                    save_name=f'UOX_tests_{DEPL_SOL}{sat_name}_{ssh_method}{correlation_name}_{isotopes_SOUHAITES[k-1]}_DRAGON5'
                    fig_name=f'HOM_CELL : UOX tests {DEPL_SOL} {SAT} {CORR} {ssh_method} - {isotopes_SOUHAITES[k-1]}'

                plt.title(fig_name)
                os.chdir(path+'/'+SAVE_DIR)
                plt.savefig(save_name+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
                os.chdir(path)
                plt.close('all')

            if Gd_tests:
                # Gd tests :
                legends=[]
                plt.figure(figsize=SIZE)
                for test_name in Gd_tests:
                    plt.plot(DRAGON_results_ALL[test_name][0], DRAGON_results_ALL[test_name][k+1],'2-',linewidth=1)
                    legends.append(f"{test_name} : D5")
                plt.xlabel('BU [MWd/t]')
                plt.grid()
                plt.legend(legends)

                if k == 0: # Comparaison des Keff
                    plt.ylabel('Keff')
                    save_name=f'Gd_tests_{DEPL_SOL}{sat_name}_{ssh_method}{correlation_name}_DRAGON5_Keff'
                    fig_name=f'HOM_CELL : Gd tests {DEPL_SOL}{SAT}{CORR}{ssh_method} - Keff'
                else : # Erreur sur isotopes
                    plt.ylabel('Isotopic density [atom/b-cm]')
                    save_name=f'Gd_tests_{DEPL_SOL}{sat_name}_{ssh_method}{correlation_name}_{isotopes_SOUHAITES[k-1]}_DRAGON5'
                    fig_name=f'HOM_CELL : Gd tests {DEPL_SOL} {SAT} {CORR} {ssh_method} - {isotopes_SOUHAITES[k-1]}'

                plt.title(fig_name)
                os.chdir(path+'/'+SAVE_DIR)
                plt.savefig(save_name+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
                os.chdir(path)
                plt.close('all')


    if visu_SERPENT==1:
        print('$$$ -------- POSTPROC.py : Serpent2 figures ')

        for k in range(len(isotopes_SOUHAITES)+1):
            if UOX_tests:
                # UOX tests :
                legends = []
                plt.figure(figsize=SIZE)
                for test_name in UOX_tests:
                    for lib_name in S2_libs:
                        for edep_id in set_edep_to_treat[lib_name]:
                            if edep_id == 0:
                                edep = "_edep0"
                            elif edep_id == 1:
                                edep = "_edep1"
                            elif edep_id == 2:
                                edep = "_edep2"
                            legends.append(f"{test_name} : {lib_name} edepmode {edep_id}")
                            plt.plot(SERPENT_results_ALL[test_name][f"{lib_name}_edep{edep_id}"][0],SERPENT_results_ALL[test_name][f"{lib_name}_edep{edep_id}"][k+1],'2-',linewidth=1)
                        #legends.append(f"{test_name} : {lib_name}")
                        #plt.plot(SERPENT_results_ALL[test_name][lib_name][0],SERPENT_results_ALL[test_name][lib_name][k+1],'2-',linewidth=1)
                plt.xlabel('BU [MWd/t]')
                plt.grid()
                    
                plt.legend(legends)
                if k == 0: # Keff
                    plt.ylabel('Keff')
                    save_name=f'UOX_tests_Serpent2_Keff'
                    fig_name=f'HOM_CELL : UOX tests - Keff'
                else : # isotopes
                    plt.ylabel('Isotopic density [atom/b-cm]')
                    save_name=f'UOX_tests_Serpent2_{isotopes_SOUHAITES[k-1]}'
                    fig_name=f'HOM_CELL : UOX tests - {isotopes_SOUHAITES[k-1]}'

                plt.title(fig_name)
                os.chdir(path+'/'+SAVE_DIR)
                plt.savefig(save_name+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
                os.chdir(path)
                plt.close('all')
            if Gd_tests:
                # Gd tests :
                legends = []
                plt.figure(figsize=SIZE)
                for test_name in Gd_tests:
                    for lib_name in S2_libs:
                        for edep_id in set_edep_to_treat[lib_name]:
                            if edep_id == 0:
                                edep = "_edep0"
                            elif edep_id == 1:
                                edep = "_edep1"
                            elif edep_id == 2:
                                edep = "_edep2"
                            legends.append(f"{test_name} : {lib_name} edepmode {edep_id}")
                            plt.plot(SERPENT_results_ALL[test_name][f"{lib_name}_edep{edep_id}"][0],SERPENT_results_ALL[test_name][f"{lib_name}_edep{edep_id}"][k+1],'2-',linewidth=1)
                        #legends.append(f"{test_name} : {lib_name}")
                        #plt.plot(SERPENT_results_ALL[test_name][lib_name][0],SERPENT_results_ALL[test_name][lib_name][k+1],'2-',linewidth=1)
                plt.xlabel('BU [MWd/t]')
                plt.grid()
                    
                plt.legend(legends)
                if k == 0: # Keff
                    plt.ylabel('Keff')
                    save_name=f'Gd_tests_Serpent2_Keff'
                    fig_name=f'HOM_CELL : Gd tests - Keff'
                else : # isotopes
                    plt.ylabel('Isotopic density [atom/b-cm]')
                    save_name=f'Gd_tests_Serpent2_{isotopes_SOUHAITES[k-1]}'
                    fig_name=f'HOM_CELL : Gd tests - {isotopes_SOUHAITES[k-1]}'

                plt.grid()
                plt.title(fig_name)
                os.chdir(path+'/'+SAVE_DIR)
                plt.savefig(save_name+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
                os.chdir(path)
                plt.close('all')





    if visu_COMP==1:
        print('$$$ -------- POSTPROC.py : Comparison DRAGON5 / Serpent2 figures ')

        for k in range(len(isotopes_SOUHAITES)+1):
            fig_U, ax_U = plt.subplots(figsize=SIZE)
            fig_Gd, ax_Gd = plt.subplots(figsize=SIZE)
            for test_name in DRAGON_results_ALL.keys():
                print(f"test name is : {test_name}")
                if test_name in UOX_tests:
                    ax_U.plot(DRAGON_results_ALL[test_name][0],DRAGON_results_ALL[test_name][k+1],'2-',linewidth=1, label=f"{test_name} : DRAGON5")
                    for lib_name in S2_libs:
                        for edep_id in set_edep_to_treat[lib_name]:
                            if edep_id == 0:
                                edep = "_edep0"
                            elif edep_id == 1:
                                edep = "_edep1"
                            elif edep_id == 2:
                                edep = "_edep2"
                            ax_U.plot(SERPENT_results_ALL[test_name][f"{lib_name}_edep{edep_id}"][0],SERPENT_results_ALL[test_name][f"{lib_name}_edep{edep_id}"][k+1],'2-',linewidth=1, label=f"{test_name} : {lib_name} edepmode {edep_id}")
                        #ax_U.plot(SERPENT_results_ALL[test_name][lib_name][0],SERPENT_results_ALL[test_name][lib_name][k+1],'2-',linewidth=1, label=f"{test_name} : {lib_name}")
                elif test_name in Gd_tests:
                    ax_Gd.plot(DRAGON_results_ALL[test_name][0],DRAGON_results_ALL[test_name][k+1],'2-',linewidth=1, label=f"{test_name} : DRAGON5")
                    for lib_name in S2_libs:
                        for edep_id in set_edep_to_treat[lib_name]:
                            if edep_id == 0:
                                edep = "_edep0"
                            elif edep_id == 1:
                                edep = "_edep1"
                            elif edep_id == 2:
                                edep = "_edep2"
                            ax_Gd.plot(SERPENT_results_ALL[test_name][f"{lib_name}_edep{edep_id}"][0],SERPENT_results_ALL[test_name][f"{lib_name}_edep{edep_id}"][k+1],'2-',linewidth=1, label=f"{test_name} : {lib_name} edepmode {edep_id}")
                        #ax_Gd.plot(SERPENT_results_ALL[test_name][lib_name][0],SERPENT_results_ALL[test_name][lib_name][k+1],'2-',linewidth=1, label=f"{test_name} : {lib_name}")
                    
            ax_U.set_xlabel('BU [MWd/t]')
            ax_Gd.set_xlabel('BU [MWd/t]')
            ax_U.grid()
            ax_Gd.grid()
            ax_U.legend(loc="best")
            ax_Gd.legend(loc="best")

            if k == 0: # Comparaison des Keff
                ax_U.set_ylabel('Keff')
                ax_Gd.set_ylabel('Keff')
                save_nameU=f'UOX_tests_{DEPL_SOL}{sat_name}_{ssh_method}{correlation_name}_COMP_Keff'
                save_nameGd=f'Gd_tests_{DEPL_SOL}{sat_name}_{ssh_method}{correlation_name}_COMP_Keff'
                fig_nameU=f'HOM_CELL : UOX tests {DEPL_SOL}{SAT} {ssh_method}{CORR}- Keff'
                fig_nameGd=f'HOM_CELL : Gd tests {DEPL_SOL}{SAT} {ssh_method}{CORR}- Keff'
            else : # comparaison isotopes
                ax_U.set_ylabel('Isotopic density [atom/b-cm]')
                ax_Gd.set_ylabel('Isotopic density [atom/b-cm]')
                save_nameU=f'UOX_tests_{DEPL_SOL}{sat_name}_{ssh_method}{correlation_name}_COMP_{isotopes_SOUHAITES[k-1]}'
                fig_nameU=f'HOM_CELL : UOX tests {DEPL_SOL}{SAT} {ssh_method}{CORR}- {isotopes_SOUHAITES[k-1]}'
                save_nameGd=f'Gd_tests_{DEPL_SOL}{sat_name}_{ssh_method}{correlation_name}_COMP_{isotopes_SOUHAITES[k-1]}'
                fig_nameGd=f'HOM_CELL : Gd tests {DEPL_SOL}{SAT} {ssh_method}{CORR}- {isotopes_SOUHAITES[k-1]}'

            ax_U.set_title(fig_nameU)
            ax_Gd.set_title(fig_nameGd)
            os.chdir(path+'/'+SAVE_DIR)
            fig_U.savefig(save_nameU+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
            fig_Gd.savefig(save_nameGd+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
            os.chdir(path)
            plt.close('all')


    if visu_DELTA==1 :
        print('$$$ -------- POSTPROC.py : ERROR DRAGON5-Serpent2 figures ')

        for k in range(len(isotopes_SOUHAITES)+1):
            fig_U,ax_U = plt.subplots(figsize=SIZE)
            fig_Gd, ax_Gd = plt.subplots(figsize=SIZE)
            for test_name in DRAGON_results_ALL.keys():
                if test_name in UOX_tests:
                    for lib_name in S2_libs:
                        for edep_id in set_edep_to_treat[lib_name]:
                            if edep_id == 0:
                                edep = "_edep0"
                            elif edep_id == 1:
                                edep = "_edep1"
                            elif edep_id == 2:
                                edep = "_edep2"
                            ax_U.plot(ERRORS_ALL[test_name][f"{lib_name}_edep{edep_id}"][0],ERRORS_ALL[test_name][f"{lib_name}_edep{edep_id}"][k+1],'2-',linewidth=1, label=f"{test_name} : {lib_name} edepmode {edep_id}")
                elif test_name in Gd_tests:
                    for lib_name in S2_libs:
                        for edep_id in set_edep_to_treat[lib_name]:
                            if edep_id == 0:
                                edep = "_edep0"
                            elif edep_id == 1:
                                edep = "_edep1"
                            elif edep_id == 2:
                                edep = "_edep2"
                            ax_Gd.plot(ERRORS_ALL[test_name][f"{lib_name}_edep{edep_id}"][0],ERRORS_ALL[test_name][f"{lib_name}_edep{edep_id}"][k+1],'2-',linewidth=1, label=f"{test_name} : {lib_name} edepmode {edep_id}")
                ax_U.set_xlabel('BU [MWd/t]')
                ax_U.grid()
                ax_U.legend(loc="best")
                ax_Gd.set_xlabel('BU [MWd/t]')
                ax_Gd.grid()
                ax_Gd.legend(loc="best")
                
            if k == 0: # Erreur sur Keff
                for step in ListeAUTOP:
                    if step <= ListeCOMPO[-1]:
                        ax_U.axvline(x=step,marker='o',color='red')
                        ax_Gd.axvline(x=step,marker='o',color='red')
                ax_U.plot([0,ListeCOMPO[-1]],[300,300],'r-.') # limite +300pcm
                ax_U.plot([0,ListeCOMPO[-1]],[-300,-300],'r-.') # limite -300pcm
                ax_U.set_ylabel('\u0394 Keff (pcm)')
                save_nameU=f'UOX_tests_{DEPL_SOL}{sat_name}_{ssh_method}{correlation_name}_ERROR_Keff'
                fig_nameU=f'HOM_CELL : UOX tests {DEPL_SOL}{SAT} {ssh_method}{CORR}- \u0394 Keff'

                ax_Gd.plot([0,ListeCOMPO[-1]],[300,300],'r-.') # limite +300pcm
                ax_Gd.plot([0,ListeCOMPO[-1]],[-300,-300],'r-.') # limite -300pcm
                ax_Gd.set_ylabel('\u0394 Keff (pcm)')
                save_nameGd=f'Gd_tests_{DEPL_SOL}{sat_name}_{ssh_method}{correlation_name}_ERROR_Keff'
                fig_nameGd=f'HOM_CELL Gd tests {DEPL_SOL}{SAT} {ssh_method}{CORR}- \u0394 Keff'
            else : # Erreur sur isotopes
                ax_U.plot([0,ListeCOMPO[-1]],[2,2],'r-.') # limite +2%
                ax_U.plot([0,ListeCOMPO[-1]],[-2,-2],'r-.') # limite -2%
                ax_Gd.plot([0,ListeCOMPO[-1]],[2,2],'r-.') # limite +2%
                ax_Gd.plot([0,ListeCOMPO[-1]],[-2,-2],'r-.') # limite -2%
                ax_U.set_ylabel('Relative error (%)')
                ax_Gd.set_ylabel('Relative error (%)')

                save_nameU=f'UOX_tests_{DEPL_SOL}{sat_name}_{ssh_method}{correlation_name}_ERROR_{isotopes_SOUHAITES[k-1]}'
                fig_nameU=f'HOM_CELL : UOX tests {DEPL_SOL}{SAT} {ssh_method}{CORR}- \u0394 {isotopes_SOUHAITES[k-1]}'
                save_nameGd=f'Gd_tests_{DEPL_SOL}{sat_name}_{ssh_method}{correlation_name}_ERROR_{isotopes_SOUHAITES[k-1]}'
                fig_nameGd=f'HOM_CELL : Gd tests {DEPL_SOL}{SAT} {ssh_method}{CORR}- \u0394 {isotopes_SOUHAITES[k-1]}'

            ax_U.grid()
            ax_Gd.grid()
            ax_U.set_title(fig_nameU)
            ax_Gd.set_title(fig_nameGd)
            os.chdir(path+'/'+SAVE_DIR)
            fig_U.savefig(save_nameU+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
            fig_Gd.savefig(save_nameGd+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
            os.chdir(path)
            plt.close('all')
        """
        # --- Print results
        for test_name in DRAGON_results_ALL.keys():
            for lib_name in S2_libs:
                for edep_id in set_edep_to_treat[lib_name]:
                    if edep_id == 0:
                        edep = ""
                    elif edep_id == 1:
                        edep = "_edep1"
                    elif edep_id == 2:
                        edep = "_edep2"
                    print(f"For test name = {test_name}, SSH : {ssh_module} {ssh_method} {CORR}, using energy deposition mode {edep_id} : The error on Keff at t=0 is {ERRORS_ALL[test_name][f"{lib_name}_edep{edep_id}"][1][0]}, Keff_DRAGON = {DRAGON_results_ALL[test_name][1][0]}, Keff_SERPENT_{lib_name} = {SERPENT_results_ALL[test_name][f"{lib_name}_edep{edep_id}"][1][0]}")
        """



























