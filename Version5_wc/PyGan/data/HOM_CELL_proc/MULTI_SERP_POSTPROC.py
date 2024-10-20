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


def MULTI_SERP_POSTPROC(pyCOMPO,ListeCOMPO,ListeAUTOP,name_geom,name_mix,suffixe,VISU_param,Nmin,S2_libs,ssh_module,ssh_method,CORR,DEPL_SOL,SAT):

    ######################################################################
    #  POST-PROCESSING OF DRAGON5 AND SERPENT2 RESULTS - ERROR PLOTTING  #
    ######################################################################

    # Retrive visualisation parameters
    visu_DRAGON=VISU_param[0]
    visu_SERPENT=VISU_param[1]
    visu_COMP=VISU_param[2]
    visu_DELTA=VISU_param[3]

    # --- Figure size
    SIZE=(6,4)

    # --- Creating SAVEDIR if not existing
    path=os.getcwd()
    a=os.path.exists('BWRresults_PyGan_'+name_geom)
    if a==False:
        os.mkdir('BWRresults_PyGan_'+name_geom)

    SAVE_DIR=f'BWRresults_PyGan_{name_geom}/{name_mix}/{suffixe}_postprocess/{ssh_module}_{ssh_method}/'
    print(SAVE_DIR)
    a=os.path.exists(SAVE_DIR)
    print(a)
    if a==False:
        os.makedirs(SAVE_DIR)


    # --- List of isotopes to be post-treated
    if name_mix == "AT10_C7_hom":
        isotopes_SOUHAITES=['U235','U236','U238','Pu239','Pu240','Pu241','Pu242','Gd154','Gd155','Gd156','Gd157','Gd158','Gd160','Xe135','Sm149']
    elif name_mix == "HOM_U5":
        isotopes_SOUHAITES=["U235", "Gd155", "Gd157", "Xe135", "Sm149"]
    elif name_mix == "HOM_U5_U8":
        isotopes_SOUHAITES=["U235", "U238", "Pu239", "Pu240", "Gd155", "Gd157", "Xe135", "Sm149"]
    elif name_mix == "HOM_UOX":
        isotopes_SOUHAITES=["U234", "U235", "U238", "Pu239", "Pu240", "Gd155", "Gd157", "Xe135", "Sm149"]
    


    if SAT == "SAT":
        SAT=" SAT "
        sat_name = "_SAT"
    elif SAT == "SATOFF":
        SAT = " SATOFF "
        sat_name = "_SATOFF"
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
         
        ISOTOPES=pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESDENS']
        print(f"Dragon isotopes = {ISOTOPES}")
        #lenISOT_DRAGON=np.shape(ISOTOPES)[0]-1
        lenISOT_DRAGON=np.shape(ISOTOPES)[0]

        DRAGON_BU=ListeCOMPO
        DRAGON_ISOTOPESDENS=np.zeros((lenISOT_DRAGON,lenBU_DRAGON))
        DRAGON_Keff=np.zeros(lenBU_DRAGON)

        #print("$$$ ---------------- DRAGON_ISOTOPESDENS shape = ",lenISOT_DRAGON,lenBU_DRAGON)

        for k in range(lenBU_DRAGON):
            DRAGON_Keff[k]=pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][k]['K-EFFECTIVE']
            #print("$$$ ---------------- ISOTOPES BU step ",k,"/",lenBU_DRAGON," = ",pyCOMPO['EDIBU_HOM']['MIXTURES'][0]['CALCULATIONS'][k]['ISOTOPESDENS'])    
            for j in range(lenISOT_DRAGON):
                #print("$$$ ---------------- ISOTOPES ",j,"/",lenISOT_DRAGON," = ",pyCOMPO['EDIBU']['MIXTURES'][0]['CALCULATIONS'][k]['ISOTOPESDENS'][j])
                DRAGON_ISOTOPESDENS[j][k]=pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][k]['ISOTOPESDENS'][j]


        # --------- List of isotopes from the DRAGON Multicompo results
        isotopes2=[]
        isotopes=[]
        for k in range(lenISOT_DRAGON):
            isotopes2=isotopes2+[pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][k]['ALIAS']]
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
        DRAGON_ALL=[
               DRAGON_BU,
               DRAGON_Keff,
           ]
        for k in range(len(isotopes_SOUHAITES)):
            DRAGON_ALL.append(DRAGON_ISOTOPESDENS[int(indices[k])])
        #print("$$$ ---------------- DRAGON_ALL",DRAGON_AL)

    # -------------------------------
    #   RETRIVE SERPENT2 RESULTS 
    # -------------------------------
    if visu_SERPENT==1 or visu_COMP==1 or visu_DELTA==1 :
        SERPENT_RESULTS = {}
        # Importing Keff and isotopes densities from Serpent2
        for lib_name in S2_libs:
            # --- Keff
            # 
            res=serpentTools.read(f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/HOM_CELL_study/{name_mix}/{name_mix}_{lib_name}_mc_res.m")

            serpent_keff=res.resdata["absKeff"]
            np.savetxt(f'serpent_keff_{lib_name}.txt',serpent_keff)
            SERPENT_keff=np.loadtxt(f'serpent_keff_{lib_name}.txt',dtype=float)
                
            # --- BU
            depl = serpentTools.read(f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/HOM_CELL_study/{name_mix}/{name_mix}_{lib_name}_mc_dep.m")
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



            SERPENT_ALL=[
                SERPENT_BU,
                serpent_keff,
            ]
            for k in range(len(isotopes_SOUHAITES)):
                SERPENT_ALL.append(SERPENT_ISOTOPESDENS[k,:])
            #print("$$$ ---------------- SERPENT_ALL",SERPENT_ALL)
            SERPENT_RESULTS[lib_name] = SERPENT_ALL
        print(f"SERPENT_RESULTS = {SERPENT_RESULTS}")    

    # -------------------------------
    #   ERROR MATRICES COMPUTATION 
    # -------------------------------
    if visu_DELTA==1 :
        ERRORS = {}
        for lib_name in S2_libs:
            ERROR=np.zeros((len(isotopes_SOUHAITES)+2,lenBU_DRAGON))
            LE=np.shape(ERROR)
            print('$$$ ------------------------ ERROR shape=',LE)
            for k in range(len(isotopes_SOUHAITES)+2):
                    
                for j in range(Nmin,lenBU_DRAGON):
                    #print('$$$ ----------------------- k=',k,'    j=',j)
                    #print('$$$ ----------------------- SERPENT_ALL[k][j]=',SERPENT_ALL[k][j])
                    if k==0: # Burmup points
                        ERROR[k][j-Nmin]=SERPENT_RESULTS[lib_name][k][j]
                    elif k==1:  # Keff values ---> compute errors in pcm 
                            ERROR[k][j-Nmin]=1.0E+5*(DRAGON_ALL[k][j]-SERPENT_RESULTS[lib_name][k][j])
                    else: # Isotopic compositions --> compute relative errors in %   
                        if SERPENT_RESULTS[lib_name][k][j]==0 :
                            ERROR[k][j-Nmin]=0
                        else:
                            ERROR[k][j-Nmin]=100*(DRAGON_ALL[k][j]-SERPENT_RESULTS[lib_name][k][j])/SERPENT_RESULTS[lib_name][k][j]
            ERRORS[lib_name] = ERROR               
        print(f"ERRORS = {ERRORS}")
        #print(f"Len iso DRAGON = {lenISOT_DRAGON}, len iso Serp = {lenISOT_SERPENT}")
    ################################################################
    #             PLOT AND SAVE FIGURES
    ################################################################
            
    colors = {"DRAGON5": "black", "oldlib": "red", "PyNjoy2016": "blue"}

    if visu_DRAGON==1:
        print('$$$ -------- POSTPROC.py : DRAGON5 figures ')

        for k in range(len(isotopes_SOUHAITES+1)):
            plt.figure(figsize=SIZE)
            plt.plot(DRAGON_ALL[0],DRAGON_ALL[k+1],'2-',linewidth=1)
            plt.xlabel('BU [MWd/t]')
            plt.grid()
            plt.legend(['DRAGON5'])

            if k == 0: # Comparaison des Keff
                plt.ylabel('Keff')
                save_name=f'{name_mix}_{DEPL_SOL}{sat_name}_{ssh_method}_DRAGON5_Keff'
                fig_name=f'{name_mix} {DEPL_SOL}{SAT}{CORR}{ssh_method} - Keff'
            else : # Erreur sur isotopes
                plt.ylabel('Isotopic density [atom/b-cm]')
                save_name=f'{name_mix}_{DEPL_SOL}{SAT}{CORR}_{ssh_method}_{isotopes_SOUHAITES[k-1]}'
                fig_name=f'{name_mix} {DEPL_SOL} {SAT} {CORR} {ssh_method} - {isotopes_SOUHAITES[k-1]}'

            plt.title(fig_name)
            os.chdir(path+'/'+SAVE_DIR)
            plt.savefig(save_name+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
            os.chdir(path)
            plt.close('all')


    if visu_SERPENT==1:
        print('$$$ -------- POSTPROC.py : Serpent2 figures ')

        for k in range(len(isotopes_SOUHAITES)+1):
            legends = []
            plt.figure(figsize=SIZE)
            for lib_name in S2_libs:
                legends.append(lib_name)
                plt.plot(SERPENT_RESULTS[lib_name][0],SERPENT_RESULTS[lib_name][k+1],'2-',linewidth=1, color=colors[lib_name])
            plt.xlabel('BU [MWd/t]')
            plt.grid()
                
            plt.legend(legends)
            if k == 0: # Keff
                plt.ylabel('Keff')
                save_name=f'{name_mix}_Serpent2_Keff'
                fig_name=f'{name_mix} - Keff'
            else : # isotopes
                plt.ylabel('Isotopic density [atom/b-cm]')
                save_name=f'{name_mix}_Serpent2_{isotopes_SOUHAITES[k-1]}'
                fig_name=f'{name_mix} - {isotopes_SOUHAITES[k-1]}'

            plt.title(fig_name)
            os.chdir(path+'/'+SAVE_DIR)
            plt.savefig(save_name+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
            os.chdir(path)
            plt.close('all')


    if visu_COMP==1:
        print('$$$ -------- POSTPROC.py : Comparison DRAGON5 / Serpent2 figures ')

        for k in range(len(isotopes_SOUHAITES)+1):
            legends = []
            plt.figure(figsize=SIZE)
            plt.plot(DRAGON_ALL[0],DRAGON_ALL[k+1],'2-',linewidth=1, color = colors["DRAGON5"])
            for lib_name in S2_libs:
                legends.append(lib_name)
                plt.plot(SERPENT_RESULTS[lib_name][0],SERPENT_RESULTS[lib_name][k+1],'2-',linewidth=1, color=colors[lib_name])
                
            plt.xlabel('BU (MWd/t)')
            plt.grid()
            plt.legend(legends)

            if k == 0: # Comparaison des Keff
                plt.ylabel('Keff')
                save_name=f'{name_mix}_{DEPL_SOL}{sat_name}_{ssh_method}{correlation_name}_COMP_Keff'
                fig_name=f'{name_mix} {DEPL_SOL}{SAT} {ssh_method}{CORR}- Keff'
            else : # comparaison isotopes
                plt.ylabel('Isotopic density [atom/b-cm]')
                save_name=f'{name_mix}_COMP_{DEPL_SOL}{sat_name}_{ssh_method}{correlation_name}_{isotopes_SOUHAITES[k-1]}'
                fig_name=f'{name_mix} {DEPL_SOL}{SAT} {ssh_method}{CORR}- {isotopes_SOUHAITES[k-1]}'

            plt.title(fig_name)
            os.chdir(path+'/'+SAVE_DIR)
            plt.savefig(save_name+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
            os.chdir(path)
            plt.close('all')
        

    if visu_DELTA==1 :
        print('$$$ -------- POSTPROC.py : ERROR DRAGON5-Serpent2 figures ')

        for k in range(len(isotopes_SOUHAITES)+1): # -1 lenISOT_DRAGON-1 ?
            legends = []
            plt.figure(figsize=SIZE)
            for lib_name in S2_libs:
                legends.append(lib_name)
                plt.plot(ERRORS[lib_name][0],ERRORS[lib_name][k+1],'2-',linewidth=1, color=colors[lib_name])
            plt.xlabel('BU [MWd/t]')
            plt.grid()
            plt.legend(legends)
            
            if k == 0: # Erreur sur Keff
                for step in ListeAUTOP:
                    if step <= ListeCOMPO[-1]:
                        plt.axvline(x=step,marker='o',color='red')
                plt.plot([0,ListeCOMPO[-1]],[300,300],'r-.') # limite +300pcm
                plt.plot([0,ListeCOMPO[-1]],[-300,-300],'r-.') # limite -300pcm
                plt.ylabel('\u0394 Keff (pcm)')
                save_name=f'{name_mix}_{DEPL_SOL}{sat_name}_{ssh_method}{correlation_name}_ERROR_Keff'
                fig_name=f'{name_mix} {DEPL_SOL}{SAT} {ssh_method}{CORR}- \u0394 Keff'
            else : # Erreur sur isotopes
                plt.plot([0,ListeCOMPO[-1]],[2,2],'r-.') # limite +2%
                plt.plot([0,ListeCOMPO[-1]],[-2,-2],'r-.') # limite -2%
                plt.ylabel('Relative error (%)')
                save_name=f'{name_mix}_{DEPL_SOL}{sat_name}_{ssh_method}{correlation_name}_ERROR_{isotopes_SOUHAITES[k-1]}'
                fig_name=f'{name_mix} {DEPL_SOL}{SAT} {ssh_method}{CORR}- \u0394 {isotopes_SOUHAITES[k-1]}'

            plt.title(fig_name)
            os.chdir(path+'/'+SAVE_DIR)
            plt.savefig(save_name+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
            os.chdir(path)
            plt.close('all')
        for lib_name in S2_libs:

            print(f"For SSH : {ssh_module} {ssh_method} {CORR} : The error on Keff at t=0 is {ERRORS[lib_name][1][0]}, Keff_DRAGON = {DRAGON_ALL[1][0]}, Keff_SERPENT_{lib_name} = {SERPENT_RESULTS[lib_name][1][0]}")
        



























