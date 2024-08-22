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


def POSTPROC_hom(pyCOMPO,ListeCOMPO,ListeAUTOP,name_geom,name_mix,suffixe,VISU_param,Nmin,GdCompo,S2_libs,ssh_module,ssh_method,CORR,SAT,DEPL_SOL):

    ########################################################
    #  Paramètres et initialisation et paramètres locaux 
    ########################################################

    # Récupération des paramètres de visualisation
    visu_DRAGON=VISU_param[0]
    visu_SERPENT=VISU_param[1]
    visu_COMP=VISU_param[2]
    visu_DELTA=VISU_param[3]

    # --- Paramètres de figures
    SIZE=(6,4)

    # --- Chemin du répertoire courant
    path=os.getcwd()
    #print("$$$ ---------- path = ",path)

    # --- Création du répertoire de stockage de résultats 
    a=os.path.exists('BWRresults_PyGan_'+name_geom)
    if a==False:
        os.mkdir('BWRresults_PyGan_'+name_geom)

    SAVE_DIR='BWRresults_PyGan_'+name_geom+'/'+name_mix+'_'+suffixe+'_postprocess/'+ssh_module+'_'+ssh_method+"/"+GdCompo+"/"
    print(SAVE_DIR)
    a=os.path.exists(SAVE_DIR)
    print(a)
    if a==False:
        os.makedirs(SAVE_DIR)

    if SAT == "SAT":
        SAT="_SAT"
    elif SAT == "SATOFF":
        SAT="_SATOFF"
    else:
        SAT=""

    # --- Numéro de MIX pour récupérer résultats homogénéisés et condensés 
    NMIX=0

    # --- Ordre des isotopoes post-processés
    if name_mix == "AT10_C7_hom":
        isotopes_SOUHAITES=['U235','U236','U238','Pu239','Pu240','Pu241','Pu242','Gd154','Gd155','Gd156','Gd157','Gd158','Gd160','Xe135','Sm149']

        if GdCompo == "Nominal":
            Gd_identifier = ""
        elif GdCompo == "Gd154":
            Gd_identifier = "_Gd154"
        elif GdCompo == "Gd155":
            Gd_identifier = "_Gd155"
        elif GdCompo == "Gd156":
            Gd_identifier = "_Gd156"
        elif GdCompo == "Gd157":
            Gd_identifier = "_Gd157"
        elif GdCompo == "Gd158":
            Gd_identifier = "_Gd158"
        elif GdCompo == "Gd160":
            Gd_identifier = "_Gd160"
        elif GdCompo == "noGd":
            Gd_identifier = "_noGd"
    elif name_mix == "HOM_U5":
        isotopes_SOUHAITES=["U235", "Gd155", "Gd157", "Xe135", "Sm149"]
        Gd_identifier = ""
    elif name_mix == "HOM_U5_U8":
        isotopes_SOUHAITES=["U235", "U238", "Pu239", "Pu240", "Gd155", "Gd157", "Xe135", "Sm149"]
        Gd_identifier = ""


    ################################################################
    #   CREACTION DES MATRICES DE DONNEES : DRAGON - SERPENT - ERROR
    ################################################################

    # ------------------------------
    #   MATRICE DES RESULTATS DRAGON 
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

        #print('$$$ ---------------- DRAGON_BU =',DRAGON_BU)
        #print('$$$ ---------------- LEN_DRAGON_BU =',len(DRAGON_BU))
        #print("$$$ ---------------- DRAGON_Keff = ",DRAGON_Keff)    
        #print("$$$ ---------------- DRAGON_ISOTOPESDENS = ",DRAGON_ISOTOPESDENS)

        # --------- Liste des isotopes recuperes dans la Multicompo
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

        #print("$$$ ---------------- DRAGON isotopes = ",isotopes)
        #print("$$$ ---------------- DRAGON isotopes souhaites = ",isotopes_SOUHAITES)
        #print("$$$ ---------------- indices correspondances = ",indices)
            
        # Store in DRAGON_ALL
        DRAGON_ALL=[
               DRAGON_BU,
               DRAGON_Keff,
           ]
        for k in range(len(isotopes_SOUHAITES)):
            DRAGON_ALL.append(DRAGON_ISOTOPESDENS[int(indices[k])])
        #print("$$$ ---------------- DRAGON_ALL",DRAGON_AL)

    # -------------------------------
    #   MATRICE DES RESULTATS SERPENT 
    # -------------------------------
    if visu_SERPENT==1 or visu_COMP==1 or visu_DELTA==1 :
        # Importing Keff and isotopes densities from Serpent2
        if "PyNjoy2016" in S2_libs:
            # 1) results from PyNjoy2016, with Predictor-Corrector BU evolution method
            # --- Keff
            # 
            res_pynjoy=serpentTools.read(f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/HOM_CELL_study/{name_mix}{Gd_identifier}_PyNjoy2016_mc_res.m")

            serpent_keff_pynjoy=res_pynjoy.resdata["absKeff"]
            np.savetxt('serpent_keff_pynjoy.txt',serpent_keff_pynjoy)
            SERPENT_keff_pynjoy=np.loadtxt('serpent_keff_pynjoy.txt',dtype=float)
                
            # --- BU
            #depFile = SERPENT_path+name_mix+"_mc_dep.m"
            depFile = f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/HOM_CELL_study/{name_mix}{Gd_identifier}_PyNjoy2016_mc_dep.m"
            dep_pynjoy = serpentTools.read(depFile)
            fuel=dep_pynjoy.materials['total']
            serpent_BU_pynjoy=fuel.burnup
            np.savetxt('serpent_BU_pynjoy.txt',serpent_BU_pynjoy)
            SERPENT_BU_pynjoy=np.loadtxt('serpent_BU_pynjoy.txt',dtype=float)
            
            # --- ISOTOPES DENSITIES
            serpent_ISOTOPESDENS=fuel.toDataFrame("adens",names=isotopes_SOUHAITES)
            np.savetxt('serpent_ISOTOPESDENS_pynjoy.txt',serpent_ISOTOPESDENS)
            SERPENT_ISOTOPESDENS=np.loadtxt('serpent_ISOTOPESDENS_pynjoy.txt',dtype=float)
            SERPENT_ISOTOPESDENS_pynjoy=np.transpose(SERPENT_ISOTOPESDENS)

            Ls1=np.shape(SERPENT_BU_pynjoy)
            print('$$$ ---------------- SERPENT_BU shape =',Ls1)

            Ls2=np.shape(SERPENT_keff_pynjoy)
            #lenISOT_SERPENT2=Ls2[0]
            #lenBU_SERPENT2=Ls2[1]
            print('$$$ ---------------- SERPENT_keff shape =',Ls2)

            Ls=np.shape(SERPENT_ISOTOPESDENS_pynjoy)
            lenISOT_SERPENT=Ls[0]
            lenBU_SERPENT=Ls[1]
            print('$$$ ---------------- SERPENT_ISOTOPESDENS shape =',Ls)
            
            # Modification unite BU pour match avec DRAGON
            serpent_keff_pynjoy=np.zeros(lenBU_SERPENT)    
            for k in range(lenBU_SERPENT):
                SERPENT_BU_pynjoy[k]=1000*SERPENT_BU_pynjoy[k]
                serpent_keff_pynjoy[k]=SERPENT_keff_pynjoy[k][0]

            #print('$$$ ---------------- SERPENT_BU =',SERPENT_BU_pynjoy)
            #print("$$$ ---------------- SERPENT_Keff = ",serpent_keff_pynjoy)    
            #print("$$$ ---------------- SERPENT_ISOTOPESDENS = ",SERPENT_ISOTOPESDENS_pynjoy)

            SERPENT_ALL_pynjoy=[
            SERPENT_BU_pynjoy,
            serpent_keff_pynjoy,
            ]
            for k in range(len(isotopes_SOUHAITES)):
                SERPENT_ALL_pynjoy.append(SERPENT_ISOTOPESDENS_pynjoy[k,:])
            #print("$$$ ---------------- SERPENT_ALL",SERPENT_ALL)

            # 1) bis : results from PyNjoy2016, with Euler BU evolution method
            if name_mix == "AT10_C7_hom":
                if Gd_identifier == "_Gd154" or Gd_identifier == "_Gd155" or Gd_identifier == "" or Gd_identifier == "_noGd":

                    # --- Keff
                    res_pynjoy_euler=serpentTools.read(f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/HOM_CELL_study/{name_mix}{Gd_identifier}_PyNjoy2016_Euler_mc_res.m")
                    serpent_keff_pynjoy_euler=res_pynjoy_euler.resdata["absKeff"]
                    np.savetxt('serpent_keff_pynjoy_euler.txt',serpent_keff_pynjoy_euler)
                    SERPENT_keff_pynjoy_euler=np.loadtxt('serpent_keff_pynjoy_euler.txt',dtype=float)
                        
                    # --- BU
                    #depFile = SERPENT_path+name_mix+"_mc_dep.m"
                    depFile_euler = f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/HOM_CELL_study/{name_mix}{Gd_identifier}_PyNjoy2016_Euler_mc_dep.m"
                    dep_pynjoy_euler = serpentTools.read(depFile_euler)
                    fuel_euler=dep_pynjoy_euler.materials['total']
                    serpent_BU_pynjoy_euler=fuel_euler.burnup
                    np.savetxt('serpent_BU_pynjoy_euler.txt',serpent_BU_pynjoy_euler)
                    SERPENT_BU_pynjoy_euler=np.loadtxt('serpent_BU_pynjoy_euler.txt',dtype=float)
                    
                    # --- ISOTOPES DENSITIES
                    serpent_ISOTOPESDENS_euler=fuel_euler.toDataFrame("adens",names=isotopes_SOUHAITES)
                    np.savetxt('serpent_ISOTOPESDENS_pynjoy_euler.txt',serpent_ISOTOPESDENS_euler)
                    SERPENT_ISOTOPESDENS_euler=np.loadtxt('serpent_ISOTOPESDENS_pynjoy_euler.txt',dtype=float)
                    SERPENT_ISOTOPESDENS_pynjoy_euler=np.transpose(SERPENT_ISOTOPESDENS_euler)

                    Ls1=np.shape(SERPENT_BU_pynjoy_euler)
                    print('$$$ ---------------- SERPENT_BU shape =',Ls1)

                    Ls2=np.shape(SERPENT_keff_pynjoy_euler)
                    #lenISOT_SERPENT2=Ls2[0]
                    #lenBU_SERPENT2=Ls2[1]
                    print('$$$ ---------------- SERPENT_keff shape =',Ls2)

                    Ls=np.shape(SERPENT_ISOTOPESDENS_pynjoy_euler)
                    lenISOT_SERPENT=Ls[0]
                    lenBU_SERPENT_euler=Ls[1]
                    print('$$$ ---------------- SERPENT_ISOTOPESDENS shape =',Ls)
                    
                    # Modification unite BU pour match avec DRAGON
                    serpent_keff_pynjoy_euler=np.zeros(lenBU_SERPENT_euler)    
                    for k in range(lenBU_SERPENT_euler):
                        SERPENT_BU_pynjoy_euler[k]=1000*SERPENT_BU_pynjoy_euler[k]
                        serpent_keff_pynjoy_euler[k]=SERPENT_keff_pynjoy_euler[k][0]

                    #print('$$$ ---------------- SERPENT_BU =',SERPENT_BU_pynjoy)
                    #print("$$$ ---------------- SERPENT_Keff = ",serpent_keff_pynjoy)    
                    #print("$$$ ---------------- SERPENT_ISOTOPESDENS = ",SERPENT_ISOTOPESDENS_pynjoy)

                    SERPENT_ALL_pynjoy_euler=[
                        SERPENT_BU_pynjoy_euler,
                        serpent_keff_pynjoy_euler,
                    ]
                    for k in range(len(isotopes_SOUHAITES)):
                        SERPENT_ALL_pynjoy_euler.append(SERPENT_ISOTOPESDENS_pynjoy_euler[k,:])


        # 2) results from Serpent2 with sss_jeff311 library
        if "oldlib" in S2_libs:
            # --- Keff
            res_oldlib=serpentTools.read(f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/HOM_CELL_study/{name_mix}{Gd_identifier}_oldlib_mc_res.m")

            serpent_keff_oldlib=res_oldlib.resdata["absKeff"]
            np.savetxt('serpent_keff_oldlib.txt',serpent_keff_oldlib)
            SERPENT_keff_oldlib=np.loadtxt('serpent_keff_oldlib.txt',dtype=float)
                
            # --- BU
            depFile = f"/home/p117902/working_dir/Serpent2_para_bateman/Linux_aarch64/HOM_CELL_study/{name_mix}{Gd_identifier}_oldlib_mc_dep.m"
            dep_oldlib = serpentTools.read(depFile)
            fuel=dep_oldlib.materials['total']
            serpent_BU_oldlib=fuel.burnup
            np.savetxt('serpent_BU_oldlib.txt',serpent_BU_oldlib)
            SERPENT_BU_oldlib=np.loadtxt('serpent_BU_oldlib.txt',dtype=float)
            
            # --- ISOTOPES DENSITIES
            serpent_ISOTOPESDENS=fuel.toDataFrame("adens",names=isotopes_SOUHAITES)
            np.savetxt('serpent_ISOTOPESDENS_oldlib.txt',serpent_ISOTOPESDENS)
            SERPENT_ISOTOPESDENS_oldlib=np.loadtxt('serpent_ISOTOPESDENS_oldlib.txt',dtype=float)
            SERPENT_ISOTOPESDENS_oldlib=np.transpose(SERPENT_ISOTOPESDENS_oldlib)

            Ls1=np.shape(SERPENT_BU_oldlib)
            print('$$$ ---------------- SERPENT_BU shape =',Ls1)

            Ls2=np.shape(SERPENT_keff_oldlib)
            print('$$$ ---------------- SERPENT_keff shape =',Ls2)

            Ls=np.shape(SERPENT_ISOTOPESDENS_oldlib)
            lenISOT_SERPENT=Ls[0]
            lenBU_SERPENT=Ls[1]
            print('$$$ ---------------- SERPENT_ISOTOPESDENS shape =',Ls)
            
            # Modification unite BU pour match avec DRAGON
            serpent_keff_oldlib=np.zeros(lenBU_SERPENT)    
            for k in range(lenBU_SERPENT):
                SERPENT_BU_oldlib[k]=1000*SERPENT_BU_oldlib[k]
                serpent_keff_oldlib[k]=SERPENT_keff_oldlib[k][0]

            #print('$$$ ---------------- SERPENT_BU =',SERPENT_BU_oldlib)
            #print("$$$ ---------------- SERPENT_Keff = ",serpent_keff_oldlib)    
            #print("$$$ ---------------- SERPENT_ISOTOPESDENS = ",SERPENT_ISOTOPESDENS_oldlib)

            SERPENT_ALL_oldlib=[
            SERPENT_BU_oldlib,
            serpent_keff_oldlib,
            ]
            for k in range(len(isotopes_SOUHAITES)):
                SERPENT_ALL_oldlib.append(SERPENT_ISOTOPESDENS_oldlib[k,:])
        #print("$$$ ---------------- SERPENT_ALL",SERPENT_ALL)

    # -------------------------------
    #   MATRICE DES ERREUR : ERROR 
    # -------------------------------
    if visu_DELTA==1 :

        ERROR_pynjoy=np.zeros((len(isotopes_SOUHAITES)+2,lenBU_DRAGON))
        ERROR_pynjoy_euler=np.zeros((len(isotopes_SOUHAITES)+2,lenBU_DRAGON))
        ERROR_oldlib=np.zeros((len(isotopes_SOUHAITES)+2,lenBU_DRAGON))
        LE=np.shape(ERROR_pynjoy)
        print('$$$ ------------------------ ERROR shape=',LE)
        for k in range(len(isotopes_SOUHAITES)+2):
                  
            for j in range(Nmin,lenBU_DRAGON):
                #print('$$$ ----------------------- k=',k,'    j=',j)
                #print('$$$ ----------------------- SERPENT_ALL[k][j]=',SERPENT_ALL[k][j])
                if k==0: # Vecteur BU
                    print(f"j = {j}")
                    if "PyNjoy2016" in S2_libs:
                        ERROR_pynjoy[k][j-Nmin]=SERPENT_ALL_pynjoy[k][j]
                        if name_mix == "AT10_C7_hom":
                            if Gd_identifier == "_Gd154" or Gd_identifier == "_Gd155" or Gd_identifier == "" or Gd_identifier == "_noGd":
                                ERROR_pynjoy_euler[k][j-Nmin]=SERPENT_ALL_pynjoy_euler[k][j]
                    if "oldlib" in S2_libs:
                        print(f"length SERPENT_ALL_oldlib = {len(SERPENT_ALL_oldlib[k])}")
                        ERROR_oldlib[k][j-Nmin]=SERPENT_ALL_oldlib[k][j]
                elif k==1:  # Vecteur Keff --> erreur en pcm
                    if "PyNjoy2016" in S2_libs:
                        print(f"length SERPENT_ALL_pynjoy = {len(SERPENT_ALL_pynjoy[k])}")
                        ERROR_pynjoy[k][j-Nmin]=1.0E+5*(DRAGON_ALL[k][j]-SERPENT_ALL_pynjoy[k][j])
                        if name_mix == "AT10_C7_hom":
                            if Gd_identifier == "_Gd154" or Gd_identifier == "_Gd155" or Gd_identifier == "" or Gd_identifier == "_noGd":
                                ERROR_pynjoy_euler[k][j-Nmin]=1.0E+5*(DRAGON_ALL[k][j]-SERPENT_ALL_pynjoy_euler[k][j])
                    if "oldlib" in S2_libs:
                        print(f"length SERPENT_ALL_oldlib = {len(SERPENT_ALL_oldlib[k])}")
                        ERROR_oldlib[k][j-Nmin]=1.0E+5*(DRAGON_ALL[k][j]-SERPENT_ALL_oldlib[k][j])
                else: # Vecteur isotopique --> erreur en %
                    if "PyNjoy2016" in S2_libs:    
                        if SERPENT_ALL_pynjoy[k][j]==0 :
                            ERROR_pynjoy[k][j-Nmin]=0
                        else:
                            ERROR_pynjoy[k][j-Nmin]=100*(DRAGON_ALL[k][j]-SERPENT_ALL_pynjoy[k][j])/SERPENT_ALL_pynjoy[k][j]
                        if name_mix == "AT10_C7_hom":
                            if Gd_identifier == "_Gd154" or Gd_identifier == "_Gd155" or Gd_identifier == "" or Gd_identifier == "_noGd":
                                if SERPENT_ALL_pynjoy_euler[k][j]==0:
                                    ERROR_pynjoy_euler[k][j-Nmin]=0
                                else:
                                    ERROR_pynjoy_euler[k][j-Nmin]=100*(DRAGON_ALL[k][j]-SERPENT_ALL_pynjoy_euler[k][j])/SERPENT_ALL_pynjoy_euler[k][j]
                    if "oldlib" in S2_libs:
                        if SERPENT_ALL_oldlib[k][j]==0:
                            ERROR_oldlib[k][j-Nmin]=0
                        else:
                            ERROR_oldlib[k][j-Nmin]=100*(DRAGON_ALL[k][j]-SERPENT_ALL_oldlib[k][j])/SERPENT_ALL_oldlib[k][j]

        #print(f"$$$ ------ SHAPE ERROR - {np.shape(ERROR)}")
        #print("$$$ ---------------- ERROR",ERROR)

        #print(f"Len iso DRAGON = {lenISOT_DRAGON}, len iso Serp = {lenISOT_SERPENT}")
    ################################################################
    #             TRACE ET SAUVEGARDE DES FIGURES
    ################################################################

    if visu_DRAGON==1:
        print('$$$ -------- POSTPROC.py : DRAGON5 figures ')

        for k in range(len(isotopes_SOUHAITES+1)):

            plt.figure()
            plt.figure(figsize=SIZE)
            plt.plot(DRAGON_ALL[0],DRAGON_ALL[k+1],'2-',linewidth=1)
            plt.xlabel('BU [MWd/t]')
            plt.grid()
            plt.legend(['DRAGON5'])

            if k == 0: # Comparaison des Keff
                plt.ylabel('Keff')
                save_name=name_mix+'_DRAGON5_Keff'
                fig_name=name_mix+f'{Gd_identifier} - Keff'
            else : # Erreur sur isotopes
                plt.ylabel('Isotopic density [atom/b-cm]')
                save_name=name_mix+f'_DRAGON5_{Gd_identifier}{SAT}_{DEPL_SOL}_{CORR}_'+isotopes_SOUHAITES[k-1]
                fig_name=name_mix+f'_{Gd_identifier} - '+isotopes_SOUHAITES[k-1]

            plt.title(fig_name)
            os.chdir(path+'/'+SAVE_DIR)
            plt.savefig(save_name+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
            os.chdir(path)
            plt.close('all')


    if visu_SERPENT==1:
        print('$$$ -------- POSTPROC.py : Serpent2 figures ')

        for k in range(len(isotopes_SOUHAITES)+1):
            legends = []
            plt.figure()
            plt.figure(figsize=SIZE)
            if "oldlib" in S2_libs:
                legends.append('Serpent2_oldlib')
                plt.plot(SERPENT_ALL_oldlib[0],SERPENT_ALL_oldlib[k+1],'2-',linewidth=1, color='red')
            if "PyNjoy2016" in S2_libs:
                legends.append('Serpent2_pynjoy')
                plt.plot(SERPENT_ALL_pynjoy[0],SERPENT_ALL_pynjoy[k+1],'2-',linewidth=1, color='blue')
                if name_mix == "AT10_C7_hom":
                    if Gd_identifier == "_Gd154" or Gd_identifier == "_Gd155" or Gd_identifier == "" or Gd_identifier == "_noGd":
                        legends.append('Serpent2_pynjoy_euler')
                        plt.plot(SERPENT_ALL_pynjoy_euler[0],SERPENT_ALL_pynjoy_euler[k+1],'2-',linewidth=1, color='green')
            plt.xlabel('BU [MWd/t]')
            plt.grid()
                
            plt.legend(legends)
            if k == 0: # Keff
                plt.ylabel('Keff')
                save_name=name_mix+f'_Serpent2_Keff{Gd_identifier}_serpent'
                fig_name=name_mix+f'{Gd_identifier} - Keff'
            else : # isotopes
                plt.ylabel('Isotopic density [atom/b-cm]')
                save_name=name_mix+f'_Serpent2_{GdCompo}_'+isotopes_SOUHAITES[k-1]
                fig_name=name_mix+f'{Gd_identifier} - '+isotopes_SOUHAITES[k-1]

            plt.title(fig_name)
            os.chdir(path+'/'+SAVE_DIR)
            plt.savefig(save_name+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
            os.chdir(path)
            plt.close('all')


    if visu_COMP==1:
        print('$$$ -------- POSTPROC.py : Comparison DRAGON5 / Serpent2 figures ')

        for k in range(len(isotopes_SOUHAITES)+1):
            legends = []
            plt.figure()
            plt.figure(figsize=SIZE)
            plt.plot(DRAGON_ALL[0],DRAGON_ALL[k+1],'2-',linewidth=1)
            if "oldlib" in S2_libs:
                legends.append('Serpent2_oldlib')
                plt.plot(SERPENT_ALL_oldlib[0],SERPENT_ALL_oldlib[k+1],'2-',linewidth=1, color='red')
                
            if "PyNjoy2016" in S2_libs:
                legends.append('Serpent2_pynjoy')
                plt.plot(SERPENT_ALL_pynjoy[0],SERPENT_ALL_pynjoy[k+1],'2-',linewidth=1, color='blue')
                if name_mix == "AT10_C7_hom":
                    if Gd_identifier == "_Gd154" or Gd_identifier == "_Gd155" or Gd_identifier == "" or Gd_identifier == "_noGd":
                        legends.append('Serpent2_pynjoy_euler')
                        plt.plot(SERPENT_ALL_pynjoy_euler[0],SERPENT_ALL_pynjoy_euler[k+1],'2-',linewidth=1, color='green')
                

            plt.xlabel('BU (MWd/t)')
            plt.grid()
            plt.legend(legends)

            if k == 0: # Comparaison des Keff
                plt.ylabel('Keff')
                save_name=name_mix+f'_COMP{Gd_identifier}_Keff'
                fig_name=name_mix+f'{Gd_identifier} - Keff'
            else : # comparaison isotopes
                plt.ylabel('Isotopic density [atom/b-cm]')
                save_name=name_mix+f'{Gd_identifier}_COMP{SAT}_{DEPL_SOL}_{CORR}_'+isotopes_SOUHAITES[k-1]
                if SAT == "_SAT":
                    fig_name=name_mix+f'{Gd_identifier} {SAT.split("_")[-1]} {CORR}- '+isotopes_SOUHAITES[k-1]
                else:
                    fig_name=name_mix+f'{Gd_identifier} {CORR} - '+isotopes_SOUHAITES[k-1]

            plt.title(fig_name)
            os.chdir(path+'/'+SAVE_DIR)
            plt.savefig(save_name+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
            os.chdir(path)
            plt.close('all')
        

    if visu_DELTA==1 :
        print('$$$ -------- POSTPROC.py : ERROR DRAGON5-Serpent2 figures ')

        for k in range(len(isotopes_SOUHAITES)+1): # -1 lenISOT_DRAGON-1 ?
            legends = []
            plt.figure()
            plt.figure(figsize=SIZE)
            print(f"Plotting for k = {k}")
            print(f"$$$ ---- ERROR_pynjoy[k+1] = {ERROR_pynjoy[k+1]}")
            if "oldlib" in S2_libs:
                legends.append('Serpent2_oldlib')
                plt.plot(ERROR_oldlib[0],ERROR_oldlib[k+1],'2-',linewidth=1, color='red')
            if "PyNjoy2016" in S2_libs:
                legends.append('Serpent2_pynjoy')
                plt.plot(ERROR_pynjoy[0],ERROR_pynjoy[k+1],'2-',linewidth=1, color='blue')
                if name_mix == "AT10_C7_hom":
                    if Gd_identifier == "_Gd154" or Gd_identifier == "_Gd155" or Gd_identifier == "" or Gd_identifier == "_noGd":
                        legends.append('Serpent2_pynjoy_euler')
                        plt.plot(ERROR_pynjoy_euler[0],ERROR_pynjoy_euler[k+1],'2-',linewidth=1, color='darkgreen')
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
                save_name=name_mix+f'{Gd_identifier}{SAT}_{DEPL_SOL}_{CORR}_ERROR_Keff'
                if SAT == "_SAT":
                    fig_name=name_mix+f'{Gd_identifier} {SAT.split("_")[-1]} {DEPL_SOL} {CORR} - \u0394 Keff'
                else:
                    fig_name=name_mix+f'{Gd_identifier} {DEPL_SOL} {CORR} - \u0394 Keff'
            else : # Erreur sur isotopes
                plt.plot([0,ListeCOMPO[-1]],[2,2],'r-.') # limite +2%
                plt.plot([0,ListeCOMPO[-1]],[-2,-2],'r-.') # limite -2%
                plt.ylabel('Relative error (%)')
                save_name=name_mix+f'{Gd_identifier}{SAT}_{DEPL_SOL}_{ssh_module}_{CORR}_ERROR_'+isotopes_SOUHAITES[k-1]
                if SAT == "_SAT":
                    fig_name=name_mix+f'{Gd_identifier} {SAT.split("_")[-1]} {DEPL_SOL} {CORR} - \u0394 '+isotopes_SOUHAITES[k-1]
                else:
                    fig_name=name_mix+f'{Gd_identifier} {DEPL_SOL} {CORR} - \u0394 '+isotopes_SOUHAITES[k-1]

            plt.title(fig_name)
            os.chdir(path+'/'+SAVE_DIR)
            plt.savefig(save_name+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
            os.chdir(path)
            plt.close('all')

        print(f"For SSH : {ssh_module} {ssh_method} {CORR} : The error on Keff at t=0 is {ERROR_pynjoy[1][0]}, Keff_DRAGON = {DRAGON_ALL[1][0]}, Keff_SERPENT = {SERPENT_ALL_pynjoy[1][0]}")
        



























