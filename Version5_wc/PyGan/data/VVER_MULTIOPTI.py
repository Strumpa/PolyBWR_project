##############################################################################
#                                                                            #
# Description : Python3 scritp for VVER pin cell and assembly optimisation   #
# Author      : L.Fede & C.Richard & M.François                              #
# Date        : 2023                                                         #
#                                                                            #
##############################################################################
#
# ----------------------------------------------------------------------------
#                                  IMPORTS                                   |
# ----------------------------------------------------------------------------
import os, shutil, sys
import pickle
import lifo
import lcm
import cle2000
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
import matplotlib.pyplot as plt
import serpentTools
from serpentTools.settings import rc
#
# ----------------------------------------------------------------------------
#                    PARAMETERS OF RESULTS VISUALISATION                     |
# ----------------------------------------------------------------------------
# Visualisation options
#     = 1 if you want to plot and save the results
#     = 0 either
# NB1 : if visu_XXX_Keff==1 or visu_XXX_BU==1 or visu_XXX_ISOTOPOESDENS==1 then visu_XXX must be equal to 1
visu_COMP=1
visu_DELTA=1
#
#
# ----------------------------------------------------------------------------
#                          PARAMETERS OF CALCULATION                         |
# ----------------------------------------------------------------------------
#

#
# use_serpentTools = 'oui' if you want to recover the results from SERPENT using serpentTools.
#                    Results will be taken from .m outputs of SERPENT.
use_serpentTools='oui'
#
# form = format of the figures generated with vverCOMP.py
#        - 'eps' (recommended)
#        - 'png'
form='eps'
#
# Nmin_XXX = minimum BU step for visualisation and comparison 
Nmin_DELTA=0
#
# name_dir = name of the directory created to store the figures
# Must be 'DRAGON_OPTI' + name
name_DIR='DRAGON_OPTI_'+'HOM_CELL'
#
# name_OPTI = name added to the figures save name
name_OPTI='_HOM_CELL'
#
# ----------------------------------------------------------------------------
#                         LISTS OF CASES TO COMPARE                          |
# ----------------------------------------------------------------------------
# DRAGON_LIST=['case1','case2','case3',...] : list of cases you wan to compare
# SERPENT_LIST=['reference'] : reference case you want to compare DRAGON to
name_geom='Gd_CELLS'

SERPENT_LIST_OVERALL=["HOM_U5", "HOM_U5_U8", "HOM_UOX", "HOM_UOX_clad", "HOM_UOX_clad_noZr", "HOM_UOX_Gd155", "HOM_UOX_Gd157", "HOM_UOX_no155157", "HOM_UOXGd"]

''' SERPENT RESULTS LIST :
--- CELL
KML2_13AU_ASS_new_new
KML2_22AU_ASS_new_new
KML2_40AU_ASS_new_new
KML2_398GO_ASS_new_new
KML2_39AWU_ASS_new_new
--- ASSEMBLY
KML2_22AU_22UOx
KML2_40AU_40UOx
KML2_439GT_44UOx
KML2_439GT_36Gd
KML2_39AWU_40UOx
KML2_39AWU_36UOx
KML2_39AWU_33Gd
'''

DRAGON_LIST_OVERALL=[
#'KML2_22AU_22UOx_LIVRABLE_TEST',
#'HOM_U5_BOC_fine_autop5_KAPS_SATOFF_USS_PT',
'2x2_UOX_UOx_autop5_PT_MOC_CALC4',
#'KML2_39AWU_33Gd_RAPPORT_final_MOC_TSPC_FULLautop_18_30_Gd'
]


''' DRAGON MULTICOMPO LIST :
--- HOM_U5 :
HOM_U5_BOC_fine_autop5_KAPS_AUTO_SUBG
HOM_U5_BOC_fine_autop5_KAPS_USS_PT
HOM_U5_BOC_fine_autop5_KAPS_USS_RSE


HOM_U5_U8_BOC_fine_autop5_KAPS_AUTO_SUBG
HOM_U5_U8_BOC_fine_autop5_KAPS_USS_PT
HOM_U5_U8_BOC_fine_autop5_KAPS_USS_RSE

'''
#
# ----------------------------------------------------------------------------
#                             LOCAL PARAMETERS                               |
# ----------------------------------------------------------------------------
#
# --- Paramètres de figures
SIZE=(6,4)
#
# --- Chemin du repertoire courant
path=os.getcwd()
print("$$$ ---------------- path = ",path)
#
# --- Creation du repertoire de stockage de resultats 
a=os.path.exists(name_DIR)
if a==False:
    os.mkdir(name_DIR)
print("$$$ ---------------- results stored in directory = ",name_DIR)
#
# Access paths to multicompo and serpent results in Resultats_serpent4PyGan
DRAGON_path='./DRAGON_MULTICOMPO/'
SERPENT_path='./Resultats_serpent4PyGan/'
#
isotopes_SOUHAITES=['U235','U238','Pu239','Pu240','Gd155','Gd157','Xe135','Sm149']
#
len_DRAGON_LIST=len(DRAGON_LIST_OVERALL)
print("$$$ ---------------- len_DRAGON_LIST = ",len_DRAGON_LIST)
DRAGON_ALL = [[] for i in range(len_DRAGON_LIST)]
SERPENT_ALL = [[] for i in range(len_DRAGON_LIST)]
ERROR_ALL = [[] for i in range(len_DRAGON_LIST)]
ls1 = [[] for i in range(len_DRAGON_LIST)]
#
LEGENDE_ERROR=DRAGON_LIST_OVERALL




###########################################################################################
#                                                                                         #
#                   RECUPERATION DES DONNEES ET CALCUL DES ERREURS                        #
#                                                                                         #
###########################################################################################

# ----------------------------------------------------------------------------
#                  BOUCLE SUR LE NOMBRE DE CAS A TRACER                      |
# ----------------------------------------------------------------------------
for i in range(len_DRAGON_LIST):

    print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
    print('$$$ ---------------- Pas n° ',i+1,'/',len_DRAGON_LIST)

    if visu_COMP==1 or visu_DELTA==1:

        # ----------------------------------------------------------------------------
        #                  CERATION DE LA MATRICE DES RESULTATS DRAGON               |
        # ----------------------------------------------------------------------------
        name_compo='COMPO_'+DRAGON_LIST_OVERALL[i]
        print("$$$ ---------------- DRAGON name_compo = ",name_compo)

        if 'ASS' in name_compo :
            DIR='EDIBU_HOM'
        else :
            DIR='EDIBU'   

        os.chdir(DRAGON_path)

        ListeCOMPO=np.loadtxt(DRAGON_LIST_OVERALL[i]+'_BUvector.txt',dtype=float)
        lenBU_DRAGON=np.shape(ListeCOMPO)[0]
        pyCOMPO=lcm.new('LCM_INP',name_compo,impx=0)

        os.chdir(path)
         
        ISOTOPES=pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESDENS']
        lenISOT_DRAGON=np.shape(ISOTOPES)[0]-1

        #print("$$$ ---------------- DRAGON_ISOTOPESDENS shape = ",lenISOT_DRAGON,lenBU_DRAGON)
            
        DRAGON_BU=ListeCOMPO
        DRAGON_ISOTOPESDENS=np.zeros((lenISOT_DRAGON,lenBU_DRAGON))
        DRAGON_Keff=np.zeros(lenBU_DRAGON)
            
        for k in range(lenBU_DRAGON):
            DRAGON_Keff[k]=pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][k]['K-EFFECTIVE']
            #print("$$$ ---------------- ISOTOPES BU step ",k,"/",lenBU_DRAGON," = ",pyCOMPO['EDIBU_HOM']['MIXTURES'][0]['CALCULATIONS'][k]['ISOTOPESDENS'])    
            for j in range(lenISOT_DRAGON):
                #print("$$$ ---------------- ISOTOPES ",j,"/",lenISOT_DRAGON," = ",pyCOMPO['EDIBU_HOM']['MIXTURES'][0]['CALCULATIONS'][k]['ISOTOPESDENS'][j])
                DRAGON_ISOTOPESDENS[j][k]=pyCOMPO[DIR]['MIXTURES'][0]['CALCULATIONS'][k]['ISOTOPESDENS'][j]

        #print('$$$ ---------------- DRAGON_BU =',DRAGON_BU)
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
        DRAGON_ALL[i]=[
               DRAGON_BU,
               DRAGON_Keff,
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
               DRAGON_ISOTOPESDENS[int(indices[10]),:],
           ]
        #print("$$$ ---------------- DRAGON_ALL[",i,"= ",DRAGON_ALL[i])


        # ----------------------------------------------------------------------------
        #                  CERATION DE LA MATRICE DES RESULTATS SERPENT               |
        # ----------------------------------------------------------------------------

        name_geo=SERPENT_LIST_OVERALL[i]
        print("$$$ ---------------- SERPENT name_compo = ",name_geo)

        # --- Keff
        res=serpentTools.read(SERPENT_path+name_geo+"_res.m")
        serpent_keff=res.resdata["absKeff"]
        np.savetxt('serpent_keff.txt',serpent_keff)
        SERPENT_keff=np.loadtxt('serpent_keff.txt',dtype=float)
            
        # --- BU
        depFile = SERPENT_path+name_geo+'_dep.m'
        dep = serpentTools.read(depFile)
        fuel=dep.materials['total']
        serpent_BU=fuel.burnup
        np.savetxt('serpent_BU.txt',serpent_BU)
        SERPENT_BU=np.loadtxt('serpent_BU.txt',dtype=float)
        
        # --- ISOTOPES DENSITIES
        serpent_ISOTOPESDENS=fuel.toDataFrame("adens",names=isotopes_SOUHAITES)
        np.savetxt('serpent_ISOTOPESDENS.txt',serpent_ISOTOPESDENS)
        SERPENT_ISOTOPESDENS=np.loadtxt('serpent_ISOTOPESDENS.txt',dtype=float)
        SERPENT_ISOTOPESDENS=np.transpose(SERPENT_ISOTOPESDENS)

        Ls1=np.shape(SERPENT_BU)
        #print('$$$ ---------------- SERPENT_BU shape =',Ls1)

        Ls2=np.shape(SERPENT_keff)
        lenISOT_SERPENT2=Ls2[0]
        lenBU_SERPENT2=Ls2[1]
        #print('$$$ ---------------- SERPENT_keff shape =',Ls2)

        Ls=np.shape(SERPENT_ISOTOPESDENS)
        lenISOT_SERPENT=Ls[0]
        lenBU_SERPENT=Ls[1]
        #print('$$$ ---------------- SERPENT_ISOTOPESDENS shape =',Ls)
        
        # Modification unite BU pour match avec DRAGON
        SERPENT_Keff=np.zeros(lenBU_SERPENT)    
        for k in range(lenBU_SERPENT):
            SERPENT_BU[k]=1000*SERPENT_BU[k]
            SERPENT_Keff[k]=SERPENT_keff[k][0]

        #print('$$$ ---------------- SERPENT_BU =',SERPENT_BU)
        #print("$$$ ---------------- SERPENT_Keff = ",SERPENT_Keff)    
        #print("$$$ ---------------- SERPENT_ISOTOPESDENS = ",SERPENT_ISOTOPESDENS)

        SERPENT_ALL[i]=[
           SERPENT_BU,
           SERPENT_Keff,
           SERPENT_ISOTOPESDENS[0,:],
           SERPENT_ISOTOPESDENS[1,:],
           SERPENT_ISOTOPESDENS[2,:],
           SERPENT_ISOTOPESDENS[3,:],
           SERPENT_ISOTOPESDENS[4,:],
           SERPENT_ISOTOPESDENS[5,:],
           SERPENT_ISOTOPESDENS[6,:],
           SERPENT_ISOTOPESDENS[7,:],
           SERPENT_ISOTOPESDENS[8,:],
           SERPENT_ISOTOPESDENS[9,:],
           SERPENT_ISOTOPESDENS[10,:],
           ]
        #print("$$$ ---------------- SERPENT_ALL[",i,"] = ",SERPENT_ALL[i])
  

if visu_DELTA==1 or visu_COMP==1 :
    #print("$$$ ---------------- DRAGON_ALL = ",DRAGON_ALL)
    LD=np.shape(DRAGON_ALL)
    print("$$$ ---------------- shape DRAGON_ALL = ",LD)

    #print("$$$ ---------------- SERPENT_ALL = ",SERPENT_ALL)
    LS=np.shape(SERPENT_ALL)    
    print("$$$ ---------------- shape SERPENT_ALL = ",LS)



# ----------------------------------------------------------------------------
#                       CERATION DE LA MATRICE DES ERREURS                    |
# ----------------------------------------------------------------------------

if visu_DELTA==1 :

    for i in range(len_DRAGON_LIST):

        ERROR=[[0] * LS[2] for _ in range (LS[1])]

        #LEr=np.shape(ERROR)
        #print('$$$ ------------------------ ERROR shape=',LEr)

        for k in range(LS[1]):
            for j in range(Nmin_DELTA,LS[2]):
                #print('$$$ ----------------------- i=',i,'    k=',k,'    j=',j)
                #print('$$$ ----------------------- SERPENT_ALL[i][k][j]=',SERPENT_ALL[i][k][j])
                if k==0: # Vecteur BU
                    ERROR[k][j-Nmin_DELTA]=SERPENT_ALL[i][k][j]
                elif k==1:  # Vecteur Keff --> erreur en pcm
                    ERROR[k][j-Nmin_DELTA]=1.0E+5*(DRAGON_ALL[i][k][j]-SERPENT_ALL[i][k][j])
                else: # Vecteur isotopique --> erreur en %
                    if SERPENT_ALL[i][k][j]==0:
                        ERROR[k][j-Nmin_DELTA]=0
                    else:
                        ERROR[k][j-Nmin_DELTA]=100*(DRAGON_ALL[i][k][j]-SERPENT_ALL[i][k][j])/SERPENT_ALL[i][k][j]
        #print("$$$ ---------------- ERROR[",l,"] = ",ERROR[l])
        ERROR_ALL[i]=ERROR

    #print("$$$ ---------------- ERROR_ALL = ",ERROR_ALL)
    LE=np.shape(ERROR_ALL)
    print("$$$ ---------------- shape ERROR_ALL = ",LE)


###########################################################################################
#                                                                                         #
#                          TRACE ET SAUVEGARDE DES FIGURES                                #
#                                                                                         #
###########################################################################################

if visu_COMP==1:
    print('$$$ ------ VVER_MULTIOPTI.py : Comparison figures')

    for k in range(LS[1]-1):
        if k == 0: # Comparaison des Keff
            plt.figure()
            plt.figure(figsize=SIZE)
            LEGENDE_COMP=[]
            for l in range(len_DRAGON_LIST):
                if l >= 1 : 
                    if DRAGON_LIST_OVERALL[l]!=DRAGON_LIST_OVERALL[l-1]:
                        plt.plot(DRAGON_ALL[l][0],DRAGON_ALL[l][k+1],'2-',linewidth=1)
                        LEGENDE_COMP=LEGENDE_COMP+[str(DRAGON_LIST_OVERALL[l])]
                    plt.plot(SERPENT_ALL[l][0],SERPENT_ALL[l][k+1],'2-',linewidth=1)
                    LEGENDE_COMP=LEGENDE_COMP+[str(SERPENT_LIST_OVERALL[l])]
                else :
                    plt.plot(DRAGON_ALL[l][0],DRAGON_ALL[l][k+1],'2-',linewidth=1)
                    plt.plot(SERPENT_ALL[l][0],SERPENT_ALL[l][k+1],'2-',linewidth=1)
                    LEGENDE_COMP=LEGENDE_COMP+[str(DRAGON_LIST_OVERALL[l])]
                    LEGENDE_COMP=LEGENDE_COMP+[str(SERPENT_LIST_OVERALL[l])]
            plt.xlabel('BU (MWj/t)')
            plt.ylabel('Keff')
            plt.grid()
            plt.legend(LEGENDE_COMP)
            save_name=name_geom+name_OPTI+'_COMP_Keff'
            fig_name=name_geom+' - Keff'
            #plt.title(fig_name)
            os.chdir(path+'/'+name_DIR)
            plt.savefig(save_name+'.'+form,bbox_inches = 'tight', format=form, dpi=1200) #enregistrement des figures dans le repertoire des resultats
            plt.savefig(save_name+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
            os.chdir(path)
            plt.close('all')
        
        else : # Erreur sur isotopes
            plt.figure()
            plt.figure(figsize=SIZE)
            for l in range(len_DRAGON_LIST):
                if l >= 1 : 
                    if DRAGON_LIST_OVERALL[l]!=DRAGON_LIST_OVERALL[l-1]:
                        plt.plot(DRAGON_ALL[l][0],DRAGON_ALL[l][k+1],'2-',linewidth=1)
                    plt.plot(SERPENT_ALL[l][0],SERPENT_ALL[l][k+1],'2-',linewidth=1)
                else :
                    plt.plot(DRAGON_ALL[l][0],DRAGON_ALL[l][k+1],'2-',linewidth=1)
                    plt.plot(SERPENT_ALL[l][0],SERPENT_ALL[l][k+1],'2-',linewidth=1)
            plt.xlabel('BU (MWj/t)')
            plt.ylabel('Concentration atomique (a/barn.cm)')
            plt.grid()
            plt.legend(LEGENDE_COMP)
            save_name=name_geom+name_OPTI+'_COMP_'+isotopes_SOUHAITES[k-1]
            fig_name=name_geom+' - '+isotopes_SOUHAITES[k-1]
            #plt.title(fig_name)
            os.chdir(path+'/'+name_DIR)
            plt.savefig(save_name+'.'+form,bbox_inches = 'tight', format=form, dpi=1200) #enregistrement des figures dans le repertoire des resultats
            plt.savefig(save_name+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats 
            os.chdir(path)
            plt.close('all')


if visu_DELTA==1 :
    print('$$$ ------ VVER_MULTIOPTI.py : DELTA figures')

    for k in range(LS[1]-1):
        if k == 0: # Erreur sur Keff
            plt.figure()
            plt.figure(figsize=SIZE)
            for l in range(len_DRAGON_LIST):
                plt.plot(ERROR_ALL[l][0],ERROR_ALL[l][k+1],'2-',linewidth=1)
            plt.plot([0,60000],[300,300],'r-.') # limite +300pcm
            plt.plot([0,60000],[-300,-300],'r-.') # limite -300pcm
            plt.xlabel('BU (MWj/t)')
            plt.ylabel('\u0394 Keff (pcm)')
            plt.grid()
            plt.legend(LEGENDE_ERROR)
            save_name=name_geom+name_OPTI+'_ERROR_Keff'
            fig_name=name_geom+' - \u0394 Keff'
            #plt.title(fig_name)
            os.chdir(path+'/'+name_DIR)
            plt.savefig(save_name+'.'+form,bbox_inches = 'tight', format=form, dpi=1200) #enregistrement des figures dans le repertoire des resultats
            plt.savefig(save_name+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats
            os.chdir(path)
            plt.close('all')
        
        else : # Erreur sur isotopes
            plt.figure()
            plt.figure(figsize=SIZE)
            for l in range(len_DRAGON_LIST):
                plt.plot(ERROR_ALL[l][0],ERROR_ALL[l][k+1],'2-',linewidth=1)
            plt.plot([0,60000],[2,2],'r-.') # limite +2%
            plt.plot([0,60000],[-2,-2],'r-.') # limite -2%
            plt.xlabel('BU (MWj/t)')
            plt.ylabel('Erreur relative (%)')
            plt.grid()
            plt.legend(LEGENDE_ERROR)
            save_name=name_geom+name_OPTI+'_ERROR_'+isotopes_SOUHAITES[k-1]
            fig_name=name_geom+' - \u0394 '+isotopes_SOUHAITES[k-1]
            #plt.title(fig_name)
            os.chdir(path+'/'+name_DIR)
            plt.savefig(save_name+'.'+form,bbox_inches = 'tight', format=form, dpi=1200) #enregistrement des figures dans le repertoire des resultats
            plt.savefig(save_name+'.png',bbox_inches = 'tight') #enregistrement des figures dans le repertoire des resultats 
            os.chdir(path)
            plt.close('all')


