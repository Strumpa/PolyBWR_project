##############################################################################
#                                                                            #
# Description : Python3 scritp for BWR test cases                            #
# Author      : R. Guasch, adapted from Lucas Fede and Mathias Francois      #
# Date        : August 2024                                                  #
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
import random as rd
import numpy as np
import matplotlib
import matplotlib.cm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from typing import List, Union
import serpentTools as st
from serpentTools.settings import rc

matplotlib.use('Agg')


##############################
########     SCHEMAS DE CALCUL
##############################

AT10_2x2_CLUSTER=['AT10_2x2_UOX', '_UOx6_autop5_PT_MOC_noCalc']

##############################
########    PARAMETRES COMMUNS
##############################

REACTION_TYPE_PARAM=['fission']
print_val_param = True
sym_allow_param = True
form_param = 'png'
n_groups_param = 2

##############################
########    NE PAS TOUCHER !!!
##############################

def BWR_CLUSTER(name_geom_param, name_case_param,reaction_type_param,n_groups_param,bu_to_print,print_val_param,sym_allow_param,form_param):
    ########## PARAMETERS ##########

    # Object to print  :
    # --- False -> display only errors
    # --- True -> display values
    print_val = print_val_param

    # Results condensated on X groups :
    # (1 or 2 available in this version)
    n_groups = n_groups_param

    # Name of the assembly among :
    # --- 'KML2_22AU_ASS'   (also work with KML2_40AU_ASS)
    # --- 'KML2_398GO_ASS'
    # --- 'KML2_39AWU_ASS'
    name_geom = name_geom_param

    # Name of the studied case (name of the DRAGON file):
    # --- example :
    # --- _case1
    # name_case = '_test_MOC_2L_CYPRIEN'
    name_case = name_case_param
    
    # Reactions to study
    # --- 'fission'
    # --- 'gamma'
    reaction_type = reaction_type_param

    # Symmetry  :
    # --- False -> display the full assembly
    # --- True -> display the underlying symmetry (lozenge for 39AWU)
    sym_allow = sym_allow_param

    # Format de sauvegarde des figures
    form = form_param

    # BU to print
    # 0 : first point
    # -1 : last point
    bu = bu_to_print

    ########## LOCAL PARAMETERS (DO NOT CHANGE) ##########

    # Graphical settings (create a scale from white to red - possible to go from bleu to red to differenciate negative differences from positive ones)
    cmap = matplotlib.cm.get_cmap('bwr')

    sym = 'diagonal'
    if name_geom == "AT10_2x2_UOX":
        sym = "none"

    ################## FILE MANAGEMENT AND COMPOS ##################

    # --- Current directory path
    path = os.getcwd()
    #print('$----------------------------'+path)
    
    # --- Creating the results storage directory
    name_DIR = 'DRAGON_ASS_' + name_geom
    a = os.path.exists(name_DIR)
    if not a:
        os.mkdir(name_DIR)
    
    # Access paths to multicompo and serpent results in Resultats_serpent4PyGan
    DRAGON_path = './DRAGON_MULTICOMPO/'
    SERPENT_path = './Resultats_serpent4PyGan/'+name_geom+'/'
    
    # Recovering multi-compos
    name_compo = 'COMPO_' + name_geom + name_case 
    print(f"name_compo = {name_compo}")
    #name_compo = "COMPO_AT10_PIN1"
    os.chdir(DRAGON_path)
    #buCOMPO = np.loadtxt(name_geom + name_case + '_BUvector.txt', dtype=float)
    lenBU_D5 = 1
    pyCOMPO = lcm.new('LCM_INP', name_compo, impx=0)
    os.chdir(path)
    #print(pyCOMPO['EDIBU_2gr']['MIXTURES'][0]['CALCULATIONS'][0]['MIXTURESVOL'])

    # Recuperation results S2
    depFile = st.read(SERPENT_path+name_geom+'_mc_dep.m')
    detFile = st.read(SERPENT_path+name_geom+'_mc_det'+str(bu)+'.m', reader = 'det').detectors

    det_bwr_2x2 = st.readDataFile(SERPENT_path+name_geom+'_mc_det'+str(bu)+'.m') 
    print(det_bwr_2x2.detectors)

    if name_geom == "AT10_2x2_UOX":
        MIXES = [0,1,2,3]
    
    # Isotopes souhaites et nb total d'isotopes
    iso_study = ['U235 ', 'U238 ', 'Pu239', 'Pu241']
    iso_study_S2 = {'U235 ': 0, 'U238 ': 2, 'Pu239': 5, 'Pu241': 7,'Gd155': 22,'Gd157': 24,'Xe135': 18,'Sm149': 30}
    reaction_name = 'NFTOT'
    n_reac = 6
    rates_names = 'fission rates'
    if reaction_type == 'gamma':
        iso_study = ['U235 ', 'U238 ', 'Pu239', 'Pu241','Xe135','Sm149']
        reaction_name = 'NG'
        n_reac = 0
        rates_names = 'radiative capture rates'
    #pyCOMPO['EDIBU'].lib()
    len_isotot = np.shape(pyCOMPO['EDIBU']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESDENS'])[0] - 1

    S2_mat_to_pin_dict = {"pin1": "24UOx_A", "pin2": "32UOx", "pin4": "45UOx_A"}
    ########## CALCULATIONS ##########

    # Creation d'une matrice correspondant à la symetrie ALAMOS
    
    # Calcul du facteur de normalisation
    prodD5 = 0.0
    prodS2 = 0.0
    ListS2 = []    
    nCell = 4
    for iso in range(len_isotot):
        isotope = pyCOMPO['EDIBU']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][iso]['ALIAS'][0:5]
        print(f"isotope = {isotope}, iso number = {iso}")
    for iso in range(len_isotot):
        print(f"iso = {iso}")
        isotope = pyCOMPO['EDIBU']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][iso]['ALIAS'][0:5]
        print(f"isotope = {isotope}")
        if isotope in ['U235 ', 'U238 ', 'Pu239', 'Pu241']:
            for mix in MIXES:
                NWT0 = pyCOMPO['EDIBU_2gr']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NWT0']
                N = pyCOMPO['EDIBU']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESDENS'][iso]
                vol = pyCOMPO['EDIBU']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESVOL'][iso]
                NFTOT = pyCOMPO['EDIBU_2gr']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NFTOT']
#                NGAMMA = pyCOMPO['EDIBU_2gr']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NG']
                for gr in range(2):
#                    prodD5 += 6*NWT0[gr]*(NFTOT[gr]+NGAMMA[gr])*N*vol
                    prodD5 += NWT0[gr]*NFTOT[gr]*N*vol
    for grS2 in range(2):
        for n_iso, isotope in enumerate(iso_study[0:4]):
            print(f"n_iso = {n_iso}, isotope = {isotope}")
            for n,name_cell in enumerate(detFile):
                print(f"HERE name_cell = {name_cell}")
                if n !=0:
                    i_pos = 10*grS2 + 6 + n_iso
                    j_pos = 10
                    #print(valSERPENT[-1])
                    num_cell = ''
                    print(name_cell[:5])
                    for i in name_cell[:5]:
                        #print(i, i=='0'and num_cell!='')
                        L = '123456789'
                        if i in L:
                            num_cell += i
                        elif i == '0' and num_cell != '':
                            num_cell += i
                    #num_cell = '_'+num_cell+'_'
                    print(f"num_cell = {num_cell}")
                    #print(num_cell)
                    n_couronnes = []
                    incr = 0
                    for nb,name_mat in enumerate(depFile):
                        print(f"nb = {nb}, name_mat = {name_mat}")
                        print(f"str(depFile[name_mat]) = {str(depFile[name_mat])}")
                        if S2_mat_to_pin_dict["pin"+num_cell] in str(depFile[name_mat]):
                            print(f"name_mat = {name_mat}")
                            n_couronnes.append(depFile[name_mat]['adens'][iso_study_S2[isotope]][bu]*depFile[name_mat]['volume'][bu])
                            incr += 1
                    n_isot = sum(n_couronnes)
                    #print(n_couronnes)
                    print(f"n_couronnes = {n_couronnes}")
                    print(f"n_isot = {n_isot}")
                    print(f"name cell = {name_cell}, i_pos = {i_pos}, j_pos = {j_pos}")
                    ListS2.append(detFile[name_cell].bins[i_pos][j_pos] * n_isot)
    prodS2 = sum(ListS2)
    print('$$$$$$$$$$$$$$$$$$$$$$$$$$', prodS2,prodD5)


    # Renormalize the results
    # TO DO : Reproduire ce qui a été fait pour 2 groupes
    
    valDRAGON_ALL_gr1=[]
    valDRAGON_ALL_gr2=[]
    valSERPENT_ALL_gr1=[]
    valSERPENT_ALL_gr2=[]
    for iso in range(len_isotot): #len_isotot
        isotope = pyCOMPO['EDIBU_2gr']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][iso]['ALIAS'][0:5]
        if isotope in iso_study:

            # Determine la position de l'isotope pour SERPENT2  
            n_iso = 0
            for n, name_iso in enumerate(iso_study):
                if name_iso == isotope:
                    n_iso = n
            for u in range(1):
                for gr in range(2):

                    # Stocke les valeurs D5
                    valDRAGON = []
                    #print(nCELL)
                    for mix in range(nCell):
                        valDRAGON.append(pyCOMPO['EDIBU_2gr']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso][reaction_name][1-gr])
                        NWT0 = pyCOMPO['EDIBU_2gr']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NWT0'][1-gr]
                        N = pyCOMPO['EDIBU']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESDENS'][iso]
                        vol = pyCOMPO['EDIBU']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESVOL'][iso]
                        valDRAGON[-1] = valDRAGON[-1]*NWT0*N*vol#mix#
                        print(valDRAGON[-1],NWT0,N,vol)
                    print("$ -------------------- isotope=",isotope,"/ gr=",gr," valDRAGON=",valDRAGON)
                    
                    # STOCKAGE VALEURS valDRAGON dans valDRAGON_ALL_gri:
                    if gr == 0 and valDRAGON_ALL_gr1 == []:
                        valDRAGON_ALL_gr1 = valDRAGON
                    elif gr == 1 and valDRAGON_ALL_gr2 == []: 
                        valDRAGON_ALL_gr2 = valDRAGON 
                    elif gr == 0 :
                        print("valDRAGON_ALL_gr1=",valDRAGON_ALL_gr1)
                        print("valDRAGON=",valDRAGON)
                        print(len(valDRAGON))
                        valDRAGON_ALL_gr1 = [valDRAGON_ALL_gr1[i]+valDRAGON[i] for i in range(len(valDRAGON))]
                    elif gr == 1 : 
                        valDRAGON_ALL_gr2 = [valDRAGON_ALL_gr2[i]+valDRAGON[i] for i in range(len(valDRAGON))]
                    print("$ -------------------- isotope=",isotope,"/ gr=",gr," valDRAGON=",valDRAGON,"/valDRAGON_ALL_gr1=",valDRAGON_ALL_gr1,"/valDRAGON_ALL_gr2=",valDRAGON_ALL_gr2)
                    # Stocke les valeurs S2
                    valSERPENT= []
                    for n,name_cell in enumerate(detFile):
                            #print(name_cell)
                            if n !=0:
                                i_pos = 10*grS2 + 6 + n_iso
                                j_pos = 10
                                valSERPENT.append(detFile[name_cell].bins[i_pos][j_pos])
                                #print(valSERPENT[-1])
                                num_cell = ''
                                for i in name_cell[:5]:
                                    L = '123456789'
                                    if i in L:
                                        num_cell += i
                                    elif i == '0' and num_cell != '':
                                        num_cell += i
                                
                                n_couronnes = []
                                incr = 0
                                for nb,name_mat in enumerate(depFile):
                                    if S2_mat_to_pin_dict["pin"+num_cell] in str(depFile[name_mat]):
                                        n_couronnes.append(depFile[name_mat]['adens'][iso_study_S2[isotope]][bu]*depFile[name_mat]['volume'][bu])
                                        incr += 1
                                n_isot = sum(n_couronnes)
                                #print(n_isot)
                                valSERPENT[-1] = valSERPENT[-1] * n_isot
                    #print("$ -------------------- isotope=",isotope,"/ gr=",gr," valSERPENT=",valSERPENT)

                    # STOCKAGE VALEURS valSERPENT dans valSERPENT_ALL_gri
                    if gr == 0 and valSERPENT_ALL_gr1 == []:
                        valSERPENT_ALL_gr1 = valSERPENT
                    elif gr == 1 and valSERPENT_ALL_gr2 == []: 
                        valSERPENT_ALL_gr2 = valSERPENT 
                    elif gr == 0 :
                        valSERPENT_ALL_gr1 = [valSERPENT_ALL_gr1[i]+valSERPENT[i] for i in range(len(valSERPENT))]
                    elif gr == 1 : 
                        valSERPENT_ALL_gr2 = [valSERPENT_ALL_gr2[i]+valSERPENT[i] for i in range(len(valSERPENT))]
                    print("$ -------------------- isotope=",isotope,"/ gr=",gr," valSERPENT=",valSERPENT,"/valSERPENT_ALL_gr1=",valSERPENT_ALL_gr1,"/valSERPENT_ALL_gr2=",valSERPENT_ALL_gr2)

    # Normalisation
    valSERPENT_ALL_gr1 = [i*prodD5/prodS2 for i in valSERPENT_ALL_gr1]
    valSERPENT_ALL_gr2 = [i*prodD5/prodS2 for i in valSERPENT_ALL_gr2]

    print("valSERPENT_ALL_gr1=",valSERPENT_ALL_gr1)
    print("valSERPENT_ALL_gr2=",valSERPENT_ALL_gr2)
    print("valDRAGON_ALL_gr1=",valDRAGON_ALL_gr1)
    print("valDRAGON_ALL_gr2=",valDRAGON_ALL_gr2)

    # Calcul des differences
    diff1 = [100*(valDRAGON_ALL_gr1[i]-valSERPENT_ALL_gr1[i])/valSERPENT_ALL_gr1[i] for i in range(len(valSERPENT_ALL_gr1))]
    diff2 = [100*(valDRAGON_ALL_gr2[i]-valSERPENT_ALL_gr2[i])/valSERPENT_ALL_gr2[i] for i in range(len(valSERPENT_ALL_gr2))]

    print("diff1=",diff1)
    print("diff2=",diff2)
##############################
########    EXECUTION FONCTION
##############################

name_geom_param = AT10_2x2_CLUSTER[0]
NAME_CASE_PARAM=[AT10_2x2_CLUSTER[1]]
BU_TO_PRINT=[0]

for name_case_param in NAME_CASE_PARAM :
    for bu_to_print in BU_TO_PRINT :
        for reaction_type_param in REACTION_TYPE_PARAM :

            print("\n\n$$$ ---------------- BEGIN :")
            print("$$$ ----------------        ASS      =",name_geom_param)
            print("$$$ ----------------        SCHEME   =",name_case_param)
            print("$$$ ----------------        REACTION =",reaction_type_param)
            print("$$$ ----------------        BU       =",bu_to_print)
            BWR_CLUSTER(name_geom_param, name_case_param,reaction_type_param,n_groups_param,bu_to_print,print_val_param,sym_allow_param,form_param)
            print("$$$ ---------------- END")






