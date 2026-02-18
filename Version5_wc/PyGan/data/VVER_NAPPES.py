##############################################################################
#                                                                            #
# Description : Python3 scritp for BWR test cases                            #
# Author      : Lucas Fede and Mathias Francois                              #
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
######## CHOIX DE L'ASSEMBLAGE
##############################

EXEC_22AU=1
EXEC_40AU=0
EXEC_398GO=0
EXEC_39AWU=0

##############################
########     SCHEMAS DE CALCUL
##############################

KML2_22AU=['KML2_22AU_ASS','_ALAMOS_evo_REGI_18_30_U2A5_inter3_noCTRA','_REL2005_evo_U2A5_inter3_noCTRA_F3_d']
KML2_40AU=['KML2_40AU_ASS','_ALAMOS_evo_REGI_18_30_U2A5_inter3_noCTRA','_REL2005_evo_REGI_18_30_U2A5_inter3_F3']
KML2_398GO=['KML2_398GO_ASS','_ALAMOS_evo_REGI_18_30_GdA3_Gd6_noCTRA','_REL2005_evo_REGI_18_30_GdA3_Gd6_F2_n2']
KML2_39AWU=['KML2_39AWU_ASS','_ALAMOS_evo_REGI_18_30_GdA3_Gd6_noCTRA','REL2005_evo_REGI_18_30_GdA3_Gd6_new']

##############################
########    PARAMETRES COMMUNS
##############################

REACTION_TYPE_PARAM=['fission','gamma']
print_val_param = True
sym_allow_param = True
form_param = 'png'
n_groups_param = 2

##############################
########    NE PAS TOUCHER !!!
##############################

def vverCOMP_A(name_geom_param, name_case_param,reaction_type_param,n_groups_param,bu_to_print,print_val_param,sym_allow_param,form_param):
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

    sym = 'triangle'
    if name_geom == 'KML2_39AWU_ASS':
        sym = 'lozenge'
    elif name_geom == 'KML2_13ARU_ASS':
        sym = 'simple'

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
    #name_compo = '_COMPO_' + name_geom + name_case
    name_compo = "COMPO_HOM_U5_BOC_fine_autop5_KAPS_SATOFF_USS_PT"
    os.chdir(DRAGON_path)
    #buCOMPO = np.loadtxt(name_geom + name_case + '_BUvector.txt', dtype=float)
    lenBU_D5 = 1
    pyCOMPO = lcm.new('LCM_INP', name_compo, impx=0)
    os.chdir(path)
    #print(pyCOMPO['EDIBU_2gr']['MIXTURES'][0]['CALCULATIONS'][0]['MIXTURESVOL'])

    # Recuperation results S2
    depFile = st.read(SERPENT_path+name_geom+'_dep.m')
    detFile = st.read(SERPENT_path+name_geom+'_det'+str(bu)+'.m', reader = 'det').detectors

    # Indices des cellules extérieures # TO DO : compléter pour les autres symétries
    loz_list = [11*i for i in range(11)]+[111+i for i in range(10)]
    loz_list_ordo = []
    for i in range(10):
        loz_list_ordo.append(11*i)
        loz_list_ordo.append(120-i)
    loz_list_ordo.append(110)
    bords = {'triangle': [0, 10, 11, 20, 21, 29, 30, 37, 38, 44, 45, 50, 51, 55, 56, 59, 60, 62, 63, 64], 'lozenge': loz_list_ordo}

    # Number of cell pins
    nCELL = 61
    if sym == 'lozenge':
        nCELL = 113
    elif sym == 'simple':
        nCELL = 6
    
    # Isotopes souhaites et nb total d'isotopes
    iso_study = ['U235 ', 'U238 ', 'Pu239', 'Pu241']
    iso_study_S2 = {'U235 ': 0, 'U238 ': 2, 'Pu239': 5, 'Pu241': 7,'Gd155': 22,'Gd157': 24,'Xe135': 18,'Sm149': 30}
    reaction_name = 'NFTOT'
    n_reac = 8
    rates_names = 'fission rates'
    if reaction_type == 'gamma':
        iso_study = ['U235 ', 'U238 ', 'Pu239', 'Pu241','Gd155','Gd157','Xe135','Sm149']
        reaction_name = 'NG'
        n_reac = 0
        rates_names = 'radiative capture rates'
    #pyCOMPO['EDIBU'].lib()
    len_isotot = np.shape(pyCOMPO['EDIBU']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESDENS'])[0] - 1

    ########## CALCULATIONS ##########

    # Creation d'une matrice correspondant à la symetrie ALAMOS
    
    # Calcul du facteur de normalisation
    prodD5 = 0.0
    prodS2 = 0.0
    ListS2 = []    
    for iso in range(len_isotot):
        isotope = pyCOMPO['EDIBU']['MIXTURES'][0]['CALCULATIONS'][0]['ISOTOPESLIST'][iso]['ALIAS'][0:5]
        if isotope in ['U235 ', 'U238 ', 'Pu239', 'Pu241']:
            for mix in range(nCELL):
                NWT0 = pyCOMPO['EDIBU_2gr']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NWT0']
                N = pyCOMPO['EDIBU']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESDENS'][iso]
                vol = pyCOMPO['EDIBU']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESVOL'][iso]
                NFTOT = pyCOMPO['EDIBU_2gr']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NFTOT']
#                NGAMMA = pyCOMPO['EDIBU_2gr']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NG']
                for gr in range(2):
#                    prodD5 += 6*NWT0[gr]*(NFTOT[gr]+NGAMMA[gr])*N*vol
                    prodD5 += 6*NWT0[gr]*NFTOT[gr]*N*vol
    for grS2 in range(2):
        for n_iso, isotope in enumerate(iso_study[0:4]):
            for n,name_cell in enumerate(detFile):
                #print(name_cell)
                if n !=0:
                    i_pos = 12*grS2 + 8 + n_iso
                    j_pos = 10
                    #print(valSERPENT[-1])
                    num_cell = ''
                    #print(name_cell[:5])
                    for i in name_cell[:5]:
                        #print(i, i=='0'and num_cell!='')
                        L = '123456789'
                        if i in L:
                            num_cell += i
                        elif i == '0' and num_cell != '':
                            num_cell += i
                    num_cell = '_'+num_cell+'_'
                    #print(num_cell)
                    n_couronnes = []
                    incr = 0
                    for nb,name_mat in enumerate(depFile):
                        if num_cell in str(depFile[name_mat]):
                            #print(name_mat)
                            n_couronnes.append(depFile[name_mat]['adens'][iso_study_S2[isotope]][bu]*depFile[name_mat]['volume'][bu]/6)
                            incr += 1
                    n_isot = sum(n_couronnes)
                    #print(n_couronnes)
                    ListS2.append(detFile[name_cell].bins[i_pos][j_pos] * n_isot)
    prodS2 = 6*sum(ListS2)
    #print('$$$$$$$$$$$$$$$$$$$$$$$$$$', prodS2,prodD5)


    # TO DO : Reproduire ce qui a été fait pour 2 groupes
    if n_groups == 2 :
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
                        for mix in range(nCELL):
                            valDRAGON.append(pyCOMPO['EDIBU_2gr']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso][reaction_name][1-gr])
                            NWT0 = pyCOMPO['EDIBU_2gr']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESLIST'][iso]['NWT0'][1-gr]
                            N = pyCOMPO['EDIBU']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESDENS'][iso]
                            vol = pyCOMPO['EDIBU']['MIXTURES'][mix]['CALCULATIONS'][bu]['ISOTOPESVOL'][iso]
                            valDRAGON[-1] = valDRAGON[-1]*NWT0*N*vol#mix#
                            #print(valDRAGON[-1],NWT0,N,vol)
                        #print("$ -------------------- isotope=",isotope,"/ gr=",gr," valDRAGON=",valDRAGON)
                        
                        # STOCKAGE VALEURS valDRAGON dans valDRAGON_ALL_gri
                        if iso==0:
                            if gr == 0 :
                                valDRAGON_ALL_gr1 = valDRAGON
                            elif gr == 1 : 
                                valDRAGON_ALL_gr2 = valDRAGON 
                        else:
                            if gr == 0 :
                                valDRAGON_ALL_gr1 = [valDRAGON_ALL_gr1[i]+valDRAGON[i] for i in range(len(valDRAGON))]
                            elif gr == 1 : 
                                valDRAGON_ALL_gr2 = [valDRAGON_ALL_gr2[i]+valDRAGON[i] for i in range(len(valDRAGON))]
                        #print("$ -------------------- isotope=",isotope,"/ gr=",gr," valDRAGON=",valDRAGON,"/valDRAGON_ALL_gr1=",valDRAGON_ALL_gr1,"/valDRAGON_ALL_gr2=",valDRAGON_ALL_gr2)
                        # Stocke les valeurs S2
                        valSERPENT= []
                        for n,name_cell in enumerate(detFile):
                             #print(name_cell)
                             if n !=0:
                                 i_pos = 12*gr + n_reac + n_iso
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
                                 num_cell = '_'+num_cell+'_'
                                 n_couronnes = []
                                 incr = 0
                                 for nb,name_mat in enumerate(depFile):
                                     if num_cell in str(depFile[name_mat]):
                                         n_couronnes.append(depFile[name_mat]['adens'][iso_study_S2[isotope]][bu]*depFile[name_mat]['volume'][bu]/6)
                                         incr += 1
                                 n_isot = sum(n_couronnes)
                                 #print(n_isot)
                                 valSERPENT[-1] = valSERPENT[-1] * n_isot
                        #print("$ -------------------- isotope=",isotope,"/ gr=",gr," valSERPENT=",valSERPENT)

                        # STOCKAGE VALEURS valSERPENT dans valSERPENT_ALL_gri
                        if iso==0:
                            if gr == 0 :
                                valSERPENT_ALL_gr1 = valSERPENT
                            elif gr == 1 : 
                                valSERPENT_ALL_gr2 = valSERPENT 
                        else:
                            if gr == 0 :
                                valSERPENT_ALL_gr1 = [valSERPENT_ALL_gr1[i]+valSERPENT[i] for i in range(len(valSERPENT))]
                            elif gr == 1 : 
                                valSERPENT_ALL_gr2 = [valSERPENT_ALL_gr2[i]+valSERPENT[i] for i in range(len(valSERPENT))]
                        #print("$ -------------------- isotope=",isotope,"/ gr=",gr," valSERPENT=",valSERPENT,"/valSERPENT_ALL_gr1=",valSERPENT_ALL_gr1,"/valSERPENT_ALL_gr2=",valSERPENT_ALL_gr2)

                       
        ###################################
        ######## valDRAGON_ALL valSERPENT_ALL
        ###################################
  
        #################
        ######## groupe 1
        #################

        # Normalisation
        valSERPENT_ALL_gr1 = [i*prodD5/prodS2 for i in valSERPENT_ALL_gr1]

        # Cellules identiques et ajout des tubes guides (SERPENT2 ONLY)
        valDRAGON_ALL_gr1_2 = valDRAGON_ALL_gr1[::-1]
        valSERPENT_ALL_gr1 = serpent_complete(valSERPENT_ALL_gr1, sym)
        valDRAGON_ALL_gr1 = valDRAGON_ALL_gr1_2

        # Ajout des tubes guides DRAGON5
        l_low, l_TG, l_TI, l_Gd = pos_pin(name_geom, True)
        nul = l_TG + l_TI
        nul.sort()
        for i in nul:
            valDRAGON_ALL_gr1.insert(i, 0.0)
        #print(len(valDRAGON))
        ind = 0
        for i in range(len(bords[sym])//2):
            somme = valDRAGON_ALL_gr1[bords[sym][ind]]+valDRAGON_ALL_gr1[bords[sym][ind+1]]
            valDRAGON_ALL_gr1[bords[sym][ind]] = somme
            valDRAGON_ALL_gr1[bords[sym][ind+1]] = somme
            ind += 2
        
        # Calcul des erreurs 
        errors = [valDRAGON_ALL_gr1[i]-valSERPENT_ALL_gr1[i] for i in range(len(valSERPENT_ALL_gr1))]
        for i in range(len(errors)):
            if abs(valSERPENT_ALL_gr1[i]) > 1e-5:
                errors[i] = 100*errors[i]/valSERPENT_ALL_gr1[i]
        vmax = max(2, max([abs(i) for i in errors]))
        vmin = -1*vmax

        # Definitions des noms de figure et fichier
        sym_name = sym
        if not sym_allow :
            sym_name = 'hexagonal'
        save_name = name_geom + '_' + sym_name + name_case + '_'+ reaction_name + '_ERROR_' + 'ALL' + '_GR_1' + '_BU_' + str(bu)
        f_name = name_geom + ' - ' + reaction_name + ' - \u0394 ' + 'ALL' + '_GR_1' + '_BU_' + str(bu)

        # Gestion des valeurs à afficher et affichage
        l_low, l_TG, l_TI, l_Gd = pos_pin(name_geom, sym_allow)
        print_values = errors
        if print_val:
             print_values = valDRAGON_ALL_gr1
        plotting(f_name, errors, print_values, name_geom, l_low, l_TI, l_TG, l_Gd, cmap, vmax, sym_allow, sym)

        # Enleve les axes
        ax = plt.gca()
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        # Gestion de la colormap
        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='bwr'), extend='both', shrink=0.85, label='\u00C9cart relatif (%)')

        # Enregistre les figures
        os.chdir(path + '/' + name_DIR)
        # plt.title(f_name)
        plt.savefig(save_name+'.'+form, bbox_inches = 'tight', format=form)
        os.chdir(path)
        plt.close('all')

        #################
        ######## groupe 2
        #################

        # Normalisation
        valSERPENT_ALL_gr2 = [i*prodD5/prodS2 for i in valSERPENT_ALL_gr2]

        # Cellules identiques et ajout des tubes guides (SERPENT2 ONLY)
        valDRAGON_ALL_gr2_2 = valDRAGON_ALL_gr2[::-1]
        valSERPENT_ALL_gr2 = serpent_complete(valSERPENT_ALL_gr2, sym)
        valDRAGON_ALL_gr2 = valDRAGON_ALL_gr2_2

        # Ajout des tubes guides DRAGON5
        l_low, l_TG, l_TI, l_Gd = pos_pin(name_geom, True)
        nul = l_TG + l_TI
        nul.sort()
        for i in nul:
            valDRAGON_ALL_gr2.insert(i, 0.0)
        ind = 0
        for i in range(len(bords[sym])//2):
            somme = valDRAGON_ALL_gr2[bords[sym][ind]]+valDRAGON_ALL_gr2[bords[sym][ind+1]]
            valDRAGON_ALL_gr2[bords[sym][ind]] = somme
            valDRAGON_ALL_gr2[bords[sym][ind+1]] = somme
            ind += 2
        
        # Calcul des erreurs 
        errors = [valDRAGON_ALL_gr2[i]-valSERPENT_ALL_gr2[i] for i in range(len(valSERPENT_ALL_gr2))]
        for i in range(len(errors)):
            if abs(valSERPENT_ALL_gr2[i]) > 1e-5:
                errors[i] = 100*errors[i]/valSERPENT_ALL_gr2[i]
        vmax = max(2, max([abs(i) for i in errors]))
        vmin = -1*vmax

        # Definitions des noms de figure et fichier
        sym_name = sym
        if not sym_allow :
            sym_name = 'hexagonal'
        save_name = name_geom + '_' + sym_name + name_case + '_'+ reaction_name + '_ERROR_' + 'ALL' + '_GR_2' + '_BU_' + str(bu)
        f_name = name_geom + ' - ' + reaction_name + ' - \u0394 ' + 'ALL' + '_GR_2' + '_BU_' + str(bu)

        # Gestion des valeurs à afficher et affichage
        l_low, l_TG, l_TI, l_Gd = pos_pin(name_geom, sym_allow)
        print_values = errors
        if print_val:
             print_values = valDRAGON_ALL_gr2
        plotting(f_name, errors, print_values, name_geom, l_low, l_TI, l_TG, l_Gd, cmap, vmax, sym_allow, sym)

        # Enleve les axes
        ax = plt.gca()
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)

        # Gestion de la colormap
        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
        plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='bwr'), extend='both', shrink=0.85, label='\u00C9cart relatif (%)')

        # Enregistre les figures
        os.chdir(path + '/' + name_DIR)
        # plt.title(f_name)
        plt.savefig(save_name+'.'+form, bbox_inches = 'tight', format=form)
        os.chdir(path)
        plt.close('all')

    else :
        print("\nNUMBERS OF GROUPS NOT SUPPORTED, PLEASE CHOOSE 1 OR 2 GROUPS\n")

# --------------------------------------------------------------------------------#
#                                                                                 #
#                                 ANCILLARY FUNCTIONS                             #
#                                                                                 #
# --------------------------------------------------------------------------------#

def plotting(fig_name, error_list, val_to_print, name_geom, lowUOx, pinTI, pinTG, pinGd, cmap, vmax, sym_allow, sym):
    plt.ion()

    errors = []
    tmp = []
    tmp_bis = []
    pos = 0
    incr = 0
    #pos_sym = [i for i in range(len(error_list))]
    #error_list = pos_sym
    #print(error_list)
    
    # Création de la liste des erreurs sur la géométrie complète à partir de la symétrie
    if sym_allow:
        errors = error_list
    else:
        if sym == 'simple':
            errors += error_list[0:3]
            errors += [error_list[1]] + error_list[3:5] + [error_list[1]]
            tmp = error_list[2:4]
            errors += tmp + [error_list[5]] + tmp[::-1]
            errors += [error_list[1]] + error_list[3:5] + [error_list[1]]
            errors += error_list[0:3]
        elif sym == 'triangle':
            # Moitié inférieure
            for i in range(0,10):
                tmp = [error_list[10-i+sum(range(12-j,12))] for j in range(i)]
                tmp_bis = [error_list[i+sum(range(11-j,11))] for j in range(i)][::-1]
                errors += tmp + [i for i in error_list[incr:incr+11-i]] + tmp_bis
                incr += 11-i
            #print(len(errors))
            # Ligne médiane
            pos = 0
            tmp = [error_list[pos]]
            for i in range(3,12)[::-1]:
                pos += i
                #print(pos)
                tmp.append(error_list[pos])
            tmp += [error_list[pos+2]] + tmp[::-1]
            pos = 0
            #print(len(tmp))
            # Moitié supérieure
            tmp_bis = errors[::-1]
            errors += tmp + tmp_bis
            #print(len(errors))
        else :
            # Moitié inférieure
            for i in range(0,10):
                incr += 10
                tmp = [error_list[20 - i + 11*j] for j in range(i+1)][::-1]              
                errors += tmp + [error_list[10-i]] + [i for i in error_list[incr:incr+11]]
            #print(len(errors))
            # Ligne médiane
            tmp = error_list[0:11]
            tmp_bis = [error_list[i*11] for i in range(11)]
            errors += tmp_bis + tmp
            # Moitié supérieure
            tmp_bis = errors[::-1]
            errors += tmp + tmp_bis
            #print(len(errors))

    #print(len(errors))

    hex_centers, h_ax = create_hex_grid(nx=25, ny=25, crop_hex=sym, sym_allow=sym_allow, do_plot=False)
    # allows to create a lattice of 321 cells depending of the geometry used

    # creating the color scale
    colors = [cmap((i+vmax)/(2*vmax)) for i in errors]
    #print(colors)

    plot_single_lattice_custom_colors(hex_centers[:, 0], hex_centers[:, 1],
                                      face_color=colors,
                                      name_geom=name_geom,
                                      min_diam=1.,
                                      plotting_gap=0.08,
                                      rotate_deg=0,
                                      line_width=0.03,
                                      lowUOx_list=lowUOx,
                                      pinGd_list=pinGd,
                                      pinTI_list=pinTI,
                                      pinTG_list=pinTG
                                      )

    # write cell n° and value of the difference between D5 and S2 results in each hexagonal cell of the lattice
    for i, hexa in enumerate(hex_centers):
        if (not i in pinTI) and (not i in pinTG):
            lab =  str(round(errors[i], 1)) #str(len(hex_centers) - i) + '\n' +
            #if val_to_print != error_list:
            #    lab += '\n'+str(round(val_to_print[i], 2))
            plt.text(hexa[0] - 0.2, hexa[1] - 0.2, lab, size=6)



def pos_pin(name_geom, sym):
    # return lists of positions of the special pins

    lowUOx = []
    pinGd = []
    pinTG = []
    pinTI = []
    if name_geom == 'KML2_39AWU_ASS':
        if not sym:  # positions of different cell types in a hexagonal lattice
            lowUOx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 22, 23, 35, 36, 49, 50, 64, 65, 80, 81, 97, 98, 115, 116,
                      134, 135,
                      154, 155, 175, 176, 195, 196, 214, 215, 232, 233, 249, 250, 265, 266, 280, 281, 294, 295, 307,
                      308, 319,
                      320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 330]
            pinTG = [57, 70, 75, 102, 107, 111, 123, 147, 160, 170, 183, 207, 219, 223, 228, 255, 260, 273]
            pinTI = [165]
            pinGd = [25, 33, 88, 109, 141, 157, 173, 189, 221, 242, 297, 305]
        else:  # positions of different cell types in a lozenge lattice
            lowUOx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 21, 32, 43, 54, 65, 76, 87, 98, 109, 120]
            pinTG = [47, 55, 60, 79, 83, 102, 115]
            pinTI = [110]
            pinGd = [22, 30, 81, 118]
    elif name_geom == 'KML2_22AU_ASS' or name_geom == 'KML2_40AU_ASS':
        if not sym:
            lowUOx = []
            pinTG = [57, 70, 75, 102, 107, 111, 123, 147, 160, 170, 183, 207, 219, 223, 228, 255, 260, 273]
            pinTI = [165]
            pinGd = []
        else:  # triangle symmetry
            lowUOx = []
            pinTG = [41, 45, 50, 58]
            pinTI = [65]
            pinGd = []
    elif name_geom == 'KML2_398GO_ASS':
        if not sym:
            lowUOx = []
            pinTG = [57, 70, 75, 102, 107, 111, 123, 147, 160, 170, 183, 207, 219, 223, 228, 255, 260, 273]
            pinTI = [165]
            pinGd = [88, 109, 141, 189, 221, 242]
        else:  # triangle symmetry
            lowUOx = []
            pinTG = [41, 45, 50, 58]
            pinTI = [65]
            pinGd = [52]
    elif name_geom == 'KML2_13ARU_ASS':
        if not sym:
            pinTI = [9]
        else:
            pinTI = [5]
    return lowUOx, pinTG, pinTI, pinGd


def serpent_complete(valSERPENT, sym):

    if sym == 'triangle':
        valSERPENT.insert(0, 0.0) # Tubes guides
        valSERPENT.insert(6, 0.0)
        valSERPENT.insert(11, 0.0)
        valSERPENT.insert(19, 0.0)
        pos = 0
        for i in range(1,11):
            val = valSERPENT[pos+1] 
            pos += i+1  
            valSERPENT.insert(pos, val)
    elif sym == 'simple':
        valSERPENT.insert(0, 0.0)
        pos = 0
        for i in range(1,3):
            val = valSERPENT[pos+1] 
            pos += i+1  
            valSERPENT.insert(pos, val)
    else :
        valSERPENT.insert(0, 0.0) # Tubes guides
        valSERPENT.insert(5, 0.0)
        valSERPENT.insert(13, 0.0)
        valSERPENT.insert(32, 0.0)
        valSERPENT.insert(36, 0.0)
        valSERPENT.insert(55, 0.0)
        valSERPENT.insert(63, 0.0)
        for i in range(1,11):
            val = valSERPENT[i]   
            valSERPENT.insert(11*i, val)

    # Renumérotation des cellules
    valSERPENT2 = []
    #print(len(valDRAGON2))
    l = len(valSERPENT)
    #print(l)
    if sym == 'triangle':
        valSERPENT2 = valSERPENT[55:66] + valSERPENT[45:55] + valSERPENT[36:45] + valSERPENT[28:36] 
        valSERPENT2 += valSERPENT[21:28] + valSERPENT[15:21] + valSERPENT[10:15] + valSERPENT[6:10]
        valSERPENT2 += valSERPENT[3:6] + valSERPENT[1:3] + [valSERPENT[0]]
    elif sym == 'simple':
        valSERPENT2 = valSERPENT[3:6] + valSERPENT[1:3] + [valSERPENT[0]]
    else:
        for i in range(0,11):
            valSERPENT2 += valSERPENT[l-11*(i+1):l-11*i]
            #print(valSERPENT2)

    return valSERPENT2

def create_hex_grid(nx: int = 4,
                    ny: int = 5,
                    min_diam: float = 1.,
                    n: int = 0,
                    align_to_origin: bool = True,
                    face_color: Union[List[float], str] = None,
                    edge_color: Union[List[float], str] = None,
                    plotting_gap: float = 0.,
                    crop_circ: float = 0.,
                    crop_hex: str = '',
                    sym_allow: str = '',
                    do_plot: bool = False,
                    rotate_deg: float = 0.,
                    keep_x_sym: bool = True,  # Don't work yet
                    h_ax: plt.Axes = None,
                    line_width: float = 0.2,
                    background_color: Union[List[float], str] = None) -> (np.ndarray, plt.Axes):
    coord_x, coord_y = make_grid(nx, ny, min_diam, n, crop_circ, crop_hex, sym_allow, rotate_deg, align_to_origin)

    if do_plot:
        h_ax = plot_single_lattice(coord_x, coord_y, face_color, edge_color, min_diam, plotting_gap, rotate_deg, h_ax,
                                   background_color, line_width)

    return np.hstack([coord_x, coord_y]), h_ax


def plot_single_lattice(coord_x, coord_y, face_color, edge_color, min_diam, plotting_gap, rotate_deg, h_ax=None,
                        background_color=None, line_width=0.2):
    """
    Adds a single lattice to the axes canvas. Multiple calls can be made to overlay few lattices.
    :return:
    """
    if face_color is None:
        face_color = (1, 1, 1, 0)  # Make the face transparent
    if edge_color is None:
        edge_color = 'k'

    if h_ax is None:
        h_fig = plt.figure(figsize=(6, 5))
        h_ax = h_fig.add_axes([0.05, 0.05, 0.9, 0.9])

    if background_color is not None:
        h_ax.set_facecolor(background_color)
    patches = []
    for curr_x, curr_y in zip(coord_x, coord_y):
        polygon = mpatches.RegularPolygon((curr_x, curr_y), numVertices=6,
                                          radius=min_diam / np.sqrt(3) * (1 - plotting_gap),
                                          orientation=np.deg2rad(-rotate_deg))
        patches.append(polygon)
    collection = PatchCollection(patches, edgecolor=edge_color, facecolor=face_color, linewidths=line_width)
    h_ax.add_collection(collection)

    h_ax.set_aspect('equal')
    h_ax.axis([coord_x.min() - 2 * min_diam, coord_x.max() + 2 * min_diam, coord_y.min() - 2 * min_diam,
               coord_y.max() + 2 * min_diam])
    # plt.plot(0, 0, 'r.', markersize=5)   # Add red point at the origin
    return h_ax


def make_grid(nx, ny, min_diam, n, crop_circ, crop_hex, sym_allow, rotate_deg, align_to_origin) -> (np.ndarray, np.ndarray):
    """
    Computes the coordinates of the hexagon centers, given the size rotation and layout specifications
    :return:
    """
    ratio = np.sqrt(3) / 2
    if n > 0:  # n variable overwrites (nx, ny) in case all three were provided
        ny = int(np.sqrt(n / ratio))
        nx = n // ny

    coord_x, coord_y = np.meshgrid(np.arange(nx), np.arange(ny), sparse=False, indexing='xy')
    coord_y = coord_y * ratio
    coord_x = coord_x.astype('float')
    coord_x[1::2, :] += 0.5
    coord_x = coord_x.reshape(-1, 1)
    coord_y = coord_y.reshape(-1, 1)

    coord_x *= min_diam  # Scale to requested size
    coord_y = coord_y.astype('float') * min_diam

    mid_x = (np.ceil(nx / 2) - 1) + 0.5 * (
            np.ceil(ny / 2) % 2 == 0)  # Pick center of some hexagon as origin for rotation or crop...
    mid_y = (np.ceil(ny / 2) - 1) * ratio  # np.median() averages center 2 values for even arrays :\
    mid_x *= min_diam
    mid_y *= min_diam

    # mid_x = (nx // 2 - (nx % 2 == 1)) * min_diam + 0.5 * (ny % 2 == 1)
    # mid_y = (ny // 2 - (ny % 2)) * min_diam * ratio

    if crop_circ > 0:
        rad = ((coord_x - mid_x) ** 2 + (coord_y - mid_y) ** 2) ** 0.5
        coord_x = coord_x[rad.flatten() <= crop_circ, :]
        coord_y = coord_y[rad.flatten() <= crop_circ, :]

    if crop_hex != '':
        # Hexagon settings
        p = 9 / 5.2
        a = 9.2
        b = 18.2
        p_bis = 9 / 5.3

        if crop_hex == 'simple':
            # Bottom hexagon equation
            eq1s = (coord_y - mid_y)
            coord_x = coord_x[eq1s.flatten() >= -2, :]
            coord_y = coord_y[eq1s.flatten() >= -2, :]
            
            if sym_allow:
                eq7 = (coord_y - mid_y) - p_bis * (coord_x - mid_x)
                coord_x = coord_x[eq7.flatten() <= 0.1, :]
                coord_y = coord_y[eq7.flatten() <= 0.1, :]
                eq8 = (coord_y - mid_y) + p_bis * (coord_x - mid_x)
                coord_x = coord_x[eq8.flatten() <= 0.1, :]
                coord_y = coord_y[eq8.flatten() <= 0.1, :]

            else:
                b = 2.4 * p
                # Bottom-right hexagon equation
                eq2 = (coord_y - mid_y) - p * (coord_x - mid_x)
                coord_x = coord_x[eq2.flatten() >= -b, :]
                coord_y = coord_y[eq2.flatten() >= -b, :]
                # Top-right hexagon equation
                eq3 = (coord_y - mid_y) + p * (coord_x - mid_x)
                coord_x = coord_x[eq3.flatten() <= b, :]
                coord_y = coord_y[eq3.flatten() <= b, :]
                # Top hexagon equation
                eq4 = (coord_y - mid_y)
                coord_x = coord_x[eq4.flatten() <= 2, :]
                coord_y = coord_y[eq4.flatten() <= 2, :]
                # Top-left hexagon equation
                eq5 = (coord_y - mid_y) - p * (coord_x - mid_x)
                coord_x = coord_x[eq5.flatten() <= b, :]
                coord_y = coord_y[eq5.flatten() <= b, :]
                # Bottom-left hexagon equation
                eq6 = (coord_y - mid_y) + p * (coord_x - mid_x)
                coord_x = coord_x[eq6.flatten() >= -b, :]
                coord_y = coord_y[eq6.flatten() >= -b, :]
                
        
        else:
            # Bottom hexagon equation
            eq4 = (coord_y - mid_y)
            coord_x = coord_x[eq4.flatten() >= -a, :]
            coord_y = coord_y[eq4.flatten() >= -a, :]
            # Bottom-right hexagon equation
            eq5 = (coord_y - mid_y) - p * (coord_x - mid_x)
            coord_x = coord_x[eq5.flatten() >= -b, :]
            coord_y = coord_y[eq5.flatten() >= -b, :]

            if not sym_allow:
                # Top hexagon equation
                eq1 = (coord_y - mid_y)
                coord_x = coord_x[eq1.flatten() <= a, :]
                coord_y = coord_y[eq1.flatten() <= a, :]
                # Top-left hexagon equation
                eq2 = (coord_y - mid_y) - p * (coord_x - mid_x)
                coord_x = coord_x[eq2.flatten() <= b, :]
                coord_y = coord_y[eq2.flatten() <= b, :]
                # Bottom-left hexagon equation
                eq3 = (coord_y - mid_y) + p * (coord_x - mid_x)
                coord_x = coord_x[eq3.flatten() >= -b, :]
                coord_y = coord_y[eq3.flatten() >= -b, :]
                # Top-right hexagon equation
                eq6 = (coord_y - mid_y) + p * (coord_x - mid_x)
                coord_x = coord_x[eq6.flatten() <= b, :]
                coord_y = coord_y[eq6.flatten() <= b, :]
                
            else:
                eq7 = (coord_y - mid_y) - p_bis * (coord_x - mid_x)
                coord_x = coord_x[eq7.flatten() <= 0.1, :]
                coord_y = coord_y[eq7.flatten() <= 0.1, :]
                if crop_hex == 'triangle':
                    eq8 = (coord_y - mid_y) + p_bis * (coord_x - mid_x)
                    coord_x = coord_x[eq8.flatten() <= 0.1, :]
                    coord_y = coord_y[eq8.flatten() <= 0.1, :]
                elif crop_hex == 'lozenge':
                    eq9 = (coord_y - mid_y)
                    coord_x = coord_x[eq9.flatten() <= 0.1, :]
                    coord_y = coord_y[eq9.flatten() <= 0.1, :]

    if not np.isclose(rotate_deg, 0):  # Check if rotation is not 0, with tolerance due to float format
        # Clockwise, 2D rotation matrix
        RotMatrix = np.array([[np.cos(np.deg2rad(rotate_deg)), np.sin(np.deg2rad(rotate_deg))],
                              [-np.sin(np.deg2rad(rotate_deg)), np.cos(np.deg2rad(rotate_deg))]])
        rot_locs = np.hstack((coord_x - mid_x, coord_y - mid_y)) @ RotMatrix.T
        # rot_locs = np.hstack((coord_x - mid_x, coord_y - mid_y))
        coord_x, coord_y = np.hsplit(rot_locs + np.array([mid_x, mid_y]), 2)

    if align_to_origin:
        coord_x -= mid_x
        coord_y -= mid_y

    return coord_x, coord_y


def plot_single_lattice_custom_colors(coord_x, coord_y, face_color, name_geom, pinGd_list, pinTG_list, pinTI_list,
                                      lowUOx_list, min_diam, plotting_gap, rotate_deg,
                                      line_width=1., h_ax=None, background_color=None, edge_color=[]):
    """
    Plot hexagonal lattice where every hexagon is colored by an individual color.
    All inputs are similar to the plot_single_lattice() except:
    :param line_width:
    :param h_ax:
    :param name_geom:
    :param pinGd_list:
    :param pinTI_list:
    :param pinTG_list:
    :param rotate_deg:
    :param plotting_gap:
    :param min_diam:
    :param coord_y:
    :param coord_x:
    :param face_color: numpy array, Nx3 or Nx4 - Color list of length |coord_x| for each hexagon face.
                                                 Each row is a RGB or RGBA values, e.g. [0.3 0.3 0.3 1]
    :param edge_color: numpy array, Nx3 or Nx4 - Color list of length |coord_x| for each hexagon edge.
                                                 Each row is a RGB or RGBA values, e.g. [0.3 0.3 0.3 1]
    :return:
    """
    if not edge_color:
        edge_color = face_color

    hatch_list = ['' for i in range(len(face_color))]
    for i in pinGd_list:
        hatch_list[i] = '/////////'
        # edge_color[i] = 'k'                        # --- Bug ? Color the full cell in black instead of edge only
    for i in pinTG_list:
        face_color[i] = [0.30, 0.30, 0.30]
    for i in pinTI_list:
        face_color[i] = [0.20, 0.20, 0.20]
    for i in range(len(hatch_list)):
        if i in lowUOx_list:
            hatch_list[i] = 'xxxxxxxxx'

    if h_ax is None:
        h_fig = plt.figure(figsize=(6, 5))
        h_ax = h_fig.add_axes([0.05, 0.05, 0.9, 0.9])
        if background_color is not None:
            h_ax.set_facecolor(background_color)

    for i, (curr_x, curr_y) in enumerate(zip(coord_x, coord_y)):
        polygon = mpatches.RegularPolygon((curr_x, curr_y), numVertices=6,
                                          radius=min_diam / np.sqrt(3) * (1 - plotting_gap),
                                          orientation=np.deg2rad(-rotate_deg),
                                          edgecolor=edge_color[i],
                                          facecolor=face_color[i],
                                          linewidth=line_width)
        if hatch_list[i] != '':
            polygon.set_hatch(hatch_list[i])
            polygon.set_facecolor('k')

        h_ax.add_artist(polygon)

    h_ax.set_aspect('equal')
    h_ax.axis([coord_x.min() - 2 * min_diam, coord_x.max() + 2 * min_diam, coord_y.min() - 2 * min_diam,
               coord_y.max() + 2 * min_diam])
    # plt.plot(0, 0, 'r.', markersize=5)   # Add red point at the origin

    return h_ax


##############################
########    EXECUTION FONCTION
##############################

    print("\n\n$$$ ------ START OF CALCULATIONS : GRAPHICAL VISUALIZATION FOR ASSEMBLIES ------ $$$")

if EXEC_22AU == 1 :

    name_geom_param = KML2_22AU[0]
    NAME_CASE_PARAM=[KML2_22AU[1],KML2_22AU[2]]
    BU_TO_PRINT=[0,48]

    for name_case_param in NAME_CASE_PARAM :
        for bu_to_print in BU_TO_PRINT :
            for reaction_type_param in REACTION_TYPE_PARAM :

                print("\n\n$$$ ---------------- BEGIN :")
                print("$$$ ----------------        ASS      =",name_geom_param)
                print("$$$ ----------------        SCHEME   =",name_case_param)
                print("$$$ ----------------        REACTION =",reaction_type_param)
                print("$$$ ----------------        BU       =",bu_to_print)
                vverCOMP_A(name_geom_param, name_case_param,reaction_type_param,n_groups_param,bu_to_print,print_val_param,sym_allow_param,form_param)
                print("$$$ ---------------- END")

elif EXEC_40AU == 1 :

    name_geom_param = KML2_40AU[0]
    NAME_CASE_PARAM=[KML2_40AU[1],KML2_40AU[2]]
    BU_TO_PRINT=[0,48]

    for name_case_param in NAME_CASE_PARAM :
        for bu_to_print in BU_TO_PRINT :
            for reaction_type_param in REACTION_TYPE_PARAM :
                vverCOMP_A(name_geom_param, name_case_param,reaction_type_param,n_groups_param,bu_to_print,print_val_param,sym_allow_param,form_param)

                print("\n\n$$$ ---------------- BEGIN :")
                print("$$$ ----------------        ASS      =",name_geom_param)
                print("$$$ ----------------        SCHEME   =",name_case_param)
                print("$$$ ----------------        REACTION =",reaction_type_param)
                print("$$$ ----------------        BU       =",bu_to_print)
                vverCOMP_A(name_geom_param, name_case_param,reaction_type_param,n_groups_param,bu_to_print,print_val_param,sym_allow_param,form_param)
                print("$$$ ---------------- END")

elif EXEC_398GO == 1 :

    name_geom_param = KML2_398GO[0]
    NAME_CASE_PARAM=[KML2_398GO[1],KML2_398GO[2]]
    BU_TO_PRINT=[0,79]

    for name_case_param in NAME_CASE_PARAM :
        for bu_to_print in BU_TO_PRINT :
            for reaction_type_param in REACTION_TYPE_PARAM :
                vverCOMP_A(name_geom_param, name_case_param,reaction_type_param,n_groups_param,bu_to_print,print_val_param,sym_allow_param,form_param)

                print("\n\n$$$ ---------------- BEGIN :")
                print("$$$ ----------------        ASS      =",name_geom_param)
                print("$$$ ----------------        SCHEME   =",name_case_param)
                print("$$$ ----------------        REACTION =",reaction_type_param)
                print("$$$ ----------------        BU       =",bu_to_print)
                vverCOMP_A(name_geom_param, name_case_param,reaction_type_param,n_groups_param,bu_to_print,print_val_param,sym_allow_param,form_param)
                print("$$$ ---------------- END")

elif EXEC_39AWU == 1 :

    name_geom_param = KML2_39AWU[0]
    NAME_CASE_PARAM=[KML2_39AWU[1],KML2_39AWU[2]]
    BU_TO_PRINT=[0,11]

    for name_case_param in NAME_CASE_PARAM :
        for bu_to_print in BU_TO_PRINT :
            for reaction_type_param in REACTION_TYPE_PARAM :
                vverCOMP_A(name_geom_param, name_case_param,reaction_type_param,n_groups_param,bu_to_print,print_val_param,sym_allow_param,form_param)

                print("\n\n$$$ ---------------- BEGIN :")
                print("$$$ ----------------        ASS      =",name_geom_param)
                print("$$$ ----------------        SCHEME   =",name_case_param)
                print("$$$ ----------------        REACTION =",reaction_type_param)
                print("$$$ ----------------        BU       =",bu_to_print)
                vverCOMP_A(name_geom_param, name_case_param,reaction_type_param,n_groups_param,bu_to_print,print_val_param,sym_allow_param,form_param)
                print("$$$ ---------------- END")

print("\n\n$$$ ------ END OF CALCULATIONS, PLEASE CHECK RESULTS ------ $$$")



