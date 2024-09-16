#####################################################################
#                                                                   #
# Description : PyGan scritp for VVER simulation with DRAGON5       #
# Author      : L. Fede                                             #
# Date        : 2023                                                #
# Purpose     : Post-processing for global values (Keff,Ni)         #
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
from matplotlib import pyplot
import matplotlib.pyplot as plt
import serpentTools
from serpentTools.settings import rc


def vverCOMP(pyCOMPO,ListeCOMPO,name_geom,VISU_param,use_serpentTools,form,Nmin):

	#########################################################################
	#  Paramètres et initialisation
	#########################################################################

	# Récupération des paramètres de visualisation
	visu_DRAGON=VISU_param[0]
	visu_SERPENT=VISU_param[1]
	visu_COMP=VISU_param[2]
	visu_DELTA=VISU_param[3]

	# --- Paramètres de figures
	SIZE=(6,4)

	# --- Options de tracé des erreurs
	LEGENDE_COMP=['DRAGON5','Serpent2']
	LEGENDE_ERROR=['DRAGON5-Serpent2']

	# --- Chemin du répertoire courant
	path=os.getcwd()
	#print("$$$ ---------- path = ",path)

	# --- Création du répertoire de stockage de résultats 
	a=os.path.exists('VVERresults_PyGan_'+name_geom)
	if a==False:
		os.mkdir('VVERresults_PyGan_'+name_geom)

        # --- Numéro de MIX pour récupérer résultats homogénéisés et condensés 
        NMIX=0

	#########################################################################
	#  Création des matrices de données DRAGON et SERPENT
	#########################################################################


	if visu_DRAGON==1 or visu_COMP==1 or visu_DELTA==1 :

		# --------- matrice DRAGON
		ISOTOPES=pyCOMPO['EDIBU']['MIXTURES'][NMIX]['CALCULATIONS'][0]['ISOTOPESDENS']
		lenISOT_DRAGON=np.shape(ISOTOPES)[0]-1
		lenBU_DRAGON=np.shape(ListeCOMPO)[0]
		print("$$$ ---------------- DRAGON_ISOTOPESDENS shape = ",lenISOT_DRAGON,lenBU_DRAGON)

		DRAGON_BU=ListeCOMPO
		DRAGON_ISOTOPESDENS=np.zeros((lenISOT_DRAGON,lenBU_DRAGON))
		DRAGON_Keff=np.zeros(lenBU_DRAGON)


		for i in range(lenBU_DRAGON):
			DRAGON_Keff[i]=pyCOMPO['EDIBU']['MIXTURES'][NMIX]['CALCULATIONS'][i]['K-EFFECTIVE']	
			for k in range(lenISOT_DRAGON):
				DRAGON_ISOTOPESDENS[k][i]=pyCOMPO['EDIBU']['MIXTURES'][NMIX]['CALCULATIONS'][i]['ISOTOPESDENS'][k]
	

		print('$$$ ---------------- DRAGON_BU =',DRAGON_BU)
		print("$$$ ---------------- DRAGON_Keff = ",DRAGON_Keff)	
		print("$$$ ---------------- DRAGON_ISOTOPESDENS = ",DRAGON_ISOTOPESDENS)


		# --------- Liste des isotopes récupérés dans la Multicompo
		isotopes2=[]
		isotopes=[]
		for k in range(lenISOT_DRAGON):
			isotopes2=isotopes2+[pyCOMPO['EDIBU']['MIXTURES'][NMIX]['CALCULATIONS'][0]['ISOTOPESLIST'][k]['ALIAS']]
		for k in range(lenISOT_DRAGON):
			if isotopes2[k][0]=='U':
				isotopes=isotopes+[isotopes2[k][0:4]]
			else:
				isotopes=isotopes+[isotopes2[k][0:5]]
		print("$$$ ---------------- DRAGON isotopes = ",isotopes)


	if visu_SERPENT==1 or visu_COMP==1 or visu_DELTA==1 :

		# --------- matrice SERPENT
		print("$$$ ---------------- vverCOMP.py : SERPENT results generated with serpentTools from .m")
		# --- Keff
		res=serpentTools.read('Resultats_serpent4PyGan/'+name_geom+"_res.m")
		serpent_keff=res.resdata["absKeff"]
		np.savetxt('serpent_keff.txt',serpent_keff)
		SERPENT_keff=np.loadtxt('serpent_keff.txt',dtype=float)

		# --- BU
		depFile = 'Resultats_serpent4PyGan/'+name_geom+'_dep.m'
		dep = serpentTools.read(depFile)
		fuel=dep.materials['total']
		serpent_BU=fuel.burnup
		np.savetxt('serpent_BU.txt',serpent_BU)
		SERPENT_BU=np.loadtxt('serpent_BU.txt',dtype=float)

		# --- ISOTOPES DENSITIES
		serpent_ISOTOPESDENS=fuel.toDataFrame("adens",names=isotopes)
		np.savetxt('serpent_ISOTOPESDENS.txt',serpent_ISOTOPESDENS)
		SERPENT_ISOTOPESDENS=np.loadtxt('serpent_ISOTOPESDENS.txt',dtype=float)

		SERPENT_ISOTOPESDENS=np.transpose(SERPENT_ISOTOPESDENS)

		Ls1=np.shape(SERPENT_BU)
		print('$$$ ---------------- SERPENT_BU shape =',Ls1)

		Ls2=np.shape(SERPENT_keff)
		lenISOT_SERPENT2=Ls2[0]
		lenBU_SERPENT2=Ls2[1]
		print('$$$ ---------------- SERPENT_keff shape =',Ls2)

		Ls=np.shape(SERPENT_ISOTOPESDENS)
		lenISOT_SERPENT=Ls[0]
		lenBU_SERPENT=Ls[1]
		print('$$$ ---------------- SERPENT_ISOTOPESDENS shape =',Ls)
		
		SERPENT_Keff=np.zeros(lenBU_SERPENT)	
		for k in range(lenBU_SERPENT):
		    SERPENT_BU[k]=1000*SERPENT_BU[k]
		    SERPENT_Keff[k]=SERPENT_keff[k][0]

		print('$$$ ---------------- SERPENT_BU =',SERPENT_BU)
		print("$$$ ---------------- SERPENT_Keff = ",SERPENT_Keff)	
		print("$$$ ---------------- SERPENT_ISOTOPESDENS = ",SERPENT_ISOTOPESDENS)

		if visu_DELTA==1 or visu_COMP==1 :
		# ---------- Vérification de la conformité des dimensions des matrices résultats
			if lenISOT_SERPENT!=lenISOT_DRAGON:
				print('$$$ ---------------- ERROR : nombre isotopes différents entre SERPENT et DRAGON')
				visu_DELTA=0
				visu_COMPO=0
			if lenBU_SERPENT!=lenBU_DRAGON:
				if lenBU_SERPENT>lenBU_DRAGON:
					print('$$$ ---------------- CAUTION : pas suffisament de pas de burn-up DRAGON --> calcul écart à SERPENT impossible')
					visu_DELTA=0
				else:
					print('$$$ ---------------- CAUTION : pas de burn-Up différents entre SERPENT et DRAGON --> utilisation de la fonction interpolation')

	#########################################################################
	#  Calcul des écarts : création matrice ERROR
	#########################################################################
	if visu_DELTA==1:

		lenBU_DELTA=lenBU_SERPENT-Nmin
		DELTA_BU=np.zeros(lenBU_DELTA)
		DELTA_Keff=np.zeros(lenBU_DELTA)
		DELTA_ISOTOPESDENS=np.zeros((lenISOT_SERPENT,lenBU_DELTA))
		
		LS2=np.shape(SERPENT_ISOTOPESDENS_interpol)
		lenISOT_SERPENT_i=LS2[0]
		lenBU_SERPENT_i=LS2[1]
		print('$$$ ---------------- SERPENT_ISOTOPESDENS_interpol shape =',LS2)

		for k in range(lenISOT_SERPENT_i):
			for i in range(Nmin,lenBU_SERPENT_i):
				DELTA_BU[i-Nmin]=SERPENT_BU[i]
				DELTA_Keff[i-Nmin]=1.0E+5*(DRAGON_Keff[i]-SERPENT_Keff[i])
				if SERPENT_ISOTOPESDENS_interpol[k][i]==0:
					DELTA_ISOTOPESDENS[k][i-Nmin]=0
				else:
					DELTA_ISOTOPESDENS[k][i-Nmin]=100*(DRAGON_ISOTOPESDENS[k][i]-SERPENT_ISOTOPESDENS_interpol[k][i])/SERPENT_ISOTOPESDENS_interpol[k][i]
		print("$$$ ---------------- DELTA_ISOTOPESDENS = ",DELTA_ISOTOPESDENS)
		print("$$$ ---------------- DELTA_Keff = ",DELTA_Keff)


	#########################################################################
	#  Tracé des figures
	#########################################################################
	# ------
	# DRAGON
	# ------
	if visu_DRAGON==1:

		if visu_DRAGON_Keff==1:
			plt.figure(figsize=SIZE)
			plt.plot(DRAGON_BU,DRAGON_Keff,'b2-')
			plt.xlabel('BU (GWj/t)')
			plt.ylabel('Keff')
			#plt.legend('DRAGON')
			plt.title('DRAGON_'+name_geom)
			plt.grid()
			os.chdir(path+'/VVERresults_PyGan_'+name_geom)
			plt.savefig(name_geom+'_DRAGON_Keff.'+form,bbox_inches = 'tight', format=form, dpi=1200) #enregistrement des figures dans le répertoire des résultats 
			plt.savefig(name_geom+'_DRAGON_Keff.png',bbox_inches = 'tight') #enregistrement des figures dans le répertoire des résultats 
			os.chdir(path)
			plt.close()

		if visu_DRAGON_ISOTOPESDENS==1:
			for k in range(lenISOT_DRAGON):
				plt.figure(figsize=SIZE)
				plt.plot(DRAGON_BU,DRAGON_ISOTOPESDENS[k],'b2-')
				plt.xlabel('BU (GWj/t)')
				plt.ylabel('Concentration atomique (a/b_cm)')
				#plt.legend('DRAGON')
				plt.title('DRAGON_'+isotopes[k]+'_'+name_geom)
				plt.grid()
				os.chdir(path+'/VVERresults_PyGan_'+name_geom)
				plt.savefig(name_geom+'_DRAGON_'+isotopes[k]+'.'+form,bbox_inches = 'tight', format=form, dpi=1200) #enregistrement des figures dans le répertoire des résultats 
				plt.savefig(name_geom+'_DRAGON_'+isotopes[k]+'.png',bbox_inches = 'tight') #enregistrement des figures dans le répertoire des résultats 
				os.chdir(path)
				plt.close()

		if visu_DRAGON_FISSIONRATES==1:
			print("$$$ --- FISSIONRATES pas encore codé")

	# ------
	# SERPENT
	# ------
	if visu_SERPENT==1:
		
		if visu_SERPENT_Keff==1:
			plt.figure(figsize=SIZE)
			plt.plot(SERPENT_BU,SERPENT_Keff,'r2-')
			plt.xlabel('BU (GWj/t)')
			plt.ylabel('Keff')
			#plt.legend('SERPENT')
			plt.title('SERPENT_'+name_geom)
			plt.grid()
			os.chdir(path+'/VVERresults_PyGan_'+name_geom)
			plt.savefig(name_geom+'_SERPENT_Keff.'+form,bbox_inches = 'tight', format=form, dpi=1200) #enregistrement des figures dans le répertoire des résultats 
			plt.savefig(name_geom+'_SERPENT_Keff.png',bbox_inches = 'tight') #enregistrement des figures dans le répertoire des résultats 
			os.chdir(path)
			plt.close()

		if visu_SERPENT_ISOTOPESDENS==1:
			for k in range(lenISOT_SERPENT):
				plt.figure(figsize=SIZE)
				plt.plot(SERPENT_BU,SERPENT_ISOTOPESDENS[k],'r2-')
				plt.xlabel('BU (GWj/t)')
				plt.ylabel('Concentration atomique (a/b_cm)')
				#plt.legend('SERPENT')
				plt.title('SERPENT_'+isotopes[k]+'_'+name_geom)
				plt.grid()
				os.chdir(path+'/VVERresults_PyGan_'+name_geom)
				plt.savefig(name_geom+'_SERPENT_'+isotopes[k]+'.'+form,bbox_inches = 'tight', format=form, dpi=1200) #enregistrement des figures dans le répertoire des résultats 
				plt.savefig(name_geom+'_SERPENT_'+isotopes[k]+'.png',bbox_inches = 'tight') #enregistrement des figures dans le répertoire des résultats 
				os.chdir(path)
				plt.close()

		if visu_SERPENT_FISSIONRATES==1:
			print("$$$ --- FISSIONRATES pas encore codé")

	# ------
	# SERPENT & DRAGON
	# ------
	if visu_COMP==1:
		
		if visu_COMP_Keff==1:
			plt.figure(figsize=SIZE)
			plt.plot(DRAGON_BU,DRAGON_Keff,'b2-')
			plt.plot(SERPENT_BU,SERPENT_Keff,'r2-')
			plt.xlabel('BU (GWj/t)')
			plt.ylabel('Keff')
			plt.legend(LEGENDES)
			plt.title('DRAGON&SERPENT_'+name_geom)
			plt.grid()
			os.chdir(path+'/VVERresults_PyGan_'+name_geom)
			plt.savefig(name_geom+'_DRAGON&SERPENT_Keff.'+form,bbox_inches = 'tight', format=form, dpi=1200) #enregistrement des figures dans le répertoire des résultats 
			plt.savefig(name_geom+'_DRAGON&SERPENT_Keff.png',bbox_inches = 'tight') #enregistrement des figures dans le répertoire des résultats 
			os.chdir(path)
			plt.close()

		if visu_COMP_ISOTOPESDENS==1 and lenISOT_DRAGON==lenISOT_SERPENT:
			for k in range(lenISOT_SERPENT):
				plt.figure(figsize=SIZE)
				plt.plot(DRAGON_BU,DRAGON_ISOTOPESDENS[k],'b2-')
				plt.plot(SERPENT_BU,SERPENT_ISOTOPESDENS[k],'r2-')
				plt.xlabel('BU (GWj/t)')
				plt.ylabel('Concentration atomique (a/b_cm)')
				plt.legend(LEGENDES)
				plt.title('DRAGON&SERPENT_'+isotopes[k]+'_'+name_geom)
				plt.grid()
				os.chdir(path+'/VVERresults_PyGan_'+name_geom)
				plt.savefig(name_geom+'_DRAGON&SERPENT_'+isotopes[k]+'.'+form,bbox_inches = 'tight', format=form, dpi=1200) #enregistrement des figures dans le répertoire des résultats 
				plt.savefig(name_geom+'_DRAGON&SERPENT_'+isotopes[k]+'.png',bbox_inches = 'tight') #enregistrement des figures dans le répertoire des résultats 
				os.chdir(path)
				plt.close()

		if visu_COMP_FISSIONRATES==1:
			print("$$$ --- FISSIONRATES pas encore codé")


	# ------
	# Ecarts DRAGON-SERPENT
	# ------
	if visu_DELTA==1:

		if visu_DELTA_Keff==1:
			plt.figure(figsize=SIZE)
			plt.plot(DELTA_BU,DELTA_Keff,'kx-')
			plt.plot([DELTA_BU[0],DELTA_BU[lenBU_DELTA-1]],[300,300],'r-.') # limite +300pcm
			plt.plot([DELTA_BU[0],DELTA_BU[lenBU_DELTA-1]],[-300,-300],'r--') # limite -300pcm
			plt.xlabel('BU (GWj/t)')
			plt.ylabel('\u0394Keff (pcm)')
			plt.legend(LEGENDES_ERROR_keff)
			plt.title('\u0394Keff - '+name_geom)
			plt.grid()
			os.chdir(path+'/VVERresults_PyGan_'+name_geom)
			plt.savefig(name_geom+'_DELTAKeff.'+form,bbox_inches = 'tight', format=form, dpi=1200) #enregistrement des figures dans le répertoire des résultats 
			plt.savefig(name_geom+'_DELTAKeff.png',bbox_inches = 'tight') #enregistrement des figures dans le répertoire des résultats
			os.chdir(path)
			plt.close()

		if visu_DELTA_ISOTOPESDENS==1 and lenISOT_DRAGON==lenISOT_SERPENT:
			for k in range(lenISOT_SERPENT):
				plt.figure(figsize=SIZE)
				plt.plot(DELTA_BU,DELTA_ISOTOPESDENS[k],'kx-')
				plt.plot([DELTA_BU[0],DELTA_BU[lenBU_DELTA-1]],[2,2],'r-.') # limite 2%
				plt.plot([DELTA_BU[0],DELTA_BU[lenBU_DELTA-1]],[-2,-2],'r--') # limite -2%
				plt.xlabel('BU (GWj/t)')
				plt.ylabel('Erreur relative (%)')
				plt.legend(LEGENDES_ERROR_iso)
				#plt.axis([0, 3, -0.05, 0.6]) # [xmin, xmax, ymin, ymax]
				plt.title('\u0394'+isotopes[k]+' - '+name_geom)
				plt.grid()
				os.chdir(path+'/VVERresults_PyGan_'+name_geom)
				plt.savefig(name_geom+'_DELTA_'+isotopes[k]+'.'+form,bbox_inches = 'tight', format=form, dpi=1200) #enregistrement des figures dans le répertoire des résultats 
				plt.savefig(name_geom+'_DELTA_'+isotopes[k]+'.png',bbox_inches = 'tight') #enregistrement des figures dans le répertoire des résultats 
				os.chdir(path)
				plt.close()

		if visu_DELTA_FISSIONRATES==1:
			print("$$$ --- FISSIONRATES pas encore codé")



















