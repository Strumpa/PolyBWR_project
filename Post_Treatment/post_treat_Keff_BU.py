import numpy as np
import matplotlib.pyplot as plt


def load_data(input_file):
    """
    input_path : str --> absolute path to data
    input_file : str --> file to be loaded : ex burnup values or keff values
    """
    data = open(input_file,"r")
    lines = data.readlines()
    return lines

def extract_Keff_Serpent(lines):
    """
    Function used to extract Keff values from Serpent _res.m file
    """
    Keff=[]
    stdv=[]
    for line in lines:
        if "ABS_KINF" in line:
            line=line.split("[")
            #print(linesplit[-1].split(" "))
            Keff.append(float(line[-1].split(" ")[2]))
            stdv.append(float(line[-1].split(" ")[3]))
    Keff,stdv = np.array(Keff),np.array(stdv)
    return Keff, stdv

def extract_Bu_Keff_Dragon(lines):
    """
    Function used to extract Bu and Keff lists in numpy array format from Dragon output file
    lines : lines read from output text file
    """
    Burnup=[]
    Keff=[]
    for line in lines:
        if "Burnup=" in line and "ECHO" not in line: #Treating the BU steps > BU = 0
            print(line)
            line=line.split(" ")
            Burnup.append(float(line[3]))
            Keff.append(float(line[7]))
    Burnup, Keff = np.array(Burnup), np.array(Keff)
    return Burnup, Keff

def Dragon5_plotter(name, title, x_data, y_data, colors, labels, markers, mfc, linestyles):
    """
    Plotter using matplotlib used to present results
    name : name of the png file to save
    title : graph title
    x_data : np.array to represent x-axis values on a 2D plot
    y_data : collection of np.arrays to represent y-axis values on a 2D plot, allowing for multiple plots
    to add : custom inupts to choose legends etc
    mec = marker edge color
    mfc = marker face color
    """
    fig, ax = plt.subplots(dpi=500)
    for i in range(len(y_data)):
        ax.plot(x_data, y_data[i], color=colors[i], label=labels[i], marker=markers[i], markersize = 3, mfc = mfc[i], mec=mfc[i],linestyle=linestyles[i], lw=1)
    ax.grid(lw=0.5, color='gray', linestyle="--")        
    ax.set_xlabel("Exposure (MWd/kg)")
    ax.set_ylabel("Keff Dragon5")
    ax.legend(loc="best")
    fig.set_size_inches([8,5])
    ax.set_title(title)
    fig.savefig(name,dpi = 500)

def Serpent2_plotter(name, title, x_data, y_data, stdv, colors, labels, markers, mfc, linestyles):
    """
    Plotter using matplotlib used to present results
    name : name of the png file to save
    title : graph title
    x_data : np.array to represent x-axis values on a 2D plot
    y_data : collection of np.arrays to represent y-axis values on a 2D plot, allowing for multiple plots
    stdv = standard deviations on stochastic outputs to Serpent2
    mec = marker edge color
    mfc = marker face color
    """
    fig, ax = plt.subplots(dpi=500)
    for i in range(len(y_data)):
        ax.errorbar(x_data, y_data[i], yerr=stdv, color=colors[i], label=labels[i], marker=markers[i], markersize = 3, mfc = mfc[i], mec=mfc[i],linestyle=linestyles[i], lw=1)
    ax.grid(lw=0.5, color='gray', linestyle="--")        
    ax.set_xlabel("Exposure (MWd/kgU)")
    ax.set_ylabel("Keff Serpent2")
    ax.legend(loc="best")
    fig.set_size_inches([8,5])
    ax.set_title(title)
    fig.savefig(name,dpi = 500)
    
def error_plotter(name, title, x_data, y_data, stdv, colors, labels, markers, mfc, linestyles):
    fig, ax = plt.subplots(dpi=500)
    for i in range(len(y_data)):
        ax.plot(x_data, y_data[i], color=colors[i], label=labels[i], marker=markers[i], markersize = 2, mfc = mfc[i], mec=mfc[i],linestyle=linestyles[i], lw=1)
    ax.grid(lw=0.5, color='gray', linestyle="--")        
    ax.set_xlabel("Exposure (MWd/kgU)")
    ax.set_ylabel("$\\Delta\\rho$ Dragon5 - Serpent2 (pcm)")
    ax.legend(loc="best")
    fig.set_size_inches([8,5])
    ax.set_title(title)
    fig.savefig(name,dpi = 500)
    
    
def compute_reactivity_diff(Dragon5_kinfs, Serpent2_kinfs, stdv):
    """
    Dragon5_kinfs : np array of Dragon5 Kinf
    Serpent2_kinfs : np array of Serpent2 Kinf
    stdv : np array containing standrad deviations for Serpent2 Kinfs
    """
    delta_rho = ((Dragon5_kinfs-1)/(Dragon5_kinfs)-(Serpent2_kinfs-1)/(Serpent2_kinfs))*100000
    print(delta_rho)
    stdv = (stdv/Serpent2_kinfs)*100000 #converting to pcm values
    print(stdv)
    return delta_rho, stdv

"""
voided_40_results = "Dragon5\\AT-10_pin_void.result"
BU,Keff_void = extract_Bu_Keff_Dragon(load_data(voided_40_results))
Bu_renorm = BU/10**3

Dragon5_plotter("Kinf_AT10_pincell_40void.png","Kinf AT-10 Pincell, 40% void", Bu_renorm, [Keff_void], ["red"], ["Kinf AT-10 Pincell 40% void"], ["D"], ["red"], ["--"])
"""

"""
Loading plottable data
"""

input_file = "Dragon5\\AT-10_pin1_SALT.result"
BU,Keff_SALT_MAV = extract_Bu_Keff_Dragon(load_data(input_file))

Bu_renorm = BU/10**3
print(Bu_renorm)


#Serpent2 reference data.
serpent_2_nom_results = "Serpent2\\AT10_pincell_mc_res.m"
Keff_Serp2,stdv = extract_Keff_Serpent(load_data(serpent_2_nom_results))

"""
Options to customize plot
"""
name = "Kinf_AT10_pincell_Compared.png"
title = "Kinf AT-10 Pincell"
colors = ["red", "magenta","blue", "purple"]
labels = ["Kinf NXT, TSPC, moderator subdivided", "Kinf NXT, TISO, moderator subdivided","Kinf SALT, TSPC, moderator subdivided", "Kinf SALT, TSPC, flux on SSH geom"]
markers = ["D","o","*","X",">"]
mfc = ["red", "magenta","blue", "purple"]
linestyles=["--", "-.","--", "-."]

error_labels = ["Error on "+label for label in labels]

Serpent2_plotter("Kinf_Serp2_AT10_pin_1", "Kinf AT10 pin 1, Serpent2 BU evolution", Bu_renorm, [Keff_Serp2], stdv, ["red"], ["Nominal Kinf, pin 1"], ["x"], ["red"], ["--"])

delta_rho,stdv = compute_reactivity_diff(Keff_SALT_MAV, Keff_Serp2, stdv)

error_plotter("Error_Kinf_AT10_pin1_SALT_MAV", "Error on Kinf of AT10 pin vs Burnup", Bu_renorm, [delta_rho], stdv, ["blue"], ["Error on Kinf, Dragon5 vs Serpent2, windmill disretization"], ["D"], ["red"], ["-."])

def compare_t0_results() :
    """
    t=0 results : checking consistency between SALT, SIBYLT and NXT, PIJ or MOC, single cell calculations

    Keff_ijk, i = SSH option, j = FLX option, k = Solver option
    i,j = [1,2,3] = [SALT, SYBILT, NXT], k = 1,2 = PIJ, MOC

    MAV = windmill discretization, subdiv = annular region in mode
    """
    print(f"Keff Seprent2 at t=0 is : {Keff_Serp2[0]}")
    print("Starting from Lucas and Mathias' work SALT+SALT")
    print("\n")
    print("------------------------- SALT+SALT MOC--------------------------")
    print("\n")
    # SALT MOC : 2 geoms
    Keff_112_MAV = 1.243044 #SALT SSH, SALT FLX, MOC Windmill discretization. 
    delta_112_MAV = (1/Keff_Serp2[0]-1/Keff_112_MAV)*1e5
    
    Keff_112_subdiv = 1.243066 #SALT SSH, SALT FLX, MOC annular discretization. 
    delta_112_subdiv = (1/Keff_Serp2[0]-1/Keff_112_subdiv)*1e5
    
    
    print(f"Error $\\Delta\\rho$ vs Serpent2 for SALT+SALT, MOC on MAV geom is {delta_112_MAV:.1f} pcm")
    print(f"Error $\\Delta\\rho$ vs Serpent2 for SALT+SALT, MOC on annular geom is {delta_112_subdiv:.1f} pcm")
    
    
    print("------------------------- SALT+SALT PIJ--------------------------")
    print("\n")
    # SALT PIJ : 2 geoms
    #MAV
    Keff_111_MAV = 1.243356 #SALT SSH, SALT FLX, PIJ annular discretization. 
    delta_111_MAV = (1/Keff_Serp2[0]-1/Keff_111_MAV)*1e5
    #Annular subdiv
    Keff_111_subdiv = 1.243379 #SALT SSH, SALT FLX, PIJ annular discretization. 
    delta_111_subdiv = (1/Keff_Serp2[0]-1/Keff_111_subdiv)*1e5
    
    
    
    print(f"Error $\\Delta\\rho$ vs Serpent2 for SALT+SALT, PIJ on MAV geom is {delta_111_MAV:.1f} pcm")
    print(f"Error $\\Delta\\rho$ vs Serpent2 for SALT+SALT, PIJ on annular geom is {delta_111_subdiv:.1f} pcm")
    
    print("\n")
    print("------------------------- Now switching to SYBILT+SALT --------------------------")
    print("\n")
    # SYBILT SSH, SALT :
    print("------------------------- SYBILT+SALT MOC --------------------------")
    print("\n")
    #SALT MOC :
    Keff_212_MAV = 1.242813 #SYBILT SSH, SALT FLX, MOC Windmill discretization. 
    delta_212_MAV = (1/Keff_Serp2[0]-1/Keff_212_MAV)*1e5
    
    Keff_212_subdiv = 1.242834 #SYBILT SSH, SALT FLX, MOC annular discretization. 
    delta_212_subdiv = (1/Keff_Serp2[0]-1/Keff_212_subdiv)*1e5
    
    print(f"Error $\\Delta\\rho$ vs Serpent2 for SYBILT+SALT, MOC on MAV geom is {delta_212_MAV:.1f} pcm")
    print(f"Error $\\Delta\\rho$ vs Serpent2 for SYBILT+SALT, MOC on annular geom is {delta_212_subdiv:.1f} pcm")
    
    #SALT PIJ :
    print("\n")
    print("------------------------- SYBILT+SALT PIJ --------------------------")
    print("\n")
    Keff_211_SS = 1.243221 ##Check this out but error greater than PIJ windmill
    Keff_211_MAV = 1.243125 #SYBILT SSH, SALT FLX, PIJ Windmill discretization. 
    delta_211_MAV = (1/Keff_Serp2[0]-1/Keff_211_MAV)*1e5
    
    Keff_211_subdiv = 1.243147 #SYBILT SSH, SALT FLX, PIJ annular discretization. 
    delta_211_subdiv = (1/Keff_Serp2[0]-1/Keff_211_subdiv)*1e5
    
    print(f"Error $\\Delta\\rho$ vs Serpent2 for SYBILT+SALT, PIJ on MAV geom is {delta_211_MAV:.1f} pcm")
    print(f"Error $\\Delta\\rho$ vs Serpent2 for SYBILT+SALT, PIJ on annular geom is {delta_211_subdiv:.1f} pcm")
    
    print("\n")
    print("------------------------- Now switching to SYBILT+SYBILT --------------------------")
    print("\n")
    # SYBILT SSH, SYBILT :
    # Only PIJ available, No MAV available <-- Invalid type of sectorization :
    print("------------------------- SYBILT+SYBILT PIJ --------------------------")
    print("\n")
    
    Keff_221_subdiv = 1.242957 #SYBILT SSH, SYBILT FLX, PIJ annular discretization. 
    delta_221_subdiv = (1/Keff_Serp2[0]-1/Keff_221_subdiv)*1e5
    
    print(f"Error $\\Delta\\rho$ vs Serpent2 for SYBILT+SYBILT, PIJ on annular geom is {delta_221_subdiv:.1f} pcm")
    
    print("\n")
    print("------------------------- Now switching to SYBILT+NXT --------------------------")
    print("\n")
    # SYBILT SSH, NXT Flux:
    # PIJ and MOC available, No MAV available <-- Invalid type of sectorization :
    print("------------------------- SYBILT+NXT PIJ --------------------------")
    print("\n")
    
    Keff_231_subdiv = 1.191973 #SYBILT SSH, NXT FLX, PIJ annular discretization. 
    delta_231_subdiv = (1/Keff_Serp2[0]-1/Keff_231_subdiv)*1e5
    
    print(f"Error $\\Delta\\rho$ vs Serpent2 for SYBILT+NXT, PIJ on annular geom is {delta_231_subdiv:.1f} pcm")
    
    print("------------------------- SYBILT+NXT MOC TSPC --------------------------")
    print("\n")
    
    Keff_232_subdiv = 1.234716 #SYBILT SSH, NXT FLX, MOC annular discretization. 
    delta_232_subdiv = (1/Keff_Serp2[0]-1/Keff_232_subdiv)*1e5
    
    print(f"Error $\\Delta\\rho$ vs Serpent2 for SYBILT+NXT, MOC TSPC on annular geom is {delta_232_subdiv:.1f} pcm")
    
    print("------------------------- SYBILT+NXT MOC TISO --------------------------")
    print("\n")
    
    Keff_232_subdiv_TISO = 1.240296 #SYBILT SSH, NXT FLX, MOC annular discretization. 
    delta_232_subdiv_TISO = (1/Keff_Serp2[0]-1/Keff_232_subdiv_TISO)*1e5
    
    print(f"Error $\\Delta\\rho$ vs Serpent2 for SYBILT+NXT, MOC TISO on annular geom is {delta_232_subdiv_TISO:.1f} pcm")
    
    print("------------------------- SYBILT+NXT PIJ TISO --------------------------")
    print("\n")
    
    Keff_231_subdiv_TISO = 1.240595 #SYBILT SSH, NXT FLX, PIJ annular discretization. 
    delta_231_subdiv_TISO = (1/Keff_Serp2[0]-1/Keff_231_subdiv_TISO)*1e5
    
    print(f"Error $\\Delta\\rho$ vs Serpent2 for SYBILT+NXT, PIJ TISO on annular geom is {delta_231_subdiv_TISO:.1f} pcm")
    
    

