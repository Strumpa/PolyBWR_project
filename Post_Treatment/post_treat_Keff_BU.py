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

input_file = "Dragon5\\AT-10_pin_NXT_subdivmode_TSPC.result"
BU,Keff_subdiv_NXT_TSPC = extract_Bu_Keff_Dragon(load_data(input_file))
input_file = "Dragon5\\AT-10_pin_NXT_subdivmode_TISO.result"
BU,Keff_subdiv_NXT_TISO = extract_Bu_Keff_Dragon(load_data(input_file))

input_file_subdiv_SALT = "Dragon5\\AT-10_pin_SALT_subdivmode.result"
BU,Keff_subdiv_SALT = extract_Bu_Keff_Dragon(load_data(input_file_subdiv_SALT))

base_input_file = "Dragon5\\AT-10_pin.result"
BU,base_Keff = extract_Bu_Keff_Dragon(load_data(base_input_file))


Bu_renorm = BU/10**3
print(Bu_renorm)


#Serpent2 reference data.
serpent_2_nom_results = "Serpent2\\AT10_pincell_mc_res.m"
Keff_Serp2,stdv = extract_Keff_Serpent(load_data(serpent_2_nom_results))

"""
Options to customize plot
"""
name = "Kinf_AT10_pincell_Compared.png"
title = "Kinf AT-10 Pincell, moderator discretiztion comparison"
colors = ["red", "magenta", "blue", "purple"]
labels = ["AT-10 pincell Kinf NXT, TSPC moderator subdivided", "AT-10 pincell Kinf NXT, TISO moderator subdivided","AT-10 pincell Kinf SALT, moderator subdivided", "AT-10 pincell Kinf SALT, flux on SSH geom"]
markers = ["D","o","X",">"]
mfc = ["red", "magenta","blue", "purple"]
linestyles=["--", "-.","--", "-."]
"""
What to plot 
"""
Dragon5_plotter(name, title, Bu_renorm, [Keff_subdiv_NXT_TSPC, Keff_subdiv_NXT_TISO, Keff_subdiv_SALT, base_Keff], colors, labels, markers, mfc, linestyles)


Serpent2_plotter("Kinf_AT10_pincell_serpent2.png", "Kinf AT-10 Pincell, Serprent2 reference results", Bu_renorm, [Keff_Serp2], stdv, ["purple"], ["Nominal Serpent2 Keff"], ["D"], ["red"], ["-."])


delta_rho_NXT_subdiv_TSPC, stdv_rho1 = compute_reactivity_diff(Keff_subdiv_NXT_TSPC, Keff_Serp2, stdv)
delta_rho_NXT_subdiv_TISO, stdv_rho4 = compute_reactivity_diff(Keff_subdiv_NXT_TISO, Keff_Serp2, stdv)

delta_rho_SALT_subdiv, stdv_rho2 = compute_reactivity_diff(Keff_subdiv_SALT, Keff_Serp2, stdv)
delta_rho_SALT_initial, stdv_rho3 = compute_reactivity_diff(base_Keff, Keff_Serp2, stdv)


error_plotter("Error_reactiv_Dragon5vsSepr2", "Evolution of $\\Delta\\rho$ (pcm) in Burnup", Bu_renorm, [delta_rho_NXT_subdiv_TSPC, delta_rho_NXT_subdiv_TISO, delta_rho_SALT_subdiv, delta_rho_SALT_initial], stdv, colors, labels, markers, mfc, linestyles)

print((1/base_Keff-1/Keff_subdiv_SALT)*1e5)