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

def plotter_function(name, title, x_data, y_data, colors, labels, markers, mfc, linestyles):
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



"""
Loading plottable data
"""

input_file = "AT-10_pin_NXT_subdivmode.result"
BU,Keff = extract_Bu_Keff_Dragon(load_data(input_file))

base_input_file = "AT-10_pin.result"
BU,base_Keff = extract_Bu_Keff_Dragon(load_data(base_input_file))

print(BU)
print(Keff)
print(base_Keff)

Bu_renorm = BU/10**3
print(Bu_renorm)

"""
Options to customize plot
"""
name = "Keff_AT10_pincell_Compared.png"
title = "Keff AT-10 Pincell, moderator discretiztion comparison"
colors = ["red", "blue"]
labels = ["AT-10 pincell Keff, moderator subdivided", "Keff AT-10 pincell, flux on SSH geom"]
markers = ["D", "X"]
mfc = ["red", "blue"]
linestyles=["--", "-."]
"""
What to plot 
"""
plotter_function(name, title, Bu_renorm, [Keff, base_Keff], colors, labels, markers, mfc, linestyles)

voided_40_results = "AT-10_pin_void"
BU,Keff = extract_Bu_Keff_Dragon(load_data(voided_40_results))
Bu_renorm = BU/10**3

plotter_function("Keff_AT10_pincell_40void.png","Keff AT-10 Pincell, 40% void",Bu_renorm, [Keff], ["red"], ["Keff AT-10 Pincell 40% void"], ["D"], ["red"], ["--"])
