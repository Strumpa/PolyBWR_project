# Purpose : python script to handle the energy mesh
# Author : R. Guasch
# Date : 12/09/2024
# Idea : open lethargy widths file from VDG: output to the energy mesh
# output a serpent2 set nfg card defining the 295 energy groups structure

import numpy as np
import matplotlib.pyplot as plt
import os
import sys


class energyMeshHandler:
    def __init__(self, lethargyWidthsFile, E0, energyUnits):

        self.lethargyWidthsFile = lethargyWidthsFile
        self.E0 = E0
        self.energyUnits = energyUnits
        self.lethargyMeshWidths = []
        self.lethargyMesh = None
        self.energyMesh = None
        self.energyMeshWidths = None

        self.readLethargyWidths()
        self.computeLethargyMesh()

        #print(f"lethargyMeshWidths = {self.lethargyMeshWidths}")
        #print(f"lethargyMesh = {self.lethargyMesh}")
        #print(f"Number of lethargy mesh = {len(self.lethargyMesh)}")

        #print(f"Number of lethargy groups = {len(self.lethargyMeshWidths)}")

        self.computeEnergyMesh()

        self.energyMesh = np.sort(self.energyMesh)
        #print(f"energyMesh = {self.energyMesh}")
        #print(f"Number of energy groups = {len(self.energyMesh)-1}")
        

    def readLethargyWidths(self):
        # read lethargy widths file
        with open(self.lethargyWidthsFile, 'r') as f:
            lines = f.readlines()
            #self.energyMesh = np.zeros(len(lines))
            #self.energyMeshWidths = np.zeros(len(lines))
            for i in range(len(lines)):
                if "CONDENSED LETHARGY WIDTHS" in lines[i]:
                    print(f"Begin reading lethargy widths from file {self.lethargyWidthsFile}")
                else:
                    line = lines[i].split()
                    for j in range(len(line)):
                        self.lethargyMeshWidths.append(float(line[j]))
        self.lethargyMeshWidths = np.array(self.lethargyMeshWidths)
        return 
    
    def computeLethargyMesh(self):
        self.lethargyMesh = np.zeros(len(self.lethargyMeshWidths))
        for i in range(len(self.lethargyMesh)):
            if i == 0:
                self.lethargyMesh[i] = -0.675
            else:
                if np.abs(self.lethargyMesh[i-1] + self.lethargyMeshWidths[i-1])<1*10**(-12):
                    #print(f"i = {i}")
                    self.lethargyMesh[i] = 0.0
                else:
                    self.lethargyMesh[i] = self.lethargyMesh[i-1] + self.lethargyMeshWidths[i-1]
        return
    
    def computeEnergyMesh(self):
        self.energyMesh = np.zeros(len(self.lethargyMesh))
        for i in range(len(self.lethargyMesh)):
            self.energyMesh[i] = self.E0*np.exp(-self.lethargyMesh[i]) # u = ln(E0/E) ==> E = E0 * exp(-u)
        return
    
    def printnfgCard(self):
        with open('SHEM295_ene.txt', 'w') as f:
            f.write('set ene 1\n')
            for i in range(len(self.energyMesh)):
                f.write(f"{self.energyMesh[i]*1E-6:.5E}  ")
                if ((i+1)%10 == 0) and (i != 0):
                    f.write('\n')
        return


if __name__ == '__main__':
    SHEM_295 = energyMeshHandler('SHEM295.txt', E0=1.0E+07, energyUnits='eV')
    SHEM_295.printnfgCard()
