# Python3 script to generate the Serpent2 ene card for XS generation
# Author: R. Guasch
# Date: 2024-11-12
# Purpose: Generate the ene card for Serpent2 XS generation
# Derived from the energyMeshHandler.py script in the PTT/energy_mesh directory

import numpy as np

class energyMeshHandler:
    def __init__(self, energy_mesh_name, D5energyMesh, E0, energyUnits):

        self.mesh_name = energy_mesh_name
        self.D5energyMesh = D5energyMesh
        self.E0 = E0
        self.energyUnits = energyUnits
        self.lethargyMeshWidths = []
        self.lethargyMesh = None
        self.energyMesh = None
        self.energyMeshWidths = None
        self.lethargyWidthsFile = None

        self.re_order_energy_mesh()
        self.printnfgCard()

    def re_order_energy_mesh(self):
        self.energyMesh = self.D5energyMesh[::-1]
        print(f"energyMesh = {self.energyMesh}")
        print(f"Number of energy groups = {len(self.energyMesh)-1}")
        return
        

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
                    print(line)
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
                    print(f"i = {i}")
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
        with open(f'{self.mesh_name}_ene.txt', 'w') as f:
            f.write(f'ene {self.mesh_name} 1\n')
            for i in range(len(self.energyMesh)):
                f.write(f"{self.energyMesh[i]*1E-6:.5E}  ")
                if ((i+1)%10 == 0) and (i != 0):
                    f.write('\n')
        return

