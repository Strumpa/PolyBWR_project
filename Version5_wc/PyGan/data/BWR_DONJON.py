# PyGan scrip to test DONJON procedures used in multiPhysics solution
# Author: R. Guasch
# Date 19/09/2024

#
# Starting point : simplePOW: a simple multiphysics example with THM: module
#
import lifo
import lcm
import cle2000
from assertS import *
import numpy as np

# construct the Lifo stack for IniDONJON
ipLifo1=lifo.new()
ipLifo1.pushEmpty("Fmap", "LCM") # Fuel Map
ipLifo1.pushEmpty("Matex", "LCM") # Material Indexation
ipLifo1.pushEmpty("Cpo", "LCM") # Compo
ipLifo1.pushEmpty("Track", "LCM") # Tracking data for FEM
ipLifo1.pushEmpty("THData", "LCM") # Thermal Hydraulic data for initialization

# call IniDONJON Cle-2000 procedure
IniDONJON = cle2000.new('IniDONJON',ipLifo1,1)
IniDONJON.exec()
print("IniDONJON execution completed")

# recover the output LCM objects
Fmap = ipLifo1.node("Fmap")
Matex = ipLifo1.node("Matex")
Cpo = ipLifo1.node("Cpo")
Track = ipLifo1.node("Track")
THData = ipLifo1.node("THData")
stateVector = Fmap["STATE-VECTOR"]
mylength = stateVector[0]*stateVector[1]
npar = stateVector[7]

print("Recovered stateVector: ", stateVector)
print("Number of parameters: ", npar)

# empty the Lifo stack
while ipLifo1.getMax() > 0:
  ipLifo1.pop();

# iteration loop
# For now : just run the first iteration with initial guess for TH params

iter = 1 # Iteration number
powi = 0.1722 # Reference at 0.1722 MW from AT10_24UOX test ?

ipLifo2 = lifo.new()
Neutronics = cle2000.new('Neutronics',ipLifo2,1)

conv = False
#while not conv:
#  iter += 1
#  if iter > 1000:
#    raise Exception("DONJON: maximum number of iterations is reached")
    
#  print("DONJON: ITERATION NUMBER:", iter)

    # construct the Lifo stack for Neutronics

ipLifo2.push(Fmap);
ipLifo2.push(Matex);
if iter == 1:
  Flux = ipLifo2.pushEmpty("Flux", "LCM")
  Power = ipLifo2.pushEmpty("Power", "LCM")
else:
  #ipLifo2.push(Flux)
  #ipLifo2.push(Power)
  print("Not iteration yet : first validate the first iteration / initialization")

ipLifo2.push(Cpo)
ipLifo2.push(Track)
ipLifo2.push(THData)
ipLifo2.push(iter)
ipLifo2.push(powi) 

# Call Neutronics component :
Neutronics.exec()
print("Neutronics.c2m execution completed")
Flux = ipLifo2.node("Flux")
Power = ipLifo2.node("Power")
Keff = Flux["K-EFFECTIVE"][0]
print(f"At initialization : Keff = {Keff}")
PowerDistribution = Power["POWER-DISTR"]
print(f"Power distribution : {PowerDistribution}")
#print(f"Uniform TH data used for initialization : {THData["THData"]}")


# empty the ipLifo2 Lifo stack
while ipLifo2.getMax() > 0:
  ipLifo2.pop();



