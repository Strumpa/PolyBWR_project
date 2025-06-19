### Take average of TH data on axial nodes to obtain Teff_fuel, T_m and rho_m for reflector calculation

import numpy as np
import os
import matplotlib.pyplot as plt

# Path to the directory containing the data files
path_to_dir = f"{os.environ['PYGAN_RESULTS']}/multiPhysics_PyGan_24UOX_cell/BiCG/EPRIvoidModel_Churchill_lockhartMartinelli/h380/mesh40_1/Data"
# Effective fuel temperature in Kelvin
Teff_per_nodes = np.loadtxt(f"{path_to_dir}/TeffFuel_24UOX_mesh40_BiCG_EPRIvoidModel_non_relaxedPOW_non_relaxedTH.txt")

# Mean coolant temperature in Kelvin
T_cool_per_nodes = np.loadtxt(f"{path_to_dir}/Twater_24UOX_mesh40_BiCG_EPRIvoidModel_non_relaxedPOW_non_relaxedTH.txt")

# Mean fuel density in kg/m^3
rho_m_per_nodes = np.loadtxt(f"{path_to_dir}/rho_24UOX_mesh40_BiCG_EPRIvoidModel_non_relaxedPOW_non_relaxedTH.txt")

z_space = np.linspace(0, 3.8, len(Teff_per_nodes))

plt.figure(figsize=(10, 6))
plt.plot(z_space, Teff_per_nodes, label='Teff_fuel', marker='o')
plt.plot(z_space, T_cool_per_nodes, label='T_m', marker='o')
plt.plot(z_space, rho_m_per_nodes, label='rho_m', marker='o')
plt.xlabel('Axial Position (m)')
plt.ylabel('Temperature (K) / Density ($kg/m^3$)')
plt.title('Axial Node Data')
plt.legend()
plt.grid()
plt.savefig('Axial_Node_Data.png')
plt.show()


## Compute the mean values for Teff_fuel, T_m and rho_m on 4 axial nodes
n_slices = 4
# Number of axial nodes per slice
n_nodes_per_slice = len(Teff_per_nodes) / n_slices
# Initialize arrays to store the mean values
Teff_mean = np.zeros(n_slices)
T_m_mean = np.zeros(n_slices)
rho_m_mean = np.zeros(n_slices)
# Loop over each slice and compute the mean values
for i in range(n_slices):
    start_index = int(i * n_nodes_per_slice)
    end_index = int((i + 1) * n_nodes_per_slice)
    Teff_mean[i] = np.mean(Teff_per_nodes[start_index:end_index])
    T_m_mean[i] = np.mean(T_cool_per_nodes[start_index:end_index])
    rho_m_mean[i] = np.mean(rho_m_per_nodes[start_index:end_index])
# Print the mean values
print("Mean values for each slice:")
print(f"Teff_fuel: {Teff_mean}")
print(f"T_m: {T_m_mean}")
print(f"rho_m: {rho_m_mean}")

print(f"Inlet T: {T_cool_per_nodes[0]}")
print(f"Inlet rho: {rho_m_per_nodes[0]}")
print(f"Outlet T: {T_cool_per_nodes[-1]}")
print(f"Outlet rho: {rho_m_per_nodes[-1]}")
