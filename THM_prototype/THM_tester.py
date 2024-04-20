# Script used to create instances of THM classes and solve het conduction in fuel rod + convection
# author : R. Guasch
# Purpose : prototyping for further developments in THM module of Donjon5

from THM_MONO import FDM_HeatConductionInFuelPin as FDM_Fuel
from THM_MONO import FVM_ConvectionInCanal as FVM_Canal
import numpy as np
from iapws import IAPWS97

#Cas 1 : base parameters
# Parameters used to create object from FDM_HeatConductioninFuelpin class
Qfiss = 0.3e9 # W/m^3
fuel_radius = 5.6e-3 # m
gap_width = 0.54e-3 # m
clad_width = 0.38e-3 # m
k_fuel = 5 # W/m/K
H_gap = 10000 # W/m^2/K
k_clad = 10 # W/m/K
I_f = 8
I_c = 3

# Paramters used to create object from FVM_ConvectioninCanal class
canal_type = "cylindrical"
canal_width = 0.5e-3 # m
Lf = 2 # m
T_in = 500 # K
Q_flow = 7000 # kg/m^2/s
P_cool = 10.8 #MPa
I_z = 10

initial_water_z0 = IAPWS97(P=P_cool,T=T_in)
h_ini = initial_water_z0.h 
Tsurf = T_in
z=0


def setup_Conduction_case(fuel_radius, I_f, gap_width, clad_width, I_c, z, Qfiss, k_fuel, k_clad, H_gap, T_surf):

    heat_conduction = FDM_Fuel(fuel_radius, I_f, gap_width, clad_width, I_c, Qfiss, k_fuel, k_clad, H_gap)

    for i in range(1,heat_conduction.N_node-1):
        print(f"i in loop= {i}")
        if i<heat_conduction.I_f-1: # setting Aij and Di values for nodes inside the fuel 
            heat_conduction.set_ADi(i, 
                                -heat_conduction.get_Di_half(i-1), 
                                heat_conduction.get_Di_half(i-1)+heat_conduction.get_Di_half(i), 
                                -heat_conduction.get_Di_half(i), 
                                heat_conduction.deltaA_f*heat_conduction.Qfiss)
        elif i==heat_conduction.I_f-1: # setting Aij and Di values for last fuel element
            heat_conduction.set_ADi(i,
                            -heat_conduction.get_Di_half(i-1),
                            (heat_conduction.get_Di_half(i-1)+heat_conduction.get_Ei_gap()),
                            -heat_conduction.get_Ei_gap(),
                            heat_conduction.deltaA_f*heat_conduction.Qfiss)
        elif i==heat_conduction.I_f: # setting Aij and Di values first fuel / gap interface
            heat_conduction.set_ADi(i, 
                                -heat_conduction.get_Ei_gap(), 
                                heat_conduction.get_Ei_gap()+heat_conduction.get_Gi(), 
                                -heat_conduction.get_Gi(), 
                                0)
        elif i==heat_conduction.I_f+1: # setting Aij and Di values second gap / clad interface
            heat_conduction.set_ADi(i, 
                                -heat_conduction.get_Gi(), 
                                heat_conduction.get_Fi_gap()+heat_conduction.get_Gi(), 
                                -heat_conduction.get_Fi_gap(), 
                                0)
        elif i>heat_conduction.I_f+1 : # setting Aij and Di for all elements in the clad, apart from the last one
            print(f"i in interior clad elements is {i}")
            heat_conduction.set_ADi(i, 
                                -heat_conduction.get_Di_half(i-1), 
                                heat_conduction.get_Di_half(i-1)+heat_conduction.get_Di_half(i), 
                                -heat_conduction.get_Di_half(i), 
                                0)
    A0,Am1 = np.zeros(heat_conduction.N_node), np.zeros(heat_conduction.N_node) 
    A0[:2] = [heat_conduction.get_Di_half(0), -heat_conduction.get_Di_half(0)]
    Am1[-2:] = [-heat_conduction.get_Di_half(heat_conduction.N_node-2), heat_conduction.get_Di_half(heat_conduction.N_node-2)+heat_conduction.get_Ei_clad()]
    D0 = heat_conduction.deltaA_f*heat_conduction.Qfiss
    Dm1 = heat_conduction.get_Ei_clad()*T_surf
    heat_conduction.set_CL(A0, Am1, D0, Dm1)

    print(heat_conduction.A_mesh_bounds)
    print(np.sqrt(heat_conduction.A_mesh_bounds/np.pi))
    print(heat_conduction.radii_at_bounds)    
    for row in heat_conduction.A:
        line = "[  "
        for elem in row:
            line+=f"{elem:.3f}   "
        line += "  ]\n"
        print(line)



#setup_Conduction_case(fuel_radius, I_f, gap_width, clad_width, I_c, z, Qfiss, k_fuel, k_clad, H_gap, Tsurf)


convection_test = FVM_Canal(Lf, T_in, Q_flow, P_cool, I_z, canal_type, rf=fuel_radius, rc=fuel_radius+gap_width+clad_width, rw=fuel_radius+gap_width+clad_width+canal_width)
convection_test.set_Fission_Power(0.3e9, variation_type="constant")
for i in range(1,convection_test.N_vol-1):
    convection_test.set_ADi(i,
                            ci=-convection_test.Q_flow/convection_test.dz,
                            ai=convection_test.Q_flow/convection_test.dz,
                            bi=0,
                            di = (np.pi*convection_test.fuel_radius**2/convection_test.A_canal)*convection_test.dz*convection_test.get_Fission_Power()[i])
A0,Am1 = np.zeros(convection_test.N_vol), np.zeros(convection_test.N_vol)
A0[0] = convection_test.Q_flow/convection_test.dz
D0 = np.pi*(convection_test.fuel_radius**2/convection_test.A_canal)*convection_test.dz/2*convection_test.get_Fission_Power()[0]+convection_test.h_z0*convection_test.Q_flow/convection_test.dz
Am1[-2:]=[-convection_test.Q_flow/convection_test.dz, convection_test.Q_flow/convection_test.dz]
Dm1 = (np.pi*convection_test.fuel_radius**2/convection_test.A_canal)*convection_test.dz*convection_test.get_Fission_Power()[-1]
convection_test.set_CL(A0,Am1,D0,Dm1)
for row in convection_test.A:
        line = "[  "
        for elem in row:
            line+=f"{elem:.3f}   "
        line += "  ]\n"
        print(line)
convection_test.h_z = convection_test.solve_h_in_canal()
print(convection_test.h_z)

Tsurf = convection_test.compute_T_surf()
print(convection_test.T_water)
print(Tsurf)
