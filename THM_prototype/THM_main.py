# Script used to create instances of THM classes and solve convection + conduction in fuel rod.
# author : R. Guasch
# Purpose : prototyping for further developments in THM module of Donjon5

from THM_MONO import FDM_HeatConductionInFuelPin as FDM_Fuel
from THM_MONO import FVM_ConvectionInCanal_MONO as FVM_Canal_MONO
import numpy as np
from iapws import IAPWS97
import matplotlib.pyplot as plt

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

# Paramters used to create object from FVM_ConvectionInCanal class
canal_type = "cylindrical"
canal_width = 0.5e-3 # m
Lf = 2 # m
T_in = 500 # K
Q_flow = 7000 # kg/m^2/s
P_cool = 10.8 #MPa
I_z = 10

initial_water_z0 = IAPWS97(P=P_cool,T=T_in)
h_ini = initial_water_z0.h 
z=0

rw=fuel_radius+gap_width+clad_width+canal_width

class Version5_THM_prototype:
    def __init__(self, case_name, canal_radius, canal_type, fuel_rod_length, T_inlet, P_inlet, Q_flow, I_z, Qfiss, Qfiss_variation_type, fuel_radius, gap_radius, clad_radius, k_fuel, H_gap, k_clad, I_f, I_c, plot_at_z):
        """
        Main constructor for THM case, first set of parameters correspond to canal properties, second set to fuel/gap/clad properties
        The structure followed is : 
        In FVM_ConvectionInCanal class : use a finite volume method to solve heat convection in the canal, then use the Dittus-Boelter correlation to obtain the convective heat transfer coef 
        between the water and the fuel rod's outer surface. This allows to solve for the temperature at this outer surface. 
        Then in the FDM_HeatConductionInFuelPin class, solve for the heat conduction using MCFD method. Compute temperature at the center of the fuel rod.
        Options to plot results can be activated giving an array of z values at which the results should be plotted.
        """
        self.name = case_name
        # canal attributes
      
        self.r_w = canal_radius # outer canal radius (m) if type is cylindrical, if type = square rw is the radius of inscribed circle in the square canal, ie half the square's side.
        self.canal_type = canal_type # cylindrical or square, used to determine the cross sectional flow area in the canal and the hydraulic diameter
        self.Lf = fuel_rod_length # fuel rod length in m
        self.T_in = T_inlet # inlet water temperature K
        self.P_cool = P_inlet # coolant pressure in MPa, assumed to be constant along the axial profile.
        self.Q_flow = Q_flow #  mass flux in kg/m^2/s, assumed to be constant along the axial profile.
        self.I_z = I_z # number of mesh elements on axial mesh

        self.Q_fiss_amp = Qfiss # amplitude of sine variation, or constant value if Qfiss_variation_type = "constant"
        self.Q_fiss_variation_type = Qfiss_variation_type # allows for a sine / cosine axial variation of the fuel power density in convection case.


        self.r_f = fuel_radius # fuel pin radius in meters
        self.gap_r = gap_radius # gap radius in meters, used to determine mesh elements for constant surface discretization
        self.clad_r = clad_radius # clad radius in meters, used to determine mesh elements for constant surface discretization
        self.k_fuel = k_fuel # thermal conductivity coefficient in fuel W/m/K
        self.H_gap = H_gap # Heat transfer coefficient through gap W/m^2/K
        self.k_clad = k_clad # thermal conductivity coefficient in clad W/m/K
        self.I_f = I_f # number of mesh elements in the fuel
        self.I_c = I_c # number of mesh elements in clad

        self.plot_results = plot_at_z

        print(f"$$$---------- THM: prototype, case treated : {self.name}.")
        print("Warning : only single phase flow treated in this implementation of heat convection in coolant canal.")

        # Prepare and solve 1D heat convection along the z direction in the canal.
        print("$$---------- Calling FVM_ConvectionInCanal_MONO class.")
        print(f"Setting up heat convection solution along the axial dimension. zmax = {self.Lf} with {self.I_z} axial elements.")
        self.convection_sol = FVM_Canal_MONO(self.Lf, self.T_in, self.Q_flow, self.P_cool, self.I_z, self.canal_type, 
                                        self.r_f, self.clad_r, self.r_w)
        self.convection_sol.set_Fission_Power(self.Q_fiss_amp, self.Q_fiss_variation_type)
        print(f"Fission power source initialized with amplitude = {self.Q_fiss_amp} and axial variation profile = {self.Q_fiss_variation_type}.")
        self.setADI_CL_convection() # setting up the linear system corresponding to the finite volume discretization + boundary condition at z=0
        print("$---------- Solving for h(z) using the Finite Volumes Method.")
        self.convection_sol.h_z = self.convection_sol.solve_h_in_canal() # solving the enthalpy evolution in the canal
        print("$---------- Solving for T_surf(z) using the Dittus-Boelter correlation. Water Properties evaluated by IAPWS97")
        self.Tsurf = self.convection_sol.compute_T_surf() # conputing and retrieving the clad surface temperatures obtained through the Dittus-Boelter correlation
    
        # Prepare and solve 1D radial heat conduction in the fuel rod, given a Clad surface temperature as a bondary condition 
        self.SetupAndSolve_Conduction_at_all_z() # creates a list of 
        
        # extend to Twater here.
         
        if self.plot_results:
            for z_val in self.plot_results:
                if z_val in self.convection_sol.z_mesh:
                    plane_indx = np.where(self.convection_sol.z_mesh==z_val)[0][0]
                    self.plot_Temperature_at_z(z,plane_indx,self.T_distrib[plane_indx])
                else: # Interpolate between nearest z values to obtain Temperature distribution at a given z.
                    second_plane_indx = np.where(self.convection_sol.z_mesh>z_val)[0][0]
                    first_plane_indx = second_plane_indx-1
                    interp_T = (self.T_distrib[first_plane_indx]+self.T_distrib[second_plane_indx])/2
                    self.plot_Temperature_at_z(z,plane_indx,interp_T)


    def SetupAndSolve_Conduction_at_all_z(self):
        self.T_distrib = []
        for axial_plane_nb in range(self.convection_sol.N_vol):
            z = self.convection_sol.z_mesh[axial_plane_nb]
            T_surf = self.convection_sol.T_surf[axial_plane_nb]
            Qfiss = self.convection_sol.get_Fission_Power()[axial_plane_nb]
            self.T_distrib.append(self.run_Conduction_In_Fuel_at_z(z,Qfiss,T_surf))


        return
    
    def setADI_CL_convection(self):
        """
        function used to set up the linear system corresponding to the finite volume discretization + boundary condition at z=0 
        """
        for i in range(1,self.convection_sol.N_vol-1):
            self.convection_sol.set_ADi_conv(i,
                                    ci=-1,
                                    ai=1,
                                    bi=0,
                                    di = self.convection_sol.q_fluid[i]*self.convection_sol.dz/(self.convection_sol.Q_flow*self.convection_sol.A_canal))
        A0,Am1 = np.zeros(self.convection_sol.N_vol), np.zeros(self.convection_sol.N_vol)
        A0[0] = 1
        D0 = self.convection_sol.h_z0 + self.convection_sol.q_fluid[i]*self.convection_sol.dz/(2*self.convection_sol.Q_flow*self.convection_sol.A_canal)
        Am1[-2:]=[-1, 1]
        Dm1 = self.convection_sol.q_fluid[i]*self.convection_sol.dz/(self.convection_sol.Q_flow*self.convection_sol.A_canal)
        self.convection_sol.set_CL_conv(A0,Am1,D0,Dm1)
        return

    def run_Conduction_In_Fuel_at_z(self,z,Qfiss_z,T_surf_z):
        print(f"$$---------- Setting up FDM_HeatConductionInFuelPin class for z = {z} m, Qfiss(z) = {Qfiss_z} W/m^3 and T_surf(z) = {T_surf_z} K")
        heat_conduction = FDM_Fuel(self.r_f, self.I_f, self.gap_r, self.clad_r, self.I_c, Qfiss_z, self.k_fuel, self.k_clad, self.H_gap, z, T_surf_z)

        for i in range(1,heat_conduction.N_node-1):

            if i<heat_conduction.I_f-1: # setting Aij and Di values for nodes inside the fuel 
                heat_conduction.set_ADi_cond(i, 
                                    -heat_conduction.get_Di_half(i-1), 
                                    heat_conduction.get_Di_half(i-1)+heat_conduction.get_Di_half(i), 
                                    -heat_conduction.get_Di_half(i), 
                                    heat_conduction.deltaA_f*heat_conduction.Qfiss)
            elif i==heat_conduction.I_f-1: # setting Aij and Di values for last fuel element
                heat_conduction.set_ADi_cond(i,
                                -heat_conduction.get_Di_half(i-1),
                                (heat_conduction.get_Di_half(i-1)+heat_conduction.get_Ei_gap()),
                                -heat_conduction.get_Ei_gap(),
                                heat_conduction.deltaA_f*heat_conduction.Qfiss)
            elif i==heat_conduction.I_f: # setting Aij and Di values first fuel / gap interface
                heat_conduction.set_ADi_cond(i, 
                                    -heat_conduction.get_Ei_gap(), 
                                    heat_conduction.get_Ei_gap()+heat_conduction.get_Gi(), 
                                    -heat_conduction.get_Gi(), 
                                    0)
            elif i==heat_conduction.I_f+1: # setting Aij and Di values second gap / clad interface
                heat_conduction.set_ADi_cond(i, 
                                    -heat_conduction.get_Gi(), 
                                    heat_conduction.get_Fi_gap()+heat_conduction.get_Gi(), 
                                    -heat_conduction.get_Fi_gap(), 
                                    0)
            elif i>heat_conduction.I_f+1 : # setting Aij and Di for all elements in the clad, apart from the last one
                heat_conduction.set_ADi_cond(i, 
                                    -heat_conduction.get_Di_half(i-1), 
                                    heat_conduction.get_Di_half(i-1)+heat_conduction.get_Di_half(i), 
                                    -heat_conduction.get_Di_half(i), 
                                    0)
        A0,Am1 = np.zeros(heat_conduction.N_node), np.zeros(heat_conduction.N_node) 
        A0[:2] = [heat_conduction.get_Di_half(0), -heat_conduction.get_Di_half(0)]
        Am1[-2:] = [-heat_conduction.get_Di_half(heat_conduction.N_node-2), heat_conduction.get_Di_half(heat_conduction.N_node-2)+heat_conduction.get_Ei_clad()]
        D0 = heat_conduction.deltaA_f*heat_conduction.Qfiss
        Dm1 = heat_conduction.get_Ei_clad()*heat_conduction.T_surf
        heat_conduction.set_CL_cond(A0, Am1, D0, Dm1)
        print(f"$---------- Solving for T(r) using the Finite Difference Method, at z = {z}.")
        heat_conduction.solve_T_in_pin()
        return heat_conduction




   
    def plot_Temperature_at_z(self, z, plane_index, Temperature_at_z):
        
        fig,ax = plt.subplots(dpi=200)
        print(f"Calculation radial mesh is = {np.sqrt(2*Temperature_at_z.A_calculation_mesh)}")
        print(f"plot mesh is {self.Temperature_at_z.plot_mesh}")
        ax.scatter(Temperature_at_z.plot_mesh, Temperature_at_z.T_distrib, marker = "x", s=5, label="Radial temperature distribution in Fuel rod.")
        #ax.scatter(np.sqrt(2*self.Temperature_at_z.A_calculation_mesh), self.Temperature_at_z.T_distrib, marker = "x", s=5, label="Radial temperature distribution in Fuel rod.")
        ax.legend(loc = "best")
        ax.grid()
        ax.set_xlabel(f"Radial position in {Temperature_at_z.plotting_units}")
        ax.set_ylabel(f"Temperature in K")
        ax.set_title(f"Temperature distribution in fuel rod at z = {z}, case 1")
        fig.savefig(f"{self.name}_Figure_plane{plane_index}")
        colors = ["red", "yellow", "green", "blue"]
        Temperature_at_z.extend_to_canal_visu(rw = self.convection_sol.wall_dist, Tw = self.convection_sol.T_water[plane_index])
        print(f"T_surf = {self.convection_sol.T_surf[plane_index]} K and T_water = {self.convection_sol.T_water[plane_index]} K")
        print(f"radii at bounds {Temperature_at_z.radii_at_bounds}")
        fig_filled, axs = plt.subplots(dpi=200)
        for i in range(len(Temperature_at_z.physical_regions_bounds)-1):
            axs.fill_between(x=Temperature_at_z.radii_at_bounds, y1=Temperature_at_z.T_distrib[1:], y2=400*np.ones(len(Temperature_at_z.radii_at_bounds)),where=(Temperature_at_z.radii_at_bounds>=self.temp_distrib_to_plot.physical_regions_bounds[i])&(self.temp_distrib_to_plot.radii_at_bounds<=self.temp_distrib_to_plot.physical_regions_bounds[i+1]), color = colors[i])
        
        #axs.fill_betweenx(y=self.temp_distrib_to_plot.T_distrib[1:], x1=self.temp_distrib_to_plot.radii_at_bounds, where=(self.temp_distrib_to_plot.radii_at_bounds<self.temp_distrib_to_plot.r_f), facecolor='red')
        #axs.fill_betweenx(self.convection_sol.T_surf[plane_index])
        ax.legend(loc = "best")
        ax.grid()
        ax.set_xlabel(f"Radial position in {Temperature_at_z.plotting_units}")
        ax.set_ylabel(f"Temperature in K")
        ax.set_title(f"Temperature distribution in fuel rod at z = {z}, case 1")
        fig_filled.savefig(f"{self.name}_Figure_plane{plane_index}_colors")

