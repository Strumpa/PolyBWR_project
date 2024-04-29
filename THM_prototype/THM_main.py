# Used to create instances of THM class and solve convection + conduction in fuel rod.
# author : R. Guasch
# Purpose : prototyping for further developments in THM module of Donjon5

from THM_MONO import FDM_HeatConductionInFuelPin as FDM_Fuel
from THM_MONO import FVM_ConvectionInCanal_MONO as FVM_Canal_MONO
import numpy as np
from iapws import IAPWS97
import matplotlib.pyplot as plt

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
        print(f"Setting up heat convection solution along the axial dimension. zmax = {self.Lf} m with {self.I_z} axial elements.")
        self.convection_sol = FVM_Canal_MONO(self.Lf, self.T_in, self.Q_flow, self.P_cool, self.I_z, self.canal_type, 
                                        self.r_f, self.clad_r, self.r_w)
        self.convection_sol.set_Fission_Power(self.Q_fiss_amp, self.Q_fiss_variation_type)
        print(f"Fission power source initialized with amplitude = {self.Q_fiss_amp} and axial variation profile = {self.Q_fiss_variation_type}.")
        self.setADI_CL_convection() # setting up the linear system corresponding to the finite volume discretization + boundary condition at z=0
        print("$---------- Solving for h(z) using the Finite Volumes Method.")
        self.convection_sol.h_z = self.convection_sol.solve_h_in_canal() # solving the enthalpy evolution in the canal
        print("$---------- Solving for T_surf(z) using the Dittus-Boelter correlation. Water Properties evaluated by IAPWS97")
        self.Tsurf = self.convection_sol.compute_T_surf() # computing and retrieving the clad surface temperatures obtained through the Dittus-Boelter correlation
    
        # Prepare and solve 1D radial heat conduction in the fuel rod, given a Clad surface temperature as a bondary condition 
        self.SetupAndSolve_Conduction_at_all_z() # creates a list of Temperature distributions in the fuel rod given a surface temperature computed by solving the conection problem
        self.get_TFuel_rowlands() # compute and store in the T_eff_fuel attribute the effective fuel temperature given by the Rowlands formula
        self.get_Tfuel_surface() # store in the T_fuel_surface attribute the fuel surface temperature computed

        # extend to Twater : adding a mesh point corresponding to the middle of the canal in the plotting array, add rw to the bounds array and add Twater to the results array
        for index_z in range(len(self.convection_sol.z_mesh)):
            self.T_distributions_axial[index_z].extend_to_canal_visu(rw = self.convection_sol.wall_dist, Tw = self.convection_sol.T_water[index_z])
         
        if self.plot_results:
            for z_val in self.plot_results:
                self.plot_Temperature_at_z(z_val)


    def SetupAndSolve_Conduction_at_all_z(self):
        self.T_distributions_axial = []
        for axial_plane_nb in range(self.convection_sol.N_vol):
            z = self.convection_sol.z_mesh[axial_plane_nb]
            T_surf = self.convection_sol.T_surf[axial_plane_nb]
            Qfiss = self.convection_sol.get_Fission_Power()[axial_plane_nb]
            self.T_distributions_axial.append(self.run_Conduction_In_Fuel_at_z(z,Qfiss,T_surf))


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
            print(self.convection_sol.q_fluid[i])
            print(self.convection_sol.dz)
            print(self.convection_sol.Q_flow)
            print(self.convection_sol.A_canal)
        A0,Am1 = np.zeros(self.convection_sol.N_vol), np.zeros(self.convection_sol.N_vol)
        A0[0] = 1
        D0 = self.convection_sol.h_z0 + self.convection_sol.q_fluid[0]*self.convection_sol.dz/(2*self.convection_sol.Q_flow*self.convection_sol.A_canal)
        Am1[-2:]=[-1, 1]
        Dm1 = self.convection_sol.q_fluid[-1]*self.convection_sol.dz/(self.convection_sol.Q_flow*self.convection_sol.A_canal)
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
                                (heat_conduction.get_Di_half(i-1)+heat_conduction.get_Ei_fuel()),
                                -heat_conduction.get_Ei_fuel(),
                                heat_conduction.deltaA_f*heat_conduction.Qfiss)
            elif i==heat_conduction.I_f: # setting Aij and Di values first fuel / gap interface
                heat_conduction.set_ADi_cond(i, 
                                    -heat_conduction.get_Ei_fuel(), 
                                    heat_conduction.get_Ei_fuel()+heat_conduction.get_Gi(), 
                                    -heat_conduction.get_Gi(), 
                                    0)
            elif i == heat_conduction.I_f+1: # setting Aij and Di values second gap / clad interface
                heat_conduction.set_ADi_cond(i, 
                                    -heat_conduction.get_Gi(), 
                                    heat_conduction.get_Fi_gap()+heat_conduction.get_Gi(), 
                                    -heat_conduction.get_Fi_gap(), 
                                    0)
            elif i == heat_conduction.I_f+2: # Treating the first clad element interface with the gap.
                heat_conduction.set_ADi_cond(i,
                                             -heat_conduction.get_Ei_gap(),
                                             (heat_conduction.get_Di_half(i)+heat_conduction.get_Ei_gap()),
                                             -heat_conduction.get_Di_half(i),
                                             0)
            elif i>heat_conduction.I_f+2 : # setting Aij and Di for all elements in the clad, apart from the last one
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
    
    def get_TFuel_rowlands(self):
        self.T_eff_fuel = np.zeros(self.convection_sol.N_vol)
        for i in range(len(self.T_distributions_axial)):
            self.T_distributions_axial[i].compute_T_eff()
            T_eff_z = self.T_distributions_axial[i].T_eff
            self.T_eff_fuel[i] = T_eff_z
        return
    def get_Tfuel_surface(self):
        self.T_fuel_surface = np.zeros(self.convection_sol.N_vol)
        for i in range(len(self.T_distributions_axial)):
            if len(self.T_distributions_axial[i].T_distrib) == self.T_distributions_axial[i].N_node:
                T_surf_fuel_z = self.T_distributions_axial[i].T_distrib[self.I_f]
            else:
                T_surf_fuel_z = self.T_distributions_axial[i].T_distrib[self.I_f+1]
            self.T_fuel_surface[i] = T_surf_fuel_z
        return




    def plot_Temperature_at_z(self, z_val):
        print(f"$$---------- Plotting Temperature distribution in rod + canal z = {z_val} m")

        print(f"z_val is {z_val}")
        print(f"z_mesh is {self.convection_sol.z_mesh}")
        if z_val in self.convection_sol.z_mesh:
            plane_index = int(np.where(self.convection_sol.z_mesh==z_val)[0][0])
            Temperature_distrib_to_plot = self.T_distributions_axial[plane_index].T_distrib
            plotting_mesh = self.T_distributions_axial[plane_index].plot_mesh
            radii_at_bounds = self.T_distributions_axial[plane_index].radii_at_bounds
            physical_regions_bounds = self.T_distributions_axial[plane_index].physical_regions_bounds
            plotting_units = self.T_distributions_axial[plane_index].plotting_units
            Tsurf = self.convection_sol.T_surf[plane_index]
            Twater = self.convection_sol.T_water[plane_index]
            Tcenter = self.T_distributions_axial[plane_index].T_center
        else: # Interpolate between nearest z values to obtain Temperature distribution at a given z.
            second_plane_index = np.where(self.convection_sol.z_mesh>z_val)[0][0]
            first_plane_index = second_plane_index-1
            plane_index = (first_plane_index+second_plane_index)/2
            print(f"plane index used is {plane_index}")
            plotting_mesh = self.T_distributions_axial[first_plane_index].plot_mesh
            radii_at_bounds = self.T_distributions_axial[first_plane_index].radii_at_bounds
            physical_regions_bounds = self.T_distributions_axial[first_plane_index].physical_regions_bounds
            Temperature_distrib_to_plot = self.T_distributions_axial[first_plane_index].T_distrib+(z_val-self.convection_sol.z_mesh[first_plane_index])*(self.T_distributions_axial[second_plane_index].T_distrib-self.T_distributions_axial[first_plane_index].T_distrib)/(self.convection_sol.z_mesh[second_plane_index]-self.convection_sol.z_mesh[first_plane_index])
            plotting_units = self.T_distributions_axial[first_plane_index].plotting_units
            Tsurf = self.convection_sol.T_surf[first_plane_index] + (z_val-self.convection_sol.z_mesh[first_plane_index])*(self.convection_sol.T_surf[second_plane_index]-self.convection_sol.T_surf[first_plane_index])/(self.convection_sol.z_mesh[second_plane_index]-self.convection_sol.z_mesh[first_plane_index])
            Twater = self.convection_sol.T_water[first_plane_index] + (z_val-self.convection_sol.z_mesh[first_plane_index])*(self.convection_sol.T_water[second_plane_index]-self.convection_sol.T_water[first_plane_index])/(self.convection_sol.z_mesh[second_plane_index]-self.convection_sol.z_mesh[first_plane_index])
            Tcenter = self.T_distributions_axial[first_plane_index].T_center + (z_val-self.convection_sol.z_mesh[first_plane_index])*(self.T_distributions_axial[second_plane_index].T_center-self.T_distributions_axial[first_plane_index].T_center)/(self.convection_sol.z_mesh[second_plane_index]-self.convection_sol.z_mesh[first_plane_index])

        
        
        if (isinstance(plane_index, int)):
            plane_index_print = plane_index
        else:
            plane_index_print = str(plane_index).split(".")[0]+str(plane_index).split(".")[1]
        print(f"at z = {z_val}, temp distrib is = {Temperature_distrib_to_plot}")
        colors = ["lime", "bisque", "chocolate", "royalblue"]
        labels = ["Fuel", "Gap", "Clad", "Water"]
        fig_filled, axs = plt.subplots(dpi=200)
        for i in range(len(physical_regions_bounds)-1):
            axs.fill_between(x=radii_at_bounds, y1=(Tcenter+50)*np.ones(len(radii_at_bounds)), y2=(Twater-50)*np.ones(len(radii_at_bounds)),where=(radii_at_bounds>=physical_regions_bounds[i])&(radii_at_bounds<=physical_regions_bounds[i+1]), color = colors[i], label = labels[i])
        axs.scatter(plotting_mesh, Temperature_distrib_to_plot, marker = "D", color="black",s=10, label="Radial temperature distribution in Fuel rod.")
        axs.legend(loc = "best")
        axs.grid()
        axs.set_xlabel(f"Radial position in {plotting_units}")
        axs.set_ylabel(f"Temperature in K")
        axs.set_title(f"Temperature distribution in fuel rod at z = {z_val}, {self.name}")
        fig_filled.savefig(f"{self.name}_Figure_plane{plane_index_print}_colors")


#Case 1 : base parameters
# Parameters used to create object from FDM_HeatConductioninFuelpin class
Qfiss1 = 0.3e9 # W/m^3
fuel_radius1 = 5.6e-3 # m
gap_width1 = 0.54e-3 # m
clad_width1 = 0.38e-3 # m
k_fuel1 = 5 # W/m/K
H_gap1 = 10000 # W/m^2/K
k_clad1 = 10 # W/m/K
I_f1 = 8
I_c1 = 3

# Paramters used to create object from FVM_ConvectionInCanal class
canal_type1 = "cylindrical" #"square" #"cylindrical"
canal_width1 = 0.5e-3 # m
Lf1 = 2 # m
T_in1 = 500 # K
Q_flow1 = 7000 # kg/m^2/s
P_cool1 = 10.8 #MPa
I_z1 = 10


rw1=fuel_radius1+gap_width1+clad_width1+canal_width1 # canal radius
gap_rad1 = fuel_radius1+gap_width1
clad_rad1 = gap_rad1+clad_width1
plot_at_z1 = [0.7,0.8,0.9]
case1 = Version5_THM_prototype("Case1_ENE6107A_project", rw1, canal_type1, Lf1, T_in1, P_cool1, Q_flow1, I_z1, Qfiss1, "constant", 
                               fuel_radius1, gap_rad1, clad_rad1, k_fuel1, H_gap1, k_clad1, I_f1, I_c1, plot_at_z1)


print(f"case 1 h_z is {case1.convection_sol.h_z} J/kg")
print(f"case 1 T_water is {case1.convection_sol.T_water} K")
print(f"case 1 Hc is {0.5*(case1.convection_sol.Hc[3]+case1.convection_sol.Hc[4])} W/m^2/K")
print(f"q_fluid1 = {case1.convection_sol.q_fluid}")

print(f"case 1 A_canal = {case1.convection_sol.A_canal} m^2")
print(f"case 1 T_surf is {case1.convection_sol.T_surf} K")


print(f"case 1 T_eff in fuel is {case1.T_eff_fuel} K")
print(f"case 1 T_surf fuel is {case1.T_fuel_surface} K")
state = IAPWS97(P=10.8,T=501.44)
print(state.h)

#Case 2 : base parameters
# Parameters used to create object from FDM_HeatConductioninFuelpin class
Qfiss2 = 0.3e9 # W/m^3
fuel_radius2 = 5.4e-3 # m
gap_width2 = 0.6e-3 # m
clad_width2 = 0.4e-3 # m
k_fuel2 = 5 # W/m/K
H_gap2 = 10000 # W/m^2/K
k_clad2 = 15 # W/m/K
I_f2 = 80
I_c2 = 3

# Paramters used to create object from FVM_ConvectionInCanal class
canal_type2 = "cylindrical"
canal_width2 = 0.35e-3 # m
Lf2 = 2.5 # m
T_in2 = 450 # K
Q_flow2 = 7200 # kg/m^2/s
P_cool2 = 10.8 #MPa
I_z2 = 100


rw2=fuel_radius2+gap_width2+clad_width2+canal_width2 # canal radius
gap_rad2 = fuel_radius2+gap_width2
clad_rad2 = gap_rad2+clad_width2
plot_at_z2 = [0.7,0.8,0.9]
case2 = Version5_THM_prototype("Case2_ENE6107A_project", rw2, "cylindrical", Lf2, T_in2, P_cool2, Q_flow2, I_z2, Qfiss2, "sine", 
                               fuel_radius2, gap_rad2, clad_rad2, k_fuel2, H_gap2, k_clad2, I_f2, I_c2, plot_at_z2)

print(f"case 2 h_z is {case2.convection_sol.h_z} J/kg")
print(f"case 2 T_water is {(case2.convection_sol.T_water[35]+case2.convection_sol.T_water[36])/2} K")
print(f"case 2 Hc is {0.5*(case2.convection_sol.Hc[35]+case2.convection_sol.Hc[36])} W/m^2/K")
print(f"q_fluid2 = {case2.convection_sol.q_fluid}")

print(f"case 2 A_canal = {case2.convection_sol.A_canal} m^2")
print(f"case 2 T_surf is {(case2.convection_sol.T_surf[35]+case2.convection_sol.T_surf[36])/2} K")


print(f"case 2 T_eff in fuel is {case2.T_eff_fuel} K")
print(f"case 2 T_surf fuel is {case2.T_fuel_surface} K")
plt.scatter(case2.convection_sol.z_mesh, case2.convection_sol.h_z)
plt.savefig("test_enthalpy")
