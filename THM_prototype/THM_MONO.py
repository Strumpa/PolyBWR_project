# Python3 class part of THM_prototype
# uses : - setting up mesh centered finite difference solutino to heat conduction in fuel rod, FDM_ConductionInFuelPin
#        - setting up finite volume discretization for solving heat convection in coolant, FVM_ConvectionInCanal
# Author : R. Guasch
# technical documentation : "Revisiting the simplified thermo-hydraulics module THM: in DONJON5 code" - A. Hébert, March 2018
# document available at http://merlin.polymtl.ca/downloads/thm.pdf

import numpy as np
from iapws import IAPWS97


class FDM_HeatConductionInFuelPin:
    def __init__(self, r_fuel, I_f, gap_width, clad_width, I_c, Qfiss, kf, kc, Hgap, z, T_surf):
        # Physical prameters
        self.r_f = r_fuel # fuel pin radius in meters
        self.I_f = I_f # number of mesh elements in the fuel
        self.gap_w = gap_width # gap width in meters, used to determine mesh elements for constant surface discretization
        self.clad_w = clad_width # clad width in meters, used to determine mesh elements for constant surface discretization
        self.I_c = I_c # number of mesh elements in clad
        self.Qfiss = Qfiss # Fission power density in W/m^3
        self.Hg = Hgap # Heat transfer coefficient through gap W/m^2/K
        self.z = z # corresponding height in m corresponding to axial discretization used in FVM_ConvectionInCanal class
        self.T_surf = T_surf # Boundary condition outer clad surface temperature computed from FVM_ConvectionInCanal class 
        
        # compute relevant quantities to initialise object
        self.N_node = I_f + I_c +2
        self.A = np.eye(self.N_node)
        self.D = np.zeros(self.N_node)
        self.compute_radii()
        self.compute_Area_meshes()
        self.radii_at_centers = np.sqrt(self.A_mesh_centers/np.pi)
        self.radii_at_bounds = np.sqrt(self.A_mesh_bounds/np.pi)
        self.kf = kf
        self.kc = kc
        # this array is probaby not needed here as one might assume that kf and kc
        # are constant in fuel/clad. I wanted to keep and option to let them vary according to temperature at node as conductive properties might differ when high temperature gradients are present.
        self.k = np.zeros(self.N_node-1) # associate a k to each mesh center element
        for i in range(len(self.k)):
            if i <self.I_f:
                self.k[i]=kf
            elif i == self.I_f:
                self.k[i]=0 # in gap!
            elif i<self.N_node:
                self.k[i]=kc
        return
    def compute_radii(self):
        self.gap_r = self.r_f+self.gap_w
        self.clad_r = self.gap_r + self.clad_w
        return
    def compute_Area_meshes(self):
        """
        building necessary meshes for solving the heat conduction equation on the constant area discretization 
        """
        self.A_mesh_bounds = []
        self.A_mesh_centers = []
        self.A_calculation_mesh = []
        A_f = np.pi*self.r_f**2
        A_gf = np.pi*self.gap_r**2
        A_cgf = np.pi*self.clad_r**2
        self.deltaA_f = A_f / self.I_f # base assumption is that delta_A is constant in each region --> delta A fuel = constant in fuel, delta A clad = constant in clad and 1 delta A gap.
        for i in range(self.I_f+1):
            self.A_mesh_bounds.append(i*self.deltaA_f)
        for i in range(self.I_f):
            self.A_mesh_centers.append(i*self.deltaA_f+self.deltaA_f/2)
            self.A_calculation_mesh.append(i*self.deltaA_f+self.deltaA_f/2)
    
        self.deltaA_g = A_gf-A_f
        self.A_mesh_bounds.append(self.A_mesh_bounds[-1]+self.deltaA_g)
        self.A_mesh_centers.append(self.A_mesh_centers[-1]+self.deltaA_f/2+self.deltaA_g/2) # last center in fuel + half of the fuel area step to get to the last fuel bound + half of the gap area step to get to the center of the gap
        self.A_calculation_mesh.append(self.A_mesh_bounds[-1])
        self.A_calculation_mesh.append(self.A_mesh_bounds[-1]+self.deltaA_g)
        self.deltaA_c = (A_cgf-A_gf)/self.I_c
        for i in range(self.I_c):
            self.A_mesh_bounds.append(self.A_mesh_bounds[-1]+self.deltaA_c)
        for i in range(self.I_c):
            self.A_mesh_centers.append(self.A_mesh_centers[-1]+self.deltaA_c)
            self.A_calculation_mesh.append(self.A_mesh_centers[-1]+self.deltaA_c)
        self.A_mesh_centers = np.array(self.A_mesh_centers)
        self.A_mesh_bounds = np.array(self.A_mesh_bounds)
        self.A_calculation_mesh = np.array(self.A_calculation_mesh)
        return
    
    def get_Di_half(self,i):
        print(f"i in get_Di_half is = {i}")
        print("getting Di for interior fuel points")
        print(f"A_mesh_bounds is = {self.A_mesh_bounds[i+1]}")
        print(f"delta Af = {self.deltaA_f}")
        print(f"k_i = {self.k}")
        print(self.I_f)
        if i > self.I_f+1:
            i=i-1
            Di_half = 4*self.A_mesh_bounds[i+1]/((self.deltaA_c/self.k[i])+(self.deltaA_c/self.k[i+1]))
        else:
            Di_half = 4*self.A_mesh_bounds[i+1]/((self.deltaA_f/self.k[i])+(self.deltaA_f/self.k[i+1]))
        print(f"D_half is = {Di_half}")
        return Di_half
    
    def get_Ei_gap(self):
        print(f"I in get_Ei_half is = {self.I_f}")
        Ei_half = 4*self.A_mesh_bounds[self.I_f]*self.k[self.I_f-1]/self.deltaA_f
        print(f"Area used in Ei_gap {self.A_mesh_bounds[self.I_f]}")
        print(f"Theoretical area: {np.pi*self.r_f**2}")
        return Ei_half
    
    def get_Ei_clad(self):
        Ei_half = 4*self.A_mesh_bounds[-1]*self.k[-1]/self.deltaA_c
        print(f"Area used in Ei_clad {self.A_mesh_bounds[-1]}")
        print(f"Theoretical area: {np.pi*self.clad_r**2}")
        return Ei_half
    
    def get_Fi_gap(self):
        print(f"k used for FI3/2 {self.k[self.I_f+1]}")
        print(f"Area used in Fi gap is {self.A_mesh_bounds[self.I_f+1]}")
        print(f"Theoretical area for Fi_gap is {np.pi*(self.gap_r)**2}")
        Fi_half = 4*self.A_mesh_bounds[self.I_f+1]*self.k[self.I_f+1]/self.deltaA_c
        print(f"Fi_half is {Fi_half}")
        return Fi_half
    
    def get_Gi(self):
        return self.Hg*self.radii_at_centers[self.I_f]
    
    def set_ADi(self, i, ci, ai, bi, di):
        # create lines for the tri-diagonal entries in
        self.A[i, i-1:i+2] = [ci, ai, bi]
        self.D[i] = di
        return
    
    def set_CL(self, A0, Am1, D0, Dm1):
        # I_and_half is the index for the I+1/2 element of th mesh which corresponds to the last point in the fuel.
        # conditions aux limites
        # A0 = A[0], Am1 = A[-1], A moins 1, 
        # D0 = D[0], Dm1 = D[-1], D moins 1.
        self.A[0], self.A[-1] = A0, Am1
        self.D[0], self.D[-1] = D0, Dm1
        return
    
    def solve_T_in_pin(self):
        self.T_distrib = np.linalg.solve(self.A, self.D)
     
    
class FVM_ConvectionInCanal:
    def __init__(self, Lf, T_in, Q_flow, P_cool, I_z, canal_type, rf, rc, rw):
        """
        Lf = fuel rod length in m
        T_in = inlet water temperature K
        Q_flow = mass flux in kg/m^2/s, assumed to be constant along the axial profile.
        P_cal = coolant pressure in MPa, assumed to be constant along the axial profile.
        I_z = number of mesh elements on axial mesh
        canal_type = cylindrical or square, used to determine the cross sectional flow area in the canal and the hydraulic diameter
        rf = fuel rod radius
        rc, rw = outer clad radius (m), outer canal radius (m) if type is cylindrical, if type = squaer rw is the radius of inscribed circle in the square canal, ie half the square's side.

        Important note : cross sectional flow area is assumed to constant along z axis. Would be interesting to expand this to treat variations in flow area as in :
        THE MODELING OF ADVANCED BWR FUEL DESIGNS WITH THE NRC FUEL DEPLETION CODES PARCS/PATHS - A. WYSOCKI et al. August 2014
        """
        # Physical parameters
        self.Lf = Lf
        self.T_in = T_in
        self.Q_flow = Q_flow
        self.P_cool = P_cool
        self.fuel_radius = rf
        self.clad_radius = rc
        self.wall_dist = rw
        initial_water_state_z0 = IAPWS97(P=P_cool,T=T_in)
        self.h_z0 = initial_water_state_z0.h*10**3 # returs enthalpy at z0 for given (P,T), this assumes 1 phase liquid water at z=0, converted to J/kg
        self.canal_type = canal_type
        # Calculating mesh parameters.
        self.N_vol = I_z
        self.dz = self.Lf/self.N_vol
        self.A, self.D = np.eye(self.N_vol), np.zeros(self.N_vol)
        self.z_mesh = np.linspace(0.5*self.dz, self.Lf-0.5*self.dz, self.N_vol) # creating z mesh from volume center to volume center. 
        
        # calculation A_canal and DH depending on the canal type : idea would be to generalize to square section "CARCEL" to use results in DONJON5.
        if self.canal_type == "cylindrical":
            self.A_canal = np.pi*self.wall_dist**2 - np.pi*self.clad_radius**2
            self.P_wetted = 2*np.pi*(self.wall_dist + self.clad_radius)
            self.DH = 4*self.A_canal / self.P_wetted # DH = 4*A/P where A = wetted cross sectional area and P is wetted perimeter = perimeter of fuel
        elif self.canal_type == "square":
            self.A_canal = (2*self.wall_dist)**2-np.pi*self.clad_radius**2
            self.DH = 4*self.A_canal / 2*(4*self.wall_dist+np.pi*self.clad_radius)


        return
    
    def set_ADi(self, i, ci, ai, bi, di):
        self.A[i, i-1:i+2] = [ci, ai, bi]
        self.D[i] = di
        return
    
    def set_CL(self, A0, Am1, D0, Dm1):
        self.A[0], self.A[-1] = A0, Am1
        self.D[0], self.D[-1] = D0, Dm1
        return
    def set_Fission_Power(self, amplitude, variation_type):
        """
        option to set fission source axial profile : 
        amblitude : fission power from fuel per unit volume W/m^3
        variation_type : string used to allow for constant, sinusoial or cosinusoidal axial profiles ("constant", "sine", "cosine" keywords allowed)
        """
        self.Power_profile = np.ones(self.N_vol)
        if variation_type == "constant":
            self.Power_profile = amplitude*self.Power_profile
        elif variation_type == "sine":
            for i in range(self.N_vol):
                self.Power_profile[i] = amplitude*np.sin(np.pi*self.z_mesh[i]/self.Lf)
        elif variation_type == "cosine":
            print("Keyword for cosine axial variation of fuel power not implemented yet")
        
        return
    
    def get_Fission_Power(self):
        """
        function to retrieve source term from axial profile used to model fission power distribution in the fuel rod
        """
        return self.Power_profile
    
    def solve_h_in_canal(self):
        self.h_z = np.linalg.solve(self.A, self.D)
        return self.h_z
    

        
    def compute_T_surf(self):
        N_MAX = 100000
        self.T_surf = np.zeros(self.N_vol)
        self.Hc = np.zeros(self.N_vol)
        self.T_water = np.zeros(self.N_vol)
        self.T_water[0] = self.T_in
        for i in range(self.N_vol):
            self.Hc[i] = (0.023)*(IAPWS97(P=self.P_cool, h=self.h_z[i]*10**-3).Liquid.Prandt)**0.4*(self.Q_flow*self.DH/(IAPWS97(P=self.P_cool, h=self.h_z[i]*10**-3).Liquid.mu))**0.8*IAPWS97(P=self.P_cool, h=self.h_z[i]*10**-3).Liquid.k/self.DH
            self.T_water[i] = IAPWS97(P=self.P_cool, h=self.h_z[i]*10**-3).T
            self.T_surf[i] = (np.pi*self.fuel_radius**2*self.get_Fission_Power()[i]/(2*np.pi*self.clad_radius)+self.Hc[i]*self.T_water[i])/self.Hc[i]
            #for iteration in range(N_MAX):
                

        return self.T_surf
    
    def set_transitoire(self, t_tot, Tini, dt):
        self.t_tot, self.dt = t_tot, dt           
        self.N_temps = round(self.t_tot / self.dt) # pas de temps (timesteps), il faut etre un nombre entier
        self.T = np.zeros((self.N_temps+1, self.N_vol)) # tableau 2D de temperature. 
        self.T[0] = Tini # Tini est une liste
        return 

    