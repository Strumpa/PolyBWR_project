# GeometricTools.py Python3 script
# Author : R. Guasch
# Purpose : Collection of tools for geometric processes

import numpy as np
import math


# NXT: formalism helper functions
def region_index_cartesian(i,j,k,nx,ny):
    return i+nx*(j-1+ny*(k-1))
def region_index_tube(r,meshr):
    print(meshr)
    for i in range(len(meshr)):
        if meshr[i]<=r and r<meshr[i+1]:
            return i+1
def region_index_CARCEL2D(g,i,j,nr,nx):
    return g+(nr+1)*(i-1+nx*(j-1))

# Angular cyclic tracking treatment :
def EQW2_angles_calculator(nangle, X, Y):
    """
    Calculate phi_k angles for  EQW2 TSPC options
    nangle : number of angles in (0,pi/2),
    for 2D specular applications, allowed values are : 8, 12, 14, 18, 20, 24, 30 = p+1 where p is prime
    X and Y : lengths of cartesian domain of size X x Y
    """
    p=nangle-1
    phi_k = []
    for k in range(p+1):
        phi_k.append(np.arctan(k*Y/((p+2-k)*X)))
    return phi_k

def phi_a(mx, my, Xp, Yp):
    """
    This function calculates the angular set for cyclic trackings
    mx,my integers, prime to one another
    Xp and Yp dimensions of the cartesian cell
    """
    return np.arctan(my*Yp/(mx*Xp))

def u_a(mx, my, Xp, Yp):
    """
    This function returns the unit vector associated with direction a(lpha)
    mx,my integers, prime to one another
    Xp and Yp dimensions of the cartesian cell
    """
    d = np.sqrt((mx*Xp)**2+(my*Yp)**2)
    return np.array([mx*Xp/d, my*Yp/d]) # (cos(phi_a), sin(phi_a))

def s_a(mx ,my, Xp, Yp):
    """
    This function returns a unit vector normal to the line with direction a(lpha)
    mx,my integers, prime to one another
    Xp and Yp dimensions of the cartesian cell
    """
    d = np.sqrt((mx*Xp)**2+(my*Yp)**2) # distance along the unfolded geometry
    return np.array([-my*Yp/d,mx*Xp/d]) # (-sin(phi_a), cos(phi_a))

def get_mxmy_pairs(nangles):
    """
    This function returns the pairs mx and my, for a given number of angles
    nangles : integer of the form p+1 where p is prime
    """
    Mx = np.zeros(nangles)
    My = np.zeros(nangles)
    for i in range(1,nangles):
        mx = i
        my = nangles-i+1
        if math.gcd(mx,my) != 1 :
            mx = mx/math.gcd(mx,my)
            my = my/math.gcd(mx,my)
        Mx[i-1] = mx
        My[i-1] = my
    
    return Mx,My

def get_equivalent_pt(rs, mx, Xp, my, Yp):
    """
    rs (np.array) = starting point on original cell
    mx,my integers, prime to one another
    Xp and Yp dimensions of the cartesian cell
    """
    return np.array([rs[0]+mx*Xp, rs[1]+my*Yp])




"""
rcomb6 = 0.4435

radii=Santamarina_radii(rcomb6, isGd=True)
print(radii)

phi_k = EQW2_angles_calculator(8, 1.295, 1.295)
print(phi_k)
"""

## Geometry definition class
Assembly_pitch = 7.62*2
Channel_box_out = 2.3975+1.1025
Channel_box_in = 2.3175+1.0225
class geom_ASSBLY:
    """
    Class used to help define and check native 3 level definition of BWR Assemly geometries
    """
    def __init__(self, pitch_A, C_list, Box_o, Chan_o, Box_i, Chan_i):
        """
        ALL dimensions in cm
        pitch_A : float, pitch for assembly 
        C_list = list of fuel pin cell object built from class geom_PIN (pins taken to be in CARCELs) assuming that the pin lattice is regular
        Box_o/i, Chan_o/i : outer and inner dimensions for outer box and coolant channel
        """
        self.pitch_A = pitch_A
        self.pitch_C = C_list[0].pitch()
        self.cell_list = C_list
        
        # Assuming square BOX and CHNL but could generalize from the class inputs
        self.Box_o = Box_o
        self.Box_i = Box_i
        self.Chan_o = Chan_o
        self.Chan_i = Chan_i
        
        # Initializing gaps and thicknesses
        self.outer_Water_gap = (self.pitch_A - self.Box_o)/2.0
        self.box_thickness = (self.Box_o - self.Box_i)/2.0
        self.chan_thickness = (self.Chan_o - self.Chan_i)/2.0


        Channel_box_xL_out = -1.1025
        Channel_box_XR_out = 2.3975


        Cell_pitch = 1.295
        L1 = Channel_box_xL_out+6.7
        L3 = 6.7-Channel_box_XR_out
        #Xtra_water = 2.0*6.7-10.0*Cell_pitch
        Xtra_water = 0.45
        print(Xtra_water)
        sum = Water_gap + Box_thickness + Xtra_water + 10*Cell_pitch + Box_thickness + Water_gap
        print(sum)
        Xtra_water_side = Xtra_water/2.0
        X1 = Water_gap+Box_thickness + Xtra_water_side
        X2 = X1 + 4*Cell_pitch 
        X3 = X2 + 3*Cell_pitch
        X4 = X3 + 3*Cell_pitch
        print(X3 + 3*Cell_pitch+Xtra_water_side)
        X_points = [0.0, X1, X2, X3, X4, 2*7.62]
        print(X_points)
        print(X4+Xtra_water_side+Box_thickness+Water_gap)


        Xmax_box=6.7-3*Cell_pitch-Xtra_water_side
        Xmin_box = -6.7+Xtra_water_side + 4*Cell_pitch

        print(Xmax_box-Xmin_box)

class geom_PIN:
    """
    Class used to help defining BWR PIN geometry
    """
    def __init__(self, pitch, fuel_radius, gap_radius, clad_radius, isGd, label):
        """
        Assuming square CARCEL pin geometry
        pitch, fuel_radius, gap_radius, clad_radius : float
        isGd : boolean value to account for presence of Gd in fuel
        """
        self.pitch = pitch
        self.outer_fuel_radius = fuel_radius
        self.gap_radius = gap_radius
        self.clad_radius = clad_radius
        self.isGd = isGd
        self.name = label
        print("Processing PIN with label "+ self.name)

    def get_Santamarina_radii(self):
        # Helper to defined fuel region radii for fuel pins --> important for evolution calculations
        # UOX : Finer discretization close to the outer radii is important for rim effects : Pu formation
        # Gd2O3 : rim effect + Gd evolution = challenge for deterministic methods, need 6 sub regions.
        """
        r_out = float, fuel radius
        isGd = bool, adapting discretization to Gd pin --> 6 radial sub regions instead.
        A. Santamarina recommandations :
        volumes for UOX pins : 50%, 80%, 95% and 100%
        volumes for Gd2O3 pins : 20%, 40%, 60%, 80%, 95% and 100%
        """
        if self.isGd==False:
            radii=[0.0,(0.5**0.5)*self.outer_fuel_radius, (0.8**0.5)*self.outer_fuel_radius, (0.95**0.5)*self.outer_fuel_radius, self.outer_fuel_radius]
            
        else :
            radii=[0.0,(0.2**0.5)*self.outer_fuel_radius, (0.4**0.5)*self.outer_fuel_radius, (0.6**0.5)*self.outer_fuel_radius, 
                      (0.8**0.5)*self.outer_fuel_radius, (0.95**0.5)*self.outer_fuel_radius, self.outer_fuel_radius]
            self.fuel_radii = radii
    
    def compute_volumes(self):
        print("compute volumes for pin not implemented yet")
        return
    
    def import_mixes(self):
        """
        idea would be to be able to import from Composition_Processor class to associate mixes with physical regions
        """
        return

