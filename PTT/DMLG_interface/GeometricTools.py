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





## Geometry definition class
# Used to check analytical volumes and compare with Serpent2/Dragon5 volumes
class geom_ASSBLY:
    """
    Class used to help define and check 3 level definition of BWR Assemly geometries
    """
    def __init__(self, pitch_A, C_list, Box_o, Chan_o, Box_i, Chan_i, controlled=False):
        """
        ALL dimensions in cm
        pitch_A : float, pitch for assembly 
        C_list = list of fuel pin cell object built from class geom_PIN (pins taken to be in CARCELs) assuming that the pin lattice is regular
        Box_o/i, Chan_o/i : outer and inner dimensions for outer box and moderator channel
        """
        self.pitch_A = pitch_A
        self.pin_names=C_list
        
        
        # Assuming square BOX and CHNL but could generalize from the class inputs
        self.Box_o = Box_o
        self.Box_i = Box_i
        self.Chan_o = Chan_o
        self.Chan_i = Chan_i
        self.isControlled = controlled
        
        # Initializing gaps and thicknesses
        self.outer_Water_gap = (self.pitch_A - self.Box_o)/2.0
        self.box_thickness = (self.Box_o - self.Box_i)/2.0
        self.chan_thickness = (self.Chan_o - self.Chan_i)/2.0


    def setPins(self, fuel_radius, gap_radius, clad_radius):
        """
        pins_names = list of pin names
        fuel_radius, gap_radius and clad_radius are radii used to define the fuel pin's geometry.
        It is assumed that all pins share the same geometry but it would be relevant to implement possibilities for different pins
        difference between Gd and UOX pins is performed through name identifier.
        """
        self.pins = []
        for pin in self.pin_names:
            if "Gd" in pin or "GD" in pin or "GADO" in pin :
                print(f"creating pin {pin}")
                self.pins.append(geom_PIN(pin, fuel_radius, gap_radius, clad_radius, isGd=True, height=1))
            else:
                self.pins.append(geom_PIN(pin, fuel_radius, gap_radius, clad_radius, isGd=False, height=1))
        return
    
    def setNumberofPins(self, pins_dict):
        """
        pins_dict = dictionnay with keys = pins_names, values = number of each pin
        """
        self.numberOfPinsperType = pins_dict
        total=0
        for type in self.numberOfPinsperType.keys():
            total+=self.numberOfPinsperType[type]
        self.Total_Nb_Pins = total
        return
    
    def computeVolumes(self):
        """
        compute total volume for each type of fuel and each subdivision in santamarina radii 
        compute volumes of structural elements and coolant/moderator
        """
        self.Volumes = {}
        for pin in self.pins:
            for node in pin.volumes.keys():
                self.Volumes[node] = pin.volumes[node]*self.numberOfPinsperType[pin.name]
        
        self.Volumes["box"] = (self.Box_o**2-self.Box_i**2 + self.Chan_o**2 - self.Chan_i**2)
        self.Volumes["clad"] = (np.pi*self.pins[0].clad_radius**2 - np.pi*self.pins[0].gap_radius**2) * self.Total_Nb_Pins
        self.Volumes["gap"] = (np.pi*self.pins[0].gap_radius**2 - np.pi*self.pins[0].outer_fuel_radius**2) * self.Total_Nb_Pins

        if self.isControlled:
            if self.ctrl_cross_symmetry:
                self.Volumes["ctrl_rod"] = self.nb_ctrl_rods*np.pi*self.ctrl_rod_radius**2/2
                self.Volumes["ctrl_cross"] = 2*(self.cross_half_length-self.cross_half_width)*self.cross_half_width + self.cross_half_width**2 - self.Volumes["ctrl_rod"]
            else:
                self.Volumes["ctrl_rod"] = self.nb_ctrl_rods*np.pi*self.ctrl_rod_radius**2
                self.Volumes["ctrl_cross"] = 4*(self.cross_half_length-self.cross_half_width)*self.cross_half_width + self.cross_half_width**2 - self.Volumes["ctrl_rod"]
            self.Volumes["moder"] = self.pitch_A**2 - self.Box_o**2 + self.Chan_i**2 - self.Volumes["ctrl_rod"] - self.Volumes["ctrl_cross"]
            # compute the sum of the existing volumes in order to compute the coolant volume
            sum_vols = 0
            for key in self.Volumes.keys():
                sum_vols+=self.Volumes[key]
        else:
            self.Volumes["moder"] = self.pitch_A**2 - self.Box_o**2 + self.Chan_i**2
            # compute the sum of the existing volumes in order to compute the coolant volume
            sum_vols = 0
            for key in self.Volumes.keys():
                sum_vols+=self.Volumes[key]
        self.Volumes["cool"] = self.pitch_A**2 - sum_vols
        return
    
    def setCtrlCross(self, cross_half_length, cross_half_width, nb_ctrl_rods, ctrl_rod_radius, ctrl_wings_symmetry):
        """
        cross_half_length, cross_half_width : float, dimensions of the control rod cross section
        nb_ctrl_rods : integer, number of control rods
        ctrl_rod_radius : float, radius of the control rod
        ctrl_wings_symmetry : boolean to account for symmetry of the control rod wings, if True : only half the wings are defined, if False : full wings are defined as a cross.
        """
        self.ctrl_rod_radius = ctrl_rod_radius
        self.cross_half_length = cross_half_length
        self.cross_half_width = cross_half_width
        self.nb_ctrl_rods = nb_ctrl_rods
        self.ctrl_cross_symmetry = ctrl_wings_symmetry
        return

    

class geom_PIN:
    """
    Class used to help defining BWR PIN geometry
    """
    def __init__(self, label, fuel_radius, gap_radius, clad_radius, isGd, height, pitch=None, sectorize=True):
        """
        Assuming square CARCEL pin geometry
        pitch, fuel_radius, gap_radius, clad_radius : float
        isGd : boolean value to account for presence of Gd in fuel
        """
        self.outer_fuel_radius = fuel_radius
        self.gap_radius = gap_radius
        self.clad_radius = clad_radius
        self.isGd = isGd
        self.name = label
        self.height = height
        self.isSectorized = sectorize
        self.pitch = pitch
        print("Processing PIN with label "+ self.name)
        if self.isSectorized:
            self.computeSantamarinaradii()
        else:
            self.pin_radii = [self.outer_fuel_radius, self.gap_radius, self.clad_radius]
        if self.pitch != None:
            self.setCellPitch(pitch)
        self.computePinVolumes()

    def computeSantamarinaradii(self):
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
            self.pin_radii=[(0.5**0.5)*self.outer_fuel_radius, (0.8**0.5)*self.outer_fuel_radius, (0.95**0.5)*self.outer_fuel_radius, self.outer_fuel_radius, self.gap_radius, self.clad_radius]
            
        else :
            self.pin_radii=[(0.2**0.5)*self.outer_fuel_radius, (0.4**0.5)*self.outer_fuel_radius, (0.6**0.5)*self.outer_fuel_radius, 
                      (0.8**0.5)*self.outer_fuel_radius, (0.95**0.5)*self.outer_fuel_radius, self.outer_fuel_radius, self.gap_radius, self.clad_radius]
        return
    
    def computePinVolumes(self):
        print("computing volumes for pins")
        volumes = []
        self.volumes = {}
        self.volume_names_to_radii = {}
        if self.isSectorized:
            if self.isGd:
                key_dict = {"1":"A","2":"B","3":"C","4":"D","5":"E","6":"F", "7":"Gap", "8":"Clad", "9":"Coolant"}
            else:
                key_dict = {"1":"A","2":"B","3":"C","4":"D", "5":"Gap", "6":"Clad", "7":"Coolant"}
        else:
            key_dict = {"1":"Fuel", "2":"Gap", "3":"Clad", "4":"Coolant"}

        self.mix_to_name_connectivity_dict = key_dict
        for radius in self.pin_radii:
            Area = np.pi*radius**2
            volumes.append(Area*self.height)
        if self.pitch != None:
            volumes.append(self.pitch**2*self.height)
        for i in range(len(volumes)):
            #if self.isGd:
            if i==0:
                self.volumes[self.name+"_"+key_dict[str(i+1)]] = volumes[i]
                self.volume_names_to_radii[self.name+"_"+key_dict[str(i+1)]] = self.pin_radii[i]
            else:
                self.volumes[self.name+"_"+key_dict[str(i+1)]] = volumes[i]-volumes[i-1]
                self.volume_names_to_radii[self.name+"_"+key_dict[str(i+1)]] = self.pin_radii[i]
        print(f"$$$ -- In computePinVolumes() volumes are {self.volumes}")
        return  
    
    def setCellPitch(self, pitch):
        self.pitch = pitch
        return
    
    
    def setMixes(self, composition_dict):
        """
        idea would be to be able to import from Composition_Processor class to associate mixes with physical regions

        composition_dict : dictionnary with keys = region names, values = dictionary with keys = isotope names, values = atomic densities
        """
        # check the number of regions / mixes.
        #print(self.volumes)
        #print(compition_dict)
        if len(self.volumes.keys()) != len(composition_dict.keys()):
            print("Error : number of regions does not match number of mixes")
            return
        else:
            self.mixes = composition_dict
        return






