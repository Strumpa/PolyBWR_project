# Python class to handle and faciliate geometry definitions in Serpent2 applications

# Author : R. Guasch 
# Date : 08/09/2025
# Purpose : assist in Validation procedure of BWR lattice calculations with DRAGON5.

class Pin_Universe_Definition:
    def __init__(self, name:str, materials:list, radii:list):
        """
        Args:
            name (str): name identifying pin universe
            materials (list): list of S2_materials associated with each pin universe cell.
            radii (list): list of floats corresponding to radii delimiting co-centric cylinders making up the inner cells of pin universe.
        """
        # Check that lenght of materials = length of radii + 1
        if (len(materials) != len(radii)+1):
            raise ValueError(f"Invalid length of materials {len(materials)} and radii {len(radii)} for Serpent2 pincell universe definition") 

        self.name = name
        self.materials = materials
        self.radii = radii
        self.pincell_card = ""


    def format_pincell_card(self):
        """
        format pincell card to faciliate Serpent2 deck generation
        """ 
        self.pincell_card += f"pin {self.name}\n"
        for i in range len(self.radii):
            self.pincell_card += f"\t {self.materials[i].name} {self.radii[i]}\n"
        # Add surrounding material
        self.pincell_card += f"\t {self.materials[-1].name} \n"

    def print_pincell_card(self):
        print(self.pincell_card)

class Lattice_Definition:
    def __init__(self, name:str, pitch:float, center_x:float, center_y:float, nx:int, ny:int, pin_universes:list):
        """class allowing to automate 2D cartesian lattice definitions in Serpent2

        Args:
            name (str): name of the lattice universe  
            pitch (float): pitch of the lattice
            center_x/y (float): x/y coordinates definiting the lattice's center on 2D x-y plane
            nx/ny (int): number of lattice elements in x and y directions 
            pin_universes (list of universes): lattice description, list of universes filling the lattice. Convention is x-increasing y-increasing ordering.
        """

        # check that nx * ny = total number of pincell universes
        if nx*ny != len(pin_universes):
            raise ValueError(f"Invalid length of pincell universes ({len(pin_universes)}) for regular cartesian {nx} by {ny} lattice.")
        
        self.name = name
        self.pitch = pitch
        self.x0 = center_x
        self.y0 = center_y
        self.nx = nx
        self.ny = ny 
        self.pin_universes_list = pin_universes

        self.lattice_card = ""

    def format_lattice_card(self):
        """
        fomatting lattice card according to Serpent2 specifications
        """

        self.lattice_card += f"lat {self.name} 1 {self.x0:.2f} {self.y0:.2f} {self.nx} {self.ny} {self.pitch:.3f}\n"
        pin_univ_counter = 1
        for pin in pin_universes_list:
            self.lattice_card += f"{pin.name} "
            pin_univ_counter += 1
            if pin_univ_counter%self.nx == 0:
                self.lattice_card += "\n" # skip a line to indicate change in y coord once each x-increasing line is full

    def print_lattice_card(self):
        print(self.lattice_card)


