## BWR_project
PhD project aimed at extending Dragon5 capabilities to treat Small Modular Reactors of BWR type, under the supervision of Prof. Alain Hébert.

# Project structure

- Version5_ev3205 : first patch after official release of Version5.0.10 March 2024  
  
- In /Dragon/data/ : - AT10_PIN.x2m allows to perform pincell calculations for the ATRIUM10 BWR assembly benchmark  
  
- AT10_ASSBLY.x2m allows to generate and track ATRIUM10 BWR assembly geometry using the G2S: -> SALT: modules    
- You will also find _TSPC and _TISO tracking studies on cartesian cell aiming at correcting a bug in NXT:
- In Version5_wc/Dragon/src/ : my development version for Dragon5, all successful modifications will be returned to Prof. Hébert and Marleau.  
- In PyGan/ : you will find the Python V5-kit wrapper functions used for running calculations and post treating Dragon5 vs Serpent2 results  

- Serpent2 : input decks for Serpent2 Monte Carlo reference Burnup Calculations
- PTT : tentative Python3 scripts to perform post treatment operations for Dragon5/Serpent2 comparisons.
   - protoype for DMLGInterface project,
   - extensive use of serpentTools package https://github.com/CORE-GATECH-GROUP/serpent-tools, Andrew Johnson, Dan Kotlyar, Stefano Terlizzi, and Gavin Ridley, “serpentTools: A Python Package for Expediting Analysis with Serpent,” Nuc. Sci. Eng, (in press) (2020). For Serpent2 post treatment.


# Compiling Version5
In order to compile Version5 with its "V5-Pykit" python interface, you just need to :
 			- cd PyGan
      - make hdf5=1
      - make tests (to execute non regression tests)
This requires that you have hdf5 installed, and that you have specified its path, along with the path to a fortran compiler in your .profile.

# Compiling Dragon5 or Donjon5
Standalone versions of Dragon5/Donjon5 can be compiled as such :
      - cd Verion5_wc/Dragon or cd Version5_wc/Donjon
      - make (hdf5=1, optional) 
      - make tests (to execute non regression tests)
It must me noted that some tests are dependant on cross sectional libraries which can be found at http://http://merlin.polymtl.ca/libraries.htm
Will update as the project evolves.
 


