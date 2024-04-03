## BWR_project
PhD project aimed at extending Dragon5 capabilities to treat Small Modular Reactors of BWR type, under the supervision of Prof. Alain Hébert. Test

# Project structure

- Version5_ev3205 : first patch after official release of Version5.0.10 March 2024  
  
- In /Dragon/data/ : - AT10_PIN.x2m allows to perform pincell calculations for the ATRIUM10 BWR assembly benchmark  
  
- AT10_ASSBLY.x2m allows to generate and track ATRIUM10 BWR assembly geometry using the G2S: -> SALT: modules    
- You will also find _TSPC and _TISO tracking studies on cartesian cell aiming at correcting a bug in NXT:
- In /Dragon/src/ : my development version for Dragon5, all successful modifications will be returned to Prof. Hébert and Marleau.  
- In PyGan/ : you will find the Python V5-kit wrapper functions used for running calculations and post treating Dragon5 vs Serpent2 results  

- Serpent2 : input decks for Serpent2 Monte Carlo reference Burnup Calculations
- Post_Treatment : tentative Python3 scripts to perform post treatment operations for Dragon5/Serpent2 comparisons.


# Compiling Version5
In order to compile Version5 with its "V5-Pykit" python interface, you just need to :
 			- cd PyGan
    			- make hdf5=1
This requires that you have hdf5 installed, and that you have specified its path, along with the path to a fortran compiler in your .profile.
Will update as the project evolves.
 


