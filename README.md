## BWR_project
PhD project aimed at extending Dragon5 capabilities to treat Small Modular Reactors of BWR type, under the supervision of Prof. Alain Hébert.

# Project structure

- Version5_wc : Working copy of Version5 envirionment used to develop BWR calculation schemes
  
- Version5 includes DRAGON : lattice code compnent, DONJON full core code component.

- Serpent2 : input decks for Serpent2 Monte Carlo reference Burnup Calculations
- PTT : tentative Python3 scripts to perform post treatment operations for Dragon5/Serpent2 comparisons.
   - extensive use of serpentTools package https://github.com/CORE-GATECH-GROUP/serpent-tools, Andrew Johnson, Dan Kotlyar, Stefano Terlizzi, and Gavin Ridley, “serpentTools: A Python Package for Expediting Analysis with Serpent,” Nuc. Sci. Eng, (in press) (2020). For Serpent2 post treatment.

- DMLG : collection of python classes and scripts allowing for easier geometry handling, DMLG_geometry_handling/ is aimed at prototyping developments to be made to efficiently support BWR assembly geometries in the Version5 environment.


# Compiling Version5
In order to compile Version5 with its "V5-Pykit" python interface, you just need to :
 			- cd PyGan
      - make hdf5=1 (openmp=1, optional)
      - make tests (to execute non regression tests)
This requires that you have hdf5 installed, and that you have specified its path, along with the path to a fortran compiler in your .profile.

# Compiling Dragon5 or Donjon5
Standalone versions of Dragon5/Donjon5 can be compiled as such :
      - cd Verion5_wc/Dragon or cd Version5_wc/Donjon
      - make (hdf5=1, optional) (openmp=1, optional)
      - make tests (to execute non regression tests)
It must me noted that some tests are dependant on cross sectional libraries which can be found at http://http://merlin.polymtl.ca/libraries.htm
Will update as the project evolves.
 


