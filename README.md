# BWR_project
PhD project aimed at extending Dragon5 capabilities to treat BWR assemblies, under the supervision of Prof. Alain Hébert.

# Project structure

- Version5_beta_ev3198 : official release of Version5.0.10 March 2024
   			- In /Dragon/data/ : - AT10_PIN.x2m allows to perform pincell calculations for the ATRIUM10 BWR assembly benchmark
					     - AT10_ASSBLY.x2m allows to generate and track ATRIUM10 BWR assembly geometry using the G2S: -> SALT: modules
					     - You will also find _TSPC and _TISO tracking studies on cartesian cell aiming at correcting a bug in NXT:
			- In /Dragon/src/ : my development version for Dragon5, all successful modifications will be returned to Prof. Hébert and Marleau.
- Serpent2 : input decks for Serpent2 Monte Carlo reference Burnup Calculations
- Post_Treatment : tentative Python3 scripts to perform post treatment operations for Dragon5/Serpent2 comparisons.

Will update as the project evolves.
 


