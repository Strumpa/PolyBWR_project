*DECK THMSCD
      REAL FUNCTION THMSCD(TEMP,FTP,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the product of the heat capacity of fuel (in J/Kg/K) times
* its density (in Kg/m^3). Version for molten salts.
*
*Copyright:
* Copyright (C) 2024 Ecole Polytechnique de Montreal.
*
*Author(s): 
* C. Garrido
*
*Parameters: input
* TEMP    fuel temperature in Kelvin.
* FTP     tpdata object with correlations to obtain properties of 
*         molten salt.
* 
*Parameters: output
* THMSCD  product of the heat capacity of fuel times its density
*         (in J/K/m^3).
*
*-----------------------------------------------------------------------
*
      USE t_saltdata
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(tpdata) FTP
      REAL TEMP
      INTEGER IMPX
*----
*  LOCAL VARIABLES
*  CP:    heat capacity in J/Kg/K
*  RHO: fuel density Kg/m^3
*----
      REAL CP,RHO,R2,R3,R4
*
      CALL THMSPT(FTP,TEMP,RHO,R2,R3,R4,CP,IMPX)
      THMSCD=RHO*CP
      RETURN
      END
