*DECK THMFRI
      SUBROUTINE THMFRI(REY,EPS,HD,FRIC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the value of the friction factor coefficient with :
* - Laminar flow correlation based on condition on Reynolds number
* - Muller Steinhagen correlation formula (single phase)
* - Churchill's correlation in two phase flows
*
*Copyright:
* Copyright (C) 2013 Ecole Polytechnique de Montreal.
*
*Author(s): 
* P. Gallet (creation)
* 07/08/2025 : Modified by M. Bellier to include Churchill
*
*Parameters: input
* REY     reynolds number 
* EPS     void fraction
* HD      hydraulic diameter
*
*Parameters: output
* FRIC    friction factor coefficient
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      REAL REY,FRIC,HD,EPS,R
*----
*  COMPUTE VALUE OF THE FRICTION FACTOR COEFFICIENT AS FUNCTION OF THE 
*  REYNOLDS NUMBER
*----

! Laminar flow
      IF (REY.LE.1187.0) THEN
            FRIC=64.0/REY
! Blasius-like correlation used by C. Huet in his python prototype
      ELSE IF (EPS.LT.0.002) THEN
            FRIC=0.3164/(REY**0.25)
! Churchill's correlation
      ELSE
            R = 0.0000004/HD !Relative roughness=Roughness/Hydraulic Diameter
            FRIC=8*(((8.0/REY)**12)+((2.475*LOG(((7/REY)**0.9)+0.27*R))
     >      **16+(37530/REY)**16)**(-1.5))**(0.0833333)
      ENDIF

      RETURN
      END
