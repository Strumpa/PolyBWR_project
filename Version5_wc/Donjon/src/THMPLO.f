*DECK THMPLO
      SUBROUTINE THMPLO(P,X,PHIL0)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the value of the corrective factor for two phase calculation
* of frictional pressure loss based on an homogeneous flow correlation
*
*Copyright:
* Copyright (C) 2014 Ecole Polytechnique de Montreal.
*
*Author(s): 
* P. Gallet
* C. Huet
* 02/2025: C. Huet - Preparation to future models
* 08/2025: M. Bellier - Implmentation of Lockhart-Martinelli correlation
*
*Parameters: input
* P       pressure (Pa)
* X       steam quality   
*
*Parameters: output
* PHIL0   corrective factor for two phase pressure loss calculation
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      REAL P,X,PHIL0
*----
*  LOCAL VARIABLES
*----
      REAL TSAT,MUL,MUG,TG,TL,R1,R2,R3, RHOL,RHOG,XLM
*----
*  COMPUTE VALUE OF THE CORRECTIVE FACTOR USING DENSITIES AND
*  VISCOSITIES OF BOTH SATURATED WATER AND DRY SATURATED STEAM
*----
*     compute the values of the thermodynamic parameters of steam and
*     liquid phases using freesteam steam tables
      CALL THMSAT(P,TSAT)
      TG=TSAT+0.01
      TL=TSAT-0.01
      CALL THMPT(P,TL,RHOL,R1,R2,MUL,R3)
      CALL THMPT(P,TG,RHOG,R1,R2,MUG,R3)
*- CORRELATION = ? 
*     PHIL0=(1+X*(RHOL/RHOG-1))/((1+X*(MUL/MUG-1))**0.25)
*- 
* - LOCKHART-MARTINELLI CORRELATION
      XLM = ((1-X)/X)**0.9*(RHOG/RHOL)**0.5*(MUG/MUL)**0.1
      PHIL0 = (1.0 + 20/XLM + 1.0/XLM**2)**0.5

      RETURN
      END