*DECK TINCHH
      SUBROUTINE TINCHH(IPMAP,NCH,IMPX,NAMCHA,TTIME,RFCHAN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute 'REF-CHAN' record in L_MAP object for history-based cases in
* hexagonal coordinates.
*
*Copyright:
* Copyright (C) 2026 Ecole Polytechnique de Montreal
*
*Author(s): 
* G. Calandrino
*
*Parameters: input
* IPMAP   pointer to fuel-map information.
* NCH     number of channels
* IMPX    print flag
* NAMCHA  channel name
* TTIME   refuelling time
*
*Parameters: output
* RFCHAN  time values at which channels are refueled inside a refueling
*         time period
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAP
      INTEGER NCH,IMPX
      CHARACTER*(*) NAMCHA
      REAL RFCHAN(NCH)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      INTEGER ISTATE(NSTATE)
      CHARACTER HNAM*8
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MIX, ICHMAP
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IHN
*
      CALL LCMSIX(IPMAP,'GEOMAP',1)
      CALL LCMGET(IPMAP,'STATE-VECTOR',ISTATE)
      IF(ISTATE(1).NE.9) CALL XABORT('TINCHH: 3-D HEXAGONAL GEOMETRY'
     +    //' REQUIRED')
      NH = ISTATE(3)
      NZ = ISTATE(5)
      NREG = ISTATE(6)
      ALLOCATE(MIX(NREG),IHN(2,NH),ICHMAP(NH))

      CALL LCMSIX(IPMAP,' ',2)
      CALL LCMGET(IPMAP,'HNAME',IHN)

      INAMH = 0
      DO 10 I=1,NH
        WRITE(HNAM,'(2A4)') IHN(1,I),IHN(2,I)
        IF(HNAM.EQ.NAMCHA) THEN
           INAMH = I
           GOTO 20
        ENDIF
  10  CONTINUE
*
  20  CALL LCMGET(IPMAP,'BMIX',MIX)
      ICHMAP(:NH)=0
      ICH=0
      DO 25 IH=1,NH
      DO 23 IZ=1,NZ
      IF(MIX((IZ-1)*NH+IH).NE.0) GO TO 24
  23  CONTINUE
      GO TO 25
  24  ICH=ICH+1
      ICHMAP(IH)=ICH
  25  CONTINUE
      IF(ICH.NE.NCH) CALL XABORT('@TINCHH: INVALID NUMBER OF CHANNELS')
  
      ICHANAM = ICHMAP(INAMH)

      IF(ICHANAM.EQ.0) CALL XABORT('TINCHH: WRONG CHANNEL NAME')
      DEALLOCATE(IHN,ICHMAP,MIX)
      RFCHAN(ICHANAM) = TTIME
      IF(IMPX.GT.0) THEN
        WRITE(6,*) 'TINCHH: REFUEL ',NAMCHA,' NUMBER ',I,' AT TIME ',
     1  TTIME
      ENDIF
      RETURN
      END
