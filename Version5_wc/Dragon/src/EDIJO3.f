*DECK EDIJO3
      SUBROUTINE EDIJO3(IPMAC2,IPTRK1,IPFLUX,IPRINT,NGCOND,IGCOND)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover ALBS information from last component of unknown array for use
* with SPH equivalence techniques. Multicell surfacic compatible
* version. It is activated with ARM keyword in ASM: module.
*
*Copyright:
* Copyright (C) 2025 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPMAC2  pointer to condensed macrolib information (L_MACROLIB
*         signature) built by EDI:.
* IPTRK1  pointer to the reference tracking object.
* IPFLUX  pointer to the reference solution (L_FLUX signature).
* IPRINT  print index.
* NGCOND  number of condensed groups.
* IGCOND  limit of condensed groups.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAC2,IPTRK1,IPFLUX
      INTEGER IPRINT,NGCOND,IGCOND(NGCOND)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      TYPE(C_PTR) JPFLUX
      INTEGER ISTATE(NSTATE)
      CHARACTER CDOOR*12
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NMC_SURF,IFR,MIX,INUM,IGEN
      REAL, ALLOCATABLE, DIMENSION(:) :: ALB,SUR,WORKD
      REAL, ALLOCATABLE, DIMENSION(:,:) :: OUTG
*----
*  RECOVER FLUX OBJECT INFORMATION
*----
      CALL LCMGET(IPFLUX,'STATE-VECTOR',ISTATE)
      NUNKNO=ISTATE(2)
      ILEAK=ISTATE(7)
*----
*  RECOVER TRACKING INFORMATION
*----
      CALL LCMGTC(IPTRK1,'TRACK-TYPE',12,CDOOR)
      CALL LCMGET(IPTRK1,'STATE-VECTOR',ISTATE)
      IF((CDOOR.NE.'EXCELL').OR.(ISTATE(7).NE.5)) THEN
        CALL XABORT('EDIJO3: MULTICELL SURFACIC OPTION NOT ACTIVATED.')
      ENDIF
      NREG=ISTATE(1)
      NUNKNO=ISTATE(2)+ISTATE(28)
      NMACRO=ISTATE(24)
      IF(NMACRO.EQ.0) CALL XABORT('EDIJO3: NO MACRO GEOMETRIES.')
      NMCEL=NMACRO
      NMERGE=NMACRO
      ALLOCATE(IGEN(NMERGE),INUM(NMCEL),NMC_SURF(NMACRO+1))
      DO IK=1,NMERGE
        IGEN(IK)=IK
      ENDDO
      DO IK=1,NMCEL
        INUM(IK)=IK
      ENDDO
      IF(NMACRO.EQ.0) CALL XABORT('EDIJO3: MACRO OPTION IS MANDATORY.')
      CALL LCMGET(IPTRK1,'NMC_SURF',NMC_SURF)
      NMIX=NMC_SURF(NMACRO+1)
      NIFR=NMC_SURF(NMACRO+1)
      ALLOCATE(IFR(NIFR),ALB(NIFR),MIX(NMIX),SUR(NMIX))
      CALL LCMGET(IPTRK1,'IFR',IFR)
      CALL LCMGET(IPTRK1,'ALB',ALB)
      CALL LCMGET(IPTRK1,'MIX',MIX)
      CALL LCMGET(IPTRK1,'SUR',SUR)
*----
*  COMPUTE THE OUTGOING CURRENT
*----
      ALLOCATE(OUTG(NGCOND,2))
      IGRFIN=0
      CALL LCMSIX(IPMAC2,'ADF',1)
      DO 70 IGRCD=1,NGCOND
      OUTG(IGRCD,:2)=0.0
      IGRDEB=IGRFIN+1
      IGRFIN=IGCOND(IGRCD)
      CALL LCMLEN(IPFLUX,'FLUX',ILON,ITYLCM)
      IF(ILON.EQ.0) CALL XABORT('EDIJO3: MISSING FLUX INFO(1).')
      JPFLUX=LCMGID(IPFLUX,'FLUX')
      DO 60 IGR=IGRDEB,IGRFIN
      CALL LCMLEL(JPFLUX,IGR,ILCMLN,ITYLCM)
      IF(ILCMLN.EQ.0) CALL XABORT('EDIJO3: MISSING FLUX INFO(2).')
      IF(ILEAK.LE.5) THEN
        IF(ILCMLN.NE.NUNKNO) CALL XABORT('EDIJO3: ARM KEYWORD MUST B'
     1  //'E SET IN ASM: MODULE(1).')
        ALLOCATE(WORKD(NUNKNO))
      ELSE IF(ILEAK.EQ.6) THEN
        IF(ILCMLN.NE.2*NUNKNO) CALL XABORT('EDIJO3: ARM KEYWORD MUST'
     1  //' BE SET IN ASM: MODULE(2).')
        ALLOCATE(WORKD(2*NUNKNO))
      ELSE
        CALL XABORT('EDIJO3: INVALID TYPE OF LEAKAGE.')
      ENDIF
      CALL LCMGDL(JPFLUX,IGR,WORKD)
      OUTC1=0.0
      OUTC2=0.0
      SURT=0.0
      DO 50 ICEL=1,NMCEL
      IKK=INUM(ICEL)
      IKG=IGEN(IKK)
      IF(IKK.EQ.0) GO TO 50
      J3=NMC_SURF(IKG+1)-NMC_SURF(IKG)
      IT=0
      DO IK=1,IKK-1
        IT=IT+(NMC_SURF(IGEN(IK)+1)-NMC_SURF(IGEN(IK)))
      ENDDO
      IS=0
      DO IK=1,ICEL-1
        IS=IS+(NMC_SURF(IGEN(INUM(IK))+1)-NMC_SURF(IGEN(INUM(IK))))
      ENDDO
      DO 40 JC=1,J3
      IF((MIX(IT+JC).EQ.IFR(IS+JC)).AND.(SUR(IS).NE.0.0)) THEN
        J1=IFR(IS+JC)
        OUTC1=OUTC1+WORKD(NREG+J1)*SUR(IS+JC)
        OUTC2=OUTC2+WORKD(NREG+J1)*SUR(IS+JC)*ALB(IS+JC)
        SURT=SURT+SUR(IS+JC)
      ENDIF
   40 CONTINUE
   50 CONTINUE
      DEALLOCATE(NMC_SURF,INUM,IGEN)
      DEALLOCATE(SUR,MIX,ALB,IFR)
      OUTG(IGRCD,1)=OUTG(IGRCD,1)+OUTC1/SURT
      OUTG(IGRCD,2)=OUTG(IGRCD,2)+OUTC2/SURT
      DEALLOCATE(WORKD)
   60 CONTINUE
   70 CONTINUE
      CALL LCMPUT(IPMAC2,'ALBS00',NGCOND*2,2,OUTG)
      IF(IPRINT.GT.3) THEN
         WRITE(6,900) (OUTG(IGR,1),IGR=1,NGCOND)
         WRITE(6,910) (OUTG(IGR,2),IGR=1,NGCOND)
         WRITE(6,'(/)')
      ENDIF
      CALL LCMSIX(IPMAC2,' ',2)
      DEALLOCATE(OUTG)
      RETURN
*
  900 FORMAT(/49H EDIJO3: OUT-CURRENTS (4J-/S) PER MACRO-GROUP ARE/
     > (1X,1P,10E13.5))
  910 FORMAT(/49H EDIJO3:  IN-CURRENTS (4J+/S) PER MACRO-GROUP ARE/
     > (1X,1P,10E13.5))
      END
