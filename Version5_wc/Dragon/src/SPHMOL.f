*DECK SPHMOL
      SUBROUTINE SPHMOL(IPMPO,ICAL,NMIL,NGROUP,NSURFD,HEDIT,VOSAP,
     & DFACT,VFLUX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover surface flux data from a MPO file generated with APOLLO3.
*
*Copyright:
* Copyright (C) 2024 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IPMPO   pointer to the MPO file.
* ICAL    index of the elementary calculation being considered.
* NMIL    number of mixtures in the elementary calculation.
* NGROUP  number of energy groups in the elementary calculation.
* NSURFD  number of surfaces in a mixture.
* HEDIT   name of output group for a (multigroup mesh, output geometry)
*         couple (generally equal to 'output_0').
* VOSAP   mixture volumes in the MPO file.
*
*Parameters: output
* DFACT   discontinuity factors.
* VFLUX   averaged volume fluxes.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMPO
      INTEGER ICAL,NMIL,NGROUP,NSURFD
      REAL VOSAP(NMIL),DFACT(NMIL,NGROUP,NSURFD),VFLUX(NMIL,NGROUP)
      CHARACTER(LEN=12) HEDIT
*----
*  LOCAL VARIABLES
*----
      INTEGER RANK,TYPE,TYPE2,DIMSR(5)
      CHARACTER RECNAM*80,HSMG*131
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: VREAL,SURF,LG
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SURFLX
*----
*  RECOVER ASSEMBLY DISCONTINUITY FACTORS
*----
      WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0,6H/flux/)')
     & TRIM(HEDIT),ICAL-1
      CALL hdf5_info(IPMPO,TRIM(RECNAM)//"NSURF",RANK,TYPE,NBYTE,DIMSR)
      CALL hdf5_info(IPMPO,TRIM(RECNAM)//"SURFFLUX",RANK,TYPE2,NBYTE,
     & DIMSR)
      IF((TYPE.EQ.99).OR.(TYPE2.EQ.99)) THEN
        CALL hdf5_list(IPMPO,TRIM(RECNAM))
        CALL XABORT('SPHMOL: UNABLE TO FIND ADF INFORMATION.')
      ENDIF
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"NSURF",NSURFD)
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"SURFFLUX",SURFLX)
      CALL hdf5_info(IPMPO,TRIM(RECNAM)//"SURF",RANK,TYPE,NBYTE,DIMSR)
      IF(TYPE.NE.99) THEN
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"SURF",SURF)
        IF(DIMSR(1).NE.NMIL*NSURFD) THEN
          WRITE(HSMG,'(24HSPHMOL: INVALID LENGTH (,I5,11H) FOR SURF ,
     &    14HGROUP. LENGTH=,I5,10H EXPECTED.)') DIMSR(1),NSURFD
          CALL XABORT(HSMG)
        ENDIF
      ELSE
*       temporary.....
        CALL hdf5_read_data(IPMPO,"/geometry/geometry_0/COORDINATE",LG)
        ALLOCATE(SURF(NSURFD))
        SURF(:NSURFD)=LG(2)
        DEALLOCATE(LG)
      ENDIF
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"TOTALFLUX",VREAL)
      DO IGR=1,NGROUP
        DO IBM=1,NMIL
          VFLUX(IBM,IGR)=VREAL(IGR)/VOSAP(IBM)
        ENDDO
      ENDDO
      print *,'NSURFD=',NSURFD,' NMIL=',NMIL,' NGROUP=',NGROUP
      DO I=1,NSURFD
        DO IGR=1,NGROUP
          DO IBM=1,NMIL
            DFACT(IBM,IGR,I)=SURFLX(I,IGR)/(VFLUX(IBM,IGR)*SURF(I))
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(VREAL,SURF,SURFLX)
      RETURN
      END
