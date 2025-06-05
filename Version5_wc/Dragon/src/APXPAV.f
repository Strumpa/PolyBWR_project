*DECK APXPAV
      SUBROUTINE APXPAV(IPAPX,IPAR,NPAR,TYPE,RVAL,IVAL,CVAL,IV,LGNEW)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To return the index of a global parameter value in the Apex file.
* Reorganize the 'paramvalues' group if required.
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
* IPAPX   pointer to the Apex file.
* IPAR    index of the global parameter.
* NPAR    total number of global parameters.
* TYPE    type of the global parameter value.
* RVAL    global parameter value if TYPE='FLOTTANT'.
* IVAL    global parameter value if TYPE='ENTIER'.
* CVAL    global parameter value if TYPE='CHAINE'.
*
*Parameters: output
* IV      index of the global parameter value.
* LGNEW   new parameter flag (=.true. if the parameter value is new).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPAPX
      INTEGER IPAR,NPAR,IV,IVAL
      REAL RVAL
      LOGICAL LGNEW
      CHARACTER TYPE*8,CVAL*(*)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (REPS=1.0E-5)
      CHARACTER RECNAM*80
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NVALUE,VINTE,VINTE_V1,
     1 DIMS_APX
      REAL, ALLOCATABLE, DIMENSION(:) :: VREAL,VREAL_V1
      CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: VCHAR,VCHAR_V1
*
      CALL hdf5_read_data(IPAPX,"/paramdescrip/NVALUE",NVALUE)
      WRITE(RECNAM,'(17H/paramvalues/PVAL,I8)') IPAR
*
      LGNEW=.TRUE.
      IF(TYPE.EQ.'FLOTTANT') THEN
         ALLOCATE(VREAL(NVALUE(IPAR)+1))
         IF(NVALUE(IPAR).EQ.0) THEN
            IV=1
            VREAL(IV)=RVAL
            NVALUE(IPAR)=1
         ELSE
            CALL hdf5_get_shape(IPAPX,RECNAM,DIMS_APX)
            ILONG=DIMS_APX(1)
            IF(ILONG.GT.NVALUE(IPAR)) CALL XABORT('APXPAV: NVALUE OVER'
     1      //'FLOW(1).')
            CALL hdf5_read_data(IPAPX,RECNAM,VREAL_V1)
            VREAL(:ILONG)=VREAL_V1(:ILONG)
            DEALLOCATE(VREAL_V1)
            DO 10 I=1,NVALUE(IPAR)
               IF(RVAL.LE.VREAL(I)*(1.+REPS))THEN
                  IV=I
                  LGNEW=RVAL.LT.VREAL(IV)*(1.-REPS)
                  GO TO 20
               ENDIF
   10       CONTINUE
            IV=NVALUE(IPAR)+1
   20       IF(LGNEW) THEN
               NVALUE(IPAR)=NVALUE(IPAR)+1
               DO 30 J=NVALUE(IPAR)-1,IV,-1
                  VREAL(J+1)=VREAL(J)
   30          CONTINUE
               VREAL(IV)=RVAL
            ENDIF
         ENDIF
         IF(LGNEW) CALL hdf5_write_data(IPAPX,TRIM(RECNAM),
     1   VREAL(:NVALUE(IPAR)))
         DEALLOCATE(VREAL)
      ELSE IF(TYPE.EQ.'ENTIER') THEN
         ALLOCATE(VINTE(NVALUE(IPAR)+1))
         IF(NVALUE(IPAR).EQ.0) THEN
            IV=1
            VINTE(IV)=IVAL
            NVALUE(IPAR)=1
         ELSE
            CALL hdf5_get_shape(IPAPX,RECNAM,DIMS_APX)
            ILONG=DIMS_APX(1)
            IF(ILONG.GT.NVALUE(IPAR)) CALL XABORT('APXPAV: NVALUE OVER'
     1      //'FLOW(2).')
            CALL hdf5_read_data(IPAPX,RECNAM,VINTE_V1)
            VINTE(:ILONG)=VINTE_V1(:ILONG)
            DEALLOCATE(VINTE_V1)
            DO 40 I=1,NVALUE(IPAR)
               IF(IVAL.LE.VINTE(I))THEN
                  IV=I
                  LGNEW=IVAL.LT.VINTE(IV)
                  GO TO 50
               ENDIF
   40       CONTINUE
            IV=NVALUE(IPAR)+1
   50       IF(LGNEW) THEN
               NVALUE(IPAR)=NVALUE(IPAR)+1
               DO 60 J=NVALUE(IPAR)-1,IV,-1
                  VINTE(J+1)=VINTE(J)
   60          CONTINUE
               VINTE(IV)=IVAL
            ENDIF
         ENDIF
         IF(LGNEW) CALL hdf5_write_data(IPAPX,TRIM(RECNAM),
     1   VINTE(:NVALUE(IPAR)))
         DEALLOCATE(VINTE)
      ELSE IF(TYPE.EQ.'CHAINE') THEN
         ALLOCATE(VCHAR(NVALUE(IPAR)+1))
         IF(NVALUE(IPAR).EQ.0) THEN
            IV=1
            VCHAR(IV)=CVAL
            NVALUE(IPAR)=1
         ELSE
            CALL hdf5_get_shape(IPAPX,RECNAM,DIMS_APX)
            ILONG=DIMS_APX(1)
            IF(ILONG.GT.NVALUE(IPAR)) CALL XABORT('APXPAV: NVALUE OVER'
     1      //'FLOW(3).')
            CALL hdf5_read_data(IPAPX,RECNAM,VCHAR_V1)
            VCHAR(:ILONG)=VCHAR_V1(:ILONG)
            DEALLOCATE(VCHAR_V1)
            DO 70 I=1,NVALUE(IPAR)
               IF(CVAL.EQ.VCHAR(I))THEN
                  IV=I
                  LGNEW=.FALSE.
                  GO TO 80
               ENDIF
   70       CONTINUE
            IV=NVALUE(IPAR)+1
   80       IF(LGNEW) THEN
               NVALUE(IPAR)=NVALUE(IPAR)+1
               VCHAR(IV)=CVAL
            ENDIF
         ENDIF
         IF(LGNEW) CALL hdf5_write_data(IPAPX,TRIM(RECNAM),
     1   VCHAR(:NVALUE(IPAR)))
         DEALLOCATE(VCHAR)
      ENDIF
*
      IF(LGNEW) THEN
        CALL hdf5_write_data(IPAPX,"/paramdescrip/NVALUE",NVALUE(:NPAR))
      ENDIF
      DEALLOCATE(NVALUE)
      RETURN
      END
