*DECK SNTSFH
      SUBROUTINE SNTSFH (IMPX,ITYPE,NHEX,LZ,MCELL,ISPLH,MAT,LOZSWP,
     >   COORDMAP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Output arrays for lozenge sweep order (direction-dependent) and
* coordinate map, both needed for resolution of the discrete ordinates
* transport equation in hexagonal geometry.
*
*Copyright:
* Copyright (C) 2025 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. A. Calloo
*
*Parameters: input
* IMPX     print parameter.
* ITYPE    geometry type (8:hexagonal 2D, 9:hexagonal 3D).
* NHEX     number of hexagons (for 3D, in one plane only).
* LZ       number of mesh elements in z-axis (including split).
* MCELL    number of macrocells to use along z-axis.
* ISPLH    mesh-splitting in 3*ISPLH**2 lozenges per hexagon.
* MAT      mixture index assigned to each element.
*
*Parameters: local
* NRINGS   number of hexagonal rings in the domain, assuming the centre 
*          hexagon counts as 1 ring.
*
*Parameters: output
* LOZSWP   lozenge sweep order depending on direction.
* COORDMAP coordinate map: mapping hexagon from DRAGON geometry indices
*          to the axial coordinate system, using p, r, s axes. The s
*          axis is redundant, which means that using p and r axes
*          effectively maps the hexagon geometry to a 2D map. Refer to 
*          the redblobgames blog for more information. 
*
*Comments:
* The lozenge under consideration is given by the position within the
* matrix. See user manual and/or data manual and/or thesis
*                                 _____
*                                /   / \
*                               / B /   \
*                         ,----(----  A  )----.
*                        /      \ C \   /      \
*                       /        \___\_/        \
*                       \   4    /     \   2    /
*                        \      /       \      /
*                         )----(   1     )----(
*                        /      \       /      \
*                       /        \_____/        \
*                       \   5    /     \   7    /
*                        \      /       \      /
*                         `----(   6     )----'
*                               \       /
*                                \_____/
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IMPX,ITYPE,NHEX,LZ,MCELL,ISPLH,MAT(ISPLH**2,3,NHEX,LZ),
     >   COORDMAP(3,NHEX)
*----
*  LOCAL VARIABLES
*----
      INTEGER, DIMENSION(3,6) :: LOZSWP,MAPCODE
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: TMPMAT
      INTEGER, ALLOCATABLE, DIMENSION(:) :: TASKSPERWAVE
*----
*  LOZENGE SWEEP ORDERING WITHIN HEXAGONS - DIRECTION DEPENDENT
*----
      LOZSWP = RESHAPE((/ 3, 2, 1, 3, 1, 2, 1, 3, 2, 1, 2, 3, 2, 1,
     >   3, 2, 3, 1 /), SHAPE(LOZSWP))
*----
*  CREATE COORDIDATE MAP FROM DRAGON INDEX TO AXIAL COORDINATES
*----
      NRINGS=INT((SQRT(  REAL((4*NHEX-1)/3)  )+1.)/2.)
      IF(NRINGS.EQ.1) CALL XABORT('NOT IMPLEMENTED FOR SINGLE HEX YET.')
      IF(NHEX.NE.1+3*NRINGS*(NRINGS-1)) CALL XABORT('SNTSFH: INVALID '
     1 //'VALUE OF NHEX(1).')
*
      MAPCODE = RESHAPE((/ -1, 0, 1, -1, 1, 0, 0, 1, -1, 1, 0, -1, 1, 
     >   -1, 0, 0, -1, 1 /), SHAPE(MAPCODE))
*
      ! It should be noted that the algorithm below effectively 
      ! reverses the y-axis. However, this should be of no consequence
      ! whatsoever as it would akin to the user defining the domain
      ! somewhat differently in the geometry. Calculations and results
      ! should be unaffected. 
      IHEX_DOM=1
      DO IRING=1,NRINGS
         ! Initialise first 'ring', i.e., centre hexagon
         IF(IRING.EQ.1) THEN
            ITMP1 = NRINGS
            ITMP2 = NRINGS
            ITMP3 = -2*(NRINGS)
            COORDMAP(1,IHEX_DOM)=ITMP1
            COORDMAP(2,IHEX_DOM)=ITMP2
            COORDMAP(3,IHEX_DOM)=ITMP3

            IHEX_DOM = IHEX_DOM+1
            ! Ignore rest of this loop and move on to next ring
            CYCLE
         ENDIF

         ! Find coordinates for hexagon when moving from ring n-1 to n
         ITMP1 = ITMP1+1
         ITMP2 = ITMP2-1
         ITMP3 = ITMP3+0
         COORDMAP(1,IHEX_DOM)=ITMP1
         COORDMAP(2,IHEX_DOM)=ITMP2
         COORDMAP(3,IHEX_DOM)=ITMP3
         IHEX_DOM = IHEX_DOM+1

         ! 'Sweep' each of the 3 axes of the hexagonal plane and their 
         ! negative directions
         DO IND=1,6
            ! Step through each hexagon per each axis
            DO IHEX=1,IRING-1
               ITMP1 = ITMP1+MAPCODE(1,IND)
               ITMP2 = ITMP2+MAPCODE(2,IND)
               ITMP3 = ITMP3+MAPCODE(3,IND)
               ! Store each of the coordinates except the last hexagon 
               ! in the last direction. This is because we already 
               ! computed that hexagon when moving from ring n-1 to n
               IF((IND.EQ.6).AND.(IHEX.EQ.IRING-1))THEN
                  CONTINUE
               ELSE
                  COORDMAP(1,IHEX_DOM)=ITMP1
                  COORDMAP(2,IHEX_DOM)=ITMP2
                  COORDMAP(3,IHEX_DOM)=ITMP3
                  IHEX_DOM = IHEX_DOM+1
               ENDIF
            ENDDO ! ihex
         ENDDO ! ind
      ENDDO ! iring

*----
*  COMPUTE NUMBER OF CONCURRENT HEXAGONS PER WAVEFRONT FOR PRINTING
*  PURPOSES ONLY
*----
      IF(MCELL > 0)THEN

      ! Build material array in axial coordinates
      NCOLS=2*NRINGS -1
      ALLOCATE(TMPMAT(ISPLH**2,3,NCOLS,NCOLS,LZ))
      TMPMAT(:,:,:,:,:)=-1
      DO IZ=1,LZ
         DO IHEX_XY=1,NHEX
            TMPMAT(:,:,COORDMAP(1,IHEX_XY),COORDMAP(2,IHEX_XY),IZ) = 
     >         MAT(:,:,IHEX_XY,IZ)
         ENDDO
      ENDDO

      ! Build TasksPerWave array
      IF(ITYPE==8)THEN
         ! 2D Hexagonal
         NWAVES=NCOLS+NCOLS-1
         ALLOCATE(TASKSPERWAVE(NWAVES))
         TASKSPERWAVE(:)=0

         DO IWAVE=1,NWAVES
            ICOUNT = 0
            DO J=MAX(1,IWAVE-NCOLS+1),MIN(NCOLS,IWAVE)
               I=IWAVE-J+1
               I=NCOLS+1-I
               IF((I.GT.NCOLS).OR.(I.LT.1)) CYCLE
               IF((J.GT.NCOLS).OR.(J.LT.1)) CYCLE
               ! If within corners of Cartesian axial coordinate map
               ! (where there are no hexagons), skip loop
               IF(TMPMAT(1,1,I,J,1).EQ.-1) CYCLE
               ICOUNT = ICOUNT + 1
            ENDDO
            TASKSPERWAVE(IWAVE) = ICOUNT
         ENDDO
      ELSEIF(ITYPE==9)THEN
         ! 3D Hexagonal
         MCELLZ = MCELL
         NWAVES=NCOLS+NCOLS+MCELLZ-2
         ALLOCATE(TASKSPERWAVE(NWAVES))
         TASKSPERWAVE(:)=0

         DO IWAVE=1,NWAVES
            ICOUNT = 0
            J_STT=MAX(1,IWAVE-NCOLS-MCELLZ+2)
            J_END=MIN(NCOLS,IWAVE)
            DO J_MC=J_STT,J_END
               J=J_MC
               I_STT=MAX(1,IWAVE-J_MC-MCELLZ+2)
               I_END=MIN(NCOLS,IWAVE-J_MC+1)
               DO I_MC=I_STT,I_END
                  I=I_MC
                  I=NCOLS+1-I
                  ! If within corners of Cartesian axial coordinate map
                  ! (where there are no hexagons), skip loop
                  IF(TMPMAT(1,1,I,J,1).EQ.-1) CYCLE
                  K_MC=IWAVE-I_MC-J_MC+2
                  ICOUNT = ICOUNT + 1
               ENDDO
            ENDDO
            TASKSPERWAVE(IWAVE) = ICOUNT
         ENDDO
      ENDIF

      DEALLOCATE(TMPMAT)

      ENDIF

*----
*  PRINT A FEW GEOMETRY CHARACTERISTICS
*----
      IF(IMPX.GT.2) THEN
         WRITE(*, 100)
         WRITE(*, 101) NCOLS
         WRITE(*, 102) NRINGS
         WRITE(*, 103)
         DO I=1,6
            WRITE(*,104) I, LOZSWP(:,I)
         ENDDO
         IF(MCELL > 0)THEN
            WRITE(*, 105) NWAVES
            WRITE(*, 106)
            DO I = 1, NWAVES
               WRITE(*, 107) TASKSPERWAVE(I)
            END DO
            DEALLOCATE(TASKSPERWAVE)
         ENDIF
      ENDIF
      IF(IMPX.GT.4) THEN
         WRITE(*, 109)
         WRITE(*, 110)
         DO I = 1, NHEX
            WRITE(*, 111) I, COORDMAP(:, I)
         END DO
      ENDIF

      RETURN
  100 FORMAT (' ')
  101 FORMAT ('NCOLS    =', I4)
  102 FORMAT ('NRINGS   =', I4)
  103 FORMAT ('LOZENGE SWEEP ORDER')
  104 FORMAT ('IND_XY:', I4, ' LOZ. ORDER:', 3I4)
  105 FORMAT ('NWAVES   =', I4)
  106 FORMAT ('TASKS PER WAVE')
  107 FORMAT (I4) 
  109 FORMAT (' ')
  110 FORMAT ('COORDINATE MAP IS GIVEN BELOW:')
  111 FORMAT ('DRAGON IND:', I4, ' AXIAL COORD:', 3I4)
      END
