!
!-----------------------------------------------------------------------
!
!Purpose:
! Generate a surfacic 2D geometry following the TDT specification.
!
!Copyright:
! Copyright (C) 2001 Ecole Polytechnique de Montreal.
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
!Author(s):
! G. Civario (CS-SI)
!
!Parameters: input/output
!  NENTRY : NUMBER OF LINKED LISTS AND FILES USED BY THE MODULE.
!  HENTRY : CHARACTER*12 NAME OF EACH LINKED LIST OR FILE.
!  IENTRY : =1 LINKED LIST; =2 XSM FILE; =3 SEQUENTIAL BINARY FILE;
!           =4 SEQUENTIAL ASCII FILE; =5 DIRECT ACCESS FILE.
!  JENTRY : =0 THE LINKED LIST OR FILE IS CREATED;
!           =1 THE LINKED LIST OR FILE IS OPEN FOR MODIFICATIONS;
!           =2 THE LINKED LIST OR FILE IS OPEN IN READ-ONLY MODE.
!  KENTRY : FILE UNIT NUMBER OR LINKED LIST ADDRESS.
!
!-----------------------------------------------------------------------
!
subroutine G2S(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
  use GANLIB
  use SALGET_FUNS_MOD
  use celluleBase
  use cellulePlaced
  use boundCond
  use ptNodes
  use pretraitement
  use derivedPSPLOT
  use track
  use segArc
  use generTabSegArc
  use generSAL

  implicit none

  integer NENTRY
  integer IENTRY,JENTRY
  type(c_ptr) KENTRY
  character*12 HENTRY
  dimension IENTRY(*),JENTRY(*),KENTRY(*),HENTRY(*)

  integer,parameter :: dimTabCelluleBase = 20000
  integer,parameter :: dimTabSegArc = 100000
  integer,parameter :: nstate=40

  type(c_ptr)  :: ipGeo,ipGeo_1
  integer      :: sizeB,sizeP,sizeSA,nbNode,nbCLP,nbFlux,nbMacro,ipSal,ipPs,ipAl, &
                  ipZa,indic,nitma,impx,drawMix,istate(nstate)
  character(len=12) :: text12
  logical      :: lmacro
  real,dimension(2) :: zoomx,zoomy
  integer,allocatable,dimension(:) :: gig,merg,imacro
  integer,dimension(10) :: datain
  real :: flott
  double precision :: dflott
  integer      :: lgMaxGig=0
  character(len=12) :: hsign
  !
  ipAl=-1 ! no Alamos file read
  ipZa=-1 ! no PropertyMap file read
  ipGeo_1=c_null_ptr ! no geometry read
  if ((NENTRY == 2).and.(IENTRY(2) == 4)) then
     !generating ps file from Salomon file
     ipPs   = FILUNIT(KENTRY(1)) ! output psfile
     g_psp_isEpsFile = (index(HENTRY(1),'.eps')/=0) !is it an eps file ?
     ipSal  = FILUNIT(KENTRY(2)) ! input Salomon/Alamos file
     ! check argument types and permissions
     if ((IENTRY(1) /= 4) .or. (JENTRY(1) /= 0)) &
          call XABORT('g2s_g2s: a new file was expected for the postscript file')
     if ((IENTRY(2) /= 4) .or. (JENTRY(2) /= 2)) &
          call XABORT('g2s_g2s: expecting Salomon file in read-only mode')
  else if ((NENTRY == 3).and.(IENTRY(3) == 4)) then
     !generating Salomon and ps files from Alamos file
     ipSal  = FILUNIT(KENTRY(1)) ! output Salomon file
     ipPs   = FILUNIT(KENTRY(2)) ! output psfile
     ipAl   = FILUNIT(KENTRY(3)) ! input Alamos file
     g_psp_isEpsFile = (index(HENTRY(2),'.eps')/=0) !is it an eps file ?
     ! check argument types and permissions
     if ((IENTRY(1) /= 4) .or. (JENTRY(1) /= 0)) &
          call XABORT('g2s_g2s: a new file was expected for the Salomon file')
     if ((IENTRY(2) /= 4) .or. (JENTRY(2) /= 0)) &
          call XABORT('g2s_g2s: a new file was expected for the postscript file')
     if ((IENTRY(3) /= 4) .or. (JENTRY(3) /= 2)) &
          call XABORT('g2s_g2s: expecting Alamos file in read-only mode')
  else if ((NENTRY == 4).and.(IENTRY(3) == 4).and.(IENTRY(4) == 4)) then
     !generating Salomon and ps files from Alamos and PropertyMap file
     ipSal  = FILUNIT(KENTRY(1)) ! output Salomon file
     ipPs   = FILUNIT(KENTRY(2)) ! output psfile
     ipAl   = FILUNIT(KENTRY(3)) ! input Alamos file
     ipZa   = FILUNIT(KENTRY(4)) ! input PropertyMap file
     g_psp_isEpsFile = (index(HENTRY(2),'.eps')/=0) !is it an eps file ?
     ! check argument types and permissions
     if ((IENTRY(1) /= 4) .or. (JENTRY(1) /= 0)) &
          call XABORT('g2s_g2s: a new file was expected for the Salomon file')
     if ((IENTRY(2) /= 4) .or. (JENTRY(2) /= 0)) &
          call XABORT('g2s_g2s: a new file was expected for the postscript file')
     if ((IENTRY(3) /= 4) .or. (JENTRY(3) /= 2)) &
          call XABORT('g2s_g2s: expecting Alamos file in read-only mode')
     if ((IENTRY(4) /= 4) .or. (JENTRY(4) /= 2)) &
          call XABORT('g2s_g2s: expecting PropertyMap file in read-only mode')
  else if ((NENTRY == 2).and.(IENTRY(2) <= 2)) then
     !generating Salomon file from LCM geometry
     ipPs   = -1        ! no postscript file
     ipSal  = FILUNIT(KENTRY(1)) ! output Salomon file
     ipGeo_1= KENTRY(2) ! geometry read
     ! check argument types and permissions
     if ((IENTRY(1) /= 4) .or. (JENTRY(1) /= 0)) &
          call XABORT('g2s_g2s: a new ASCII file was expected for writing geometry')
     if ((IENTRY(2) > 2) .or. (JENTRY(2) /= 2)) &
          call XABORT('g2s_g2s: expecting LCM geometry in read-only mode(1)')
  else if ((NENTRY == 3).and.(IENTRY(3) <= 2)) then
     !generating Salomon file and ps file from LCM geometry
     ipSal  = FILUNIT(KENTRY(1)) ! output Salomon file
     ipPs   = FILUNIT(KENTRY(2)) ! output psfile
     g_psp_isEpsFile = (index(HENTRY(2),'.eps')/=0) !is it an eps file ?
     ipGeo_1= KENTRY(3) ! geometry read
     ! check argument types and permissions
     if ((IENTRY(1) /= 4) .or. (JENTRY(1) /= 0)) &
          call XABORT('g2s_g2s: a new file was expected for writing geometry')
     if ((IENTRY(2) /= 4) .or. (JENTRY(2) /= 0)) &
          call XABORT('g2s_g2s: a new file was expected for the postscript file')
     if ((IENTRY(3) > 2) .or. (JENTRY(3) /= 2)) &
          call XABORT('g2s_g2s: expecting LCM geometry in read-only mode(2)')
  else
     call XABORT('g2s_g2s: you must provide 2, 3 or 4 arguments')
  end if

  impx=1
  drawMix = 0
  zoomx = (/ 0.0, 1.0 /)
  zoomy = (/ 0.0, 1.0 /)
  typgeo=0
  lmacro=.false.
  10 call REDGET(indic,nitma,flott,text12,dflott)
  if (indic == 10) go to 20
  if (indic /= 3) call XABORT('g2s_g2s: character data expected.')
  if (text12 == 'EDIT') then
    ! read the print index.
    call REDGET(indic,impx,flott,text12,dflott)
    if (indic /= 1) call XABORT('g2s_g2s: integer data expected(1).')
  else if (text12 == 'DRAWNOD') then
     drawMix=1
  else if (text12 == 'DRAWMIX') then
     drawMix=2
  else if (text12 == 'DRAWELEM') then
     drawMix=3
  else if (text12 == 'DRAWMACRO') then
     drawMix=4
  else if (text12 == 'ZOOMX') then
    call REDGET(indic,nitma,zoomx(1),text12,dflott)
    if (indic /= 2) call XABORT('g2s_g2s: real data expected(1).')
    call REDGET(indic,nitma,zoomx(2),text12,dflott)
    if (indic /= 2) call XABORT('g2s_g2s: real data expected(2).')
    if ((zoomx(1).lt.0.0).or.(zoomx(2).le.zoomx(1)).or.(zoomx(2).gt.1.0)) then
      call XABORT('g2s_g2s: invalid zoom factors in x.')
    endif
  else if (text12 == 'ZOOMY') then
    call REDGET(indic,nitma,zoomy(1),text12,dflott)
    if (indic /= 2) call XABORT('g2s_g2s: real data expected(3).')
    call REDGET(indic,nitma,zoomy(2),text12,dflott)
    if (indic /= 2) call XABORT('g2s_g2s: real data expected(4).')
    if ((zoomy(1).lt.0.0).or.(zoomy(2).le.zoomy(1)).or.(zoomy(2).gt.1.0)) then
      call XABORT('g2s_g2s: invalid zoom factors in y.')
    endif
  else if (text12 == 'ALAMOS') then
    if (NENTRY == 2) call XABORT('g2s_g2s: three entries required.')
    if (ipAl == -1) call XABORT('g2s_g2s: no RHS Salomon file.')
    call REDGET(indic,typgeo,flott,text12,dflott)
    if (indic /= 1) call XABORT('g2s_g2s: integer data expected(2).')
  else if (text12 == 'MACRO') then
    lmacro=.true.
  else if (text12 == ';') then
     go to 20
  else
     call XABORT('g2s_g2s: '//text12//' is an invalid keyword.')
  end if
  go to 10

  !conversion of Alamos file into a Salomon file
  20 if (ipAl /= -1) call g2s_convert(impx,ipAl,ipZa,ipSal)

  sizeB = 0   !cellules de base
  sizeP = 0   !cellules placees
  sizeSA = 0  !elements geometriques
  if (c_associated(ipGeo_1)) then
     !initialisation des differents tableaux
     call initializeData(dimTabCelluleBase,dimTabSegArc)

     ! check if ipGeo_1 contains a macro-geometry
     call LCMGTC(ipGeo_1,'SIGNATURE',12,hsign)
     if(hsign /='L_GEOM') call XABORT('g2s_g2s: L_GEOM object expected')
     call LCMGET(ipGeo_1,'STATE-VECTOR',istate)
     if(istate(1)==31) THEN
       ! process a macro-geometry
       ipGeo=c_null_ptr
       allocate(merg(dimTabCelluleBase),imacro(dimTabCelluleBase),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("g2s_g2s: allocation pb(1)")
       call g2s_macrogeom(impx,ipGeo_1,dimTabCelluleBase,dimTabSegArc,sizeSA, &
       & nbNode,nbFlux,merg)
       imacro(:nbFlux)=1
       nbCLP=0
     else
       ! copy the input geometric object
       call lcmop(ipGeo,'geom_copy',0,1,0)
       call lcmequ(ipGeo_1,ipGeo)

       !unfold the geometry
       call g2s_unfold(ipGeo,impx)

       !pretraitement des donnees lues (remplace la partie python)
       !+completion des cellules de base et remplissage du tableau
       !des cellules placees
       call prepareData(ipGeo,sizeB,sizeP,lgMaxGig)
       if (impx > 0) write(*,*) 'fin   : prepareData lgMaxGig=',lgMaxGig

       !en sortie, toutes les cellules de base ont tous leurs
       !champs remplis, et le tableau des cellules placees est pret

       !eclatement des cellules
       call splitCells(sizeP,sizeSA)
 
       !creation de nouveaux segements aux interfaces des cellules
       !et elimination des doublons
       call addSegsAndClean(sizeSA)

       !prise en compte des conditions aux limites
       call appliBoundariConditions(ipGeo,sizeSA,nbCLP)

       !calcul des nodes delimites par les elements
       allocate(merg(dimTabCelluleBase),imacro(dimTabCelluleBase),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("g2s_g2s: allocation pb(2)")
       call createNodes(sizeSA,dimTabCelluleBase,lmacro,nbNode,merg,imacro)
       if (sizeSA > dimTabSegArc) call XABORT('g2s_g2s: sizeSA overflow')
     
       !calcul des arrays gig et merg
       allocate(gig(nbNode*lgMaxGig),stat=alloc_ok)
       if (alloc_ok /= 0) call XABORT("g2s_g2s: allocation pb(3)")
       call generateTrack(sizeP,sizeSA,nbNode,lgMaxGig,gig,merg)
       nbFlux=maxval(merg(:nbNode))
       nbMacro=maxval(imacro(:nbNode))
       deallocate(gig)
     endif
  else
     if (JENTRY(NENTRY) == 0) call XABORT('g2s_g2s: a RHS ascii file is expected')
     !initialisation de TabSegArc
     call SALGET(datain,4,ipSal,0,'dimensions for geometry')
     nbNode=datain(3)
     sizeSA=datain(4)
     rewind(ipSal)
     allocate(TabSegArc(sizeSA),stat=alloc_ok)
     if (alloc_ok /= 0) call XABORT("g2s_g2s: generateTabSegArc => allocation pb")
     call initializebCData()  
     allocate(merg(nbNode),imacro(nbNode),stat=alloc_ok)
     if (alloc_ok /= 0) call XABORT("g2s_g2s: allocation pb(4)")
     call generateTabSegArc(ipSal,sizeSA,nbNode,nbCLP,nbFlux,nbMacro,merg,imacro,impx)
  endif

  !impression des segArc charges
  if (ipPs /= -1) call drawSegArc(ipPs,sizeSA,drawMix,imacro,zoomx,zoomy)

  if (c_associated(ipGeo_1)) then
     !creation du fichier de commande SAL
     call generateSALFile(ipSal,sizeSA,nbNode,nbCLP,nbFlux,nbMacro,merg,imacro)
     if (c_associated(ipGeo)) call LCMCL(ipGeo,2)
  endif
  deallocate(imacro,merg)

  write(6,*) "  At end of G2S:"
  write(6,*) "    ",sizeSA,"segs or arcs"
  write(6,*) "    ",nbNode,"nodes"
  write(6,*) "    ",nbFlux,"fluxes"
  write(6,*) "    ",nbMacro,"macros"
  write(6,*) "    ",nbCLP,"boundary conditions other than default"

  !liberation de la memoire allouee
  call destroyData(sizeB,sizeP)
end subroutine G2S

subroutine initializeData(dimTabCelluleBase,dimTabSegArc)
  use celluleBase
  use cellulePlaced
  use boundCond
  use segArc
  implicit none
  integer,intent(in) :: dimTabCelluleBase,dimTabSegArc

  call initializeTabCelluleBase(dimTabCelluleBase)
  call initializeTabCellulePlaced()
  allocate(tabSegArc(dimTabSegArc))
  call initializebCData()  
end subroutine initializeData

subroutine destroyData(szB,szP)
  use celluleBase
  use cellulePlaced
  use boundCond
  use segArc
  implicit none
  integer,intent(in) :: szB,szP

  call destroyTabCelluleBase(szB)
  deallocate(tabSegArc)
  call destroyTabCellulePlaced(szP)
  call destroybCData()
end subroutine destroyData
