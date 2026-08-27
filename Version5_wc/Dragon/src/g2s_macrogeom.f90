!
!-----------------------------------------------------------------------
!
!Purpose:
! Convert a MACRO-GEOM geometry format towards Salomon format.
!
!Copyright:
! Copyright (C) 2026 Ecole Polytechnique de Montreal.
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
!Author(s):
! A. Hebert
!
!Parameters: input
! impx              print flag (impx=0 for no print).
! ipgeo             pointer to the geometry LCM object (L_GEOM signature).
! dimTabCelluleBase allocated space for merg
! dimTabSegArc      allocated space for TabSegArc
!
!Parameters: output
! szSA              number of elements
! nbNode            number of nodes
! nbFlux            number of fluxes
! merg              merging indices per node
!
!-----------------------------------------------------------------------
!
subroutine g2s_macrogeom(impx,ipgeo,dimTabCelluleBase,dimTabSegArc,iszSA, &
& nbNode,nbFlux,merg)
    use constUtiles
    use segArc
    use GANLIB
    type(C_PTR) ipgeo
    integer,intent(in) :: impx,dimTabCelluleBase,dimTabSegArc
    integer,intent(out) :: iszSA,nbNode,nbFlux
    integer,intent(inout) :: merg(dimTabCelluleBase)

    ! local variables
    integer, parameter :: nstate=40
    integer :: istate(nstate)
    double precision :: sArea,point1(2),point2(2)
    type(t_segArc) :: sa

    ! allocatable arrays
    integer, dimension(:), allocatable :: itnode
    real, dimension(:,:), allocatable :: node

    ! ----------------------------------
    ! Recover the macro-geometry from ipgeo
    ! ----------------------------------
    if(impx>0) write(6,'(/40h g2s_macrogeom: process a macro-geometry)')
    call LCMGET(ipgeo,'STATE-VECTOR',istate)
    if(istate(1).ne.31) call XABORT('g2s_macrogeom: istate(1)=40 expected')
    call LCMSIX(ipgeo,'MACRO-GEOM',1)
      call LCMGET(ipgeo,'STATE-VECTOR',istate)
      idim=istate(1)
      nbNode=istate(2)
      if(nbNode>dimTabCelluleBase) call XABORT('g2s_macrogeom: tabSegArc overflow')
      nmerg=istate(4)
      allocate(node(idim,nbNode),itnode(nbNode))
      call LCMGET(ipgeo,'NODE',node)
      call LCMGET(ipgeo,'ITNODE',itnode)
      if(impx>0) write(6,200) istate(:5)
      do im=1,nbNode
        if(itnode(im)==1) then
          if(idim.lt.4) call XABORT('g2s_macrogeom: idim overflow(1).')
        else if(itnode(im)==2) then
          call XABORT('g2s_macrogeom: invalid itnode value.')
        else if(itnode(im)>=3) then
          if(idim.lt.2*itnode(im)) call XABORT('g2s_macrogeom: idim overflow(3).')
        endif
      enddo
      if(istate(5)==3) call XABORT('g2s_macrogeom: 3d not supported.')
    call LCMSIX(ipgeo,' ',2)

    ! ----------------------------------
    ! Accumulate surfacic elements
    ! ----------------------------------
    iszSA=0
    do im=1,nbNode
      if(itnode(im)==1) then
        point1(:2)=(/node(1,im), node(3,im)/) ; point2(:2)=(/node(2,im), node(3,im)/)
        call g2s_adelem(im,point1,point2,iszSA,tabSegArc,dimTabSegArc)
        point1(:2)=(/node(2,im), node(3,im)/) ; point2(:2)=(/node(2,im), node(4,im)/)
        call g2s_adelem(im,point1,point2,iszSA,tabSegArc,dimTabSegArc)
        point1(:2)=(/node(2,im), node(4,im)/) ; point2(:2)=(/node(1,im), node(4,im)/)
        call g2s_adelem(im,point1,point2,iszSA,tabSegArc,dimTabSegArc)
        point1(:2)=(/node(1,im), node(4,im)/) ; point2(:2)=(/node(1,im), node(3,im)/)
        call g2s_adelem(im,point1,point2,iszSA,tabSegArc,dimTabSegArc)
      else
        ! Shoelace formula
        sArea = 0.0d0
        do ic=1,2*itnode(im),2
          point1(1) = node(ic,im)
          point1(2) = node(ic+1,im)
          if(ic+1==2*itnode(im)) then
            point2(1) = node(1,im)
            point2(2) = node(2,im)
          else
            point2(1) = node(ic+2,im)
            point2(2) = node(ic+3,im)
          endif
          sArea=sArea+(point2(1)-point1(1))*(point2(2)+point1(2))
        enddo
        do ic=1,2*itnode(im),2
          point1(1) = node(ic,im)
          point1(2) = node(ic+1,im)
          if(ic+1==2*itnode(im)) then
            point2(1) = node(1,im)
            point2(2) = node(2,im)
          else
            point2(1) = node(ic+2,im)
            point2(2) = node(ic+3,im)
          endif
          if(sArea<0.0d0) then
            call g2s_adelem(im,point1,point2,iszSA,tabSegArc,dimTabSegArc)
          else if(sArea>0.0d0) then
            call g2s_adelem(im,point2,point1,iszSA,tabSegArc,dimTabSegArc)
          else
            call XABORT('g2s_macrogeom: Shoelace formula failure')
          endif
        enddo
      endif
    enddo

    !----
    !  Remix homogenized indices
    !----
    if(nmerg.gt.0) then
      call LCMSIX(ipgeo,'MACRO-GEOM',1)
        call LCMGET(ipgeo,'IREMIX',merg)
        nbFlux=0
        do i=1,nbNode
          nbFlux=max(nbFlux,merg(i))
        enddo
      call LCMSIX(ipgeo,' ',2)
    else
      nbFlux=nbNode
      do i=1,nbNode
        merg(i)=i
      enddo
    endif
    return
    !
    200 FORMAT(/22H macrogeometry options/1X,21(1H-)/ &
    & 7H IDIM  ,I8,34H   (FIRST DIMENSION OF ARRAY NODE)/ &
    & 7H NZS   ,I8,23H   (NUMBER OF MIXTURES)/ &
    & 7H ISYM  ,I8,39H   (=0/1/2/3: UNDEFINED/NONE/REFL/ROTA)/ &
    & 7H NMERGE,I8,33H   (=0/NUMBER OF MERGED MIXTURES)/ &
    & 7H IDIM  ,I8,16H   (=2/3: 2D/3D))
    !
    contains
      subroutine g2s_adelem(im,point1,point2,iszSA,tabSegArc,dimTabSegArc)
        integer im,dimTabSegArc
        double precision point1(2),point2(2)
        type(t_segArc) :: tabSegArc(dimTabSegArc)
        !
        do i=1,iszSA
          sa = tabSegArc(i)
          if(sa%typ /= 1) call XABORT('g2s_adelem: only segments are supported')
          if((abs(point2(1)-sa%x)<gSALeps).and.(abs(point2(2)-sa%y)<gSALeps).and. &
          & ((abs(point1(1)-sa%dx)<gSALeps).and.(abs(point1(2)-sa%dy)<gSALeps))) then
            tabSegArc(i)%noded=im
            return
          endif
        enddo
        iszSA=iszSA+1
        if(iszSA>dimTabSegArc) call XABORT('g2s_adelem: tabSegArc overflow')
        tabSegArc(iszSA)%typ=1
        tabSegArc(iszSA)%x=point1(1)
        tabSegArc(iszSA)%y=point1(2)
        tabSegArc(iszSA)%dx=point2(1)
        tabSegArc(iszSA)%dy=point2(2)
        tabSegArc(iszSA)%nodeg=im
        tabSegArc(iszSA)%noded=0
        tabSegArc(iszSA)%neutronicMixg=im
        tabSegArc(iszSA)%neutronicMixd=0
      end subroutine g2s_adelem
 end subroutine g2s_macrogeom
