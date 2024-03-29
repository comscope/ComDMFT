      subroutine solid
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer :: maxbv,maxbj,maxw2,ind
c     **************   GROUP **************************************
      call memory_solid(maxbv,maxbj,maxw2)
      call get_bz
      call group(maxbv,maxbj,maxw2)
!     ------- Real mesh checkings ------------------------------------
      call r_mesh_checking(mdiv)
      if(maswrk) write(iun,'(a35,3i5)')
     &  ' R_MESH_CHECKING: MDIV is set to be',(mdiv(ind),ind=1,3)
      call r_mesh_checking(nrdiv)
      if(maswrk) write(iun,'(a36,3i5)')
     &  ' R_MESH_CHECKING: NRDIV is set to be',(nrdiv(ind),ind=1,3)
!     -----------------------------------------------------------------     
      call bravay
      call symlmm
      if(maxnsymb.eq.0) maxnsymb=1
      call get_tau
      call getmts
      call radmsh
      maxnrad=maxval(nrad)
      call start_par_solid
      call getind
c     -------------------------------------------------------------------      
      call init_vh      
      call memory1_solid
      call scf_0
      call neighbor_atm
c     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      allocate(sovr(numplw))
      if(complex_ro) allocate(sovi(numplw))
      call warp
c     ------- Here for pair correlation function -------------------	
      call g_pair_print_prepare
c$$$      call print_xsf
c     ----------------------------------------------------------------
      call memory_flapw
      call rkmaxtest
      call solid_scf
      call memory_close

      end


      subroutine rkmaxtest
      use parallel_mod
      use manager_mod
      use solid_mod
      use atom_mod
      
      implicit none      
      integer :: isort
      double precision :: pi

      pi=4.0d0*datan2(1.0d0,1.0d0)
      
      if (maswrk) then
        write(iun,*) 'rkmax==l test'
        do isort=1, nsort
          write(iun,'(a, i5,3x,a, f12.6)')
     $      'isort', isort, 'rk/l',
     $      smt(isort)*cutoff*2.0d0*pi/par/lmb(isort)
          write(cun,'(a, i5,3x,a, f12.6)')
     $      'isort', isort, 'rk/l',
     $      smt(isort)*cutoff*2.0d0*pi/par/lmb(isort)          
        enddo
        call flush(cun)
      endif
      end
