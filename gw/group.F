      subroutine group(lmm,jmm,maxw2)
***************************************************************
*     Program for preparation of transformation matrices       *
*     for spherical harmonics according to the point group.    *
*     To find spherical harmonic which after applying the      *
*     declared group operation will reduce to y(l) multiply    *
*     vector y(m) to matrice u(m,m1). ip(j) points to the atom *
*     which goes to the atom #j after the same operation.      *
*     input:                                                  *
*     rbas0 - translational vectors in units of lat.par.    *
*     tau - positions of atoms in units of lat.par.         *
*     lmm - lmax i.e 0:s,1:p,2:d,3:f,etc                    *
*     natom  - number of atoms in primitive cell            *
*     ngroup - number of group elements founded             *
*     u - wigner's matrices                                 *
*     ip  (for comments see above)                          *
*     shift - shifts if group is non-symmorphic one         *
*     inv = 0 no  inversion                                 *
*     inv = 1 yes inversion, ifia=1                         *
*     inv = 2 yes inversion, ifia=-1                        *
*     note that:                                                 *
*     u    (l,g) = u(m'm,l,g) = u(iwig,g)               *
*     m'm                                              *
*     where iwig=l*(2*l-1)*(2*l+1)/3+(2*l+1)*(l+m)+l+m'+1        *
***************************************************************
      use manager_mod
      use parallel_mod
      use solid_mod
      use input_mod
      implicit none
      integer, intent(in) :: lmm,jmm,maxw2
      integer :: nt,ngen,ir,l,m,m1,iwig,ig,jatom,iatom,ia,ib,ic,ii,jj
      real*8 :: x(3),x1(3),t(3),sd(3),dshift,tempmat(3,3)
      integer, allocatable :: figen(:)
      real*8, allocatable :: fac(:),su1(:,:),gen(:,:),sgen(:,:),ih(:),
     &  ihh(:),sig(:)
      complex*16, allocatable :: su1_j(:,:),genj(:,:),ihj(:),ihhj(:),
     &  sigj(:)
      real*8 u1(3,3)

      integer :: ii_m,jj_m,i_index,jjval,jjdim,llval,iival
      real*8 :: jjval_dble

      integer, external :: iwig_j_ind

c$$$  character*10 :: ch_g
      inv=0
      allocate(fac(0:8*maxw2))
      call facc(fac,8*maxw2)
      nt=index(symgen,'_',.true.)
      allocate(genj(maxwig_j,48))
      allocate(gen(maxwig,48))
      allocate(sgen(3,144))
      allocate(figen(48))
      allocate(ih(maxwig))
      allocate(ihj(maxwig_j))
      allocate(su1(-lmm:lmm,-lmm:lmm))
      allocate(su1_j(-jmm:jmm,-jmm:jmm))
      if (symgen.ne."input") then
c     The symmetry operations were read in subroutine readini and the
c     symmop_rotation and symmop_shift data structures initialized
c     there.
        call parsgn(symgen,nt,gen,genj,sgen,figen,lmm,jmm,ngen,ih,ihj,
     &    su1,su1_j,fac)
      endif
      allocate(ihh(maxwig))
      allocate(sig(maxwig))
      allocate(ihhj(maxwig_j))
      allocate(sigj(maxwig_j))
      if (symgen.eq."input") then
c     
        call fractional2cartesian(ngroup,symmop_rotations,
     &    symmop_translations,rbas)

        ifia = 1
        u=0.0d0
        uj=0.0d0
        do ig = 1, ngroup
c$$$          if (maswrk) write(iun,*) 'sym:',ig
          call rotmat_mt_basis(symmop_rotations(:,:,ig),
     $      u(1,ig),uj(1,ig),lmm,jmm)
        enddo
        call dcopy(3*ngroup,symmop_translations,1,shift,1)
      else
        u=0.0d0
        uj=0.0d0
        call sgroup(gen,genj,sgen,ngen,u,uj,shift,ngroup,48,lmm,jmm,ih,
     &    ihj,ihh,ihhj,sig,sigj,su1,su1_j,figen)
      endif

      deallocate(symmop_rotations,symmop_translations)
      deallocate(genj,gen,sgen,figen,ih,ihj,su1,su1_j,ihh,sig,ihhj,sigj)
      symmorph=.true.

      do ir=1,ngroup
        call zone0_cart(shift(1,ir),rbas,shift(1,ir))
        dshift=shift(1,ir)**2+shift(2,ir)**2+shift(3,ir)**2
        if(dshift.gt.1.d-5) symmorph=.false.
        dshift=0.d0
        do l=0,lmm
          do m=-l,l
            do m1=-l,l
              iwig=l*(2*l-1)*(2*l+1)/3+(2*l+1)*(l+m1)+l+m+1
              if(m1.eq.m) then
                dshift=dshift+abs(u(iwig,ir)-(-1.d0)**l)
              else
                dshift=dshift+abs(u(iwig,ir))
              endif
            enddo
          enddo
        enddo
        if(dshift.lt.1.d-5) then
          inv=1+(1-ifia(ir))/2
          inv_num=ir
        endif
      enddo
c     ----- Check of symmetry operations relative to atomic positions --------------

c$$$      if (maswrk) then
c$$$        
c$$$        open(unit=112, file='rot_matrix_mt_basis' )
c$$$        do ig=1,ngroup
c$$$          
c$$$          tempmat(1,1)=u(10,ig)
c$$$          tempmat(2,1)=u(8,ig)
c$$$          tempmat(3,1)=u(9,ig)
c$$$          tempmat(1,2)=u(4,ig)
c$$$          tempmat(2,2)=u(2,ig)
c$$$          tempmat(3,2)=u(3,ig)
c$$$          tempmat(1,3)=u(7,ig)
c$$$          tempmat(2,3)=u(5,ig)
c$$$          tempmat(3,3)=u(6,ig)          
c$$$          write(112, '(3(f12.6,2x))')
c$$$     $      shift(1,ig),shift(2,ig), shift(3,ig)
c$$$          write(112, '(3(f12.6,2x))') tempmat(:,1)
c$$$          write(112, '(3(f12.6,2x))') tempmat(:,2)
c$$$          write(112, '(3(f12.6,2x))') tempmat(:,3)
c$$$          write(112, '(1000000000(f12.6,2x))') u(:,ig)
c$$$          write(112, '(1000000000(2f12.6,2x))') uj(:,ig)
c$$$        enddo
c$$$        close(112)
c$$$      endif


c$$$      if (maswrk) then

c$$$        write(iun,*) 'iwig_j'
c$$$        do llval=0, lmm
c$$$          do iival=-1,1, 2
c$$$            
c$$$            if ((llval .eq. 0) .and. (iival .eq. -1)) cycle
c$$$            write(iun,*) 'll, 2ii', llval, iival              
c$$$            i_index=(iival+1)/2+1
c$$$            jjval=llval*2+iival ! 2j
c$$$            jjval_dble=llval+iival/2.0d0 ! j
c$$$            jjdim=jjval+1       ! j dimension
c$$$            
c$$$            do ii_m=-jjval, jjval, 2
c$$$              write(iun,'(20000(i5))')
c$$$     $          (iwig_j(llval,iival,ii_m,jj_m),
c$$$     $          jj_m=-jjval, jjval, 2)
c$$$            enddo
c$$$          enddo
c$$$        enddo
        
        
c$$$        open(unit=112, file='rot_matrix_mt_basis_uj' )
c$$$        do ig=1,ngroup
c$$$          
c$$$          tempmat(1,1)=u(10,ig)
c$$$          tempmat(2,1)=u(8,ig)
c$$$          tempmat(3,1)=u(9,ig)
c$$$          tempmat(1,2)=u(4,ig)
c$$$          tempmat(2,2)=u(2,ig)
c$$$          tempmat(3,2)=u(3,ig)
c$$$          tempmat(1,3)=u(7,ig)
c$$$          tempmat(2,3)=u(5,ig)
c$$$          tempmat(3,3)=u(6,ig)
c$$$          write(112,*) 'sym: ', ig
c$$$          write(112, '(3(f12.6,2x))')
c$$$     $      shift(1,ig),shift(2,ig), shift(3,ig)
c$$$          write(112, '(3(f12.6,2x))') tempmat(:,1)
c$$$          write(112, '(3(f12.6,2x))') tempmat(:,2)
c$$$          write(112, '(3(f12.6,2x))') tempmat(:,3)
c$$$
c$$$          do llval=0, lmm
c$$$            do iival=-1,1, 2
c$$$
c$$$              if ((llval .eq. 0) .and. (iival .eq. -1)) cycle
c$$$              write(112,*) 'll, 2ii', llval, iival              
c$$$              i_index=(iival+1)/2+1
c$$$              jjval=llval*2+iival ! 2j
c$$$              jjval_dble=llval+iival/2.0d0 ! j
c$$$              jjdim=jjval+1     ! j dimension
c$$$
c$$$              do ii_m=-jjval, jjval, 2
c$$$                write(112,'(20000(i5))')
c$$$     $            (iwig_j_ind(llval,iival,ii_m,jj_m),
c$$$     $            jj_m=-jjval, jjval, 2)
c$$$              enddo              
c$$$
c$$$              do ii_m=-jjval, jjval, 2
c$$$                write(112,'(20000(2f20.12,2x))')
c$$$     $            (uj(iwig_j_ind(llval,iival,ii_m,jj_m),ig),
c$$$     $            jj_m=-jjval, jjval, 2)
c$$$              enddo
c$$$            enddo
c$$$          enddo
c$$$          write(112,*) 
c$$$              
c$$$        enddo
c$$$        close(112)
c$$$      endif      
      

      if (maswrk) open(unit=112, file='sym_operation_atom' )     
      if(.not.allocated(tshift)) allocate(tshift(3,natom,ngroup))
      do ig=1,ngroup

        if(maswrk) write(112,'(i5,a,i5)') ig, ' operation', ifia(ig)
        if(maswrk) write(112,'(a)') '   rotation'
        tempmat(1,1)=u(10,ig)
        tempmat(2,1)=u(8,ig)
        tempmat(3,1)=u(9,ig)
        tempmat(1,2)=u(4,ig)
        tempmat(2,2)=u(2,ig)
        tempmat(3,2)=u(3,ig)
        tempmat(1,3)=u(7,ig)
        tempmat(2,3)=u(5,ig)
        tempmat(3,3)=u(6,ig)
        do ii=1,3
          if(maswrk) write(112,'(3f12.6)') (tempmat(ii,jj), jj=1,3)
        enddo
        if(maswrk) write(112,'(a)') '   shift'
        if(maswrk) write(112,'(3f12.6)') shift(1:3,ig)

        if(maswrk) write(112,'(a)')
     $    'jatom, iatom, r_{iatom}-r_{R tau_jatom +tshift} cartesian '
     $    //'r_{iatom}-r_{R tau_jatom +tshift} latticevec'
c$$$  endif
        
        jat:do jatom=1,natom
        call rotate(tau(1,jatom),tau(2,jatom),tau(3,jatom),x(1),
     &    x(2),x(3),u(2,ig),2)  !!! = A*t^(-1)_a
        x=x+shift(:,ig)
        do iatom=1,natom
          x1=tau(:,iatom)-x
          do ia=-3,3
            do ib=-3,3
              do ic=-3,3
                t(1)=ia*rbas(1,1)+ib*rbas(1,2)+ic*rbas(1,3)
                t(2)=ia*rbas(2,1)+ib*rbas(2,2)+ic*rbas(2,3)
                t(3)=ia*rbas(3,1)+ib*rbas(3,2)+ic*rbas(3,3)
                sd=t-x1
                if(abs(sd(1))+abs(sd(2))+abs(sd(3)).lt.1.d-6) then
                  ip(iatom,ig)=jatom
                  tshift(:,iatom,ig)=t
                  if(maswrk) write(112,
     $              '(i5,i5,3f12.6,3i5)')
     $              jatom,iatom,
     $              t(1),t(2),t(3),ia,ib,ic
                  
                  cycle jat
                endif
              enddo
            enddo
          enddo
        enddo
        if(maswrk) then
          write(iun,*)'Not find jatom for ig =',ig,' iatom =',jatom
          write(iun,*)tau(1,jatom),' => ',x(1)
          write(iun,*)tau(2,jatom),' => ',x(2)
          write(iun,*)tau(3,jatom),' => ',x(3)
        endif
        call ending
      enddo jat

      if(maswrk) write(112,*)
      if(maswrk) write(112,*)
      if(maswrk) write(112,*)        
      enddo
      if (maswrk) close(112)
c     ----- Check of symmetry operations relative to external m. field ---------
      if(irel==2.and.magn==2) call check_b
c     --------------------------------------------------------------------------
      invers=.false.
      if(inv.eq.1) invers=.true.
      if(maswrk) then
        write(iun,*)'Number of elements in symmetry group is ',ngroup
        if(inv.eq.1) write(iun,*)'Inversional operation is presented'
        if(inv.eq.2) write(iun,*)
     $    'Inversional operation is not presented'
        if(symmorph) write(iun,*)'The group is symmorphic'
        if(.not.symmorph) write(iun,*)'The group is not symmorphic'
      endif
      deallocate(fac)
      real_things=.false.
      complex_ro=.true.
      if(invers) then
        if(maxval(abs(shift(:,inv_num)))<1.d-6) complex_ro=.false.
      endif
c     -------- Cartesian rotation matrices for optimization --------------
      do ig=1,ngroup
        u_opt(1,1,ig)=u(10,ig)
        u_opt(2,1,ig)=u(8,ig)
        u_opt(3,1,ig)=u(9,ig)
        u_opt(1,2,ig)=u(4,ig)
        u_opt(2,2,ig)=u(2,ig)
        u_opt(3,2,ig)=u(3,ig)
        u_opt(1,3,ig)=u(7,ig)
        u_opt(2,3,ig)=u(5,ig)
        u_opt(3,3,ig)=u(6,ig)
      enddo	
c     ------ Multiplication table -----------------------------------------	
      call group_table
c$$$  
c$$$  do ig=1,ngroup
c$$$  c$$$        if (maswrk) then
c$$$  if(maswrk) write(993,'(i5,a,i5)') ig, ' operation', ifia(ig)
c$$$  if(maswrk) write(993,'(a)') '   rotation'
c$$$  tempmat(1,1)=u(10,ig)
c$$$  tempmat(2,1)=u(8,ig)
c$$$  tempmat(3,1)=u(9,ig)
c$$$  tempmat(1,2)=u(4,ig)
c$$$  tempmat(2,2)=u(2,ig)
c$$$  tempmat(3,2)=u(3,ig)
c$$$  tempmat(1,3)=u(7,ig)
c$$$  tempmat(2,3)=u(5,ig)
c$$$  tempmat(3,3)=u(6,ig)
c$$$  do ii=1,3
c$$$  if(maswrk) write(993,'(3f12.6)') (tempmat(ii,jj), jj=1,3)
c$$$  enddo
c$$$  if(maswrk) write(993,'(a)') '   shift'
c$$$  if(maswrk) write(993,'(3f12.6)') shift(1:3,ig)
c$$$  enddo
c$$$  
c$$$  if (maswrk) then
c$$$  do ig=1,ngroup
c$$$  if (ig .lt. 10) then
c$$$  write(ch_g, '(i1)') ig
c$$$  else
c$$$  write(ch_g, '(i2)') ig
c$$$  endif
c$$$  open(unit=11 ,file='u_'//trim(adjustl(ch_g))//'.dat')
c$$$  open(unit=12 ,file='uj_'//trim(adjustl(ch_g))//'.dat')          
c$$$  write(11,'(100000f12.6)') u(:,ig)
c$$$  write(12,'(100000f12.6)') uj(:,ig)
c$$$  close(11)
c$$$  close(12)
c$$$  enddo
c$$$  endif
c$$$  
      call timel('**** GROUP finished ****************')					
      end
