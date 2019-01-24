      subroutine green_0(correlation)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use lda_dmft_mod
      implicit none
      logical, intent(in) :: correlation
c$$$      character*1 :: intt
      integer :: k,ispin,len,k0,i,j
      real*8 :: de,green0
      real*8 :: maxerror
      integer :: fail_ind,ibnd,jbnd,aa,bb,cc, intloc(2)
      complex*16, allocatable :: umat(:,:),proj_mat(:,:),
     $  tmp1(:,:),tmp2(:,:),deltarhomat1(:,:),deltarhomat2(:,:)
c$$$      character*300 :: ch
c     ------------- Forming G_FULL(0-) --> G_FULL_0 ----------------------


      
      if(.not.correlation) g_full_00=0.d0
      if(correlation) g_full_0=(0.d0,0.d0)
      if(me3_tau==0) then
        do ispin=1,nspin
          do k=1,ndim3_k(me3_k+1)
            k0=n3_mpi_k(me3_k+1)+k
c$$$            
c$$$            if (k0 .lt. 10) then
c$$$              write(ch, '(a, i1, a)') "check_", k0, "_"
c$$$            elseif (k0 .lt. 100) then
c$$$              write(ch, '(a, i2, a)') "check_", k0, "_"
c$$$            elseif (k0 .lt. 1000) then
c$$$              write(ch, '(a, i3, a)') "check_", k0, "_"
c$$$            elseif (k0 .lt. 10000) then
c$$$              write(ch, '(a, i4, a)') "check_", k0, "_"
c$$$            endif
      
      
            if(.not.correlation) then
              g_full_00(:,k,ispin)=0.d0
              do i=1,n_bnd(k0,ispin)
                de=e_bnd(i,k0,ispin)-chem_pot
                g_full_00(i,k,ispin)=-green0(de,tau_mesh(n_tau))
              enddo
              if (is_lda_dmft) then
                allocate(umat(nbndf,nbndf))
                allocate(proj_mat(nbndf,n_dmft_wan_max))
                allocate(tmp1(nbndf,n_dmft_wan_max))
                allocate(tmp2(nbndf,nbndf))
                allocate(deltarhomat1(nbndf,nbndf))
                allocate(deltarhomat2(nbndf,nbndf))                                
                
                umat=0.0d0
                call dft_basis_rotmat_k(k,k0,ispin,umat)
c$$$  if (dmft_rho_cal_count .eq. 1) then
c$$$  call is_identity_dcmplx(umat,n_bnd(k0,ispin),nbndf)
c$$$  endif

c$$$                write(179,*) 'umat 0'
c$$$                do i=1,n_bnd(k0,ispin)                
c$$$                  write(179,'(10000(2f12.6, 2x))')
c$$$     $              (umat(i,j), j=1,n_bnd(k0,ispin))
c$$$                enddo
                
                call force_mat_unitary(n_bnd(k0,ispin),nbndf,umat)

c$$$                write(179,*) 'umat'
c$$$                do i=1,n_bnd(k0,ispin)                
c$$$                  write(179,'(10000(2f12.6, 2x))')
c$$$     $              (umat(i,j), j=1,n_bnd(k0,ispin))
c$$$                enddo                
                
                call check_if_unitary(n_bnd(k0,ispin),nbndf,umat,
     $            fail_ind, maxerror)
c$$$  write(179,'(a,3i5,f20.12)') 'unitary_check max error',
c$$$  $            ispin,k,k0,maxerror                
                if (fail_ind .ne. 0) then
                  write(iun,*) 'unitary_check fail', maxerror
                  call ending
                endif
                
                do ibnd=1,nbndf
                  do jbnd=1, nbndf
                    g_full_0(ibnd,jbnd,k,ispin)=0.0d0
                  enddo
                enddo
!     rho
                do i=1,n_bnd(k0,ispin)
                  de=e_bnd(i,k0,ispin)-chem_pot
                  g_full_0(i,i,k,ispin)=-green0(de,tau_mesh(n_tau))
c$$$                  if (maswrk) then
c$$$                    write(iun,'(a,3x,i5, f12.6, 5x, 10000(2f12.6, 2x))')
c$$$     $                intt, k0,de*13.605,g_full_0(i,i,k,ispin)
c$$$                  endif
                enddo

c$$$                open(unit=189,file=trim(ch)//'rho'//intt//'_0.dat')
c$$$
c$$$                do i=1,n_bnd(k0,ispin)                
c$$$                  write(189,'(10000(2f12.6,2x))')
c$$$     $              (g_full_0(i,j,k,ispin),j=1,n_bnd(k0,ispin))
c$$$                enddo
c$$$                close(189)
!     add rho_dmft

                deltarhomat1=0.0d0
                deltarhomat2=0.0d0                                
                
                proj_mat=0.0d0
                do aa=1,n_bnd(k0,ispin)
                  do bb=1, n_dmft_bnd(k0)
                    do cc=1, n_dmft_wan(k0)
                      proj_mat(aa,cc)=proj_mat(aa,cc)
     $                  +umat(aa,dmft_bnd_index(bb,k0))
     $                  *disentangle_mat(bb,cc,k0)
                    enddo
                  enddo
                enddo
                
                tmp1=0.0d0
                call zgemm('n','n',n_bnd(k0,ispin),n_dmft_wan(k0),
     $            n_dmft_wan(k0),(1.d0,0.d0),proj_mat,nbndf,
     $            dmft_den_mat(1,1,k0),
     $            n_dmft_wan_max,(0.d0,0.d0),tmp1,nbndf)
                call zgemm('n','c',n_bnd(k0,ispin),n_bnd(k0,ispin),
     $            n_dmft_wan(k0),(1.d0,0.d0),tmp1,nbndf,proj_mat,nbndf,
     &            (0.d0,0.d0),deltarhomat1,nbndf)


c$$$                open(unit=189,file=trim(ch)//'rho'//intt//'_1.dat')                
c$$$
c$$$                do i=1,n_bnd(k0,ispin)
c$$$                  write(189,'(10000(2f12.6,2x))')                  
c$$$     $              (deltarhomat1(i,j),j=1,n_bnd(k0,ispin))
c$$$                enddo
c$$$                close(189)
                
                tmp2=0.0d0
                
                call zgemm('n','c',n_bnd(k0,ispin),n_bnd(k0,ispin),
     $            n_dmft_wan(k0),(1.d0,0.d0),proj_mat,nbndf,
     $            proj_mat,nbndf,
     &            (0.d0,0.d0),tmp2,nbndf)
                
                call check_hermitian(tmp2,n_bnd(k0,ispin),nbndf,
     $            'subtracting_double_counted_rho')
                
                do aa=1,n_bnd(k0,ispin)
                  de=e_bnd(aa,k0,ispin)-chem_pot                                          
                  do bb=1,n_bnd(k0,ispin)                  
                    do cc=1,n_bnd(k0,ispin)
                      deltarhomat2(bb,cc)
     $                  =deltarhomat2(bb,cc)
     $                  -tmp2(bb,aa)*green0(de,tau_mesh(n_tau))
     $                  *tmp2(aa,cc)
                    enddo
                  enddo
                enddo

c$$$                open(unit=189,file=trim(ch)//'rho'//intt//'_2.dat')                
c$$$
c$$$
c$$$                do i=1,n_bnd(k0,ispin)
c$$$                  write(189,'(10000(2f12.6,2x))')                                    
c$$$     $              (deltarhomat2(i,j),j=1,n_bnd(k0,ispin))
c$$$                enddo                                
c$$$                close(189)

                do bb=1,n_bnd(k0,ispin)                  
                  do cc=1,n_bnd(k0,ispin)
                    g_full_0(bb,cc,k,ispin)
     $                =g_full_0(bb,cc,k,ispin)
     $                +deltarhomat1(bb,cc)
     $                -deltarhomat2(bb,cc)                    
                  enddo
                enddo
c$$$                if (intt .eq. '0') then
c$$$                  intloc=maxloc(cdabs(deltarhomat1-deltarhomat2))                  
c$$$                  write(*,'(a,i5,f12.6, 2i5,2f12.6)') 'deltarho k0',k0,
c$$$     $              maxval(cdabs(deltarhomat1-deltarhomat2)),intloc,
c$$$     $              (e_bnd(intloc(1),k0,ispin)-chem_pot)*13.605d0,
c$$$     $              (e_bnd(intloc(2),k0,ispin)-chem_pot)*13.605d0
c$$$                endif
                
                deallocate(umat)
                deallocate(proj_mat)
                deallocate(tmp1)
                deallocate(tmp2)
                deallocate(deltarhomat1)
                deallocate(deltarhomat2)                                

              endif
            else
              call ferm_unpack_tau(g_full_0(1,1,k,ispin),
     &          g_full(1,1,1,1,k,ispin),
     &          n_bnd(k0,ispin),nbndf,nbndf,2)
              g_full_0(:,:,k,ispin)=-g_full_0(:,:,k,ispin)
            endif
          enddo                 !! over k
        enddo                   !! over ispin
      endif
      if(nproc_tau/=1) then
        if(.not.correlation) then
          len=nbndf*maxval(ndim3_k)*nspin_0
          if(len/=0) call dgop(g_full_00,len,'  +',comm_pnk)
          if (is_lda_dmft) then
            len=nbndf**2*ndim3_k(me3_k+1)*nspin_0
            if (len/=0) call dgop(g_full_0,2*len,'  +',comm_pnk)
          endif
        else
          len=nbndf**2*maxval(ndim3_k)*nspin_0
          if(len/=0) call dgop(g_full_0,2*len,'  +',comm_pnk)
        endif
      endif
      end
