      subroutine sig_omega0_solid(n,ind_k,ispin)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,ind_k,ispin
      integer :: i,j,ind_omega,i_omega,n0
      complex*16, allocatable :: tmp(:,:,:),tmp0(:,:)
      character*30 :: chk
      integer :: ii, jj, kk,k
      n0=nbndf
      allocate(tmp(n0,n0,2))
      allocate(tmp0(n0,0:n_omega))
      tmp=(0.d0,0.d0)
      tmp0=(0.d0,0.d0)
      do ind_omega=1,ndim3_omega
        i_omega=me3_tau*ndim3_omega+ind_omega-1
        do j=1,n
          tmp0(j,i_omega)=dcmplx(sig_c_omega(j,j,1,ind_omega,ind_k,
     &      ispin),
     &      sig_c_omega(j,j,2,ind_omega,ind_k,
     &      ispin))
        enddo
        if(i_omega<=1) then
          i=i_omega+1
          call ferm_unpack_omega(tmp(1,1,i),
     &      sig_c_omega(1,1,1,ind_omega,ind_k,
     &      ispin),n,n0,n0)
        endif
      enddo                     !! over ind_omega

c$$$      k=n3_mpi_k(me3_k+1)+ind_k
c$$$      if (k .lt. 10) then
c$$$        write(chk, '(i1)') k
c$$$      elseif (k .lt. 100) then
c$$$        write(chk, '(i2)') k
c$$$      elseif (k .lt. 1000) then
c$$$        write(chk, '(i3)') k
c$$$      elseif (k .lt. 10000) then
c$$$        write(chk, '(i4)') k        
c$$$      endif
c$$$
c$$$      open(unit=55, file='sig_'//adjustl(trim(chk))//'.dat')
c$$$      do kk=1, 2
c$$$        do ii=1, n
c$$$          write(55, '(100000(2f20.12, 4x))')
c$$$     $      (tmp(ii,jj,kk), jj=1, n)
c$$$        enddo
c$$$      enddo
c$$$      close(55)
      
      if(nproc_tau/=1) then
        call dgop(tmp0,2*n0*(n_omega+1),'  +',comm_pnk)
        call dgop(tmp,2*n0**2*2,'  +',comm_pnk)
      endif
      call sigma_0_z_ren(n,n0,tmp,tmp0)
      call pack_hermit(tmp,sig_c_0(1,1,ind_k,ispin),n,n0,n0,
     &  (0.d0,0.d0),(1.d0,0.d0))
      call pack_hermit(tmp(1,1,2),z_ren(1,1,ind_k,ispin),n,n0,
     &  n0,(0.d0,0.d0),(1.d0,0.d0))
      deallocate(tmp,tmp0)
      end
