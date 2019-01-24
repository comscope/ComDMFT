
      subroutine z_wan_rr_mt(z0,zt)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use wannier90_mod      
      implicit none
      include 'mpif.h'
      complex*16, intent(in) :: z0(nfun,num_orb_proj,nspin,
     &  ndim_kk(me+1))      
      complex*16, intent(out) :: zt(nfun,num_orb_proj,nspin,
     &  ndim_kk(me+1))
      integer :: i,j,k,k0,me1,ispin,ind_k,ind_i,ind,ir
      integer, allocatable :: ndim_nf(:),n_mpi_nf(:)
      complex*16, allocatable :: z1(:,:),zz(:,:)
      me1=me+1
      allocate(ndim_nf(nproc))
      allocate(n_mpi_nf(nproc))
      call size_shift_par(nfun,nproc,ndim_nf,n_mpi_nf)
      allocate(z1(nqdiv,nfun))
      allocate(zz(nqdiv,nfun))
      zt=(0.d0,0.d0)
      do ispin=1,nspin
        do j=1,num_orb_proj
          zz=(0.d0,0.d0)
          do i=1,nfun
            do ind_k=1,ndim_kk(me1)
              k=n_mpi_kk(me1)+ind_k
              call zone1_number(pnt(:,k),rb0,ndiv,ind)
              zz(ind,i)=z0(i,j,ispin,ind_k)
            enddo
          enddo
          if(nproc/=1) call mpi_allreduce_dcmplx
     $      (zz,nqdiv*nfun,mpi_sum,MPI_COMM_WORLD)
          z1=(0.d0,0.d0)
          do ind_i=1,ndim_nf(me1)
            i=n_mpi_nf(me1)+ind_i
            call fft3(ndiv(1),ndiv(2),ndiv(3),1,zz(1,i),1)
            do ir=1,nqdiv
              z1(ir,i)=zz(ir,i)
            enddo
          enddo
          if(nproc/=1) call mpi_allreduce_dcmplx
     $      (z1,nqdiv*nfun,mpi_sum,MPI_COMM_WORLD)
          do i=1,nfun
            do ind_k=1,ndim_kk(me1)
              ir=n_mpi_kk(me1)+ind_k
              zt(i,j,ispin,ind_k)=z1(ir,i)/dble(nqdiv)! ind_k: index in fftgrid
            enddo
          enddo
        enddo
      enddo
      deallocate(ndim_nf,n_mpi_nf,zz,z1)
      end
