      subroutine z_wan_rr_int(z0,zi)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use wannier90_mod      
      implicit none
      include 'mpif.h'
      complex*16, intent(in) :: z0(nbasmpw,num_orb_proj,nspin,
     &  ndim_kk(me+1))      
      complex*16, intent(out) :: zi(nr_full,num_orb_proj,nrel,nspin,
     &  ndim_kk(me+1))
      integer :: i,j,k,k0,i0,me1,ispin,ind_k,ind,ir,irl,nbas,ibas0
      double precision :: pi2,fac,v(3),fc
      complex*16 :: cc
      complex*16, allocatable :: z1(:,:,:),zz(:)
!      write(100,*)'z0',z0,'amega',amega,'gbs',gbs,'i_kref',i_kref,
!     & 'indgb',indgb,'me',me,'n_mpi_kk',n_mpi_kk,'nbask',nbask,
!     & 'nbasmpw',nbasmpw,'ndim_kk',ndim_kk,'ndiv',ndiv,'nproc',nproc,
!     & 'nqdiv',nqdiv,'nr_full',nr_full,'nrdiv',nrdiv,'nrel',nrel,
!     & 'nspin',nspin,'num_orb_proj',num_orb_proj,'pi',pi,'pnt',pnt,
!     & 'r0b0',r0b0,'r_pnt',r_pnt,'rb0',rb0,'rbas',rbas
!      call mpi_bcast(nproc_k,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)      


      me1=me+1
      pi2=pi+pi
      allocate(z1(nqdiv,nr_full,nrel))
      allocate(zz(nr_full))
      fac=1.d0/sqrt(amega)/dble(nqdiv)
      zi=(0.d0,0.d0)
      do ispin=1,nspin
        do j=1,num_orb_proj
          z1=(0.d0,0.d0)
          do ind_k=1,ndim_kk(me1)
            k=n_mpi_kk(me1)+ind_k
            k0=i_kref(k)
            nbas=nbask(k0)/nrel                        

            do irl=1,nrel
              ibas0=(irl-1)*nbas
              zz=(0.d0,0.d0)              
              do i=1,nbas
                i0=indgb(i,k)
                call zone1_number(gbs(:,i0),rbas,nrdiv,ind)
                zz(ind)=zz(ind)+z0(ibas0+i,j,ispin,ind_k)
              enddo
              
              call fft3(nrdiv(1),nrdiv(2),nrdiv(3),1,zz,1)
              call zone1_number(pnt(:,k),rb0,ndiv,ind)
              do i=1,nr_full
c     ----- Multiplication with the phase factors ------------------
                v(:)=r_pnt(1,i)*r0b0(:,1)+r_pnt(2,i)*r0b0(:,2)
     &            +r_pnt(3,i)*r0b0(:,3)
                fc=pi2*dot_product(pnt(:,k),v)
                cc=cdexp(dcmplx(0.d0,fc))
                z1(ind,i,irl)=zz(i)*cc
              enddo
            enddo            
          enddo
            
          if(nproc/=1) call mpi_allreduce_dcmplx
     $      (z1,nqdiv*nr_full*nrel,mpi_sum,MPI_COMM_WORLD)
          call fft3(ndiv(1),ndiv(2),ndiv(3),nr_full*nrel,z1,1)

          do i=1,nr_full
            do irl=1,nrel
              do ind_k=1,ndim_kk(me1)
                ir=n_mpi_kk(me1)+ind_k
                zi(i,j,irl,ispin,ind_k)=z1(ir,i,irl)*fac! i: fft grid, ind_k: fft grid
              enddo
            enddo
          enddo
        enddo
      enddo
      deallocate(zz,z1)
      end
      
