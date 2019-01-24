      subroutine v_mm_r_from_k
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      include "mpif.h"
      integer :: k,k0,ind,iatom,jatom,isort,jsort,i0,j0,l1,l,m,m1,l1m1,
     &  lm,ir0,ir,ind_k,k1,ind_ir0
      double precision :: pi2
      complex*16, allocatable :: sc(:,:,:)
      allocate(sc(nlb2,nlb2,nqdiv))
      pi2=pi+pi
      sc=(0.d0,0.d0)
      do ind_k=1,ndim3_k(me3_k+1)
        k0=n3_mpi_k(me3_k+1)+ind_k
        do k1=1,k_star(k0)
          k=k_list(k1,k0)
          call zone1_number(pnt(:,k),rb0,ndiv,ind)
          call sym_vmt0(k,k0,ind_k,sc(1,1,ind))
        enddo
      enddo                     !! over ind_k
      if(nproc_k/=1) call mpi_allreduce_dcmplx
     $  (sc,nlb2**2*nqdiv,mpi_sum,comm_tau)
      do j0=1,nlb2
        do i0=1,nlb2
          call fft3(ndiv(1),ndiv(2),ndiv(3),1,sc(i0,j0,:),-1)
        enddo
      enddo
      sc=sc/dfloat(nqdiv)
      if(allocated(s_tt_r0)) deallocate(s_tt_r0)
      allocate(s_tt_r0(limlpb,limlpb,ndim_nrr(me3_k+1)))
      s_tt_r0=0.d0
      do ind_ir0=1,ndim_nrr(me3_k+1)
        ir0=n_mpi_nrr(me3_k+1)+ind_ir0
        ir=nrr_red(1,ir0)
        iatom=nrr_red(2,ir0)
        jatom=nrr_red(3,ir0)
        jsort=is(jatom)
        j0=iob2(jatom)-1
        isort=is(iatom)
        i0=iob2(iatom)-1
        do l1=0,lmpb(jsort)
          do l=0,lmpb(isort)
            do m1=-l1,l1
              l1m1=l1*(l1+1)+m1+1
              do m=-l,l
                lm=l*(l+1)+m+1
                s_tt_r0(lm,l1m1,ind_ir0)=sc(i0+lm,j0+l1m1,ir)
              enddo
            enddo
          enddo
        enddo
      enddo
      deallocate(sc)
      end
