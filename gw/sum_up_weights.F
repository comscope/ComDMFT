      subroutine sum_up_weights(weigd,e,n0)
	use atom_mod
	use etot_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: n0
      real*8, intent(in) :: e
      real*8, intent(out) :: weigd(n0,npnt,nspin)
	integer :: j,i,ispin,ind_k,k,n
	real*8 :: gg
	complex*16 :: om0,omi
      weigd=0.d0
      om0=dcmplx(e,0.d0)
      omi=dcmplx(0.d0,e_small)
	do ispin=1,nspin
	  do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
	    n=n_bnd(k,ispin)
          if(ubi=='dft'.or.ubi==' hf') then
            do i=1,n
              gg=imag((1.d0,0.d0)/(om0+omi+chem_pot-e_bnd(i,k,ispin)))
	        weigd(i,k,ispin)=-wgt(k)*gg/pi
	      enddo
          else if(ubi==' qp') then
            do i=1,n
              gg=imag((1.d0,0.d0)/(om0+omi+chem_pot_qp-e_qp(i,k,ispin)))
              do j=1,n
	          weigd(j,k,ispin)=weigd(j,k,ispin)
     &	                        -wgt(k)*gg*abs(q_qp(i,j,ind_k,ispin))**2
     &                          /pi
              enddo
	      enddo
	    endif
        enddo
      enddo
	if(nproc_k/=1) call dgop(weigd,n0*npnt*nspin,'  +',comm_pnt)
	weigd=weigd*2.d0/dfloat(nspin)/nrel
      end
