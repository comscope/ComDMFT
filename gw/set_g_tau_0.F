      subroutine set_g_tau_0(ispin,ind_k,key,g_tau)
c     --- Construct G_full(tau) from non-interacting data ------
c     key = 0  - from LDA(HF) - E_BND
c     key = 1  - from QP      - E_QP      
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      integer, intent(in) :: ispin,ind_k,key
      real*8, intent(out) :: g_tau(nbndf,nbndf,2,ndim3_tau)
	integer :: k,n,i,ind_tau,i_tau,it,jb,ib
	real*8 :: g_diag,de,green0
	complex*16, allocatable :: g_cmpl(:,:,:)
	allocate(g_cmpl(nbndf,nbndf,2))
	k=n3_mpi_k(me3_k+1)+ind_k
      n=n_bnd(k,ispin)
      g_tau=0.d0
	do ind_tau=1,ndim3_tau
	  i_tau=me3_tau*ndim3_tau+ind_tau-1
	  g_cmpl=(0.d0,0.d0)
	  do it=1,2
	    if(it==2) i_tau=n_tau-i_tau
	    if(key==0) then
	      do i=1,n
	        de=e_bnd(i,k,ispin)-chem_pot
	        g_cmpl(i,i,it)=dcmplx(green0(de,tau_mesh(i_tau)),0.d0)
            enddo
          else if(key==1) then
     	      do i=1,n
              de=e_qp(i,k,ispin)-chem_pot_qp
	        g_diag=green0(de,tau_mesh(i_tau))
	        do jb=1,n
	          do ib=1,n
	            g_cmpl(ib,jb,it)=g_cmpl(ib,jb,it)
     &	                          +q_qp(ib,i,ind_k,ispin)
     &				              *g_diag
     &                            *conjg(q_qp(jb,i,ind_k,ispin))
	          enddo
	        enddo
            enddo 
	    endif
	  enddo
		call ferm_pack_tau(g_cmpl,g_tau(1,1,1,ind_tau),n,nbndf,nbndf)
	enddo
	deallocate(g_cmpl)
      end
