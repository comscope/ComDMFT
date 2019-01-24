      subroutine qp_green(ispin,ind_k)
	use atom_mod
	use etot_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      integer, intent(in) :: ispin,ind_k
	integer :: k,n,i,ind_tau,i_tau,ic,ib,jb,ind_omega,i_omega
	real*8 :: tu,de,green0,w_n,eps
	complex*16 :: g_diag
	complex*16, allocatable :: g_cmpl(:,:,:)
	eps=1.d-16
	allocate(g_cmpl(nbndf,nbndf,2))
	gc_tau=0.d0
	gc_omega=0.d0
	k=n3_mpi_k(me3_k+1)+ind_k
      n=n_bnd(k,ispin)
	do ind_tau=1,ndim3_tau
	  i_tau=me3_tau*ndim3_tau+ind_tau-1
	  g_cmpl=(0.d0,0.d0)
	  do ic=1,2
	    if(ic==1) tu=tau_mesh(i_tau)
	    if(ic==2) tu=tau_mesh(n_tau-i_tau)
c ----- Full QP Green's function (TAU) -----------------------------
     	    do i=1,n
            de=e_qp(i,k,ispin)-chem_pot_qp
	      g_diag=green0(de,tu)
	      do jb=1,n
	        do ib=1,n
	          g_cmpl(ib,jb,ic)=g_cmpl(ib,jb,ic)
     & 			                +q_qp(ib,i,ind_k,ispin)
     &				            *g_diag
     &                          *conjg(q_qp(jb,i,ind_k,ispin))
	        enddo
	      enddo
          enddo 
c ----- Minus exchange part -----------------------------
     	    do i=1,n
            de=e_bnd(i,k,ispin)-chem_pot
	      g_diag=green0(de,tu)
	      g_cmpl(i,i,ic)=g_cmpl(i,i,ic)-g_diag
          enddo
c          call equal(2,g_cmpl(1,1,ic),eps,nbndf**2)
        enddo
c -----------------------------------------------------------      
		call ferm_pack_tau(g_cmpl,gc_tau(1,1,1,ind_tau),n,nbndf,nbndf)
	enddo
	do ind_omega=1,ndim3_omega
	  i_omega=me3_tau*ndim3_omega+ind_omega-1
	  w_n=w_omega(i_omega)
	  g_cmpl=(0.d0,0.d0)
c ----- Full QP Green's function (OMEGA) -----------------------------
     	  do i=1,n
          de=e_qp(i,k,ispin)-chem_pot_qp
	    g_diag=(1.d0,0.d0)/(dcmplx(0.d0,w_n)-de)
	    do jb=1,n
	      do ib=1,n
	        g_cmpl(ib,jb,1)=g_cmpl(ib,jb,1)
     & 	 	                 +q_qp(ib,i,ind_k,ispin)
     &				         *g_diag
     &                       *conjg(q_qp(jb,i,ind_k,ispin))
	      enddo
	    enddo
        enddo 
c ----- Minus exchange part -----------------------------
     	  do i=1,n
          de=e_bnd(i,k,ispin)-chem_pot
	    g_diag=(1.d0,0.d0)/(dcmplx(0.d0,w_n)-de)
	    g_cmpl(i,i,1)=g_cmpl(i,i,1)-g_diag
        enddo
c        call equal(2,g_cmpl(1,1,1),eps,nbndf**2)
c -----------------------------------------------------------
	  call ferm_pack_omega(g_cmpl,gc_omega(1,1,1,ind_omega),n,nbndf,
     &	                     nbndf)
	enddo  !! over ind_omega
	deallocate(g_cmpl)
      end
