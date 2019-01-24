      subroutine w_cheb_prepare_red(w_n,w_n_cheb,w_tau)
      use atom_mod
      use heg_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
	real*8, intent(out) :: w_n(n_pbtot_red,n_pbtot_red,0:n_nu,nqdiv),
     &  	                 w_n_cheb(n_pbtot_red,n_pbtot_red,
     &                                0:n_nu_asy-2,nqdiv),
     &  	                 w_tau(n_pbtot_red,n_pbtot_red,ndim3_tau,
     &                             nqdiv)
	integer :: k,ig,ind_q,iq,n,kst,ind_nu,i_nu,nn
	real*8, allocatable :: wnu(:,:,:)
	complex*16, allocatable :: tmp(:,:)
      n=n_pbtot_red
      nn=n*n
      allocate(wnu(n,n,ndim3_nu))
c ----------------------------------------------------------------
      allocate(tmp(n,n))
      w_n=0.d0
      w_n_cheb=0.d0
      w_tau=0.d0
      do ind_q=1,ndim3_k
        k=ndim3_k*me3_k+ind_q
	  do kst=1,k_star(k)
	    iq=k_list(kst,k)
	    ig=k_group(iq)
	    wnu=0.d0
          do ind_nu=1,ndim3_nu
            i_nu=me3_tau*ndim3_nu+ind_nu-1
            tmp=w_red_q(:,:,ind_nu,ind_q)
	      call sym_w_red(iq,k,tmp)
	      call pack_hermit(tmp,w_n(1,1,i_nu,iq),n,n,n,0.d0,1.d0)
            wnu(:,:,ind_nu)=w_n(:,:,i_nu,iq)
          enddo  !! over ind_nu
	    if(nproc_tau/=1) call dgop(w_n(1,1,0,iq),nn*(n_nu+1),'  +',
     &	                             comm_pnk)
c --- Prepare Chebyshev's representation for W-interpolation -----
	    call dgemm('n','t',nn,n_nu_asy-1,n_nu_asy-1,1.d0,
     &	             w_n(1,1,n_nu_exa+1,iq),nn,cheb_nu(0,1),
     &               n_nu_asy-1,0.d0,w_n_cheb(1,1,0,iq),nn)
          call from_nu_to_tau_ba(wnu,n,w_tau(1,1,1,iq),n,n)
        enddo   !! over kst
      enddo   !! over ind_q
      deallocate(tmp,wnu)
	if(nproc_k/=1) then
	  call dgop(w_n,nn*nqdiv*(n_nu+1),'  +',comm_pnt)
	  call dgop(w_n_cheb,nn*nqdiv*(n_nu_asy-1),'  +',comm_pnt)
	  call dgop(w_tau,nn*nqdiv*ndim3_tau,'  +',comm_pnt)
	endif
      end