      subroutine output_gx_band_tau(ispin,k,n0)
	use atom_mod
	use etot_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: ispin,k,n0
	integer :: i_len,m,n,l1,l2,m0,i_tau,i1_tau,l,ll,i,ii
	real*8 :: fc,t1,t2,de,green0
	real*8, allocatable :: g_out(:,:)
	m=5
	n=n0/m
	if(n*m<n0) n=n+1
	fc=evolt/2.d0
	allocate(g_out(0:n_tau,m))
	if(maswrk) then
	  i_len=len_trim(allfile)
	  open(3,file=allfile(1:i_len)//'_Gx_tau_band_k')
	  if(ispin/=1.or.k/=1) call goto_end(3)
	  write(3,*)' Ispin = ',ispin,' K = ',k
	endif
	do i=1,n
        l1=m*(i-1)+1
        l2=min(m*i,n_bnd(k,ispin))
	  if(maswrk) write(3,*)' Bands from ',l1,' to ',l2
	  m0=l2-l1+1
	  if(m0<=0) cycle
	  g_out=0.d0
	  do i_tau=0,n_tau/2
	    t1=tau_mesh(i_tau)
	    i1_tau=n_tau-i_tau
	    t2=tau_mesh(i1_tau)
	    do l=1,m0
	      ll=l1+l-1
            de=e_bnd(ll,k,ispin)-chem_pot
	      g_out(i_tau,l)=green0(de,t1)
	      g_out(i1_tau,l)=green0(de,t2)
	    enddo
	  enddo
	  if(maswrk) then
	    do i_tau=0,n_tau
	      write(3,'(f9.3,5(1x,f12.7))')tau_mesh(i_tau),
     &		               (g_out(i_tau,ii),ii=1,m0)
	    enddo
	  endif
	enddo
	if(maswrk) close(3)
	deallocate(g_out)
      end