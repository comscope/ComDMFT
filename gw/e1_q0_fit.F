      subroutine e1_q0_fit(e1)
c	 Deals with (1-VP)^(-1)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	real*8, intent(inout) :: e1(n_line,0:n_nu)
	integer :: k,i_len,i_nu,nd
	real*8 :: q2,hh,a0,an,ai,fc,q,heg_e1_small_q,pi2a,ee1,ee2,f,aa,bb
	real*8, allocatable :: e0(:),e2(:)
	if(nqdiv==1) return
	pi2a=(pi+pi)/par
c -------- Real Axis Values ------------------------------------
	i_len=len_trim(allfile)
	if(maswrk) open(3,file=allfile(1:i_len)//'_em1_gg0_k')
	do k=2,n_line
	  q2=dot_product(q_line(:,k),q_line(:,k))
	  q=sqrt(q2)*pi2a
	  if(maswrk) write(3,'(8(1x,f8.5))')q,(e1(k,i_nu),i_nu=0,6)
	enddo
	if(maswrk) close(3)
	if(maswrk) then
	  open(3,file=allfile(1:i_len)//'_em1_gg0_nu')
	  do i_nu=0,n_nu
	    write(3,'(e9.2,7(1x,f8.5))')w_nu(i_nu),(e1(k,i_nu),k=1,n_line)
	  enddo
	  close(3)
	endif
c ------------- Initial fitting ---------------------------------------
	allocate(e0(0:n_nu))
	if(n_line<3) then
c --- Non-interacting HEG values -------------------------
        do i_nu=0,n_nu
          e0(i_nu)=heg_e1_small_q(w_nu(i_nu),0.d0,ro_el_gas)
        enddo
      else if(n_line==3) then
c --- Extrapolated (linear) values -------------------------
	  do i_nu=0,n_nu
	    hh=e1(2,i_nu)-e1(3,i_nu)
	    e0(i_nu)=e1(2,i_nu)+hh
	  enddo
      else if(n_line>3) then
c --- Extrapolated (quadratic) values -------------------------
	  do i_nu=0,n_nu
	    ee1=e1(2,i_nu)-e1(3,i_nu)
	    ee2=e1(3,i_nu)-e1(4,i_nu)
	    hh=ee1-ee2
	    e0(i_nu)=e1(2,i_nu)+ee1+hh
	  enddo
	endif
	e1(1,:)=e0
c --------------------------------------------------------------------
	do i_nu=1,n_nu
	  if(e1(1,i_nu)<e1(1,i_nu-1)) then
          e1(1,i_nu)=e1(1,i_nu-1)
	  endif
	enddo
c --------------------------------------------------------
	a0=e1(1,0)
	if(a0<0.d0) then
	  an=e1(1,n_nu)
	  if(abs(an-a0)>1.d-10) then
	  do i_nu=0,n_nu
	    ai=e1(1,i_nu)
	    e1(1,i_nu)=an-(ai-an)*an/(a0-an)
	  enddo
	  endif
	endif
	if(metal) e1(1,0)=0.d0
c --------------------------------------------------------
	an=e1(1,n_nu)
	if(an>1.d0) then
	  a0=e1(1,0)
	  if(abs(an-a0)>1.d-10) then
	    fc=(1.d0-a0)/(an-a0)
	    do i_nu=0,n_nu
	      ai=e1(1,i_nu)
	      e1(1,i_nu)=a0+fc*(ai-a0)
	    enddo
	  endif
	endif
c ----- Plasmon pole variant -----------------------------------
      if(metal) then
	  allocate(e2(0:n_nu))
        e2=e1(1,:)
        f=e2(0)+0.9d0*(e2(n_nu)-e2(0))
        nd=0
        a0=abs(f-e2(0))
        do i_nu=1,n_nu
          an=abs(f-e2(i_nu))
          if(an<a0) then
            a0=an
            nd=i_nu
          endif
        enddo
        a0=e2(0)/(e2(0)-1.d0)
        an=e2(nd)/(e2(nd)-1.d0)
        aa=(an-a0)/w_nu(nd)**2
        do i_nu=1,n_nu
          bb=w_nu(i_nu)**2*aa+a0
          bb=1.d0/bb
          bb=1.d0-bb
          if(abs(bb)>1.d-9) e2(i_nu)=1.d0/bb
        enddo
        e1(1,:)=e2
        deallocate(e2)
      endif
c -----------------------------------------------------
	if(maswrk) then
	  open(3,file=allfile(1:i_len)//'_em1_gg0_fit_nu')
	  do i_nu=0,n_nu
	    write(3,'(e12.5,2(1x,f18.10))')w_nu(i_nu),e0(i_nu),e1(1,i_nu)
	  enddo
	  close(3)
	endif
	deallocate(e0)
      end
      
      
      subroutine e2_q0_fit(e1)
c	 Deals with (1-VP)^(-1)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	real*8, intent(inout) :: e1(n_line,0:n_nu)
	integer :: k,i_len,i_nu,nd
	real*8 :: q2,a0,an,ai,fc,q,heg_e1_small_q,pi2a,f,aa,bb
	real*8, allocatable :: e0(:),e2(:)
	if(nqdiv==1) return
	pi2a=(pi+pi)/par
c -------- Real Axis Values ------------------------------------
	i_len=len_trim(allfile)
	if(maswrk) open(3,file=allfile(1:i_len)//'_em1_gg0_k')
	do k=2,n_line
	  q2=dot_product(q_line(:,k),q_line(:,k))
	  q=sqrt(q2)*pi2a
	  if(maswrk) write(3,'(8(1x,f8.5))')q,(e1(k,i_nu),i_nu=0,6)
	enddo
	if(maswrk) close(3)
	if(maswrk) then
	  open(3,file=allfile(1:i_len)//'_em1_gg0_nu')
	  do i_nu=0,n_nu
	    write(3,'(e9.2,7(1x,f8.5))')w_nu(i_nu),(e1(k,i_nu),k=1,n_line)
	  enddo
	  close(3)
	endif
c ------------- Initial fitting ---------------------------------------
	allocate(e0(0:n_nu))
	if(n_line<3) then
c --- Non-interacting HEG values -------------------------
        do i_nu=0,n_nu
          e0(i_nu)=heg_e1_small_q(w_nu(i_nu),0.d0,ro_el_gas)
        enddo
      else
c -----------  E1(q) = e1(0) + a q^2 ----------------------------     
	  do i_nu=0,n_nu
	    e0(i_nu)=(4.d0*e1(2,i_nu)-e1(3,i_nu))/3.d0
	  enddo
	endif
	e1(1,:)=e0
c --------------------------------------------------------------------
	do i_nu=1,n_nu
	  if(e1(1,i_nu)<e1(1,i_nu-1)) then
          e1(1,i_nu)=e1(1,i_nu-1)
	  endif
	enddo
c --------------------------------------------------------
	a0=e1(1,0)
	if(a0<0.d0) then
	  an=e1(1,n_nu)
	  if(abs(an-a0)>1.d-10) then
	  do i_nu=0,n_nu
	    ai=e1(1,i_nu)
	    e1(1,i_nu)=an-(ai-an)*an/(a0-an)
	  enddo
	  endif
	endif
	if(metal) e1(1,0)=0.d0
c --------------------------------------------------------
	an=e1(1,n_nu)
	if(an>1.d0) then
	  a0=e1(1,0)
	  if(abs(an-a0)>1.d-10) then
	  fc=(1.d0-a0)/(an-a0)
	  do i_nu=0,n_nu
	    ai=e1(1,i_nu)
	    e1(1,i_nu)=a0+fc*(ai-a0)
	  enddo
	  endif
	endif
c ----- Plasmon pole variant -----------------------------------
      if(metal) then
	  allocate(e2(0:n_nu))
        e2=e1(1,:)
        f=e2(0)+0.9d0*(e2(n_nu)-e2(0))
        nd=0
        a0=abs(f-e2(0))
        do i_nu=1,n_nu
          an=abs(f-e2(i_nu))
          if(an<a0) then
            a0=an
            nd=i_nu
          endif
        enddo
        a0=e2(0)/(e2(0)-1.d0)
        an=e2(nd)/(e2(nd)-1.d0)
        aa=(an-a0)/w_nu(nd)**2
        do i_nu=1,n_nu
          bb=w_nu(i_nu)**2*aa+a0
          bb=1.d0/bb
          bb=1.d0-bb
          if(abs(bb)>1.d-9) e2(i_nu)=1.d0/bb
        enddo
        e1(1,:)=e2
        deallocate(e2)
      endif
c -----------------------------------------------------
	if(maswrk) then
	  open(3,file=allfile(1:i_len)//'_em1_gg0_fit_nu')
	  do i_nu=0,n_nu
	    write(3,'(e12.5,2(1x,f18.10))')w_nu(i_nu),e0(i_nu),e1(1,i_nu)
	  enddo
	  close(3)
	endif
	deallocate(e0)
      end
      