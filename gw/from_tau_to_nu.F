      subroutine from_tau_to_nu_aa_g(f_nu,f_tau,n)
c   ------- Transform TAU ---> NU representation --------------------------
c	F_tau - real symmetric, but is stored in square form
c     F(tau) == F(Beta-tau)
c     F_nu is also real symmetric and it is required in square form
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: f_tau(n,ndim3_tau)
	real*8, intent(out) :: f_nu(n,ndim3_nu)
      integer :: ind_tau,ind_nu,i,i_nu,i1,m,m0
      real*8 :: yy1
	real*8, allocatable :: at(:,:),an(:,:)
	allocate(at(n,ndim3_tau))
	allocate(an(n,0:n_nu_d))
	do ind_tau=1,ndim3_tau
	  at(:,ind_tau)=2.d0*f_tau(:,ind_tau)
	enddo
	call from_tau_to_nu_cos_d(an,at,n)
	f_nu=0.d0
      m=interp_nu_d
      m0=m-1
	do i=1,n
	  do ind_nu=1,ndim3_nu
	    i_nu=me3_tau*ndim3_nu+ind_nu-1
	    if(i_nu<=n_nu_exa) then
	      f_nu(i,ind_nu)=an(i,i_nu)
	    else if(i_nu<n_nu) then
	      i1=ind_nu_int(i_nu)
	      call fint(an(i,i1:i1+m0),w_nu_d(i1),m,w_nu(i_nu),yy1)
	      f_nu(i,ind_nu)=yy1
	    else if(i_nu==n_nu) then
	      f_nu(i,ind_nu)=an(i,n_nu_d)
	    endif
	  enddo
	enddo
	deallocate(at,an)
      end
      
      subroutine from_tau_to_nu_xi(f_nu,f_tau,n0,n)
c   ------- Transform TAU ---> NU representation --------------------------
c   ---- We transform totally general real F_TAU with N:N elements ------
c	F_tau - is stored in (t,beta-t) form
c     F(tau) ne F(Beta-tau)
c     F_nu is required in complex*16 [0:n_nu] form
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n,n0
	real*8, intent(in) :: f_tau(n0,n0,2,ndim3_tau)
	complex*16, intent(out) :: f_nu(n,n,ndim3_nu)
      integer :: ind_tau,ind_nu,i,nn,it,j
	real*8, allocatable :: at(:,:,:,:)
	complex*16, allocatable :: an(:,:)
	nn=n*n
	allocate(at(n,n,2,ndim3_tau))
	allocate(an(n,n))
	do ind_tau=1,ndim3_tau
	  do it=1,2
	    at(1:n,1:n,it,ind_tau)=f_tau(1:n,1:n,it,ind_tau)
	  enddo
	enddo
	call from_tau_to_nu_re(f_nu,at,nn)
	do ind_nu=1,ndim3_nu
	  an=f_nu(:,:,ind_nu)
	  do j=1,n
	    do i=1,j-1
	      f_nu(i,j,ind_nu)=an(i,j)+(0.d0,1.d0)*an(j,i)
	    enddo
	    f_nu(j,j,ind_nu)=an(j,j)
	    do i=j+1,n
	      f_nu(i,j,ind_nu)=an(j,i)-(0.d0,1.d0)*an(i,j)
	    enddo
	  enddo
	enddo
	deallocate(at,an)
c	allocate(bt(n,n,2,ndim3_tau))
c	do ind_tau=1,ndim3_tau
c	  do it=1,2
c	      call unpack_hermit(bt(1,1,it,ind_tau),f_tau(1,1,it,ind_tau),
c     &	                       n,n0,n,(0.d0,0.d0),(1.d0,0.d0))
c        enddo
c      enddo
c      call from_tau_to_nu_gen(f_nu,bt,nn)
c	deallocate(bt)
      end
      
      subroutine from_tau_to_nu_sin_d(f_nu,f_tau,n)
c   ------- Transform TAU ---> NU representation ----------------------
c   This is a basic operation B(v) = Int_0^beta/2 B(tau) sin(v*tau)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: f_tau(n,ndim3_tau)
	real*8, intent(out) :: f_nu(n,0:n_nu_d)
	call dgemm('n','n',n,n_nu_d+1,ndim3_tau,1.d0,f_tau,n,
     &	         ab_nu_from_tau(1,0,2),ndim3_tau,0.d0,f_nu,n)
	if(nproc_tau/=1) call dgop(f_nu,n*(n_nu_d+1),'  +',comm_pnk)
      end
      
      subroutine from_tau_to_nu_re_d(f_nu,f_tau,n)
c   ------- Transform TAU ---> NU representation --------------------------
c   ---- We transform totally general real F_TAU with N elements ------
c	F_tau - is stored in (t,beta-t) form
c     F(tau) ne F(Beta-tau)
c     F_nu is required in complex*16 [0:n_nu_d] form
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: f_tau(n,2,ndim3_tau)
	complex*16, intent(out) :: f_nu(n,ndim3_nu_d(me3_tau+1))
      integer :: i_tau,ind_tau,ind_nu,i_nu
	real*8, allocatable :: at(:,:),bt(:,:),an(:,:),bn(:,:)
	allocate(at(n,ndim3_tau))
	allocate(bt(n,ndim3_tau))
	allocate(an(n,0:n_nu_d))
	allocate(bn(n,0:n_nu_d))
      do ind_tau=1,ndim3_tau
        i_tau=me3_tau*ndim3_tau+ind_tau-1
	  at(:,ind_tau)=f_tau(:,1,ind_tau)+f_tau(:,2,ind_tau)
	  bt(:,ind_tau)=f_tau(:,1,ind_tau)-f_tau(:,2,ind_tau)
	enddo
	call from_tau_to_nu_cos_d(an,at,n)
	call from_tau_to_nu_sin_d(bn,bt,n)
      do ind_nu=1,ndim3_nu_d(me3_tau+1)
	  i_nu=n3_mpi_nu_d(me3_tau+1)+ind_nu-1
	  f_nu(:,ind_nu)=dcmplx(an(:,i_nu),bn(:,i_nu))
	enddo
	deallocate(at,bt,an,bn)
      end
      
      subroutine from_tau_to_nu_re(f_nu,f_tau,n)
c   ------- Transform TAU ---> NU representation --------------------------
c   ---- We transform totally general real F_TAU with N elements ------
c	F_tau - is stored in (t,beta-t) form
c     F(tau) ne F(Beta-tau)
c     F_nu is required in complex*16 [0:n_nu] form
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: f_tau(n,2,ndim3_tau)
	complex*16, intent(out) :: f_nu(n,ndim3_nu)
      integer :: i_tau,ind_tau,ind_nu,i_nu,i1,i,m,m0
      real*8 :: yy1,yy2
	real*8, allocatable :: at(:,:),bt(:,:),an(:,:),bn(:,:)
	allocate(at(n,ndim3_tau))
	allocate(bt(n,ndim3_tau))
	allocate(an(n,0:n_nu_d))
	allocate(bn(n,0:n_nu_d))
      do ind_tau=1,ndim3_tau
        i_tau=me3_tau*ndim3_tau+ind_tau-1
	  at(:,ind_tau)=f_tau(:,1,ind_tau)+f_tau(:,2,ind_tau)
	  bt(:,ind_tau)=f_tau(:,1,ind_tau)-f_tau(:,2,ind_tau)
	enddo
	call from_tau_to_nu_cos_d(an,at,n)
	call from_tau_to_nu_sin_d(bn,bt,n)
      m=interp_nu_d
      m0=m-1
	do ind_nu=1,ndim3_nu
	  i_nu=me3_tau*ndim3_nu+ind_nu-1
	  if(i_nu<=n_nu_exa) then
	    f_nu(:,ind_nu)=dcmplx(an(:,i_nu),bn(:,i_nu))
	  else if(i_nu<n_nu) then
	    i1=ind_nu_int(i_nu)
	    do i=1,n
	      call fint(an(i,i1:i1+m0),w_nu_d(i1),m,w_nu(i_nu),yy1)
	      call fint(bn(i,i1:i1+m0),w_nu_d(i1),m,w_nu(i_nu),yy2)
	      f_nu(i,ind_nu)=dcmplx(yy1,yy2)
	    enddo
	  else if(i_nu==n_nu) then
	    f_nu(:,ind_nu)=dcmplx(an(:,n_nu_d),bn(:,n_nu_d))
	  endif
	enddo
	deallocate(at,bt,an,bn)
      end
      
      subroutine from_tau_to_nu_gen(f_nu,f_tau,n)
c   ------- Transform TAU ---> NU representation ----------------------
c   ---- We transform totally general complex F_TAU with N elements ---
c	F_tau - is stored in (t,beta-t) form
c     F(tau) ne F(Beta-tau)
c     F_nu is required in complex*16 [0:n_nu] form
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	complex*16, intent(in) :: f_tau(n,2,ndim3_tau)
	complex*16, intent(out) :: f_nu(n,ndim3_nu)
      integer :: i_tau,ind_tau,ind_nu,i,i_nu,i1,m,m0
      real*8 :: yy1,yy2
	real*8, allocatable :: at(:,:),bt(:,:),an(:,:),bn(:,:)
	allocate(at(n,ndim3_tau))
	allocate(bt(n,ndim3_tau))
	allocate(an(n,0:n_nu_d))
	allocate(bn(n,0:n_nu_d))
c -------- Real part transformation -----------------------------
      do ind_tau=1,ndim3_tau
        i_tau=me3_tau*ndim3_tau+ind_tau-1
	  at(:,ind_tau)=real(f_tau(:,1,ind_tau))+real(f_tau(:,2,ind_tau))
	  bt(:,ind_tau)=real(f_tau(:,1,ind_tau))-real(f_tau(:,2,ind_tau))
	enddo
	call from_tau_to_nu_cos_d(an,at,n)
	call from_tau_to_nu_sin_d(bn,bt,n)
	f_nu=0.d0
      m=interp_nu_d
      m0=m-1
	do i=1,n
	  do ind_nu=1,ndim3_nu
	    i_nu=me3_tau*ndim3_nu+ind_nu-1
	    if(i_nu<=n_nu_exa) then
	      f_nu(i,ind_nu)=dcmplx(an(i,i_nu),bn(i,i_nu))
	    else if(i_nu<n_nu) then
	      i1=ind_nu_int(i_nu)
	      call fint(an(i,i1:i1+m0),w_nu_d(i1),m,w_nu(i_nu),yy1)
	      call fint(bn(i,i1:i1+m0),w_nu_d(i1),m,w_nu(i_nu),yy2)
	      f_nu(i,ind_nu)=dcmplx(yy1,yy2)
	    else if(i_nu==n_nu) then
	      f_nu(i,ind_nu)=dcmplx(an(i,n_nu_d),bn(i,n_nu_d))
	    endif
	  enddo
	enddo
c -------- Imaginary part transformation -----------------------------
      do ind_tau=1,ndim3_tau
        i_tau=me3_tau*ndim3_tau+ind_tau-1
	  at(:,ind_tau)=imag(f_tau(:,1,ind_tau))+imag(f_tau(:,2,ind_tau))
	  bt(:,ind_tau)=imag(f_tau(:,1,ind_tau))-imag(f_tau(:,2,ind_tau))
	enddo
	call from_tau_to_nu_cos_d(an,at,n)
	call from_tau_to_nu_sin_d(bn,bt,n)
	do i=1,n
	  do ind_nu=1,ndim3_nu
	    i_nu=me3_tau*ndim3_nu+ind_nu-1
	    if(i_nu<=n_nu_exa) then
	      f_nu(i,ind_nu)=f_nu(i,ind_nu)+(0.d0,1.d0)
     &	                             *dcmplx(an(i,i_nu),bn(i,i_nu))
	    else if(i_nu<n_nu) then
	      i1=ind_nu_int(i_nu)
	      call fint(an(i,i1:i1+m0),w_nu_d(i1),m,w_nu(i_nu),yy1)
	      call fint(bn(i,i1:i1+m0),w_nu_d(i1),m,w_nu(i_nu),yy2)
	      f_nu(i,ind_nu)=f_nu(i,ind_nu)+(0.d0,1.d0)
     &	                             *dcmplx(yy1,yy2)
	    else if(i_nu==n_nu) then
	      f_nu(i,ind_nu)=f_nu(i,ind_nu)+(0.d0,1.d0)
     &	                             *dcmplx(an(i,n_nu_d),bn(i,n_nu_d))
	    endif
	  enddo
	enddo
	deallocate(at,bt,an,bn)
      end
      
      subroutine from_tau_to_nu_full(f_nu,f_tau,n)
c   ------- Transform TAU ---> NU representation --------------------------
c   ---- We transform totally general complex F_TAU with N elements ------
c	F_tau - is stored in (t,beta-t) form
c     F(tau) ne F(Beta-tau)
c     F_nu is required in complex*16 [-n_nu(1) : n_nu(2)] form
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	complex*16, intent(in) :: f_tau(n,2,ndim3_tau)
	complex*16, intent(out) :: f_nu(n,ndim3_nu,2)
      integer :: i_tau,ind_tau,ind_nu,i,i_nu,i1,m,m0
      real*8 :: yy1,yy2,yy3,yy4
	complex*16, allocatable :: at(:,:),bt(:,:),an(:,:),bn(:,:)
	allocate(at(n,ndim3_tau))
	allocate(bt(n,ndim3_tau))
	allocate(an(n,0:n_nu_d))
	allocate(bn(n,0:n_nu_d))
      do ind_tau=1,ndim3_tau
        i_tau=me3_tau*ndim3_tau+ind_tau-1
	  at(:,ind_tau)=f_tau(:,1,ind_tau)+f_tau(:,2,ind_tau)
	  bt(:,ind_tau)=f_tau(:,1,ind_tau)-f_tau(:,2,ind_tau)
	enddo
	call from_tau_to_nu_cos_d(an,at,2*n)
	call from_tau_to_nu_sin_d(bn,bt,2*n)
	f_nu=(0.d0,0.d0)
      m=interp_nu_d
      m0=m-1
	do i=1,n
	  do ind_nu=1,ndim3_nu
	    i_nu=me3_tau*ndim3_nu+ind_nu-1
	    if(i_nu<=n_nu_exa) then
	      f_nu(i,ind_nu,1)=an(i,i_nu)-(0.d0,1.d0)*bn(i,i_nu)
	      f_nu(i,ind_nu,2)=an(i,i_nu)+(0.d0,1.d0)*bn(i,i_nu)
	    else if(i_nu<n_nu) then
	      i1=ind_nu_int(i_nu)
	      call fint(real(an(i,i1:i1+m0)),w_nu_d(i1),m,w_nu(i_nu),yy1)
	      call fint(real(bn(i,i1:i1+m0)),w_nu_d(i1),m,w_nu(i_nu),yy2)
	      call fint(imag(an(i,i1:i1+m0)),w_nu_d(i1),m,w_nu(i_nu),yy3)
	      call fint(imag(bn(i,i1:i1+m0)),w_nu_d(i1),m,w_nu(i_nu),yy4)
	      f_nu(i,ind_nu,1)=dcmplx(yy1,yy3)-(0.d0,1.d0)*dcmplx(yy2,yy4)
	      f_nu(i,ind_nu,2)=dcmplx(yy1,yy3)+(0.d0,1.d0)*dcmplx(yy2,yy4)
	    else if(i_nu==n_nu) then
	      f_nu(i,ind_nu,1)=an(i,n_nu_d)-(0.d0,1.d0)*bn(i,n_nu_d)
	      f_nu(i,ind_nu,2)=an(i,n_nu_d)+(0.d0,1.d0)*bn(i,n_nu_d)
	    endif
	  enddo
	enddo
	deallocate(at,bt,an,bn)
      end
      
      subroutine from_tau_to_nu_cos_d(f_nu,f_tau,n)
c   ------- Transform TAU ---> NU representation ----------------------
c   This is a basic operation A(v) = Int_0^beta/2 A(tau) cos(v*tau)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: f_tau(n,ndim3_tau)
	real*8, intent(out) :: f_nu(n,0:n_nu_d)
	call dgemm('n','n',n,n_nu_d+1,ndim3_tau,1.d0,f_tau,n,
     &	         ab_nu_from_tau(1,0,1),ndim3_tau,0.d0,f_nu,n)
	if(nproc_tau/=1) call dgop(f_nu,n*(n_nu_d+1),'  +',comm_pnk)
      end
      
      subroutine from_tau_to_nu_bb(f_nu,n0,f_tau,n1,n)
c   ------- Transform TAU ---> NU  --------------------------
c     Input:
c	F_tau - hermitian, but is stored in square R*8 form = B
c     Output:
c     F_nu is hermitian and stored in B square R*8 form
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n,n0,n1
	real*8, intent(in) :: f_tau(n1,n1,ndim3_tau)
	real*8, intent(out) :: f_nu(n0,n0,ndim3_nu)
      integer :: nn,i,j,ij,ind_tau,nblock,nstep,ind_nu,j0,jm,i_nu,i1,m,
     &           m0
      real*8 :: yy1
	real*8, allocatable :: at(:,:),an(:,:)
	nblock=1000
	nstep=max(1,nblock/n)
	nn=n*nstep
	allocate(at(nn,ndim3_tau))
	allocate(an(nn,0:n_nu_d))
	j0=0
      m=interp_nu_d
      m0=m-1
	do while(j0<n)
	  jm=min(nstep,n-j0)
c -------------------------- SIN transformation ------------------
	  at=0.d0
	  do ind_tau=1,ndim3_tau
	    ij=0
	    do j=j0+1,j0+jm
	      do i=1,n
	        ij=ij+1
	        at(ij,ind_tau)=f_tau(i,j,ind_tau)
	      enddo
	    enddo
	  enddo
	  call from_tau_to_nu_sin_d(an,at,nn)
	  ij=0
	  do j=j0+1,j0+jm
	    do i=1,n
	      ij=ij+1
	      do ind_nu=1,ndim3_nu
	        i_nu=me3_tau*ndim3_nu+ind_nu-1
	        if(i_nu<=n_nu_exa) then
	          f_nu(i,j,ind_nu)=an(ij,i_nu)
	        else if(i_nu<n_nu) then
	          i1=ind_nu_int(i_nu)
	          call fint(an(ij,i1:i1+m0),w_nu_d(i1),m,w_nu(i_nu),yy1)
	          f_nu(i,j,ind_nu)=yy1
	        else if(i_nu==n_nu) then
	          f_nu(i,j,ind_nu)=an(ij,n_nu_d)
	        endif
	      enddo
	    enddo
	  enddo
	  j0=j0+jm
	enddo
	deallocate(at,an)
      end
      
      subroutine from_tau_to_nu_ba(f_nu,n0,f_tau,n1,n)
c   ------- Transform TAU ---> NU  --------------------------
c     Input:
c	F_tau - hermitian, but is stored in square R*8 form = A
c     Output:
c     F_nu is hermitian and stored in A square R*8 form
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n,n0,n1
	real*8, intent(in) :: f_tau(n1,n1,ndim3_tau)
	real*8, intent(out) :: f_nu(n0,n0,ndim3_nu)
      integer :: nn,i,j,ij,ind_tau,nblock,nstep,ind_nu,j0,jm,i_nu,i1,m,
     &           m0
      real*8 :: yy1
	real*8, allocatable :: at(:,:),an(:,:)
	nblock=1000
	nstep=max(1,nblock/n)
	nn=n*nstep
	allocate(at(nn,ndim3_tau))
	allocate(an(nn,0:n_nu_d))
	j0=0
      m=interp_nu_d
      m0=m-1
	do while(j0<n)
	  jm=min(nstep,n-j0)
c -------------------------- COS transformation ------------------
	  at=0.d0
	  do ind_tau=1,ndim3_tau
	    ij=0
	    do j=j0+1,j0+jm
	      do i=1,n
	        ij=ij+1
	        at(ij,ind_tau)=f_tau(i,j,ind_tau)
	      enddo
	    enddo
	  enddo
	  call from_tau_to_nu_cos_d(an,at,nn)
	  ij=0
	  do j=j0+1,j0+jm
	    do i=1,n
	      ij=ij+1
	      do ind_nu=1,ndim3_nu
	        i_nu=me3_tau*ndim3_nu+ind_nu-1
	        if(i_nu<=n_nu_exa) then
	          f_nu(i,j,ind_nu)=an(ij,i_nu)
	        else if(i_nu<n_nu) then
	          i1=ind_nu_int(i_nu)
	          call fint(an(ij,i1:i1+m0),w_nu_d(i1),m,w_nu(i_nu),yy1)
	          f_nu(i,j,ind_nu)=yy1
	        else if(i_nu==n_nu) then
	          f_nu(i,j,ind_nu)=an(ij,n_nu_d)
	        endif
	      enddo
	    enddo
	  enddo
	  j0=j0+jm
	enddo
	deallocate(at,an)
      end
      
      subroutine from_tau_to_nu_ab_c(f_nu,n0,f_tau,n1,n)
c   ------- Transform TAU ---> NU representation --------------------------
c	F_tau - real symmetric, but is stored in square form
c     F(tau) ne F(Beta-tau)
c     F_nu is required in square complex*16 form
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n,n0,n1
	real*8, intent(in) :: f_tau(n1,n1,0:n_tau)
	complex*16, intent(out) :: f_nu(n0,n0,0:n_nu)
      integer :: nn,i,j,ij,i_tau,i1_tau,i_nu,ind_tau,i1,m,m0
      real*8 :: yy1,yy2
	real*8, allocatable :: at(:,:),bt(:,:),an(:,:),bn(:,:)
	nn=n*(n+1)/2
	allocate(at(nn,ndim3_tau))
	allocate(an(nn,0:n_nu_d))
	allocate(bt(nn,ndim3_tau))
	allocate(bn(nn,0:n_nu_d))
	do ind_tau=1,ndim3_tau
	  i_tau=me3_tau*ndim3_tau+ind_tau-1
	  i1_tau=n_tau-i_tau
	  ij=0
	  do j=1,n
	    do i=1,j
	      ij=ij+1
	      at(ij,ind_tau)=f_tau(i,j,i_tau)+f_tau(i,j,i1_tau)
	      bt(ij,ind_tau)=f_tau(i,j,i_tau)-f_tau(i,j,i1_tau)
	    enddo
	  enddo
	enddo
	call from_tau_to_nu_cos_d(an,at,nn)
	call from_tau_to_nu_sin_d(bn,bt,nn)
	f_nu=(0.d0,0.d0)
      m=interp_nu_d
      m0=m-1
	ij=0
	do j=1,n
	  do i=1,j
	    ij=ij+1
	    do i_nu=0,n_nu
	      if(i_nu<=n_nu_exa) then
	        f_nu(i,j,i_nu)=dcmplx(an(ij,i_nu),bn(ij,i_nu))
	      else if(i_nu<n_nu) then
	        i1=ind_nu_int(i_nu)
	        call fint(an(ij,i1:i1+m0),w_nu_d(i1),m,w_nu(i_nu),yy1)
	        call fint(bn(ij,i1:i1+m0),w_nu_d(i1),m,w_nu(i_nu),yy2)
	        f_nu(i,j,i_nu)=dcmplx(yy1,yy2)
	      else if(i_nu==n_nu) then
	        f_nu(i,j,i_nu)=dcmplx(an(ij,n_nu_d),bn(ij,n_nu_d))
	      endif
	      f_nu(j,i,i_nu)=f_nu(i,j,i_nu)
	    enddo
	  enddo
	enddo
	deallocate(at,bt,an,bn)
      end
      
      subroutine from_tau_to_nu_aa_r(f_nu,n0,f_tau,n1,n)
c   ------- Transform TAU ---> NU representation --------------------------
c	F_tau - real symmetric, but is stored in square form
c     F(tau) == F(Beta-tau)
c     F_nu is also real symmetric and it is required in square form
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n,n0,n1
	real*8, intent(in) :: f_tau(n1,n1,0:n_tau/2)
	real*8, intent(out) :: f_nu(n0,n0,0:n_nu)
      integer :: nn,i,j,ij,i_tau,i_nu,ind_nu,ind_tau,i1,m,m0
      real*8 :: yy1
	real*8, allocatable :: at(:,:),an(:,:)
	nn=n*(n+1)/2
	allocate(at(nn,ndim3_tau))
	allocate(an(nn,0:n_nu_d))
	do ind_tau=1,ndim3_tau
	  i_tau=me3_tau*ndim3_tau+ind_tau-1
	  ij=0
	  do j=1,n
	    do i=1,j
	      ij=ij+1
	      at(ij,ind_tau)=2.d0*f_tau(i,j,i_tau)
	    enddo
	  enddo
	enddo
	call from_tau_to_nu_cos_d(an,at,nn)
	f_nu=0.d0
      m=interp_nu_d
      m0=m-1
	ij=0
	do j=1,n
	  do i=1,j
	    ij=ij+1
	    do ind_nu=1,ndim3_nu
	      i_nu=me3_tau*ndim3_nu+ind_nu-1
	      if(i_nu<=n_nu_exa) then
	        f_nu(i,j,i_nu)=an(ij,i_nu)
	      else if(i_nu<n_nu) then
	        i1=ind_nu_int(i_nu)
	        call fint(an(ij,i1:i1+m0),w_nu_d(i1),m,w_nu(i_nu),yy1)
	        f_nu(i,j,i_nu)=yy1
	      else if(i_nu==n_nu) then
	        f_nu(i,j,i_nu)=an(ij,n_nu_d)
	      endif
	      f_nu(j,i,i_nu)=f_nu(i,j,i_nu)
	    enddo
	  enddo
	enddo
	deallocate(at,an)
	if(nproc_tau/=1) call dgop(f_nu,n0*n0*(n_nu+1),'  +',comm_pnk)
      end
      
      subroutine from_tau_to_nu_aa_s(f_nu,f_tau,n)
c   ------- Transform TAU ---> NU representation (sinus) -----------
c     F(tau) = - F(Beta-tau)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: f_tau(n,0:n_tau/2)
	real*8, intent(out) :: f_nu(n,0:n_nu)
      integer :: i,i_tau,i_nu,ind_nu,ind_tau,i1,m,m0
      real*8 :: yy1
	real*8, allocatable :: at(:,:),an(:,:)
	allocate(at(n,ndim3_tau))
	allocate(an(n,0:n_nu_d))
	do ind_tau=1,ndim3_tau
	  i_tau=me3_tau*ndim3_tau+ind_tau-1
	  at(:,ind_tau)=2.d0*f_tau(:,i_tau)
	enddo
	call from_tau_to_nu_sin_d(an,at,n)
	f_nu=0.d0
      m=interp_nu_d
      m0=m-1
	do i=1,n
	  do ind_nu=1,ndim3_nu
	    i_nu=me3_tau*ndim3_nu+ind_nu-1
	    if(i_nu<=n_nu_exa) then
	      f_nu(i,i_nu)=an(i,i_nu)
	    else if(i_nu<n_nu) then
	      i1=ind_nu_int(i_nu)
	      call fint(an(i,i1:i1+m0),w_nu_d(i1),m,w_nu(i_nu),yy1)
	      f_nu(i,i_nu)=yy1
	    else if(i_nu==n_nu) then
	      f_nu(i,i_nu)=an(i,n_nu_d)
	    endif
	  enddo
	enddo
	deallocate(at,an)
	if(nproc_tau/=1) call dgop(f_nu,n*(n_nu+1),'  +',comm_pnk)
      end
      
       subroutine from_tau_to_nu_full_seq(f_nu,f_tau,n,tr)
c   ------- Transform TAU ---> NU representation --------------------------
c   ---- We transform totally general complex F_TAU with N elements ------
c	F_tau - is stored in (t,beta-t) form
c     F(tau) ne F(Beta-tau)
c     F_nu is required in complex*16 [-n_nu(1) : n_nu(2)] form
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
      real*8, intent(in) :: tr(n_tau/2,0:n_nu_d,2)
	complex*16, intent(in) :: f_tau(n,0:n_tau)
	complex*16, intent(out) :: f_nu(n,0:n_nu,2)
      integer :: i_tau,i,i_nu,i1,m,m0,i1_tau
      real*8 :: yy1,yy2,yy3,yy4
	complex*16, allocatable :: at(:,:),bt(:,:),an(:,:),bn(:,:)
	allocate(at(n,n_tau/2))
	allocate(bt(n,n_tau/2))
	allocate(an(n,0:n_nu_d))
	allocate(bn(n,0:n_nu_d))
      do i_tau=1,n_tau/2
        i1_tau=n_tau-i_tau
	  at(:,i_tau)=f_tau(:,i_tau)+f_tau(:,i1_tau)
	  bt(:,i_tau)=f_tau(:,i_tau)-f_tau(:,i1_tau)
	enddo
	call dgemm('n','n',2*n,n_nu_d+1,n_tau/2,1.d0,at,2*n,tr(1,0,1),
     &	         n_tau/2,0.d0,an,2*n)
	call dgemm('n','n',2*n,n_nu_d+1,n_tau/2,1.d0,bt,2*n,tr(1,0,2),
     &	         n_tau/2,0.d0,bn,2*n)
	f_nu=(0.d0,0.d0)
      m=interp_nu_d
      m0=m-1
	do i=1,n
	  do i_nu=0,n_nu
	    if(i_nu<=n_nu_exa) then
	      f_nu(i,i_nu,1)=an(i,i_nu)-(0.d0,1.d0)*bn(i,i_nu)
	      f_nu(i,i_nu,2)=an(i,i_nu)+(0.d0,1.d0)*bn(i,i_nu)
	    else if(i_nu<n_nu) then
	      i1=ind_nu_int(i_nu)
	      call fint(real(an(i,i1:i1+m0)),w_nu_d(i1),m,w_nu(i_nu),yy1)
	      call fint(real(bn(i,i1:i1+m0)),w_nu_d(i1),m,w_nu(i_nu),yy2)
	      call fint(imag(an(i,i1:i1+m0)),w_nu_d(i1),m,w_nu(i_nu),yy3)
	      call fint(imag(bn(i,i1:i1+m0)),w_nu_d(i1),m,w_nu(i_nu),yy4)
	      f_nu(i,i_nu,1)=dcmplx(yy1,yy3)-(0.d0,1.d0)*dcmplx(yy2,yy4)
	      f_nu(i,i_nu,2)=dcmplx(yy1,yy3)+(0.d0,1.d0)*dcmplx(yy2,yy4)
	    else if(i_nu==n_nu) then
	      f_nu(i,i_nu,1)=an(i,n_nu_d)-(0.d0,1.d0)*bn(i,n_nu_d)
	      f_nu(i,i_nu,2)=an(i,n_nu_d)+(0.d0,1.d0)*bn(i,n_nu_d)
	    endif
	  enddo
	enddo
	deallocate(at,bt,an,bn)
      end