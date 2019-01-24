      subroutine from_nu_to_tau0(n,f_nu,f_tau)
c   ------- Transform F(NU) ---> F(TAU=0)  --------------------------
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: f_nu(n,0:n_nu)
	real*8, intent(out) :: f_tau(n)
	call dgemm('n','n',n,1,n_nu+1,1.d0,f_nu(1,0),n,
     &	         tau_from_nu(0,1,1),n_nu+1,0.d0,f_tau,n)
      f_tau=f_tau/2.d0
      end
      
      
      subroutine from_nu_to_tau_aa_g(f_nu,f_tau,n)
c   ------- Transform NU ---> TAU  --------------------------
c     Input:
c     F_nu is real nonsymmetric; stored in rectangular real*8 form
c     Output:
c	F_tau - real nonsymmetric; stored in rectangular form
c     F(tau) = F(Beta-tau)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: f_nu(n,ndim3_nu)
	real*8, intent(out) :: f_tau(n,ndim3_tau)
      integer :: ind_nu,ind_tau,i_nu,n0,nsh,nn,ii,i0,n1
	real*8, allocatable :: at(:,:),an(:,:)
	nsh=8192
	n0=min(n,nsh)
	nn=n/nsh
	if(mod(n,nsh)/=0) nn=nn+1
	allocate(at(nsh,ndim3_tau))
	allocate(an(nsh,0:n_nu))
	i0=0
	do ii=1,nn
	  n1=min(n-i0,nsh)
	  an=0.d0
	  do ind_nu=1,ndim3_nu
	    i_nu=me3_tau*ndim3_nu+ind_nu-1
	    an(1:n1,i_nu)=f_nu(i0+1:i0+n1,ind_nu)
	  enddo
	  if(nproc_tau/=1) call dgop(an,nsh*(n_nu+1),'  +',comm_pnk)
	  call from_nu_to_tau_cos(nsh,an,at)
	  do ind_tau=1,ndim3_tau
	    f_tau(i0+1:i0+n1,ind_tau)=0.5d0*at(1:n1,ind_tau)
	  enddo
	  i0=i0+n1
	enddo
	deallocate(at,an)
      end
      
      subroutine from_nu_to_tau_sin(n,b_nu,b_t)
c   ------- Transform B(NU) ---> B(TAU)  --------------------------
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: b_nu(n,0:n_nu)
	real*8, intent(out) :: b_t(n,ndim3_tau)
	call dgemm('n','n',n,ndim3_tau,n_nu+1,1.d0,b_nu(1,0),n,
     &	         tau_from_nu(0,1,2),n_nu+1,0.d0,b_t,n)
      end
      
      
      subroutine from_nu_to_tau_cos(n,a_nu,a_tau)
c   ------- Transform A(NU) ---> A(TAU)  --------------------------
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: a_nu(n,0:n_nu)
	real*8, intent(out) :: a_tau(n,ndim3_tau)
	call dgemm('n','n',n,ndim3_tau,n_nu+1,1.d0,a_nu(1,0),n,
     &	         tau_from_nu(0,1,1),n_nu+1,0.d0,a_tau,n)
      end
      
      
      subroutine from_nu_to_tau_cc(f_nu,f_tau,n)
c   ------- Transform NU ---> TAU  --------------------------
c     Input:
c     F_nu is given by its F(-v) and F(v) components (N elements)
c     Output:
c	F_tau - C*16 nonsymmetric; stored in N elements (t;beta-t)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	complex*16, intent(in) :: f_nu(n,ndim3_nu,2)
	complex*16, intent(out) :: f_tau(n,2,ndim3_tau)
      integer :: ind_nu,ind_tau,i_nu
	complex*16, allocatable :: at(:,:),an(:,:),bt(:,:),bn(:,:)
	allocate(at(n,ndim3_tau))
	allocate(an(n,0:n_nu))
	allocate(bt(n,ndim3_tau))
	allocate(bn(n,0:n_nu))
	an=(0.d0,0.d0)
	bn=(0.d0,0.d0)
	do ind_nu=1,ndim3_nu
	  i_nu=me3_tau*ndim3_nu+ind_nu-1
	  an(:,i_nu)=0.5d0*(f_nu(:,ind_nu,1)+f_nu(:,ind_nu,2))
	  bn(:,i_nu)=(0.d0,0.5d0)*(f_nu(:,ind_nu,1)-f_nu(:,ind_nu,2))
	enddo
	if(nproc_tau/=1) then
	  call dgop(an,2*n*(n_nu+1),'  +',comm_pnk)
	  call dgop(bn,2*n*(n_nu+1),'  +',comm_pnk)
	endif
	call from_nu_to_tau_cos(2*n,an,at)
	call from_nu_to_tau_sin(2*n,bn,bt)
	do ind_tau=1,ndim3_tau
	  f_tau(:,1,ind_tau)=0.5d0*(at(:,ind_tau)+bt(:,ind_tau))
	  f_tau(:,2,ind_tau)=0.5d0*(at(:,ind_tau)-bt(:,ind_tau))
	enddo
	deallocate(at,an,bt,bn)
      end
      
      
      subroutine from_nu_to_tau_c(f_nu,f_tau,n)
c   ------- Transform NU ---> TAU  --------------------------
c     Input:
c     F_nu is C*16 nonsymmetric; stored in N elements [0:n_nu]
c     Output:
c	F_tau - R*8 nonsymmetric; stored in N elements (t;beta-t)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	complex*16, intent(in) :: f_nu(n,ndim3_nu)
	real*8, intent(out) :: f_tau(n,2,ndim3_tau)
      integer :: ind_nu,ind_tau,i_nu
	real*8, allocatable :: at(:,:),an(:,:),bt(:,:),bn(:,:)
	allocate(at(n,ndim3_tau))
	allocate(an(n,0:n_nu))
	allocate(bt(n,ndim3_tau))
	allocate(bn(n,0:n_nu))
	an=0.d0
	bn=0.d0
	do ind_nu=1,ndim3_nu
	  i_nu=me3_tau*ndim3_nu+ind_nu-1
	  an(:,i_nu)=real(f_nu(:,ind_nu))
	  bn(:,i_nu)=imag(f_nu(:,ind_nu))
	enddo
	if(nproc_tau/=1) then
	  call dgop(an,n*(n_nu+1),'  +',comm_pnk)
	  call dgop(bn,n*(n_nu+1),'  +',comm_pnk)
	endif
	call from_nu_to_tau_cos(n,an,at)
	call from_nu_to_tau_sin(n,bn,bt)
	do ind_tau=1,ndim3_tau
	  f_tau(:,1,ind_tau)=0.5d0*(at(:,ind_tau)+bt(:,ind_tau))
	  f_tau(:,2,ind_tau)=0.5d0*(at(:,ind_tau)-bt(:,ind_tau))
	enddo
	deallocate(at,an,bt,bn)
      end
      
      subroutine from_nu_to_tau_ba(f_nu,n0,f_tau,n1,n)
c   ------- Transform NU ---> TAU  --------------------------
c     Input:
c     F_nu is Hermitian and stored in R*8 square A form
c     Output:
c	F_tau - hermitian, but is stored in square form = A
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n,n0,n1
	real*8, intent(in) :: f_nu(n0,n0,ndim3_nu)
	real*8, intent(out) :: f_tau(n1,n1,ndim3_tau)
      integer :: nn,i,j,ij,ind_tau,i_nu,nblock,nstep,ind_nu,j0,jm
	real*8, allocatable :: at(:,:),an(:,:)
	nblock=1000
	nstep=max(1,nblock/n)
	nn=n*nstep
	allocate(at(nn,ndim3_tau))
	allocate(an(nn,0:n_nu))
	j0=0
	do while(j0<n)
	  jm=min(nstep,n-j0)
c -------------------------- COS transformation ------------------
	  an=0.d0
	  do ind_nu=1,ndim3_nu
	    i_nu=me3_tau*ndim3_nu+ind_nu-1
	    ij=0
	    do j=j0+1,j0+jm
	      do i=1,n
	        ij=ij+1
	        an(ij,i_nu)=f_nu(i,j,ind_nu)
	      enddo
	    enddo
	  enddo
	  if(nproc_tau/=1) call dgop(an,nn*(n_nu+1),'  +',comm_pnk)
	  call from_nu_to_tau_cos(nn,an,at)
	  do ind_tau=1,ndim3_tau
	    ij=0
	    do j=j0+1,j0+jm
	      do i=1,n
	        ij=ij+1
	        f_tau(i,j,ind_tau)=at(ij,ind_tau)
	      enddo
	    enddo
	  enddo
	  j0=j0+jm
	enddo
	deallocate(at,an)
      end
      
       subroutine from_nu_to_tau_aa_r(f_nu,n0,f_tau,n1,n)
c   ------- Transform NU ---> TAU  --------------------------
c     Input:
c     F_nu is real symmetric but stored in square real*8 form
c     Output:
c	F_tau - real symmetric, but is stored in square form
c     F(tau) = F(Beta-tau)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n,n0,n1
	real*8, intent(in) :: f_nu(n0,n0,0:n_nu)
	real*8, intent(out) :: f_tau(n1,n1,0:n_tau/2)
      integer :: nn,i,j,ij,ind_tau,i_nu,i_tau
	real*8, allocatable :: at(:,:),an(:,:)
	nn=n*(n+1)/2
	allocate(at(nn,ndim3_tau))
	allocate(an(nn,0:n_nu))
	do i_nu=0,n_nu
	  ij=0
	  do j=1,n
	    do i=1,j
	      ij=ij+1
	      an(ij,i_nu)=f_nu(i,j,i_nu)
	    enddo
	  enddo
	enddo
	call from_nu_to_tau_cos(nn,an,at)
	f_tau=0.d0
	do ind_tau=1,ndim3_tau
	  i_tau=me3_tau*ndim3_tau+ind_tau-1
	  ij=0
	  do j=1,n
	    do i=1,j
	      ij=ij+1
	      f_tau(i,j,i_tau)=0.5d0*at(ij,ind_tau)
	      f_tau(j,i,i_tau)=f_tau(i,j,i_tau)
	    enddo
	  enddo
	enddo
	if(nproc_tau/=1) call dgop(f_tau,n1**2*(n_tau/2+1),'  +',
     &	                         comm_pnk)
	deallocate(at,an)
      end
      
      
      subroutine from_nu_to_tau_cc_seq(f_nu,f_tau,n,trn)
c   ------- Transform NU ---> TAU  --------------------------
c     Input:
c     F_nu is given by its F(-v) and F(v) components (N elements)
c     Output:
c	F_tau - C*16 nonsymmetric; stored in N elements (t;beta-t)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
      real*8, intent(in) :: trn(0:n_nu,0:n_tau/2,2)
	complex*16, intent(in) :: f_nu(n,0:n_nu,2)
	complex*16, intent(out) :: f_tau(n,0:n_tau)
      integer :: i_tau,i_nu,i1_tau
	complex*16, allocatable :: at(:,:),an(:,:),bt(:,:),bn(:,:)
	allocate(at(n,0:n_tau/2))
	allocate(an(n,0:n_nu))
	allocate(bt(n,0:n_tau/2))
	allocate(bn(n,0:n_nu))
	an=(0.d0,0.d0)
	bn=(0.d0,0.d0)
	do i_nu=0,n_nu
	  an(:,i_nu)=0.5d0*(f_nu(:,i_nu,1)+f_nu(:,i_nu,2))
	  bn(:,i_nu)=(0.d0,0.5d0)*(f_nu(:,i_nu,1)-f_nu(:,i_nu,2))
	enddo
	call dgemm('n','n',2*n,n_tau/2+1,n_nu+1,1.d0,an,2*n,trn(0,0,1),
     &	         n_nu+1,0.d0,at,2*n)
	call dgemm('n','n',2*n,n_tau/2+1,n_nu+1,1.d0,bn,2*n,trn(0,0,2),
     &	         n_nu+1,0.d0,bt,2*n)
	do i_tau=0,n_tau/2
	  i1_tau=n_tau-i_tau
	  f_tau(:,i_tau)=0.5d0*(at(:,i_tau)+bt(:,i_tau))
	  f_tau(:,i1_tau)=0.5d0*(at(:,i_tau)-bt(:,i_tau))
	enddo
	deallocate(at,an,bt,bn)
      end