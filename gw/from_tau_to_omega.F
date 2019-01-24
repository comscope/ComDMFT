      subroutine from_tau_to_omega(f_omega,f_tau,n)
c   ------- Transform TAU ---> OMEGA  --------------------------
c     Input:
c	F_tau - Real General, and is stored in N elements
c     Output:
c     F_omega is complex*16 and stored in C*16 N elements
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: f_tau(n,0:n_tau)
	complex*16, intent(out) :: f_omega(n,0:n_omega)
      integer :: i1_tau,i_tau,i_omega,ind_tau,i1,i,m,m0
      real*8 :: yy1,yy2
	real*8, allocatable :: at(:,:),ao(:,:),bt(:,:),bo(:,:)
c -------------------------- COS + SIN transformations ------------------
	allocate(at(n,ndim3_tau))
	allocate(ao(n,0:n_omega_d))
	allocate(bt(n,ndim3_tau))
	allocate(bo(n,0:n_omega_d))
	at=0.d0
	bt=0.d0
	do ind_tau=1,ndim3_tau
	  i_tau=me3_tau*ndim3_tau+ind_tau-1
	  i1_tau=n_tau-i_tau
	  at(:,ind_tau)=f_tau(:,i_tau)-f_tau(:,i1_tau)
	  bt(:,ind_tau)=f_tau(:,i_tau)+f_tau(:,i1_tau)
	enddo
	call from_tau_to_omega_cos_d(ao,at,n)
	call from_tau_to_omega_sin_d(bo,bt,n)
      m=interp_omega_d
      m0=m-1
	do i_omega=0,n_omega
	  if(i_omega<=n_omega_exa) then
	    f_omega(:,i_omega)=dcmplx(ao(:,i_omega),bo(:,i_omega))
	  else if(i_omega<n_omega) then
	    i1=ind_omega_int(i_omega)
	    do i=1,n
	      call fint(ao(i,i1:i1+m0),w_omega_d(i1),m,w_omega(i_omega),
     &                yy1)
	      call fint(bo(i,i1:i1+m0),w_omega_d(i1),m,w_omega(i_omega),
     &                yy2)
	      f_omega(i,i_omega)=dcmplx(yy1,yy2)
	    enddo
	  else if(i_omega==n_omega) then
	    f_omega(:,i_omega)=dcmplx(ao(:,n_omega_d),bo(:,n_omega_d))
	  endif
	enddo
	deallocate(at,ao,bt,bo)
      end
      
      subroutine from_tau_to_omega_vrt_ferm(f_omega,f_tau,n)
c   ------- Transform TAU ---> OMEGA  --------------------------
c     Input:
c	F_tau - C*16, N elements in (t;b-t) form
c     Output:
c     F_omega - C*16, N elemnents
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	complex*16, intent(in) :: f_tau(n,2,ndim3_tau)
	complex*16, intent(out) :: f_omega(n,ndim3_omega,2)
      integer :: ind_tau,i_omega,i1,i,m,m0,ind_omega
      real*8 :: yy1,yy2,yy1i,yy2i
      complex*16 :: a00,b00
	complex*16, allocatable :: at(:,:),ao(:,:),bt(:,:),bo(:,:)
	allocate(at(n,ndim3_tau))
	allocate(ao(n,0:n_omega_d))
	allocate(bt(n,ndim3_tau))
	allocate(bo(n,0:n_omega_d))
	do ind_tau=1,ndim3_tau
	  at(:,ind_tau)=f_tau(:,1,ind_tau)-f_tau(:,2,ind_tau)
	  bt(:,ind_tau)=(0.d0,1.d0)*(f_tau(:,1,ind_tau)
     &	                          +f_tau(:,2,ind_tau))
	enddo
	call from_tau_to_omega_cos_d(ao,at,2*n)
	call from_tau_to_omega_sin_d(bo,bt,2*n)
      m=interp_omega_d
      m0=m-1
	do ind_omega=1,ndim3_omega
	  i_omega=me3_tau*ndim3_omega+ind_omega-1
	  if(i_omega<=n_omega_exa) then
	    f_omega(:,ind_omega,1)=ao(:,i_omega)-bo(:,i_omega)
	    f_omega(:,ind_omega,2)=ao(:,i_omega)+bo(:,i_omega)
	  else if(i_omega<n_omega) then
	    i1=ind_omega_int(i_omega)
	    do i=1,n
	      call fint(real(ao(i,i1:i1+m0)),w_omega_d(i1),m,
     &	              w_omega(i_omega),yy1)
	      call fint(real(bo(i,i1:i1+m0)),w_omega_d(i1),m,
     &	              w_omega(i_omega),yy2)
	      call fint(imag(ao(i,i1:i1+m0)),w_omega_d(i1),m,
     &	              w_omega(i_omega),yy1i)
	      call fint(imag(bo(i,i1:i1+m0)),w_omega_d(i1),m,
     &	              w_omega(i_omega),yy2i)
	      a00=dcmplx(yy1,yy1i)
	      b00=dcmplx(yy2,yy2i)
	      f_omega(i,ind_omega,1)=a00-b00
	      f_omega(i,ind_omega,2)=a00+b00
	    enddo
	  else if(i_omega==n_omega) then
	    f_omega(:,ind_omega,1)=ao(:,n_omega_d)-bo(:,n_omega_d)
	    f_omega(:,ind_omega,2)=ao(:,n_omega_d)+bo(:,n_omega_d)
	  endif
	enddo
	deallocate(at,ao,bt,bo)
      end
      
      subroutine from_tau_to_omega_sin_d(f_omega,f_tau,n)
c   ------- Transform TAU ---> OMEGA representation -------------------
c   This is a besic operation B(w) = Int_0^beta/2 B(tau) sin(w*tau)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: f_tau(n,ndim3_tau)
	real*8, intent(out) :: f_omega(n,0:n_omega_d)
	call dgemm('n','n',n,n_omega_d+1,ndim3_tau,1.d0,f_tau,n,
     &	         ab_omega_from_tau(1,0,2),ndim3_tau,0.d0,f_omega,n)
	if(nproc_tau/=1) call dgop(f_omega,n*(n_omega_d+1),'  +',comm_pnk)
      end
      
      subroutine from_tau_to_omega_sin(f_omega,f_tau,n)
c   ------- Transform TAU ---> OMEGA representation -------------------
c   This is a besic operation B(w) = Int_0^beta/2 B(tau) sin(w*tau)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: f_tau(n,0:n_tau/2)
	real*8, intent(out) :: f_omega(n,ndim3_omega)
	call dgemm('n','n',n,ndim3_omega,n_tau/2+1,1.d0,f_tau(1,0),n,
     &	         ab_omega_from_tau(0,1,2),n_tau/2+1,0.d0,f_omega,n)
      end
      
      subroutine from_tau_to_omega_par(f_omega,f_tau,n)
c   ------- Transform TAU ---> OMEGA  --------------------------
c     Input:
c	F_tau - Real General, and is stored in N elements in (t;b-t)
c     Output:
c     F_omega is complex*16 and stored in C*16 N elements
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: f_tau(n,2,ndim3_tau)
	complex*16, intent(out) :: f_omega(n,ndim3_omega)
      integer :: ind_omega,ind_tau,i_omega,i1,i,m,m0
      real*8 :: yy1,yy2
	real*8, allocatable :: at(:,:),ao(:,:),bt(:,:),bo(:,:)
	allocate(at(n,ndim3_tau))
	allocate(ao(n,0:n_omega_d))
	allocate(bt(n,ndim3_tau))
	allocate(bo(n,0:n_omega_d))
c -------------------------- COS + SIN transformations ------------------
	at=0.d0
	bt=0.d0
	do ind_tau=1,ndim3_tau
	  at(:,ind_tau)=f_tau(:,1,ind_tau)-f_tau(:,2,ind_tau)
	  bt(:,ind_tau)=f_tau(:,1,ind_tau)+f_tau(:,2,ind_tau)
	enddo
	call from_tau_to_omega_cos_d(ao,at,n)
	call from_tau_to_omega_sin_d(bo,bt,n)
      m=interp_omega_d
      m0=m-1
	do ind_omega=1,ndim3_omega
	  i_omega=me3_tau*ndim3_omega+ind_omega-1
	  if(i_omega<=n_omega_exa) then
	    f_omega(:,ind_omega)=dcmplx(ao(:,i_omega),bo(:,i_omega))
	  else if(i_omega<n_omega) then
	    i1=ind_omega_int(i_omega)
	    do i=1,n
	      call fint(ao(i,i1:i1+m0),w_omega_d(i1),m,w_omega(i_omega),
     &                yy1)
	      call fint(bo(i,i1:i1+m0),w_omega_d(i1),m,w_omega(i_omega),
     &                yy2)
	      f_omega(i,ind_omega)=dcmplx(yy1,yy2)
	    enddo
	  else if(i_omega==n_omega) then
	    f_omega(:,ind_omega)=dcmplx(ao(:,n_omega_d),bo(:,n_omega_d))
	  endif
	enddo
	deallocate(at,ao,bt,bo)
      end
      
      subroutine from_tau_to_omega_d(i_omega,f_omega,f_tau,n)
c   ------- Transform TAU ---> specific OMEGA_D  ----------------------
c     Input:
c	F_tau - Real General, and is stored in N elements (t;b-t)
c     Output:
c     F_omega is complex*16 and stored in C*16 N elements
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n,i_omega
	real*8, intent(in) :: f_tau(n,2,ndim3_tau)
	complex*16, intent(out) :: f_omega(n)
      integer :: i_tau,ind_tau
	real*8, allocatable :: at(:,:),ao(:),bt(:,:),bo(:)
c -------------------------- COS + SIN transformations ------------------
	allocate(at(n,ndim3_tau))
	allocate(ao(n))
	allocate(bt(n,ndim3_tau))
	allocate(bo(n))
	at=0.d0
	bt=0.d0
	do ind_tau=1,ndim3_tau
	  i_tau=me3_tau*ndim3_tau+ind_tau-1
	  at(:,ind_tau)=f_tau(:,1,ind_tau)-f_tau(:,2,ind_tau)
	  bt(:,ind_tau)=f_tau(:,1,ind_tau)+f_tau(:,2,ind_tau)
	enddo
	call dgemm('n','n',n,1,ndim3_tau,1.d0,at,n,
     &	         ab_omega_from_tau(1,i_omega,1),ndim3_tau,0.d0,ao,n)
	call dgemm('n','n',n,1,ndim3_tau,1.d0,bt,n,
     &	         ab_omega_from_tau(1,i_omega,2),ndim3_tau,0.d0,bo,n)
	f_omega(:)=dcmplx(ao(:),bo(:))
	if(nproc_tau/=1) call dgop(f_omega,2*n,'  +',comm_pnk)
	deallocate(at,ao,bt,bo)
      end
      
      subroutine from_tau_to_omega_cos_d(f_omega,f_tau,n)
c   ------- Transform TAU ---> OMEGA representation -------------------
c   This is a besic operation A(w) = Int_0^beta/2 A(tau) cos(w*tau)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: f_tau(n,ndim3_tau)
	real*8, intent(out) :: f_omega(n,0:n_omega_d)
	call dgemm('n','n',n,n_omega_d+1,ndim3_tau,1.d0,f_tau,n,
     &	         ab_omega_from_tau(1,0,1),ndim3_tau,0.d0,f_omega,n)
	if(nproc_tau/=1) call dgop(f_omega,n*(n_omega_d+1),'  +',comm_pnk)
      end
      
      subroutine from_tau_to_omega_baa(f_omega,n0,f_tau,n1,n)
c   ------- Transform TAU ---> OMEGA  --------------------------
c     Input:
c	F_tau - hermitian, but is stored in square form = A;B
c     Output:
c     F_omega is complex*16 and stored in A,B form = (Re,Im)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n,n0,n1
	real*8, intent(in) :: f_tau(n1,n1,2,ndim3_tau)
	real*8, intent(out) :: f_omega(n0,n0,2,ndim3_omega)
      integer :: nn,i,j,ij,ind_tau,nblock,nstep,ind_omega,j0,jm,it,
     &           i_omega,i1,m,m0
      real*8 :: yy1
	real*8, allocatable :: at(:,:),ao(:,:)
	nblock=1000
	nstep=max(1,nblock/n)
	nn=n*nstep
	allocate(at(nn,ndim3_tau))
	allocate(ao(nn,0:n_omega_d))
      m=interp_omega_d
      m0=m-1
	j0=0
	do while(j0<n)
	  jm=min(nstep,n-j0)
c -------------------------- COS + SIN transformations ------------------
        do it=1,2
	    at=0.d0
	    do ind_tau=1,ndim3_tau
	      ij=0
	      do j=j0+1,j0+jm
	        do i=1,n
	          ij=ij+1
	          at(ij,ind_tau)=f_tau(i,j,it,ind_tau)
	        enddo
	      enddo
	    enddo
	    if(it==1) call from_tau_to_omega_cos_d(ao,at,nn)
	    if(it==2) call from_tau_to_omega_sin_d(ao,at,nn)
	    ij=0
	    do j=j0+1,j0+jm
	      do i=1,n
	        ij=ij+1
	        do ind_omega=1,ndim3_omega
	          i_omega=me3_tau*ndim3_omega+ind_omega-1
	          if(i_omega<=n_omega_exa) then
	            f_omega(i,j,it,ind_omega)=ao(ij,i_omega)
	          else if(i_omega<n_omega) then
	            i1=ind_omega_int(i_omega)
	            call fint(ao(ij,i1:i1+m0),w_omega_d(i1),m,
     &                      w_omega(i_omega),yy1)
	            f_omega(i,j,it,ind_omega)=yy1
	          else if(i_omega==n_omega) then
	            f_omega(i,j,it,ind_omega)=ao(ij,n_omega_d)
	          endif
	        enddo
	      enddo
	    enddo
	  enddo   !! over it
	  j0=j0+jm
	enddo
	deallocate(at,ao)
      end
      
      subroutine from_tau_to_omega_abb(f_omega,n0,f_tau,n1,n)
c   ------- Transform TAU ---> OMEGA  --------------------------
c     Input:
c	F_tau - NonSymmetric, and is stored in square R*8 (t;beta-t) form
c     Output:
c     F_omega is complex*16 and stored in C*16 intact form
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n,n0,n1
	real*8, intent(in) :: f_tau(n1,n1,2,ndim3_tau)
	complex*16, intent(out) :: f_omega(n0,n0,ndim3_omega)
      integer :: nn,i,j,ij,ind_tau,ind_omega,i_omega,i1,m,m0
      real*8 :: yy1,yy2
	real*8, allocatable :: at(:,:),ao(:,:),bt(:,:),bo(:,:)
	nn=n*n
	allocate(at(nn,ndim3_tau))
	allocate(ao(nn,0:n_omega_d))
	allocate(bt(nn,ndim3_tau))
	allocate(bo(nn,0:n_omega_d))
	do ind_tau=1,ndim3_tau
	  ij=0
	  do j=1,n
	    do i=1,n
	      ij=ij+1
	      at(ij,ind_tau)=f_tau(i,j,1,ind_tau)-f_tau(i,j,2,ind_tau)
	      bt(ij,ind_tau)=f_tau(i,j,1,ind_tau)+f_tau(i,j,2,ind_tau)
	    enddo
	  enddo
	enddo
	call from_tau_to_omega_cos_d(ao,at,nn)
	call from_tau_to_omega_sin_d(bo,bt,nn)
      m=interp_omega_d
      m0=m-1
	ij=0
	do j=1,n
	  do i=1,n
	    ij=ij+1
	    do ind_omega=1,ndim3_omega
	      i_omega=me3_tau*ndim3_omega+ind_omega-1
	      if(i_omega<=n_omega_exa) then
	        f_omega(i,j,ind_omega)=dcmplx(ao(ij,i_omega),
     &                                      bo(ij,i_omega))
	      else if(i_omega<n_omega) then
	        i1=ind_omega_int(i_omega)
	        call fint(ao(ij,i1:i1+m0),w_omega_d(i1),m,
     &                  w_omega(i_omega),yy1)
	        call fint(bo(ij,i1:i1+m0),w_omega_d(i1),m,
     &                  w_omega(i_omega),yy2)
	        f_omega(i,j,ind_omega)=dcmplx(yy1,yy2)
	      else if(i_omega==n_omega) then
	        f_omega(i,j,ind_omega)=dcmplx(ao(ij,n_omega_d),
     &                                      bo(ij,n_omega_d))
	      endif
	    enddo
	  enddo
	enddo
	deallocate(at,ao,bt,bo)
      end
      
      subroutine from_tau_to_omega_vrt(f_omega,f_tau,n)
c   ------- Transform TAU ---> OMEGA  --------------------------
c     Input:
c	F_tau - C*16, N elements in (t;b-t) form
c     Output:
c     F_omega - C*16, N elemnents
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	complex*16, intent(in) :: f_tau(n,2,ndim3_tau)
	complex*16, intent(out) :: f_omega(n,0:n_omega)
      integer :: ind_tau,i_omega,i1,i,m,m0
      real*8 :: yy1,yy2,yy1i,yy2i
      complex*16 :: a00,b00
	complex*16, allocatable :: at(:,:),ao(:,:),bt(:,:),bo(:,:)
	allocate(at(n,ndim3_tau))
	allocate(ao(n,0:n_omega_d))
	allocate(bt(n,ndim3_tau))
	allocate(bo(n,0:n_omega_d))
	do ind_tau=1,ndim3_tau
	  at(:,ind_tau)=f_tau(:,1,ind_tau)-f_tau(:,2,ind_tau)
	  bt(:,ind_tau)=(0.d0,1.d0)*(f_tau(:,1,ind_tau)
     &	                          +f_tau(:,2,ind_tau))
	enddo
	call from_tau_to_omega_cos_d(ao,at,2*n)
	call from_tau_to_omega_sin_d(bo,bt,2*n)
      m=interp_omega_d
      m0=m-1
	do i_omega=0,n_omega
	  if(i_omega<=n_omega_exa) then
	    f_omega(:,i_omega)=ao(:,i_omega)+bo(:,i_omega)
	  else if(i_omega<n_omega) then
	    i1=ind_omega_int(i_omega)
	    do i=1,n
	      call fint(real(ao(i,i1:i1+m0)),w_omega_d(i1),m,
     &	              w_omega(i_omega),yy1)
	      call fint(real(bo(i,i1:i1+m0)),w_omega_d(i1),m,
     &	              w_omega(i_omega),yy2)
	      call fint(imag(ao(i,i1:i1+m0)),w_omega_d(i1),m,
     &	              w_omega(i_omega),yy1i)
	      call fint(imag(bo(i,i1:i1+m0)),w_omega_d(i1),m,
     &	              w_omega(i_omega),yy2i)
	      a00=dcmplx(yy1,yy1i)
	      b00=dcmplx(yy2,yy2i)
	      f_omega(i,i_omega)=a00+b00
	    enddo
	  else if(i_omega==n_omega) then
	    f_omega(:,i_omega)=ao(:,n_omega_d)+bo(:,n_omega_d)
	  endif
	enddo
	deallocate(at,ao,bt,bo)
      end
      
      subroutine from_tau_to_omega_seq(f_omega,f_tau,n,tran)
c   ------- Transform TAU ---> OMEGA  --------------------------
c     Input:
c	F_tau - C*16, N elements in (0:n_tau) form
c     Output:
c     F_omega - C*16, N elemnents (-+ omega)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: tran(n_tau/2,0:n_omega_d,2)
	complex*16, intent(in) :: f_tau(n,0:n_tau)
	complex*16, intent(out) :: f_omega(n,0:n_omega,2)
      integer :: i_tau,i_omega,i1,i,m,m0,i1_tau
      real*8 :: yy1,yy2,yy1i,yy2i
      complex*16 :: a00,b00
	complex*16, allocatable :: at(:,:),ao(:,:),bt(:,:),bo(:,:)
	allocate(at(n,n_tau/2))
	allocate(ao(n,0:n_omega_d))
	allocate(bt(n,n_tau/2))
	allocate(bo(n,0:n_omega_d))
	do i_tau=1,n_tau/2
	  i1_tau=n_tau-i_tau
	  at(:,i_tau)=f_tau(:,i_tau)-f_tau(:,i1_tau)
	  bt(:,i_tau)=(0.d0,1.d0)*(f_tau(:,i_tau)+f_tau(:,i1_tau))
	enddo
	call dgemm('n','n',2*n,n_omega_d+1,n_tau/2,1.d0,at,2*n,
     &	         tran(1,0,1),n_tau/2,0.d0,ao,2*n)
	call dgemm('n','n',2*n,n_omega_d+1,n_tau/2,1.d0,bt,2*n,
     &	         tran(1,0,2),n_tau/2,0.d0,bo,2*n)
      m=interp_omega_d
      m0=m-1
	do i_omega=0,n_omega
	  if(i_omega<=n_omega_exa) then
	    f_omega(:,i_omega,1)=ao(:,i_omega)-bo(:,i_omega)
	    f_omega(:,i_omega,2)=ao(:,i_omega)+bo(:,i_omega)
	  else if(i_omega<n_omega) then
	    i1=ind_omega_int(i_omega)
	    do i=1,n
	      call fint(real(ao(i,i1:i1+m0)),w_omega_d(i1),m,
     &	              w_omega(i_omega),yy1)
	      call fint(real(bo(i,i1:i1+m0)),w_omega_d(i1),m,
     &	              w_omega(i_omega),yy2)
	      call fint(imag(ao(i,i1:i1+m0)),w_omega_d(i1),m,
     &	              w_omega(i_omega),yy1i)
	      call fint(imag(bo(i,i1:i1+m0)),w_omega_d(i1),m,
     &	              w_omega(i_omega),yy2i)
	      a00=dcmplx(yy1,yy1i)
	      b00=dcmplx(yy2,yy2i)
	      f_omega(i,i_omega,1)=a00-b00
	      f_omega(i,i_omega,2)=a00+b00
	    enddo
	  else if(i_omega==n_omega) then
	    f_omega(:,i_omega,1)=ao(:,n_omega_d)-bo(:,n_omega_d)
	    f_omega(:,i_omega,2)=ao(:,n_omega_d)+bo(:,n_omega_d)
	  endif
	enddo
	deallocate(at,ao,bt,bo)
      end