      subroutine from_omega_to_tau(n,f_omega,f_tau)
c   ------- Transform OMEGA ---> TAU representation -------------
c     F_OMEGA - C*16 with N*N elements
c     F_TAU   - R*8 with N*N elements in (t;beta-t) form
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	complex*16, intent(in) :: f_omega(n,n,ndim3_omega)
	real*8, intent(out) :: f_tau(n,n,2,ndim3_tau)
      integer :: ind_tau,ind_omega
	real*8, allocatable :: at(:,:,:,:),ao(:,:,:,:)
	allocate(at(n,n,2,ndim3_tau))
	allocate(ao(n,n,2,ndim3_omega))
	do ind_omega=1,ndim3_omega
	  ao(:,:,1,ind_omega)=real(f_omega(:,:,ind_omega))
	  ao(:,:,2,ind_omega)=imag(f_omega(:,:,ind_omega))
	enddo
	call from_omega_to_tau_baa(ao,n,at,n,n)
	do ind_tau=1,ndim3_tau
	  f_tau(:,:,1,ind_tau)=0.5d0*(at(:,:,1,ind_tau)+at(:,:,2,ind_tau))
	  f_tau(:,:,2,ind_tau)=0.5d0*(at(:,:,2,ind_tau)-at(:,:,1,ind_tau))
	enddo
      deallocate(ao,at)
      end
      
      subroutine from_omega_to_tau_v(n,f_omega,f_tau)
c   ------- Transform OMEGA ---> TAU representation ----------------
c     F(-w)=F*(w)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	complex*16, intent(in) :: f_omega(n,0:n_omega)
	real*8, intent(out) :: f_tau(n,0:n_tau)
      integer :: ind_tau,i_tau
	real*8, allocatable :: at(:,:),bt(:,:),ao(:,:)
	allocate(at(n,ndim3_tau))
	allocate(ao(n,0:n_omega))
	allocate(bt(n,ndim3_tau))
      ao=real(f_omega)
	call from_omega_to_tau_cos(n,ao,at)
      ao=imag(f_omega)
	call from_omega_to_tau_sin(n,ao,bt)
	f_tau=0.d0
	do ind_tau=1,ndim3_tau
	  i_tau=me3_tau*ndim3_tau+ind_tau-1
	  f_tau(:,i_tau)=0.5d0*(at(:,ind_tau)+bt(:,ind_tau))
	  f_tau(:,n_tau-i_tau)=0.5d0*(bt(:,ind_tau)-at(:,ind_tau))
	enddo
	if(nproc_tau/=1) call dgop(f_tau,n*(n_tau+1),'  +',comm_pnk)
      deallocate(ao,at,bt)
      end
      
      subroutine from_omega_to_tau_sin(n,b_omega,b_t)
c   ------- Transform B(OMEGA) ---> B(TAU)  --------------------------
c     This is the basic transform :
c     A(t) = 4/b Sum_w>0 sin(wt) A(w)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: b_omega(n,0:n_omega)
	real*8, intent(out) :: b_t(n,ndim3_tau)
	call dgemm('n','n',n,ndim3_tau,n_omega+1,1.d0,b_omega(1,0),
     &	         n,tau_from_omega(0,1,2),n_omega+1,0.d0,b_t,n)
      end
      
      subroutine from_omega_to_tau_cos(n,a_omega,a_tau)
c   ------- Transform A(OMEGA) ---> A(TAU)  --------------------------
c     This is the basic transform :
c     A(t) = 4/b Sum_w>0 cos(wt) A(w)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: a_omega(n,0:n_omega)
	real*8, intent(out) :: a_tau(n,ndim3_tau)
	call dgemm('n','n',n,ndim3_tau,n_omega+1,1.d0,a_omega(1,0),
     &	         n,tau_from_omega(0,1,1),n_omega+1,0.d0,a_tau,n)
      end
      
      subroutine from_omega_to_tau_baa_seq(f_omega,n0,f_tau,n1,n)
c   ------- Transform OMEGA ---> TAU  --------------------------
c     Input:
c     F_omega is complex*16 and stored in A,B form = (Re,Im)
c     Output:
c	F_tau - hermitian, but is stored in square form = A;B
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n,n0,n1
	real*8, intent(in) :: f_omega(n0,n0,2,0:n_omega)
	real*8, intent(out) :: f_tau(n1,n1,2,0:n_tau/2)
      integer :: nn,i,j,ij,ind_tau,i_omega,nblock,nstep,j0,jm,it,i_tau
	real*8, allocatable :: at(:,:),ao(:,:)
	nblock=1000
	nstep=max(1,nblock/n)
	nn=n*nstep
	allocate(at(nn,ndim3_tau))
	allocate(ao(nn,0:n_omega))
	f_tau=0.d0
	j0=0
	do while(j0<n)
	  jm=min(nstep,n-j0)
c -------------------------- COS + SIN transformations ------------------
        do it=1,2
	    ao=0.d0
	    do i_omega=0,n_omega
	      ij=0
	      do j=j0+1,j0+jm
	        do i=1,n
	          ij=ij+1
	          ao(ij,i_omega)=f_omega(i,j,it,i_omega)
	        enddo
	      enddo
	    enddo
	    if(it==1) call from_omega_to_tau_cos(nn,ao,at)
	    if(it==2) call from_omega_to_tau_sin(nn,ao,at)
	    do ind_tau=1,ndim3_tau
	      i_tau=me3_tau*ndim3_tau+ind_tau-1
	      ij=0
	      do j=j0+1,j0+jm
	        do i=1,n
	          ij=ij+1
	          f_tau(i,j,it,i_tau)=at(ij,ind_tau)
	        enddo
	      enddo
	    enddo
	  enddo   !! over it
	  j0=j0+jm
	enddo
	if(nproc_tau/=1) call dgop(f_tau,n1*n1*2*(n_tau/2+1),'  +',
     &	                         comm_pnk)
	deallocate(at,ao)
      end
      
      subroutine from_omega_to_tau_baa(f_omega,n0,f_tau,n1,n)
c   ------- Transform OMEGA ---> TAU  --------------------------
c     Input:
c     F_omega is complex*16 and stored in A,B form = (Re,Im)
c     Output:
c	F_tau - hermitian, but is stored in square form = A;B
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n,n0,n1
	real*8, intent(in) :: f_omega(n0,n0,2,ndim3_omega)
	real*8, intent(out) :: f_tau(n1,n1,2,ndim3_tau)
      integer :: nn,i,j,ij,ind_tau,i_omega,nblock,nstep,ind_omega,j0,jm,
     &           it
	real*8, allocatable :: at(:,:),ao(:,:)
	nblock=1000
	nstep=max(1,nblock/n)
c	nstep=min(n,nblock)
	nn=n*nstep
	allocate(at(nn,ndim3_tau))
	allocate(ao(nn,0:n_omega))
	j0=0
	do while(j0<n)
	  jm=min(nstep,n-j0)
c -------------------------- COS + SIN transformations ------------------
        do it=1,2
	    ao=0.d0
	    do ind_omega=1,ndim3_omega
	      i_omega=me3_tau*ndim3_omega+ind_omega-1
	      ij=0
	      do j=j0+1,j0+jm
	        do i=1,n
	          ij=ij+1
	          ao(ij,i_omega)=f_omega(i,j,it,ind_omega)
	        enddo
	      enddo
	    enddo
	    if(nproc_tau/=1) call dgop(ao,nn*(n_omega+1),'  +',comm_pnk)
	    if(it==1) call from_omega_to_tau_cos(nn,ao,at)
	    if(it==2) call from_omega_to_tau_sin(nn,ao,at)
	    do ind_tau=1,ndim3_tau
	      ij=0
	      do j=j0+1,j0+jm
	        do i=1,n
	          ij=ij+1
	          f_tau(i,j,it,ind_tau)=at(ij,ind_tau)
	        enddo
	      enddo
	    enddo
	  enddo   !! over it
	  j0=j0+jm
	enddo
	deallocate(at,ao)
      end
      
      subroutine from_omega_to_tau_all(n,f_omega,f_tau)
c   ------- Transform OMEGA ---> TAU representation -------------
c     F_OMEGA - C*16 with N elements
c     F_TAU   - R*8 with N elements in (t;beta-t) form
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	complex*16, intent(in) :: f_omega(n,ndim3_omega)
	real*8, intent(out) :: f_tau(n,2,ndim3_tau)
      integer :: ind_tau,ind_omega,i_omega
	real*8, allocatable :: at(:,:),ao(:,:)
	allocate(at(n,ndim3_tau))
	allocate(ao(n,0:n_omega))
	ao=0.d0
	do ind_omega=1,ndim3_omega
	  i_omega=me3_tau*ndim3_omega+ind_omega-1
	  ao(:,i_omega)=real(f_omega(:,ind_omega))
	enddo
	if(nproc_tau/=1) call dgop(ao,n*(n_omega+1),'  +',comm_pnk)
	call dgemm('n','n',n,ndim3_tau,n_omega+1,1.d0,ao(1,0),n,
     &	         tau_from_omega(0,1,1),n_omega+1,0.d0,at,n)
	do ind_tau=1,ndim3_tau
	  f_tau(:,1,ind_tau)=at(:,ind_tau)
	  f_tau(:,2,ind_tau)=-at(:,ind_tau)
	enddo
	ao=0.d0
	do ind_omega=1,ndim3_omega
	  i_omega=me3_tau*ndim3_omega+ind_omega-1
	  ao(:,i_omega)=imag(f_omega(:,ind_omega))
	enddo
	if(nproc_tau/=1) call dgop(ao,n*(n_omega+1),'  +',comm_pnk)
	call dgemm('n','n',n,ndim3_tau,n_omega+1,1.d0,ao(1,0),n,
     &	         tau_from_omega(0,1,2),n_omega+1,0.d0,at,n)
	do ind_tau=1,ndim3_tau
	  f_tau(:,1,ind_tau)=f_tau(:,1,ind_tau)+at(:,ind_tau)
	  f_tau(:,2,ind_tau)=f_tau(:,2,ind_tau)+at(:,ind_tau)
	enddo
      deallocate(ao,at)
      f_tau=f_tau/2.d0
      end
      
       subroutine from_omega_to_tau0(n,f_omega,f_tau)
c   ------- Transform F(OMEGA) ---> F(TAU=0)  --------------------------
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: f_omega(n,0:n_omega)
	real*8, intent(out) :: f_tau(n)
	call dgemm('n','n',n,1,n_omega+1,1.d0,f_omega(1,0),n,
     &	         tau_from_omega(0,1,1),n_omega+1,0.d0,f_tau,n)
      f_tau=f_tau/2.d0
      end
      
      subroutine from_omega_nu0_to_tau(n,f_omega,f_tau,omega_tau)
c   ------- Transform OMEGA ---> TAU representation ----------------
c     F(-w)=F*(w)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: omega_tau(0:n_omega,0:n_tau,2)
	complex*16, intent(in) :: f_omega(n,0:n_omega)
	complex*16, intent(out) :: f_tau(n,0:n_tau)
      integer :: i_tau
	real*8, allocatable :: at(:,:),bt(:,:),ao(:,:)
	allocate(at(n,0:n_tau))
	allocate(bt(n,0:n_tau))
	allocate(ao(n,0:n_omega))
      ao=real(f_omega)
	call dgemm('n','n',n,n_tau+1,n_omega+1,1.d0,ao,n,omega_tau(0,0,1),
     &	         n_omega+1,0.d0,at,n)
      ao=imag(f_omega)
	call dgemm('n','n',n,n_tau+1,n_omega+1,1.d0,ao,n,omega_tau(0,0,2),
     &	         n_omega+1,0.d0,bt,n)
	f_tau=0.d0
	do i_tau=0,n_tau/2
	  f_tau(:,i_tau)=0.5d0*(at(:,i_tau)+bt(:,i_tau))
	  if(i_tau==n_tau/2) cycle
	  f_tau(:,n_tau-i_tau)=0.5d0*(bt(:,i_tau)-at(:,i_tau))
	enddo
      deallocate(ao,at,bt)
      end