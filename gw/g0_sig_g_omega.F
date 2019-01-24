      subroutine g0_sig_g_omega(g0,s,g,key,n0,n)
c	key=0 - for G0 : (1 + G S) G0=G
c	key=1 - for G  : (1 - G0 S) G=G0
c	key=2 - for Sig:  G G0 S = G - G G0 G^(-1)
	use manager_mod
	use parallel_mod
      implicit none
	integer, intent(in) :: n0,n
	complex*16, intent(inout) :: g0(n0,n0,ndim3_omega),
     &	                           g(n0,n0,ndim3_omega),
     &                             s(n0,n0,ndim3_omega)
	integer :: info,i,key,j,ind_omega,lwork
	integer, allocatable :: ipiv(:)
	complex*16, allocatable :: work(:),a(:,:),b(:,:)
	lwork=3*n
	allocate(a(n,n),b(n,n),ipiv(n),work(lwork))
	if(key==0) then
c ---------   We are solving for G0 -------------------------
	  do ind_omega=1,ndim3_omega
	    call zgemm('n','n',n,n,n,(1.d0,0.d0),g(1,1,ind_omega),n0,
     &		       s(1,1,ind_omega),n0,(0.d0,0.d0),a,n)
	    do i=1,n
	      a(i,i)=a(i,i)+(1.d0,0.d0)
	    enddo
	    g0(1:n,1:n,ind_omega)=g(1:n,1:n,ind_omega)
	    call zgetrf(n,n,a,n,ipiv,info)
	    call zgetrs('n',n,n,a,n,ipiv,g0(1,1,ind_omega),n0,info)
	  enddo
	else if(key==1) then
c ---------   We are solving for G -------------------------
	  do ind_omega=1,ndim3_omega
	    call zgemm('n','n',n,n,n,(-1.d0,0.d0),g0(1,1,ind_omega),n0,
     &		       s(1,1,ind_omega),n0,(0.d0,0.d0),a,n)
	    do i=1,n
	      a(i,i)=a(i,i)+(1.d0,0.d0)
	    enddo
	    g(1:n,1:n,ind_omega)=g0(1:n,1:n,ind_omega)
	    call zgetrf(n,n,a,n,ipiv,info)
	    call zgetrs('n',n,n,a,n,ipiv,g(1,1,ind_omega),n0,info)
	  enddo
	else if(key==2) then
c ---------   We are solving for Sigma -------------------------
	  do ind_omega=1,ndim3_omega
	    call zgemm('n','n',n,n,n,(1.d0,0.d0),g(1,1,ind_omega),n0,
     &		       g0(1,1,ind_omega),n0,(0.d0,0.d0),a,n)
	    do i=1,n
	      do j=1,n
	        b(i,j)=a(j,i)
			s(i,j,ind_omega)=g(j,i,ind_omega)
	      enddo
	    enddo
	    call zgetrf(n,n,s(1,1,ind_omega),n0,ipiv,info)
	    call zgetrs('n',n,n,s(1,1,ind_omega),n0,ipiv,b,n,info)
	    do i=1,n
	      do j=1,n
			s(i,j,ind_omega)=g(i,j,ind_omega)-b(j,i)
	      enddo
	    enddo
	    call zgetrf(n,n,a,n,ipiv,info)
	    call zgetrs('n',n,n,a,n,ipiv,s(1,1,ind_omega),n0,info)
	  enddo
	endif
	deallocate(work,a,b,ipiv)
      end