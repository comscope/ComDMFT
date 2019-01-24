      subroutine tos_c(st,sd,ndim,gc_diag,dgc_diag)
	use manager_mod
	use units_mod
      implicit none
	integer, intent(in) :: ndim
	real*8, intent(in) :: gc_diag(ndim,0:n_omega),
     &	                  dgc_diag(ndim,0:n_omega)
	real*8, intent(inout) :: st,sd
	integer :: i,idim,ii,k
	real*8 :: w_n,dx,gcd,dgcd,b2
	real*8, allocatable :: aa(:,:),bb(:,:),cc(:,:),dd(:,:),h(:),y(:)
	allocate(aa(n_omega,2),bb(n_omega,2),cc(n_omega,2),dd(n_omega,2),
     &	     h(n_omega),y(0:n_omega))
	b2=2.d0/betta_t
	do i=1,n_omega
	  h(i)=w_omega(i)-w_omega(i-1)
	enddo
      do idim=1,ndim
	  y=gc_diag(idim,:)
	  call spline_start(h,y,aa(:,1),bb(:,1),cc(:,1),dd(:,1),n_omega)
	  y=dgc_diag(idim,:)
	  call spline_start(h,y,aa(:,2),bb(:,2),cc(:,2),dd(:,2),n_omega)
	  w_n=w_omega(0)
	  ii=1
	  do while(w_n<=omega_max)
	    do k=ii,n_nu
	      if(w_n>=w_omega(k-1).and.w_n<=w_omega(k)) then
	        ii=k
	        exit
	      endif
	    enddo
	    dx=w_n-w_omega(ii-1)
	    gcd=aa(ii,1)+dx*(bb(ii,1)+dx*(cc(ii,1)+dx*dd(ii,1)))
	    dgcd=aa(ii,2)+dx*(bb(ii,2)+dx*(cc(ii,2)+dx*dd(ii,2)))
	    st=st+gcd
	    sd=sd+dgcd
	    w_n=w_n+b2*pi
	  enddo
	enddo
	st=b2*st
	sd=b2*sd
      end