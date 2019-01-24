      subroutine f_om_nu_maz(f_maz,f_omr,f_omi,n_maz,n_0,w_0,npn,nsp)
c	Calculates functions for n_maz Matsubara's freuencies from splines
c	Same for Nu-values
	use manager_mod
	use parallel_mod
	use units_mod
      implicit none
	integer :: npn,nsp,i,ii,kk,ispin,n_maz,n_0,im
	real*8 :: f_omr(0:n_0,npn,nsp),f_omi(0:n_0,npn,nsp),w_n,
     &	      dx,wcnu_r,wcnu_i,w_0(0:n_0)
	complex*16 :: f_maz(0:n_maz,npn,nsp)
	real*8, allocatable :: h(:),aa(:,:),bb(:,:),cc(:,:),dd(:,:)
	allocate(aa(n_0,2),bb(n_0,2),cc(n_0,2),dd(n_0,2),h(n_0))
	do i=1,n_0
	  h(i)=w_0(i)-w_0(i-1)
	enddo
	do ispin=1,nsp
	  do i=1,npn
	    call spline_start(h,f_omr(:,i,ispin),aa(:,1),bb(:,1),
     & 	                  cc(:,1),dd(:,1),n_0)
	    call spline_start(h,f_omi(:,i,ispin),aa(:,2),bb(:,2),
     &	                  cc(:,2),dd(:,2),n_0)
	    ii=1
	    do im=0,n_maz
	      w_n=w_0(0)+2.d0*pi/betta_t*im
	      do kk=ii,n_0
	        if(w_n>=w_0(kk-1).and.w_n<=w_0(kk)) then
	          ii=kk
	          exit
	        endif
	      enddo
	      dx=w_n-w_0(ii-1)
	      wcnu_r=aa(ii,1)+dx*(bb(ii,1)+dx*(cc(ii,1)+dx*dd(ii,1)))
	      wcnu_i=aa(ii,2)+dx*(bb(ii,2)+dx*(cc(ii,2)+dx*dd(ii,2)))
	      f_maz(im,i,ispin)=dcmplx(wcnu_r,wcnu_i)
	    enddo
	  enddo
	enddo
	deallocate(aa,bb,cc,dd,h)
      end