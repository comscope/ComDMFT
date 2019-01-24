      subroutine interp_w_nu_shell(f_in,chb,f_out,x_out,n_out,n)
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: n_out,n
	real*8, intent(in) :: f_in(n,0:n_nu),chb(n,0:n_nu_asy-2),
     &	                    x_out(0:n_out)
	real*8, intent(out) :: f_out(n,0:n_out)
	integer :: io
	do io=0,n_out
	  call boson_interp_nu_v(n,f_in,chb,x_out(io),f_out(1,io))
	enddo   !! over io
      end
      
      
      subroutine interp_w_nu_c(f_in,chb,f_out,x_out,n_out,n)
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: n_out,n
	real*8, intent(in) :: x_out(0:n_out)
	complex*16, intent(in) :: f_in(n,0:n_nu),chb(n,0:n_nu_asy-2)
	complex*16, intent(out) :: f_out(n,0:n_out)
	real*8, allocatable :: t1(:,:),t2(:,:),t(:,:)
	integer :: io
	allocate(t1(n,0:n_nu))
	allocate(t2(n,0:n_nu_asy-2))
	allocate(t(n,2))
	do io=0,n_out
	  t1=real(f_in)
	  t2=real(chb)
	  call boson_interp_nu_v(n,t1,t2,x_out(io),t(1,1))
	  t1=imag(f_in)
	  t2=imag(chb)
	  call boson_interp_nu_v(n,t1,t2,x_out(io),t(1,2))
	  f_out(:,io)=dcmplx(t(:,1),t(:,2))
	enddo   !! over io
	deallocate(t1,t2,t)
      end