      subroutine interp_array_opt(f_in,f_out,x_in,x_out,n_in,n_out,m0,
     &                            n0,m1,n1,m,n,order)
c ------- Lagrange's interpolation ----------
c     order = 1 linear two-point interpolation
c     order = 3 four-point interpolation
c     order = 5 six-point interpolation
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      integer, intent(in) :: n_in,n_out,m0,n0,m1,n1,m,n,order
	real*8, intent(in) :: f_in(m0,n0,0:n_in),x_in(0:n_in),
     &	                    x_out(0:n_out)
	real*8, intent(out) :: f_out(m1,n1,0:n_out)
	integer :: j,k,nprev,j0,ii,io,j1,ord
	real*8 :: s
	f_out=0.d0
	ord=min(order,n_in)
	ii=1   !! ascending order
	if(n_out>0) then
	  if(x_out(1)<x_out(0)) ii=2  !! descending order
	endif
	if(ii==1) nprev=0
	if(ii==2) nprev=n_in
	do io=0,n_out
	  if(x_in(0)-x_out(io)>1.d-9.or.x_out(io)-x_in(n_in)>1.d-9) cycle
c --- Renew NPREV -------------------
        k=nprev
        if(ii==1) then
          do j=nprev,n_in-1
            if(x_in(j)>x_out(io)) exit
            k=j
          enddo
        else if(ii==2) then
          do j=nprev,0,-1
            if(x_in(j)<x_out(io)) exit
            k=j-1
          enddo
        endif
        nprev=k
        j0=nprev-ord/2
        j0=max(j0,0)
        j0=min(j0,n_in-ord)
c ----- Interpolation itself ----------------------- 
        f_out(1:m,1:n,io)=0.d0
        do j=j0,j0+ord
          s=1.d0
          do j1=j0,j-1
            s=s*(x_out(io)-x_in(j1))/(x_in(j)-x_in(j1))
          enddo
          do j1=j+1,j0+ord
            s=s*(x_out(io)-x_in(j1))/(x_in(j)-x_in(j1))
          enddo
          f_out(1:m,1:n,io)=f_out(1:m,1:n,io)+f_in(1:m,1:n,j)*s
        enddo 
      enddo
      end