      subroutine interp_g_tau_bnd(lagr,n0,n,n_in,n_out,tau_in,tau_out,
     &                            gx0_in,g_in,gx_in,g_out)
      implicit none
      integer, intent(in) :: lagr,n0,n,n_in,n_out
      real*8, intent(in) :: tau_in(0:n_in),tau_out(0:n_out)
      complex*16, intent(in) :: gx0_in(n0,n0,0:n_in),
     &                          gx_in(n0,n0,0:n_out)
	complex*16, intent(inout) :: g_in(n0,n0,0:n_in)
	complex*16, intent(out) :: g_out(n0,n0,0:n_out)
	integer :: i_tau,i,j
c --------- Subract G_X on old mesh ----------------------------
      do i_tau=0,n_in
        do j=1,n
          do i=1,n
            g_in(i,j,i_tau)=g_in(i,j,i_tau)-gx0_in(i,j,i_tau)
          enddo
        enddo
      enddo
c ------ Interpolate correlated part ----------------------------- 
      call interp_array_opt(g_in,g_out,tau_in,tau_out,n_in,n_out,2*n0,
     &                      n0,2*n0,n0,2*n,n,lagr)
c -------- Add exchange part on the new mesh ----------------------
      do i_tau=0,n_out
        do j=1,n
          do i=1,n
            g_out(i,j,i_tau)=g_out(i,j,i_tau)+gx_in(i,j,i_tau)
          enddo
        enddo
      enddo       
c --------- Return G_X on old mesh ----------------------------
      do i_tau=0,n_in
        do j=1,n
          do i=1,n
            g_in(i,j,i_tau)=g_in(i,j,i_tau)+gx0_in(i,j,i_tau)
          enddo
        enddo
      enddo  
      end