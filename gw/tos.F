      subroutine tos(x,ttos,tdos,key)
c     key=0 for free ( or LDA) electrons
c     key=1 for exchange only
c     key=2 for GW or GWG
c     key=3 for QP
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer, intent(in) :: key
	real*8, intent(in) :: x
	real*8, intent(out) :: ttos,tdos
	ttos=0.d0
	tdos=0.d0
        if(key==0) then
          call tos_lda(x,ttos,tdos)
        else if(key==1) then
	    call tos_x_solid(x,ttos,tdos)
        else if(key==2) then
		  call tos_x_solid(x,ttos,tdos)
          call tos_c_solid(x,ttos)
        else if(key==3) then
	    call tos_qp_solid(x,ttos)
	  endif
      end