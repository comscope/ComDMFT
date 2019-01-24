      subroutine output_boson_q_nu(bos,ch)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	real*8, intent(in) :: bos(0:n_nu,n_line)
	character*3, intent(in) :: ch
	real*8 :: qq,pi2a
	integer :: k,i_len,i_nu,m,m1,k0
	character*4 :: ch0
	pi2a=(pi+pi)/par
	if(ubi=='dft') ch0='.dft'
	if(ubi==' hf') ch0='.hf '
	if(ubi==' qp') ch0='.qp '
	if(ubi==' gw') ch0='.gw '
	if(ubi=='gwg') ch0='.gwg'
	m=min(n_line,5)
	m1=min(n_nu,4)
	if(maswrk) then
	  i_len=len_trim(allfile)
	  open(3,file=allfile(1:i_len)//'_NU_Q_'//ch//ch0)
	  do i_nu=0,n_nu
 	    write(3,'(f11.5,5(1x,e12.5))')w_nu(i_nu)*evolt/2,
     &	                               (bos(i_nu,k),k=1,m)
	  enddo
	  close(3)
	  open(3,file=allfile(1:i_len)//'_Q_NU_'//ch//ch0)
	  do k=1,n_line
          k0=kline_in_npnt(k)
          qq=sqrt(dot_product(pnt(:,k0),pnt(:,k0)))*pi2a
 	    write(3,'(f11.5,5(1x,e12.5))')qq,
     &	                               (bos(i_nu,k),i_nu=0,m1)
	  enddo
	  close(3)
	endif
      end