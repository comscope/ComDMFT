      subroutine set_adm(ad)
	use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
	use units_mod
	implicit none
	real*8, intent(out) :: ad
	integer :: k
	real*8 :: hh
      k=0
      hh=0.0d0
	if(ubi=='dft'.and.iter==1) then
	  ad=admix
	else if(ubi==' hf'.and.iter_hf==1) then
	  ad=1.d0
	else if(ubi==' gw'.and.iter_gw==1) then
	  ad=1.d0
	else if(ubi==' qp'.and.iter_qp==1) then
	  ad=1.d0
	else
	  if(ubi=='dft') then
	    hh=(acc_it_gw-admix)/(iter_dft-1)
	    k=0
	    k=iter-k-1
	  else if(ubi==' hf') then
	    hh=(acc_it_gw-adm_gw)/iter_hf
	    k=iter_dft
	    k=iter-k-1
	  else if(ubi==' gw') then
	    hh=(acc_it_gw-adm_gw)/iter_gw
	    k=iter_dft+iter_hf
	    k=iter-k-1
	  else if(ubi==' qp') then
	    hh=(acc_it_gw-adm_gw)/iter_qp
	    k=iter_dft+iter_hf+iter_gw
	    k=iter-k-1
	  endif
	  if(ubi=='dft') ad=admix+hh*k
	  if(ubi/='dft') ad=adm_gw+hh*k
	endif
      end
