      real*8 function free_energy_x(de)
	use manager_mod
      implicit none
	real*8 :: de,x,ffd
	x=de*betta_t
	if(x.lt.-100.d0) then
	  ffd=de
	else if(x<0.d0) then
	  ffd=de-log(1.d0+exp(x))/betta_t
	else if(x<100.d0) then
	  ffd=-log(1.d0+exp(-x))/betta_t
	else
	  ffd=0.d0
	endif
	free_energy_x=ffd
      end