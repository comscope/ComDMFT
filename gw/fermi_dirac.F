      real*8 function fermi_dirac(de)
	use manager_mod
      implicit none
	real*8 :: de,x,ffd
	x=de*betta_t
	if(x.lt.-100.d0) then
	  ffd=1.d0
	else if(x.gt.100.d0) then
	  ffd=0.d0
	else
	  ffd=1.d0/(1.d0+exp(x))
	endif
	fermi_dirac=ffd
      end