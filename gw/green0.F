      double precision function green0(de,tau)
	use manager_mod
      implicit none
	real*8, intent(in) :: de,tau
	real*8 :: gg0,gg1,gg,dem
	dem=150.d0
	if(de.lt.0.d0) then
	  if(-de*betta_t.gt.dem) then
	    gg0=0.d0
	  else
	    gg0=exp(betta_t*de)
	  endif
	  if(-de*(betta_t-tau).gt.dem) then
	    gg1=0.d0
	  else
	    gg1=exp((betta_t-tau)*de)
	  endif
	  gg=-gg1/(1.d0+gg0)
	else
	  if(de*betta_t.gt.dem) then
	    gg0=0.d0
	  else
	    gg0=exp(-betta_t*de)
	  endif
	  if(de*tau.gt.dem) then
	    gg1=0.d0
	  else
	    gg1=exp(-de*tau)
	  endif
	  gg=-gg1*(1.d0-gg0/(1.d0+gg0))
	endif
	green0=gg
      end