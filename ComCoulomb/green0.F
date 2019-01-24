      double precision function green0(de,tau)
      ! eigenstate green's function at tau (0<tau<\beta)
      use manager_mod
      implicit none
      double precision, intent(in) :: de,tau
      double precision :: gg0,gg1,gg,dem
      dem=150.d0 
      if(de.lt.0.d0) then       ! for hole
!     fermi occupation=1/(1.0d0+gg0)
        if(-de*betta_t.gt.dem) then
          gg0=0.d0
        else
          gg0=dexp(betta_t*de)
        endif
        if(-de*(betta_t-tau).gt.dem) then
          gg1=0.d0
        else
          gg1=dexp((betta_t-tau)*de)
        endif
        gg=-gg1/(1.d0+gg0)
      else                      ! for electron
        if(de*betta_t.gt.dem) then
          gg0=0.d0
        else
          gg0=dexp(-betta_t*de)
        endif
        if(de*tau.gt.dem) then
          gg1=0.d0
        else
          gg1=dexp(-de*tau)
        endif
        gg=-gg1*(1.d0-gg0/(1.d0+gg0))
      endif
      green0=gg
      end
