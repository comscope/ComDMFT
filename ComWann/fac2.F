      double precision function fac2(l)
      integer :: l,i,ll
      if(l.eq.-1) then
        fac2=1
        return
      endif
      if(l.eq.1) then
        fac2=1
        return
      endif
      fac2=1
      ll=(l-1)/2
      do i=1,ll
        fac2=fac2*(2*i+1)
      enddo
      end
