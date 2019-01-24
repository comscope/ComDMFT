      integer function lget(ilm)
c$$$  return l quantum number neglecting spin degree of freedom
      integer, intent(in) :: ilm
      lget=floor(dsqrt(dble(ilm)-0.1d0))
      end
