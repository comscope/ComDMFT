      integer function lget(ilm)
	integer, intent(in) :: ilm
      lget=dsqrt(dfloat(ilm)-0.1d0)
      end