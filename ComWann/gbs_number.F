      integer function gbs_number(v)
c     Gives the position of vector V in the list of GBS-vectors -----      
      use comwann_mod
      implicit none
      double precision, intent(in) :: v(3)
      integer :: i1,i2,i3
      double precision :: det,det3
      det=det3(gbas(1,1),gbas(1,2),gbas(1,3))
      i1=idnint(det3(v,gbas(1,2),gbas(1,3))/det)
      i2=idnint(det3(gbas(1,1),v,gbas(1,3))/det)
      i3=idnint(det3(gbas(1,1),gbas(1,2),v)/det)
      gbs_number=indplw(i1,i2,i3)
      end
