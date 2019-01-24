      subroutine g_x_tau(ispin,k,gx_tau,t,ndim,x)
c$$$     ------calculate eq. 17
      

      use manager_mod
      use solid_mod
      implicit none 
      integer, intent(in) :: ispin,k,ndim
      double precision, intent(in) :: t,x
      double precision, intent(out) :: gx_tau(ndim)
      integer :: i
      double precision :: de,green0
c     ------------------------------------------------------------
      do i=1,ndim
        de=e_bnd(i,k,ispin)-x
        gx_tau(i)=green0(de,t)
      enddo
      end
