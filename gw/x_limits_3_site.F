      subroutine x_limits_3_site(y,z,ti,tj,tk,si,sj,sk,zmin,zmax,info)
      implicit none
      real*8, intent(in) :: y,z,ti(3),tj(3),tk(3),si,sj,sk
      integer, intent(out) :: info
      real*8, intent(out) :: zmin,zmax
      real*8 :: zi2,zj2,zk2,zi_min,zi_max,zj_min,zj_max,zk_min,zk_max
      info=0
      zi2=si*si-(y-ti(2))**2-(z-ti(3))**2
      if(zi2<=0.d0) return
      zj2=sj*sj-(y-tj(2))**2-(z-tj(3))**2
      if(zj2<=0.d0) return
      zk2=sk*sk-(y-tk(2))**2-(z-tk(3))**2
      if(zk2<=0.d0) return
      zi2=sqrt(zi2)
      zj2=sqrt(zj2)
      zk2=sqrt(zk2)
      zi_min=ti(1)-zi2
      zi_max=ti(1)+zi2
      zj_min=tj(1)-zj2
      zj_max=tj(1)+zj2
      zk_min=tk(1)-zk2
      zk_max=tk(1)+zk2
      zmin=max(zi_min,zj_min,zk_min)
      zmax=min(zi_max,zj_max,zk_max)
      if(zmax<zmin) return
      info=1
      end