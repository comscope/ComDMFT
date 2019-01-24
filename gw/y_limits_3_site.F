      subroutine y_limits_3_site(x,z,ti,tj,tk,si,sj,sk,key,zout,info)
c     key = 1 To find a minimum
c     key = 2 To find a maximum
      implicit none
      integer, intent(in) :: key
      real*8, intent(in) :: x,z,ti(3),tj(3),tk(3),si,sj,sk
      integer, intent(out) :: info
      real*8, intent(out) :: zout
      real*8 :: zi2,zj2,zk2,zi_min,zi_max,zj_min,zj_max,zk_min,zk_max,
     &          z_min,z_max
      info=0
      zi2=si*si-(x-ti(1))**2-(z-ti(3))**2
      if(zi2<=0.d0) return
      zj2=sj*sj-(x-tj(1))**2-(z-tj(3))**2
      if(zj2<=0.d0) return
      zk2=sk*sk-(x-tk(1))**2-(z-tk(3))**2
      if(zk2<=0.d0) return
      zi2=sqrt(zi2)
      zj2=sqrt(zj2)
      zk2=sqrt(zk2)
      zi_min=ti(2)-zi2
      zi_max=ti(2)+zi2
      zj_min=tj(2)-zj2
      zj_max=tj(2)+zj2
      zk_min=tk(2)-zk2
      zk_max=tk(2)+zk2
      z_min=max(zi_min,zj_min,zk_min)
      z_max=min(zi_max,zj_max,zk_max)
      if(z_max<z_min) return
      if(key==1) zout=z_min
      if(key==2) zout=z_max
      info=1
      end