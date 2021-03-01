C     File: flapwmbpt_input_sub.F
      
C     calculate real grid mesh for density and wavefunctions

      
      subroutine real_grid(avec, bvec, cvec, kmax,kratio,
     $  ii_max,jj_max, kk_max)
      implicit none

      double precision, intent(in) :: avec(3),bvec(3),cvec(3), kmax,
     $  kratio
      integer, intent(out) :: ii_max, jj_max, kk_max
      
CF2py intent(in) avec(3),bvec(3),cvec(3), kmax, kratio    
CF2py intent(out) ii_max, jj_max, kk_max

      integer :: ii, jj, kk, maxdiv
      double precision :: tempmat(3,3), tempvec(3), vol_tmp,krad_tmp,
     $  krad, vol
      double precision, external :: rad_inscribed_sphere       


      maxdiv=100

      vol=1.0d10
      

c$$$      print *, avec
c$$$      print *, bvec
c$$$      print *, cvec
c$$$      print *, kmax
c$$$      print *, kratio
c$$$      print *, maxdiv

      ii=-20
      jj=-20
      kk=-20
      tempmat=0.0d0
      tempmat(:,1)=avec*ii
      tempmat(:,2)=bvec*jj
      tempmat(:,3)=cvec*kk
      call cross(tempvec, bvec*jj,cvec*kk)
      vol_tmp=dabs(sum(avec*tempvec)*ii)
      krad_tmp=rad_inscribed_sphere(tempmat)*kratio      
c$$$      print *, 'here', krad_tmp, vol_tmp
 

      do ii =-maxdiv, maxdiv
        do jj =-maxdiv, maxdiv
          do kk =-maxdiv, maxdiv          
            if (ii*jj*kk  .ne. 0) then
              tempmat=0.0d0
              tempmat(:,1)=avec*ii
              tempmat(:,2)=bvec*jj
              tempmat(:,3)=cvec*kk
              call cross(tempvec, bvec*jj,cvec*kk)
              vol_tmp=dabs(sum(avec*tempvec)*ii)
              krad_tmp=rad_inscribed_sphere(tempmat)*kratio
              
              if ((krad_tmp>kmax) .and. (vol_tmp .lt. vol)) then
c$$$                print *, ii, jj, kk, krad_tmp, kmax, vol_tmp, vol
                ii_max=abs(ii)
                jj_max=abs(jj)
                kk_max=abs(kk)
                vol=vol_tmp
              endif
            endif
          enddo
        enddo
      enddo

      end
      


      double precision function rad_inscribed_sphere(a)
      implicit none
      double precision, intent(in) :: a(3,3)
      double precision :: r,v(3),v_abs
      call cross(v,a(1,2),a(1,3))
      v_abs=sqrt(dot_product(v,v))
      r=abs(dot_product(v,a(:,1)))/v_abs
      call cross(v,a(1,1),a(1,3))
      v_abs=sqrt(dot_product(v,v))
      r=min(r,abs(dot_product(v,a(:,2)))/v_abs)
      call cross(v,a(1,1),a(1,2))
      v_abs=sqrt(dot_product(v,v))
      r=min(r,abs(dot_product(v,a(:,3)))/v_abs)
      rad_inscribed_sphere=r/2.d0
      end function rad_inscribed_sphere
      
      
      
      subroutine cross(a,b,c)
      implicit none
      double precision,  intent(in) ::  b(3),c(3)
      double precision, intent(out) :: a(3)
      a(1)=b(2)*c(3)-b(3)*c(2)
      a(2)=b(3)*c(1)-b(1)*c(3)
      a(3)=b(1)*c(2)-b(2)*c(1)
      end subroutine cross
      
C End File flapwmbpt_input_sub.F
