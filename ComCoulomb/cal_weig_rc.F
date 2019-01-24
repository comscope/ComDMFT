      subroutine cal_weig_rc

      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      use wannier90_mod
      use units_mod
      implicit none
      include "mpif.h"

      integer :: k1,k2,k3,k0,k,i_star,ig,tmparr(3),ind,m
      double precision, allocatable :: tmppnt(:,:),dist_r(:)
      double precision :: tmp1, tmp2,tmp3,tempvec(3),tempvec2(3),
     $  tempvec3(3),weig_rc_tmp(ndiv(1),ndiv(2),ndiv(3)),dist_rc_tmp,
     $  maxdev,gtild(3),v(3)
      double precision, external :: dcrossdot

      do m=1,3
        qb0(m,1) = gbas(m,1)/ndiv(1)
        qb0(m,2) = gbas(m,2)/ndiv(2) 
        qb0(m,3) = gbas(m,3)/ndiv(3)
      enddo

      allocate(weig_rc(npnt))
      weig_rc=0.0d0
      allocate(tmppnt(3,npnt))
      tmppnt=0.0d0
      allocate(dist_r(npnt))
      dist_r=0.0d0            
      
      if(we_need_u .and. p_local_k) then
        weig_rc_tmp=0.0d0
        
        dist_rc_tmp=
     $    dist_rc*
     $    (dcrossdot(gbas(:,1),gbas(:,2),gbas(:,3))*3.0d0/4.0d0/pi)
     $    **(1.0d0/3.0d0)

c$$$        dist_rc_tmp=dist_rc*(2.0d0*pi/par)
        
        call cal_kpoint_weight_distance1(gbas,qb0,
     $    ndiv,dist_rc_tmp,weig_rc_tmp,maswrk,iun)

        if (maswrk) write(iun,*) 'weig_cal result'
        do k1=1, ndiv(1)
          do k2=1, ndiv(2)
            do k3=1, ndiv(3)
              if (maswrk) then
                write(iun,'(3i5,f12.6)') k1,k2,k3,
     $            weig_rc_tmp(k1,k2,k3)

              endif
            enddo
          enddo
        enddo        

        do k0=1, npnt
          do i_star=1, k_star(k0)
            k=k_list(i_star,k0)
            ig=k_group(k)
            call zone1_number(pnt(:,k),rb0,ndiv,ind)
            call indexing(nqdiv,3,ndiv,-1,ind,tmparr)          
            weig_rc(k0)=weig_rc(k0)
     $        +weig_rc_tmp(tmparr(1),tmparr(2),tmparr(3))
     $        /dble(k_star(k0))
          enddo
        enddo
        
        if (maswrk) then
          write(iun,*) 'weig_rc', sum(weig_rc*wgt(1:npnt))        
          do k=1,npnt      
            write(iun,'(1i5,4f12.6)')
     $        k,pnt(:,k),weig_rc(k)
          enddo
        endif
      endif

      end


      subroutine cal_kpoint_weight_distance1(gbas,qb0,
     $  ndiv,dist_rc,weig,maswrk,iun)

      implicit none

      integer, intent(in) :: ndiv(3),iun
      double precision, intent(in) :: gbas(3,3), qb0(3,3),
     $  dist_rc
      logical :: maswrk
      double precision, intent(out) :: weig(ndiv(1),ndiv(2),ndiv(3))


      integer :: i1,i2,i3,j1,j2,j3,kj1,kj2,kj3,kjvec(3),kvec(3),
     $  ii,g1,g2,g3,nn,jj,cnt,bin(ndiv(1),ndiv(2),ndiv(3)),p1,p2,p3

      double precision :: dist2k,dist2k_tmp,dist2g,dist2g_tmp,a1,a2,a3

      integer, external :: lcm


C$$$      call random_seed()

c$$$      write(*,*) lcm(ndiv(2),ndiv(3))
c$$$      write(*,*) lcm(ndiv(1),lcm(ndiv(2),ndiv(3)))
      nn=lcm(ndiv(1), lcm(ndiv(2),ndiv(3)))
c$$$      write(*,*) nn
      if (nn .lt. 100) then
        nn=(100/nn+1)*nn
      endif
      if (maswrk) write(iun,*) 'fine grid for weig_rc', nn
C$$$  if ((nn-(nn/2)*2) .ne. 0) then
C$$$        nn=nn+1
C$$$      endif
C$$$      pi=datan2(1.0d0,1.0d0)*4.0d0
C$$$      allocate(kshift(3,nn,nn,nn))
      

      weig=0.0d0

c$$$      write(*,*) 'nn', nn
      bin=0

      do i1=1, nn
C$$$        call random_number(a1)        
        a1=(i1-0.5d0)/dble(nn)
        j1=nint(a1*ndiv(1))
        p1=j1+1
        if (p1 .gt. ndiv(1)) then
          p1=p1-ndiv(1)
        endif

        do i2=1, nn
C$$$          call random_number(a2)        
          a2=(i2-0.5d0)/dble(nn)
          j2=nint(a2*ndiv(2))
          p2=j2+1
          if (p2 .gt. ndiv(2)) then
            p2=p2-ndiv(2)
          endif
          

          do i3=1, nn
C$$$            call random_number(a3)        
            a3=(i3-0.5d0)/dble(nn)
            j3=nint(a3*ndiv(3))

            p3=j3+1
            if (p3 .gt. ndiv(3)) then
              p3=p3-ndiv(3)
            endif

C$$$  jvec: nearest k point without metric
C$$$  k-a=k-j+j-a

c$$$            write(*,'(3(i5,f12.6))')
c$$$     $        j1, a1*ndiv(1),j2, a2*ndiv(2),j3, a3*ndiv(3)
C$$$  distance from closest gamma point
            

            bin(p1,p2,p3)=bin(p1,p2,p3)+1
            
            dist2g=1000000.0d0
            do g1=-2, 2
              do g2=-2, 2
                do g3=-2,2
                  dist2g_tmp
     $              =dsqrt(sum((matmul(gbas,
     $              (/a1-g1,a2-g2,a3-g3/)
     $              ))**2))                  
                  if (dist2g_tmp .lt. dist2g) then                  
                    dist2g=dist2g_tmp
                  endif
                enddo
              enddo
            enddo
            
            if (dist2g .gt. dist_rc) then
              dist2k=1000000.0d0
              kjvec=0
C$$$              cnt=0
              do kj1=-2, 2
                do kj2=-2, 2
                  do kj3=-2, 2
                    dist2k_tmp
     $                =dsqrt(sum((matmul(qb0,
     $                (/kj1+j1-a1*ndiv(1),
     $                kj2+j2-a2*ndiv(2),
     $                kj3+j3-a3*ndiv(3)/)
     $                ))**2))

                    if (dist2k_tmp .lt. dist2k-0.0001d0) then
C$$$                      cnt=1
                      dist2k=dist2k_tmp
                      kjvec=(/kj1,kj2,kj3/)
C$$$                    elseif (dist2k_tmp .lt. dist2k+0.0001d0) then
C$$$                      cnt=cnt+1
C$$$                      kjvec(:,cnt)=(/kj1,kj2,kj3/)
                    endif
                  enddo
                enddo
              enddo
C$$$              write(*,'(4i5)') i1,i2,i3,cnt
C$$$  kvec: nearest k point with metric            
C$$$  do jj=1, cnt
              kvec=kjvec+(/j1,j2,j3/)
              do ii=1, 3
                if (kvec(ii) .ge. ndiv(ii)) then
                  kvec(ii)=kvec(ii)-ndiv(ii)
                elseif (kvec(ii) .lt. 0) then
                  kvec(ii)=kvec(ii)+ndiv(ii)                  
                endif
              enddo
              
              weig(kvec(1)+1,kvec(2)+1,kvec(3)+1)
     $          =weig(kvec(1)+1,kvec(2)+1,kvec(3)+1)
     $          +1.0d0/dble(nn)**3*ndiv(1)*ndiv(2)*ndiv(3)
C$$$  enddo
C$$$              if ((kvec(1) .eq. 0) .and. (kvec(2) .eq. 1)
C$$$     $          .and. (kvec(3) .eq. 4)) then
C$$$                write(*,*) i1, i2, i3
C$$$              endif
            endif
          enddo
        enddo
      enddo

      do i1=1, ndiv(1)
        do i2=1, ndiv(2)
          do i3=1, ndiv(3)
            if (bin(i1,i2,i3) .ne. bin(1,1,1)) then
              if (maswrk) then
                write(iun,*) 'weig rc error'
                call flush(iun)
              endif
              call ending
            endif
          enddo
        enddo
      enddo

      end

      integer function lcm(a,b)
      implicit none
      integer, intent(in) :: a,b
      integer, external :: gcd
      lcm = a*b / gcd(a,b)
      end function lcm

      integer function gcd(a,b)
      implicit none
      integer,intent(in) :: a,b
      integer :: t,bb,aa
      bb=b
      aa=a
      do while (bb/=0)
        t = bb
        bb = mod(aa,bb)
        aa = t
      end do
      gcd = abs(aa)
      end function gcd




      subroutine dcross(a,b,c)
**********************************************************
*     cross product (ax,ay,az)=(bx,by,bz)*(cx,cy,cz)      *
**********************************************************
      implicit none
      double precision, intent(in) :: b(3),c(3)
      double precision, intent(out) :: a(3)      
      a(1)=b(2)*c(3)-b(3)*c(2)
      a(2)=b(3)*c(1)-b(1)*c(3)
      a(3)=b(1)*c(2)-b(2)*c(1)
      end



      double precision function dcrossdot(a, b, c)
      implicit none
      double precision, intent(in) :: a(3),b(3),c(3)
      double precision :: d(3)
      
      call dcross(d, b, c)
      dcrossdot=sum(d*a)
      end
      
