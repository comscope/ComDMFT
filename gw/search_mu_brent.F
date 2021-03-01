      subroutine search_mu_brent(x1,y0,x0)
      use manager_mod
      use parallel_mod
      implicit none
      real*8, intent(in) :: x1,y0
      real*8, intent(out) :: x0

      integer :: key
      double precision :: aa, bb, machep,t, tos_aa,tos_bb,
     $  mu_target, tos_ef
      double precision, external :: cal_tos,find_fermi_level

      if(ubi=='dft') key=0
      if(ubi==' hf') key=1
      if(ubi==' qp') key=3
      if(ubi==' gw') key=2
      
      if(maswrk) then
        write(iun,*)' Chemical potential search'

        if(key==0) then
          write(iun,*)' Chem_pot_LDA Tos Dtos:'
        else if(key==1) then
          write(iun,*)' Chem_pot_x Tos Dtos:'
        else if(key==2) then
          write(iun,*)' Chem_pot_GW   Tos_X   Tos_C   Tos_SUM:'
        else if(key==3) then
          write(iun,*)' Chem_pot_QP Tos:'
        endif
      endif

      aa=-10.0d0
      bb=10.0d0
      machep = epsilon ( machep )
      t = 1.0d-10
      x0=x1

c$$$  check if we have different sign at aa and bb
      tos_aa=cal_tos(aa,y0)
      tos_bb=cal_tos(bb,y0)

      if ((tos_aa .gt. 0.0d0) .or. (tos_bb .lt. 0.0d0)) then
        if (maswrk) then
          write(iun,*) "wrong initial setup in search fermi level"
          call flush(iun)
        endif
        call ending
      endif
      
      mu_target=find_fermi_level(aa, bb, machep,t,
     $  cal_tos,y0)

      x0=mu_target
      
      tos_ef=cal_tos(x0,y0)

      if (maswrk) then
        write(iun,'(a5, 3f12.6)') 'Ef', x0, tos_ef, y0
      endif                  

      end

      double precision function cal_tos(mu_in, nele)

      use manager_mod      

      implicit none
      double precision,intent(in) :: mu_in, nele

      integer :: key
      double precision :: ttos,tdos

      if(ubi=='dft') key=0
      if(ubi==' hf') key=1
      if(ubi==' qp') key=3
      if(ubi==' gw') key=2
      

      call tos(mu_in,ttos,tdos,key)

      cal_tos=ttos-nele
      end

      double precision function find_fermi_level (a,b,machep,t,f,nele)


      implicit none
      include "mpif.h"      

      double precision, intent(in) :: a,b,machep,t,nele
      double precision, external :: f

      integer :: ii, maxcnt
      double precision ::  c,d,e,fa,fb,fc,m,p,q,r,s,sa,sb,tol
!     
!     Make local copies of A and B.
!     
      find_fermi_level=0.0d0
      sa = a
      sb = b
      fa = f (sa,nele)
      fb = f (sb,nele)

      c = sa
      fc = fa
      e = sb - sa
      d = e

      maxcnt=300
      do ii=1, maxcnt
        if ( abs ( fc ) < abs ( fb ) ) then
          sa = sb
          sb = c
          c = sa
          fa = fb
          fb = fc
          fc = fa
        end if
        tol = 2.0D+00 * machep * abs ( sb ) + t
        m = 0.5D+00 * ( c - sb )
        if ( abs ( m ) <= tol .or. fb == 0.0D+00 ) then
          exit
        end if
        if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then
          e = m
          d = e
        else
          s = fb / fa
          if ( sa == c ) then
            p = 2.0D+00 * m * s
            q = 1.0D+00 - s
          else
            q = fa / fc
            r = fb / fc
            p = s * ( 2.0D+00 * m * q * ( q - r ) -
     $        ( sb - sa ) * ( r - 1.0D+00 ) )
            q = ( q - 1.0D+00 ) * ( r - 1.0D+00 ) * ( s - 1.0D+00 )
          end if
          if ( 0.0D+00 < p ) then
            q = - q
          else
            p = - p
          end if

          s = e
          e = d
          
          if ( 2.0D+00 * p < 3.0D+00 * m * q - abs ( tol * q ) .and. 
     &      p < abs ( 0.5D+00 * s * q ) ) then
            d = p / q
          else
            e = m
            d = e
          end if
          
        end if
        
        sa = sb
        fa = fb
        
        if ( tol < abs ( d ) ) then
          sb = sb + d
        else if ( 0.0D+00 < m ) then
          sb = sb + tol
        else
          sb = sb - tol
        end if
        
        fb = f ( sb,nele)
        
        if ( ( 0.0D+00 < fb .and. 0.0D+00 < fc ) .or.
     $    ( fb <= 0.0D+00 .and. fc <= 0.0D+00 ) ) then
          c = sa
          fc = fa
          e = sb - sa
          d = e
        end if
        
      end do
      
      if (ii .eq. maxcnt) then
        print *, 'fermi level search fail'
        call ending                
      endif      
      
      find_fermi_level = sb

      end      
