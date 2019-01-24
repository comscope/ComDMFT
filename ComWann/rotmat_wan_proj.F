      subroutine rotmat_wan_proj(proj_x_dir,proj_z_dir,jval,jdim,cubic,
     $  rotmat)
      
      implicit none
      integer, intent(in) :: jdim,cubic
      double precision, intent(in) :: jval,
     $  proj_x_dir(3), proj_z_dir(3)
      complex*16, intent(out) :: rotmat(jdim, jdim)

      integer :: lval,ii,jj,ii_m,jj_m
      double precision :: beta,gamma,
     $  proj_x_unit(3), proj_z_unit(3),proj_y_unit(3),
!     $  rot1(jdim, jdim),
!     $  rot2(jdim, jdim),rot3(jdim, jdim),
     $  alpha,
     $  rotmat_ordinary(3,3)
      double precision, external :: ddot
      complex*16 :: ai
!     $  tempmat3(jdim,jdim)
      complex*16, allocatable :: to_real_harmonics(:,:)


      if (jdim .ne. nint(2*jval+1)) then
        write(*,*) 'something wrong in jdim'
      endif


!      pi=datan2(1.0d0,1.0d0)*4.0d0
      ai=dcmplx(0.0d0,1.0d0)

      proj_z_unit=proj_z_dir
      proj_x_unit=proj_x_dir
!     beta evaluation
      proj_z_unit=proj_z_unit
     $  /dsqrt(ddot(3,proj_z_unit,1,proj_z_unit,1))
      proj_x_unit=proj_x_unit
     $  /dsqrt(ddot(3,proj_x_unit,1,proj_x_unit,1))

      if (dsqrt(sum(proj_x_dir*proj_x_dir)) .lt. 1.0d-6) then
        write(*,*) proj_x_dir, 'x is too small'
      endif
      if (dsqrt(sum(proj_z_dir*proj_z_dir)) .lt. 1.0d-6) then
        write(*,*) proj_z_dir, 'z is too small'
      endif



      if (sum(proj_x_dir*proj_z_dir)
     $  /dsqrt(sum(proj_z_dir*proj_z_dir))
     $  /dsqrt(sum(proj_x_dir*proj_x_dir))
     $  .gt. 1.0d-6) then
        write(*,*) 'x and z is not orthogonal'
      endif

      call dcross(proj_y_unit, proj_z_unit, proj_x_unit)

      rotmat_ordinary(1,:)=proj_x_unit
      rotmat_ordinary(2,:)=proj_y_unit
      rotmat_ordinary(3,:)=proj_z_unit


      beta=acos(rotmat_ordinary(3,3))

      if (dabs(dsin(beta)) .gt. 1.0d-6) then
        alpha=datan2(rotmat_ordinary(2,3),rotmat_ordinary(1,3))
        gamma=datan2(rotmat_ordinary(3,2),-rotmat_ordinary(3,1))
      else
        alpha=0.0d0
        gamma=datan2(rotmat_ordinary(2,1),rotmat_ordinary(2,2))
      endif

      rotmat=0.0d0

      call Wigner_D_mat(jval,jdim,-alpha,-beta,-gamma,rotmat)

      if (cubic.eq.1) then
        
        if (dabs(modulo(jval, 0.5d0)) .lt. 1.0d-6) then
          
          lval=nint(jval)
          
          allocate(to_real_harmonics(2*lval+1, 2*lval+1))
          to_real_harmonics=0.0d0
          do ii=1, 2*lval+1
            ii_m=ii-lval-1
            if (ii_m .ne. 0) then
              do jj=-abs(ii_m)+lval+1, abs(ii_m)+lval+1, 2*abs(ii_m)
                jj_m=jj-lval-1
                if ((ii_m .lt. 0) .and. (jj_m .lt. 0)) then
                  to_real_harmonics(ii,jj)
     $              =to_real_harmonics(ii,jj)
     $              +ai/dsqrt(2.0d0)
                elseif ((ii_m .lt. 0) .and. (jj_m .gt. 0)) then
                  to_real_harmonics(ii,jj)
     $              =to_real_harmonics(ii,jj)
     $              -ai/dsqrt(2.0d0)*(-1)**abs(ii_m)
                elseif ((ii_m .gt. 0) .and. (jj_m .lt. 0)) then
                  to_real_harmonics(ii,jj)
     $              =to_real_harmonics(ii,jj)
     $              +1.0d0/dsqrt(2.0d0)
                elseif ((ii_m .gt. 0) .and. (jj_m .gt. 0)) then
                  to_real_harmonics(ii,jj)
     $              =to_real_harmonics(ii,jj)
     $              +1.0d0/dsqrt(2.0d0)*(-1)**abs(ii_m)
                endif
              enddo
            else
              to_real_harmonics(ii,ii)
     $          =to_real_harmonics(ii,ii)+1.0d0
            endif
          enddo
          
          
          
          rotmat
     $      =matmul(to_real_harmonics,
     $      matmul(rotmat,
     $      transpose(dconjg(to_real_harmonics))))
          deallocate(to_real_harmonics)
          
        endif
      endif      
      
      end
      




      subroutine Wigner_D_mat(jval,jdim,alpha,beta,gamma,Wigner_D)
      implicit none
      integer,intent(in) :: jdim
      double precision, intent(in) :: jval,alpha, beta, gamma
      complex*16, intent(out) :: Wigner_D(jdim,jdim)
      integer :: ii,jj
      double precision :: ii_m, jj_m
      double precision,external :: djnm
      
      complex*16 :: ai

      ai=dcmplx(0.0d0,1.0d0)

      Wigner_D=0.0d0

      do ii=1, jdim
        ii_m=-jval+dble(ii-1)
        do jj=1, jdim
          jj_m=-jval+dble(jj-1)
          Wigner_D(ii,jj)
     $      =Wigner_D(ii,jj)
     $      +djnm(jval,ii_m,jj_m,beta)
     $      *cdexp(-ai*alpha*ii_m)*cdexp(-ai*gamma*jj_m)
        enddo
      enddo

      end





      function djnm (j, n, m, beta_rad)
!     d_{nm}^j(beta)
      implicit none
!     -------------------------------------------------------------------
!     input: angular momentum quantum numbers j, m, n (all real)
!     beta_rad = Euler angle beta (radian)
!     -------------------------------------------------------------------
!--------------------------------------------------------------------
!     formal arguments
!--------------------------------------------------------------------
      double precision, intent(in) :: j, m, n, beta_rad
!--------------------------------------------------------------------
!     local variables
!--------------------------------------------------------------------
      integer :: itmin1, itmin2, itmin, itmax1, itmax2, itmax, ij1, ij2,
     $  it, iphase, ia, ib, ic
      double precision :: djnm, cosb2, sinb2, sqrt_fac, sumt, denom,
     $  term
!--------------------------------------------------------------------
!     external functions
!--------------------------------------------------------------------
      double precision, external :: fac10 ! computes factorial(n)/10**n
!--------------------------------------------------------------------
!     program starts here
!--------------------------------------------------------------------
      cosb2 = cos(beta_rad/2.0d0)
      sinb2 = sin(beta_rad/2.0d0)
!--------------------------------------------------------------------
!     determine lower and upper limits for summation index it; these
!     are derived from the requirement that all factorials n! in the
!     denominator are restricted to values with n >=0.
!--------------------------------------------------------------------
      itmin1 = 0
      itmin2 = nint(m-n)
      itmin = max(itmin1,itmin2)
      itmax1 = nint(j+m)
      itmax2 = nint(j-n)
      itmax = min(itmax1,itmax2)
!     write (6,'(10X,A,2I6)') ' itmin, itmax = ', itmin, itmax
      ij1 = nint(j-m)
      ij2 = nint(j+n)
      sqrt_fac = sqrt(fac10(itmax1)*fac10(ij1)*fac10(ij2)*fac10(itmax2))
!     
      sumt = 0.0d0
      do it = itmin, itmax
        iphase = (-1)**it
        ia = itmax1 - it
        ib = itmax2 - it
        ic = it + nint(n-m)
!     write (6,'(10X,A,5I6)') ' it, iphase, ia, ib, ic  = ', it, iphase, ia, ib, ic
        denom = fac10(ia) * fac10(ib) * fac10(it) * fac10(ic)
        term = iphase * cosb2**(ia+ib) * sinb2**(it+ic) / denom
        sumt = sumt + term
      end do
      djnm = sqrt_fac * sumt
!     
      return
      end function djnm


      function fac10 (n)
      implicit none
!     -----------------------------------------------
!     function fac10(n) calculates factorial(n)/10**n
!     -----------------------------------------------
!     input: integer n >= 0 (you may want to check this
!     in the program calling this function)
!     -----------------------------------------------
!------------------------------------------------
!     formal arguments
!------------------------------------------------
      integer, intent(in) :: n
!------------------------------------------------
!     local variables
!------------------------------------------------
      integer :: i
      double precision :: fac10, q
!     -----------------------------------------------
      if (n == 0) then
        fac10 = 1.0d0
      else
        fac10 = 1.0d0
        q= 1.0d0
        do i = 1, n
          fac10 = fac10 * q / 10.0d0
          q = q + 1.0d0
        end do
      endif
!     
      return
      end function fac10

