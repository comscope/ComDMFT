      subroutine rotmat_mt_basis(r,g,g_j,lmax,jmax)

      use solid_mod
      use parallel_mod
      use manager_mod            
      
      implicit none

      integer, intent(in):: lmax,jmax
      real*8,  intent(in):: r(3,3)
      real*8,  intent(out):: g(maxwig)
      complex*16, intent(out):: g_j(maxwig_j)


c$$$  integer :: l,m,m1,iwig,ist,ii,iwig0,jj,im,im1,k,ik, inverse_sign
c$$$  real*8 :: su,ss,teta,fi,om,alpha,beta,gamma
c$$$  complex*16 :: suj,cu, ai
c$$$  real*8 :: su1(-lmax:lmax,-lmax:lmax),rotmat_ordinary(3,3)
c$$$  complex*16 :: su1_j(jmax+1,jmax+1)
c$$$  real*8, external :: wign, det3
c$$$  complex*16, external :: wign_j
      integer :: inverse_sign,llval,mdim,ii,ii_m,jj,jj_m,iival,i_index,
     $  jjval,jjdim,lival
      real*8 :: rotmat_ordinary(3,3), alpha,beta,gamma,
     $  jjval_dble,llval_dble,pi


      complex*16 :: ai
      complex*16, allocatable :: rotmat(:,:),
     $  to_real_harmonics(:,:)

      integer, external :: iwig_ind, iwig_j_ind
      real*8, external :: det3

      ai=dcmplx(0.0d0,1.0d0)
      pi=4.0d0*datan2(1.0d0,1.0d0)

      rotmat_ordinary(:,1)=r(:,1)/dsqrt(sum(r(:,1)*r(:,1)))
      rotmat_ordinary(:,2)=r(:,2)/dsqrt(sum(r(:,2)*r(:,2)))
      rotmat_ordinary(:,3)=r(:,3)/dsqrt(sum(r(:,3)*r(:,3)))


      if (det3(rotmat_ordinary(:,1),
     $  rotmat_ordinary(:,2),rotmat_ordinary(:,3)) <0.0) then
        inverse_sign=-1
        rotmat_ordinary=rotmat_ordinary*-1.0d0
      else
        inverse_sign=1
      endif
      
      beta=acos(rotmat_ordinary(3,3))          
      
      if (dabs(dsin(beta)) .gt. 1.0d-6) then
        alpha=datan2(rotmat_ordinary(2,3),rotmat_ordinary(1,3))
        gamma=datan2(rotmat_ordinary(3,2),-rotmat_ordinary(3,1))
      else
        alpha=0.0d0
        gamma=datan2(rotmat_ordinary(2,1),rotmat_ordinary(2,2))
      endif
      
      g=0.0d0
      g_j=0.0d0

      if (alpha <-1.0d-6) alpha=alpha+2*pi
      if (gamma <-1.0d-6) gamma=gamma+2*pi      
      
c$$$      if (maswrk) then
c$$$        write(iun,'(3f12.6)') rotmat_ordinary(1,:)*inverse_sign
c$$$        write(iun,'(3f12.6)') rotmat_ordinary(2,:)*inverse_sign
c$$$        write(iun,'(3f12.6)') rotmat_ordinary(3,:)*inverse_sign                
c$$$        write(iun,'(a, 3f12.6, i5)')
c$$$     $    'alpha, beta, gamma, inversion', alpha, beta, gamma,
c$$$     $    inverse_sign
c$$$        call flush(iun)
c$$$      endif


      do llval=0, lmax
        
        llval_dble=llval*1.0d0
        mdim=2*llval+1
        allocate(rotmat(mdim,mdim))
        allocate(to_real_harmonics(mdim,mdim))          
        rotmat=0.0d0
        to_real_harmonics=0.0d0
        
        call Wigner_D_mat(llval_dble,mdim,-alpha,-beta,-gamma,rotmat)
        
        if (inverse_sign <0) then
          rotmat=rotmat*(-1)**llval
        endif
        
        to_real_harmonics=0.0d0
        do ii=1, 2*llval+1
          ii_m=ii-llval-1
          if (ii_m .ne. 0) then
            do jj=-abs(ii_m)+llval+1, abs(ii_m)+llval+1, 2*abs(ii_m)
              jj_m=jj-llval-1
              if ((ii_m .lt. 0) .and. (jj_m .lt. 0)) then
                to_real_harmonics(ii,jj)
     $            =to_real_harmonics(ii,jj)
     $            +ai/dsqrt(2.0d0)
              elseif ((ii_m .lt. 0) .and. (jj_m .gt. 0)) then
                to_real_harmonics(ii,jj)
     $            =to_real_harmonics(ii,jj)
     $            -ai/dsqrt(2.0d0)*(-1)**abs(ii_m)
              elseif ((ii_m .gt. 0) .and. (jj_m .lt. 0)) then
                to_real_harmonics(ii,jj)
     $            =to_real_harmonics(ii,jj)
     $            +1.0d0/dsqrt(2.0d0)
              elseif ((ii_m .gt. 0) .and. (jj_m .gt. 0)) then
                to_real_harmonics(ii,jj)
     $            =to_real_harmonics(ii,jj)
     $            +1.0d0/dsqrt(2.0d0)*(-1)**abs(ii_m)
              endif
            enddo
          else
            to_real_harmonics(ii,ii)
     $        =to_real_harmonics(ii,ii)+1.0d0
          endif
        enddo
        
        rotmat
     $    =matmul(to_real_harmonics,
     $    matmul(rotmat,
     $    transpose(dconjg(to_real_harmonics))))
        
        
        do ii=1, 2*llval+1
          ii_m=ii-llval-1
          do jj=1, 2*llval+1
            jj_m=jj-llval-1          
            g(iwig_ind(llval,ii_m,jj_m))=dble(rotmat(ii,jj))
          enddo
        enddo
        deallocate(rotmat)
        deallocate(to_real_harmonics)
        
        
        do iival=-1,1, 2
          if ((llval .eq. 0) .and. (iival .eq. -1)) cycle
          i_index=(iival+1)/2+1
          jjval=llval*2+iival   ! 2j
          jjval_dble=llval+iival/2.0d0 ! j
          jjdim=jjval+1         ! j dimension
c$$$  lival=llval+llval+(iival+1)/2
          
          allocate(rotmat(jjdim,jjdim))
          rotmat=0.0d0
          
          call Wigner_D_mat(jjval_dble,jjdim,-alpha,-beta,-gamma,rotmat)
          
          if (inverse_sign <0) then
            rotmat=rotmat*(-1)**llval
          endif
          
          do ii_m=-jjval, jjval, 2
            ii=(ii_m+jjval)/2+1
            do jj_m=-jjval, jjval, 2
              jj=(jj_m+jjval)/2+1
              g_j(iwig_j_ind(llval,iival,ii_m,jj_m))=rotmat(ii,jj)
            enddo
          enddo
          deallocate(rotmat)
        enddo
      enddo
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
  
  
c$$$  c$$$      double precision function djnm(jval,mp,mm,beta_rad)
c$$$  c$$$c$$$      mp=m1
c$$$  c$$$c$$$      mm=m2
c$$$  c$$$      implicit none
c$$$  c$$$      double precision, intent(in) :: jval, mp, mm, beta_rad
c$$$  c$$$
c$$$  c$$$      integer :: kk,kmin,kmax
c$$$  c$$$
c$$$  c$$$      double precision, external :: ifactorial
c$$$  c$$$
c$$$  c$$$      kmin=min(0, nint(mm-mp))
c$$$  c$$$      kmax=max(0, nint(jval+mm), nint(jval-mp))
c$$$  c$$$      
c$$$  c$$$      djnm=0.0d0
c$$$  c$$$
c$$$  c$$$      do kk=kmin, kmax
c$$$  c$$$        if (
c$$$  c$$$     $    (jval+mm-kk.ge. 0) .and.
c$$$  c$$$     $    (jval-kk-mp.ge. 0) .and.
c$$$  c$$$     $    (kk-mm+mp.ge. 0) .and.
c$$$  c$$$     $    (kk .ge. 0)
c$$$  c$$$     $    ) then
c$$$  c$$$          djnm=djnm+(-1)**(kk-mm+mp)
c$$$  c$$$     $      *dsqrt(
c$$$  c$$$     $      ifactorial(jval+mm)*ifactorial(jval-mm)
c$$$  c$$$     $      *ifactorial(jval+mp)*ifactorial(jval-mp)
c$$$  c$$$     $      )
c$$$  c$$$     $      /(
c$$$  c$$$     $      ifactorial(jval+mm-kk)*ifactorial(jval-kk-mp)
c$$$  c$$$     $      *ifactorial(kk-mm+mp)*ifactorial(kk)
c$$$  c$$$     $      )
c$$$  c$$$     $      *dcos(beta_rad/2.0d0)**(2*jval-2*kk+mm-mp)
c$$$  c$$$     $      *dsin(beta_rad/2.0d0)**(2*kk-mm+mp)
c$$$  c$$$        endif
c$$$  c$$$      enddo
c$$$  c$$$      
c$$$  c$$$      end
c$$$  c$$$
c$$$  c$$$
c$$$  c$$$      DOUBLE PRECISION function iFactorial(j)
c$$$  c$$$      IMPLICIT NONE
c$$$  c$$$      INTEGER, intent(in) :: j
c$$$  c$$$      INTEGER :: i
c$$$  c$$$      DOUBLE PRECISION :: x
c$$$  c$$$      if (j<0) print *,
c$$$  c$$$     $  "iFactorial defined only for non-negative numbers!"
c$$$  c$$$      x=1
c$$$  c$$$      iFactorial = x
c$$$  c$$$      if (j.eq.1) return
c$$$  c$$$      DO i=2,j
c$$$  c$$$        x = x*i
c$$$  c$$$      END DO
c$$$  c$$$      iFactorial = x
c$$$  c$$$      return
c$$$  c$$$      end function iFactorial
c$$$  c$$$
      double precision function djnm (j, n, m, beta_rad)
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
      double precision :: cosb2, sinb2, sqrt_fac, sumt, denom,
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


c$$$  
c$$$  
c     
      
c-----------------------------------------------------------------------
c     
      integer function iwig_ind(l,m,m1)
c     
c     The vector g stores the transformation matrices for all spherical
c     harmonics up to l=lmax. The resulting indexing becomes non-trivial
c     therefore it seemed reasonable to encapsulate it in a function of
c     its own. The expression is essentially stolen from subroutine
c     parsop [parsbl.F]. 
c     
      implicit none
      integer, intent(in) :: l  ! l quantum number
      integer, intent(in) :: m  ! m quantum number 
      integer, intent(in) :: m1 ! m' quantum number
c     
      iwig_ind = 1
      if ((m.gt.l).or.(-m.gt.l).or.(m1.gt.l).or.(-m1.gt.l)) return
      iwig_ind = l*(2*l-1)*(2*l+1)/3+(2*l+1)*(l+m1)+l+m+1
c     
      return
      end


      integer function iwig_j_ind(l,ii,m,m1)
c     
c     The vector g stores the transformation matrices for all spherical
c     harmonics up to l=lmax. The resulting indexing becomes non-trivial
c     therefore it seemed reasonable to encapsulate it in a function of
c     its own. The expression is essentially stolen from subroutine
c     parsop [parsbl.F]. 
c     
      implicit none
      integer, intent(in) :: l  ! l quantum number
      integer, intent(in) :: ii ! i quantum number*2      
      integer, intent(in) :: m  ! m quantum number*2
      integer, intent(in) :: m1 ! m' quantum number*2
c     
      iwig_j_ind = 1
      if ((m.gt.2*l+ii).or.(-m.gt.2*l+ii)
     $  .or.(m1.gt.2*l+ii).or.(-m1.gt.2*l+ii)) return

      
      iwig_j_ind = (8*l**3+12*l**2+10*l+3)/3+l*(2*l+1)*ii
     $  +(l+l+ii+1)*(m1+ii)/2+(ii+m)/2
c     
      return
      end      
c     
c-----------------------------------------------------------------------
c     
c$$$  real*8 function almm1(l,m,m1)
c$$$  c     
c$$$  c     The a^l_{mm'} coefficient of Eq. (6.4) of [1].
c$$$  c     
c$$$  c     [1] J. Ivanic, K. Ruedenberg, "Rotation matrices for real spherical
c$$$  c     harmonics. Direct determination by recursion", J. Phys. Chem. 100
c$$$  c     (1996) 6342-6347, doi: 10.1021/jp953350u
c$$$  c     
c$$$  implicit none
c$$$  integer, intent(in) :: l  ! l quantum number
c$$$  integer, intent(in) :: m  ! m quantum number 
c$$$  integer, intent(in) :: m1 ! m' quantum number
c$$$  c     
c$$$  real*8 e                  ! enumerator
c$$$  real*8 d                  ! denominator
c$$$  c     
c$$$  e = (l+m)*(l-m)
c$$$  d = (l+m1)*(l-m1)
c$$$  almm1 = sqrt(e/d)
c$$$  c     
c$$$  return
c$$$  end
c$$$  c     
c$$$  c-----------------------------------------------------------------------
c$$$  c     
c$$$  real*8 function blmm1(l,m,m1)
c$$$  c     
c$$$  c     The b^l_{mm'} coefficient of Eq. (6.5) of [1].
c$$$  c     
c$$$  c     [1] J. Ivanic, K. Ruedenberg, "Rotation matrices for real spherical
c$$$  c     harmonics. Direct determination by recursion", J. Phys. Chem. 100
c$$$  c     (1996) 6342-6347, doi: 10.1021/jp953350u
c$$$  c     
c$$$  implicit none
c$$$  integer, intent(in) :: l  ! l quantum number
c$$$  integer, intent(in) :: m  ! m quantum number 
c$$$  integer, intent(in) :: m1 ! m' quantum number
c$$$  c     
c$$$  real*8 e                  ! enumerator
c$$$  real*8 d                  ! denominator
c$$$  c     
c$$$  e = (l+m)*(l+m-1)
c$$$  d = (l+m1)*(l-m1)
c$$$  blmm1 = sqrt(e/d)
c$$$  c     
c$$$  return
c$$$  end
c$$$  c     
c$$$  c-----------------------------------------------------------------------
c$$$  c     
c$$$  real*8 function clmm1(l,m,m1)
c$$$  c     
c$$$  c     The c^l_{mm'} coefficient of Eq. (7.6) of [1].
c$$$  c     
c$$$  c     [1] J. Ivanic, K. Ruedenberg, "Rotation matrices for real spherical
c$$$  c     harmonics. Direct determination by recursion", J. Phys. Chem. 100
c$$$  c     (1996) 6342-6347, doi: 10.1021/jp953350u
c$$$  c     
c$$$  implicit none
c$$$  integer, intent(in) :: l  ! l quantum number
c$$$  integer, intent(in) :: m  ! m quantum number 
c$$$  integer, intent(in) :: m1 ! m' quantum number
c$$$  c     
c$$$  real*8 e                  ! enumerator
c$$$  real*8 d                  ! denominator
c$$$  c     
c$$$  e = (l+m)*(l-m)
c$$$  d = (l+m1)*(l+m1-1)
c$$$  clmm1 = sqrt(e/d)
c$$$  c     
c$$$  return
c$$$  end
c$$$  c     
c$$$  c-----------------------------------------------------------------------
c$$$  c     
c$$$  real*8 function dlmm1(l,m,m1)
c$$$  c     
c$$$  c     The d^l_{mm'} coefficient of Eq. (7.7) of [1].
c$$$  c     
c$$$  c     [1] J. Ivanic, K. Ruedenberg, "Rotation matrices for real spherical
c$$$  c     harmonics. Direct determination by recursion", J. Phys. Chem. 100
c$$$  c     (1996) 6342-6347, doi: 10.1021/jp953350u
c$$$  c     
c$$$  implicit none
c$$$  integer, intent(in) :: l  ! l quantum number
c$$$  integer, intent(in) :: m  ! m quantum number 
c$$$  integer, intent(in) :: m1 ! m' quantum number
c$$$  c     
c$$$  real*8 e                  ! enumerator
c$$$  real*8 d                  ! denominator
c$$$  c     
c$$$  e = (l+m)*(l+m-1)
c$$$  d = (l+m1)*(l+m1-1)
c$$$  dlmm1 = sqrt(e/d)
c$$$  c     
c$$$  return
c$$$  end
c     
c-----------------------------------------------------------------------


c-----------------------------------------------------------------------
c     
!     > \brief Construct transformation of spherical harmonics
!     >
!     > The GW code needs to know how spherical harmonic functions transform
!     > under the symmetry operations associated with the spacegroup of the
!     > material. From libraries such a pymatgen the transformation matrices
!     > of the symmetry operations can be obtained easily. However, only the
!     > \f$ Y_{lm} \f$ with \f$ l=1 \f$ transform according to them. This
!     > routine uses recursion relations by Ivanic et al. [1,2,3] to express
!     > the transformations for \f$ l>1 \f$ in terms of the transformations
!     > for \f$ l-1 \f$ and \f$ l=1 \f$.
!     >
!     > Note that the input rotation matrix specifies the transformation
!     >
!     > \f{eqnarray}{
!     >   \left(\begin{array}{c}
!     >      y' \\
!     >      z' \\
!     >      x' 
!     >   \end{array}\right) = \mathbf{R}
!     >   \left(\begin{array}{c}
!     >      y \\
!     >      z \\
!     >      x 
!     >   \end{array}\right)
!     > \f}
!     >
!     > For the complex case a matrix \f$\mathbf{D}\f$ is needed that
!     > transforms total momentum functions. Currently those are
!     > not available.
!     >
!     > ### References ###
!     >
!     > [1] J. Ivanic, K. Ruedenberg, "Rotation matrices for real spherical
!     >     harmonics. Direct determination by recursion", J. Phys. Chem.
!     >     <b>100</b> (1996) 6342-6347, doi:
!     >     <a href="https://dx.doi.org/10.1021/jp953350u">
!     >     10.1021/jp953350u</a>.
!     >
!     > [2] J. Ivanic, K. Ruedenberg, "Errata: Rotation matrices for real
!     >     spherical harmonics. Direct determination by recursion", J. Phys.
!     >     Chem. <b>102</b> (1998) 9099-9100, doi: 
!     >     <a href="https://dx.doi.org/10.1021/jp9833350">
!     >     10.1021/jp9833350</a>.
!     >
!     > [3] C.H. Choi, J. Ivanic, M.S. Gordon, K. Ruedenberg, "Rapid and
!     >     stable determination of rotation matrices between spherical 
!     >     harmonics by direct recursion", J. Chem. Phys. <b>111</b> (1999)
!     >     8825-8831, doi:
!     >     <a href="https://dx.doi.org/10.1063/1.480229">
!     >     10.1063/1.480229</a>.
!     >
c$$$  subroutine bld_sh_trans(r,g,g_j,lmax)
c$$$  !     >
c$$$  use solid_mod             ! for maxwig
c$$$  use manager_mod           ! for irel
c$$$  implicit none
c$$$  !     >
c$$$  !     >    Maximum \f$ l \f$ quantum number
c$$$  integer, intent(in):: lmax
c$$$  !     >
c$$$  !     >    The transformation for \f$ l=1 \f$
c$$$  real*8,  intent(in):: r(-1:1,-1:1)
c$$$  !     >
c$$$  !     >    The transformations for
c$$$  !     >    \f$ l=0,...,l_{max}\f$ stored as a vector
c$$$  !     >    \f$ g=\{R(l=0),R(l=1),...,R(l=l_{max})\} \f$
c$$$  real*8,  intent(out)::g(maxwig)
c$$$  !     >
c$$$  !     >    The transformations for \f$ l=0,...,l_{max} \f$ stored as a vector
c$$$  !     >    \f$ g_j=\{R(l=0),R(l=1),...,R(l=l_{max})\} \f$. These
c$$$  !     >    transformations are complex valued for the relativistic case where
c$$$  !     >    the basis functions are complex.
c$$$  complex*16, intent(out)::g_j(maxwig)
c$$$  c     
c$$$  integer l                 ! l quantum number
c$$$  integer m,m1              ! m and m' quantum numbers
c$$$  c     
c$$$  integer iwig              ! indexing function
c$$$  real*8  almm1             ! a^l_{m,m1} coefficient of Eq.(6.4) of [1]
c$$$  real*8  blmm1             ! b^l_{m,m1} coefficient of Eq.(6.5) of [1]
c$$$  real*8  clmm1             ! c^l_{m,m1} coefficient of Eq.(7.6) of [1]
c$$$  real*8  dlmm1             ! d^l_{m,m1} coefficient of Eq.(7.7) of [1]
c$$$  c     
c$$$  g_j=0.0d0
c$$$  g(1) = 1.0d0
c$$$  if (lmax.ge.1) then
c$$$  call dcopy(3*3,r,1,g(2:),1)
c$$$  !     Clear the rest of the array to avoid trouble with funky reads
c$$$  !     of positions l,m=l+1,m1. Such coefficients are multiplied
c$$$  !     by 0 anyway but if the memory locations contain garbage we
c$$$  !     might still get floating point exceptions. 
c$$$  call dcopy(maxwig-10,0.0d0,0,g(11:),1)
c$$$  endif
c$$$  do l = 1, lmax-1
c$$$  c     
c$$$  c     First do m1 = -l-1
c$$$  c     
c$$$  c     Eq.(7.9c) of [2].
c$$$  c     
c$$$  m1 = l
c$$$  do m = l+1, 2, -1
c$$$  g(iwig(l+1,-m,-m1-1))
c$$$  +      = clmm1(l+1,m,m1+1)*(r(0,-1)*g(iwig(l,-m,m1))
c$$$  +      +r(0,1)*g(iwig(l,-m,-m1)))
c$$$  +      + 0.5d0*dlmm1(l+1,m,m1+1)*(r(-1,-1)*g(iwig(l,m-1,m1))
c$$$  +      +r(-1,1)*g(iwig(l,m-1,-m1))
c$$$  +      +r(1,-1)*g(iwig(l,-m+1,m1))
c$$$  +      +r(1,1)*g(iwig(l,-m+1,-m1)))
c$$$  +      + 0.5d0*dlmm1(l+1,-m,m1+1)*(r(-1,-1)*g(iwig(l,m+1,m1))
c$$$  +      +r(-1,1)*g(iwig(l,m+1,-m1))
c$$$  +      -r(1,-1)*g(iwig(l,-m-1,m1))
c$$$  +      -r(1,1)*g(iwig(l,-m-1,-m1)))
c$$$  enddo
c$$$  m = 1
c$$$  g(iwig(l+1,-m,-m1-1))
c$$$  +    = clmm1(l+1,m,m1+1)*(r(0,-1)*g(iwig(l,-m,m1))
c$$$  +    +r(0,1)*g(iwig(l,-m,-m1)))
c$$$  +    + sqrt(0.5d0)*dlmm1(l+1,m,m1+1)*(r(-1,-1)*g(iwig(l,m-1,m1))
c$$$  +    +r(-1,1)*g(iwig(l,m-1,-m1)))
c$$$  +    + 0.5d0*dlmm1(l+1,-m,m1+1)*(r(-1,-1)*g(iwig(l,m+1,m1))
c$$$  +    +r(-1,1)*g(iwig(l,m+1,-m1))
c$$$  +    -r(1,-1)*g(iwig(l,-m-1,m1))
c$$$  +    -r(1,1)*g(iwig(l,-m-1,-m1)))
c$$$  c     
c$$$  c     Eq.(7.9a) of [2].
c$$$  c     
c$$$  g(iwig(l+1,0,-m1-1))
c$$$  +    = clmm1(l+1,0,m1+1)*(r(0,-1)*g(iwig(l,0,m1))
c$$$  +    +r(0,1)*g(iwig(l,0,-m1)))
c$$$  +    - sqrt(0.5d0)*dlmm1(l+1,0,m1+1)*(r(1,-1)*g(iwig(l,1,m1))
c$$$  +    +r(1,1)*g(iwig(l,1,-m1))
c$$$  +    +r(-1,-1)*g(iwig(l,-1,m1))
c$$$  +    +r(-1,1)*g(iwig(l,-1,-m1)))
c$$$  c     
c$$$  c     Eq.(7.9b) of [2].
c$$$  c     
c$$$  m = 1
c$$$  g(iwig(l+1,m,-m1-1))
c$$$  +    = clmm1(l+1,m,m1+1)*(r(0,-1)*g(iwig(l,m,m1))
c$$$  +    +r(0,1)*g(iwig(l,m,-m1)))
c$$$  +    + sqrt(0.5d0)*dlmm1(l+1,m,m1+1)*(r(1,-1)*g(iwig(l,m-1,m1))
c$$$  +    +r(1,1)*g(iwig(l,m-1,-m1)))
c$$$  +    - 0.5d0*dlmm1(l+1,-m,m1+1)*(r(1,-1)*g(iwig(l,m+1,m1))
c$$$  +    +r(1,1)*g(iwig(l,m+1,-m1))
c$$$  +    +r(-1,-1)*g(iwig(l,-m-1,m1))
c$$$  +    +r(-1,1)*g(iwig(l,-m-1,-m1)))
c$$$  do m = 2, l+1
c$$$  g(iwig(l+1,m,-m1-1))
c$$$  +      = clmm1(l+1,m,m1+1)*(r(0,-1)*g(iwig(l,m,m1))
c$$$  +      +r(0,1)*g(iwig(l,m,-m1)))
c$$$  +      + 0.5d0*dlmm1(l+1,m,m1+1)*(r(1,-1)*g(iwig(l,m-1,m1))
c$$$  +      +r(1,1)*g(iwig(l,m-1,-m1))
c$$$  +      -r(-1,-1)*g(iwig(l,-m+1,m1))
c$$$  +      -r(-1,1)*g(iwig(l,-m+1,-m1)))
c$$$  +      - 0.5d0*dlmm1(l+1,-m,m1+1)*(r(1,-1)*g(iwig(l,m+1,m1))
c$$$  +      +r(1,1)*g(iwig(l,m+1,-m1))
c$$$  +      +r(-1,-1)*g(iwig(l,-m-1,m1))
c$$$  +      +r(-1,1)*g(iwig(l,-m-1,-m1)))
c$$$  enddo
c$$$  c     
c$$$  c     Now do m1 = -l, ..., l
c$$$  c     
c$$$  do m1 = -l, l
c$$$  c     
c$$$  c     Eq.(6.7) of [2]. Note that when m=l blmm1(l+1,-m,m1)=0.
c$$$  c     
c$$$  do m = l+1, 2, -1
c$$$  g(iwig(l+1,-m,m1)) 
c$$$  +        = almm1(l+1,m,m1)*r(0,0)*g(iwig(l,-m,m1))
c$$$  +        + 0.5d0*blmm1(l+1,m,m1)*(r(-1,0)*g(iwig(l,m-1,m1))
c$$$  +        +r(1,0)*g(iwig(l,-m+1,m1)))
c$$$  +        + 0.5d0*blmm1(l+1,-m,m1)*(r(-1,0)*g(iwig(l,m+1,m1))
c$$$  +        -r(1,0)*g(iwig(l,-m-1,m1)))
c$$$  enddo
c$$$  m = 1
c$$$  g(iwig(l+1,-m,m1)) 
c$$$  +      = almm1(l+1,m,m1)*r(0,0)*g(iwig(l,-m,m1))
c$$$  +      + sqrt(0.5d0)*blmm1(l+1,m,m1)*r(-1,0)*g(iwig(l,m-1,m1))
c$$$  +      + 0.5d0*blmm1(l+1,-m,m1)*(r(-1,0)*g(iwig(l,m+1,m1))
c$$$  +      -r(1,0)*g(iwig(l,-m-1,m1)))
c$$$  c     
c$$$  c     Eq.(6.3) of [1]
c$$$  c     
c$$$  g(iwig(l+1,0,m1))
c$$$  +      = almm1(l+1,0,m1)*r(0,0)*g(iwig(l,0,m1))
c$$$  +      - blmm1(l+1,0,m1)*(r(1,0)*g(iwig(l,1,m1))
c$$$  +      +r(-1,0)*g(iwig(l,-1,m1)))/sqrt(2.0d0)
c$$$  c     
c$$$  c     Eq.(6.6) of [1]. Note that when m=l blmm1(l+1,-m,m1)=0.
c$$$  c     
c$$$  m = 1
c$$$  g(iwig(l+1,m,m1))
c$$$  +      = almm1(l+1,m,m1)*r(0,0)*g(iwig(l,m,m1))
c$$$  +      + sqrt(0.5d0)*blmm1(l+1,m,m1)*r(1,0)*g(iwig(l,m-1,m1))
c$$$  +      - 0.5d0*blmm1(l+1,-m,m1)*(r(1,0)*g(iwig(l,m+1,m1))
c$$$  +      +r(-1,0)*g(iwig(l,-m-1,m1)))
c$$$  do m = 2, l+1
c$$$  g(iwig(l+1,m,m1))
c$$$  +        = almm1(l+1,m,m1)*r(0,0)*g(iwig(l,m,m1))
c$$$  +        + 0.5d0*blmm1(l+1,m,m1)*(r(1,0)*g(iwig(l,m-1,m1))
c$$$  +        -r(-1,0)*g(iwig(l,-m+1,m1)))
c$$$  +        - 0.5d0*blmm1(l+1,-m,m1)*(r(1,0)*g(iwig(l,m+1,m1))
c$$$  +        +r(-1,0)*g(iwig(l,-m-1,m1)))
c$$$  enddo
c$$$  enddo
c$$$  c     
c$$$  c     Finally do m1 = l+1
c$$$  c     
c$$$  c     Eq.(7.8b) of [1]
c$$$  c     
c$$$  m1 = l
c$$$  do m = l+1, 2, -1
c$$$  g(iwig(l+1,-m,m1+1))
c$$$  +      = clmm1(l+1,m,m1+1)*(r(0,1)*g(iwig(l,-m,m1))
c$$$  +      -r(0,-1)*g(iwig(l,-m,-m1)))
c$$$  +      + 0.5d0*dlmm1(l+1,m,m1+1)*(r(-1,1)*g(iwig(l,m-1,m1))
c$$$  +      -r(-1,-1)*g(iwig(l,m-1,-m1))
c$$$  +      +r(1,1)*g(iwig(l,-m+1,m1))
c$$$  +      -r(1,-1)*g(iwig(l,-m+1,-m1)))
c$$$  +      + 0.5d0*dlmm1(l+1,-m,m1+1)*(r(-1,1)*g(iwig(l,m+1,m1))
c$$$  +      -r(-1,-1)*g(iwig(l,m+1,-m1))
c$$$  +      -r(1,1)*g(iwig(l,-m-1,m1))
c$$$  +      +r(1,-1)*g(iwig(l,-m-1,-m1)))
c$$$  enddo
c$$$  m = 1
c$$$  g(iwig(l+1,-m,m1+1))
c$$$  +    = clmm1(l+1,m,m1+1)*(r(0,1)*g(iwig(l,-m,m1))
c$$$  +    -r(0,-1)*g(iwig(l,-m,-m1)))
c$$$  +    + sqrt(0.5d0)*dlmm1(l+1,m,m1+1)*(r(-1,1)*g(iwig(l,m-1,m1))
c$$$  +    -r(-1,-1)*g(iwig(l,m-1,-m1)))
c$$$  +    + 0.5d0*dlmm1(l+1,-m,m1+1)*(r(-1,1)*g(iwig(l,m+1,m1))
c$$$  +    -r(-1,-1)*g(iwig(l,m+1,-m1))
c$$$  +    -r(1,1)*g(iwig(l,-m-1,m1))
c$$$  +    +r(1,-1)*g(iwig(l,-m-1,-m1)))
c$$$  c     
c$$$  c     Eq.(7.5) of [1]
c$$$  c     
c$$$  g(iwig(l+1,0,m1+1))
c$$$  +    = clmm1(l+1,0,m1+1)*(r(0,1)*g(iwig(l,0,m1))
c$$$  +    -r(0,-1)*g(iwig(l,0,-m1)))
c$$$  +    - sqrt(0.5d0)*dlmm1(l+1,0,m1+1)*(r(1,1)*g(iwig(l,1,m1))
c$$$  +    -r(1,-1)*g(iwig(l,1,-m1))
c$$$  +    +r(-1,1)*g(iwig(l,-1,m1))
c$$$  +    -r(-1,-1)*g(iwig(l,-1,-m1)))
c$$$  c     
c$$$  c     Eq.(7.8a) of [1]
c$$$  c     
c$$$  m = 1
c$$$  g(iwig(l+1,m,m1+1))
c$$$  +    = clmm1(l+1,m,m1+1)*(r(0,1)*g(iwig(l,m,m1))
c$$$  +    -r(0,-1)*g(iwig(l,m,-m1)))
c$$$  +    + sqrt(0.5d0)*dlmm1(l+1,m,m1+1)*(r(1,1)*g(iwig(l,m-1,m1))
c$$$  +    -r(1,-1)*g(iwig(l,m-1,-m1)))
c$$$  +    - 0.5d0*dlmm1(l+1,-m,m1+1)*(r(1,1)*g(iwig(l,m+1,m1))
c$$$  +    -r(1,-1)*g(iwig(l,m+1,-m1))
c$$$  +    +r(-1,1)*g(iwig(l,-m-1,m1))
c$$$  +    -r(-1,-1)*g(iwig(l,-m-1,-m1)))
c$$$  do m = 2, l+1
c$$$  g(iwig(l+1,m,m1+1))
c$$$  +      = clmm1(l+1,m,m1+1)*(r(0,1)*g(iwig(l,m,m1))
c$$$  +      -r(0,-1)*g(iwig(l,m,-m1)))
c$$$  +      + 0.5d0*dlmm1(l+1,m,m1+1)*(r(1,1)*g(iwig(l,m-1,m1))
c$$$  +      -r(1,-1)*g(iwig(l,m-1,-m1))
c$$$  +      -r(-1,1)*g(iwig(l,-m+1,m1))
c$$$  +      +r(-1,-1)*g(iwig(l,-m+1,-m1)))
c$$$  +      - 0.5d0*dlmm1(l+1,-m,m1+1)*(r(1,1)*g(iwig(l,m+1,m1))
c$$$  +      -r(1,-1)*g(iwig(l,m+1,-m1))
c$$$  +      +r(-1,1)*g(iwig(l,-m-1,m1))
c$$$  +      -r(-1,-1)*g(iwig(l,-m-1,-m1)))
c$$$  enddo
c$$$  enddo
c$$$  if (irel.eq.2) then
c$$$  !     I guess the recursion relations for the transformation
c$$$  !     of total momentum functions should be feasible with what
c$$$  !     is known for spherical harmonics. However, I have not 
c$$$  !     found those recursion relations yet.
c$$$  c$$$  if (maswrk) then
c$$$  c$$$  write(iun,*)'Subroutine bld_sh_trans:'
c$$$  c$$$  write(iun,*)'Symmetry operations for total momentum functions'
c$$$  c$$$  write(iun,*)'not implemented yet.'
c$$$  c$$$  endif
c$$$  call ending
c$$$  endif
c$$$  end
      
