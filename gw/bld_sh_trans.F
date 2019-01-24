c-----------------------------------------------------------------------
c
!> \brief Construct transformation of spherical harmonics
!>
!> The GW code needs to know how spherical harmonic functions transform
!> under the symmetry operations associated with the spacegroup of the
!> material. From libraries such a pymatgen the transformation matrices
!> of the symmetry operations can be obtained easily. However, only the
!> \f$ Y_{lm} \f$ with \f$ l=1 \f$ transform according to them. This
!> routine uses recursion relations by Ivanic et al. [1,2,3] to express
!> the transformations for \f$ l>1 \f$ in terms of the transformations
!> for \f$ l-1 \f$ and \f$ l=1 \f$.
!>
!> Note that the input rotation matrix specifies the transformation
!>
!> \f{eqnarray}{
!>   \left(\begin{array}{c}
!>      y' \\
!>      z' \\
!>      x' 
!>   \end{array}\right) = \mathbf{R}
!>   \left(\begin{array}{c}
!>      y \\
!>      z \\
!>      x 
!>   \end{array}\right)
!> \f}
!>
!> For the complex case a matrix \f$\mathbf{D}\f$ is needed that
!> transforms total momentum functions. Currently those are
!> not available.
!>
!> ### References ###
!>
!> [1] J. Ivanic, K. Ruedenberg, "Rotation matrices for real spherical
!>     harmonics. Direct determination by recursion", J. Phys. Chem.
!>     <b>100</b> (1996) 6342-6347, doi:
!>     <a href="https://dx.doi.org/10.1021/jp953350u">
!>     10.1021/jp953350u</a>.
!>
!> [2] J. Ivanic, K. Ruedenberg, "Errata: Rotation matrices for real
!>     spherical harmonics. Direct determination by recursion", J. Phys.
!>     Chem. <b>102</b> (1998) 9099-9100, doi: 
!>     <a href="https://dx.doi.org/10.1021/jp9833350">
!>     10.1021/jp9833350</a>.
!>
!> [3] C.H. Choi, J. Ivanic, M.S. Gordon, K. Ruedenberg, "Rapid and
!>     stable determination of rotation matrices between spherical 
!>     harmonics by direct recursion", J. Chem. Phys. <b>111</b> (1999)
!>     8825-8831, doi:
!>     <a href="https://dx.doi.org/10.1063/1.480229">
!>     10.1063/1.480229</a>.
!>
      subroutine bld_sh_trans(r,g,g_j,lmax)
!>
      use solid_mod        ! for maxwig
      use manager_mod      ! for irel
      implicit none
!>
!>    Maximum \f$ l \f$ quantum number
      integer, intent(in):: lmax
!>
!>    The transformation for \f$ l=1 \f$
      real*8,  intent(in):: r(-1:1,-1:1)
!>
!>    The transformations for
!>    \f$ l=0,...,l_{max}\f$ stored as a vector
!>    \f$ g=\{R(l=0),R(l=1),...,R(l=l_{max})\} \f$
      real*8,  intent(out)::g(maxwig)
!>
!>    The transformations for \f$ l=0,...,l_{max} \f$ stored as a vector
!>    \f$ g_j=\{R(l=0),R(l=1),...,R(l=l_{max})\} \f$. These
!>    transformations are complex valued for the relativistic case where
!>    the basis functions are complex.
      complex*16, intent(out)::g_j(maxwig)
c
      integer l            ! l quantum number
      integer m,m1         ! m and m' quantum numbers
c
      integer iwig         ! indexing function
      real*8  almm1        ! a^l_{m,m1} coefficient of Eq.(6.4) of [1]
      real*8  blmm1        ! b^l_{m,m1} coefficient of Eq.(6.5) of [1]
      real*8  clmm1        ! c^l_{m,m1} coefficient of Eq.(7.6) of [1]
      real*8  dlmm1        ! d^l_{m,m1} coefficient of Eq.(7.7) of [1]
c
      g_j=0.0d0
      g(1) = 1.0d0
      if (lmax.ge.1) then
        call dcopy(3*3,r,1,g(2:),1)
        ! Clear the rest of the array to avoid trouble with funky reads
        ! of positions l,m=l+1,m1. Such coefficients are multiplied
        ! by 0 anyway but if the memory locations contain garbage we
        ! might still get floating point exceptions. 
        call dcopy(maxwig-10,0.0d0,0,g(11:),1)
      endif
      do l = 1, lmax-1
c
c       First do m1 = -l-1
c
c       Eq.(7.9c) of [2].
c
        m1 = l
        do m = l+1, 2, -1
          g(iwig(l+1,-m,-m1-1))
     +    = clmm1(l+1,m,m1+1)*(r(0,-1)*g(iwig(l,-m,m1))
     +                        +r(0,1)*g(iwig(l,-m,-m1)))
     +    + 0.5d0*dlmm1(l+1,m,m1+1)*(r(-1,-1)*g(iwig(l,m-1,m1))
     +                              +r(-1,1)*g(iwig(l,m-1,-m1))
     +                              +r(1,-1)*g(iwig(l,-m+1,m1))
     +                              +r(1,1)*g(iwig(l,-m+1,-m1)))
     +    + 0.5d0*dlmm1(l+1,-m,m1+1)*(r(-1,-1)*g(iwig(l,m+1,m1))
     +                               +r(-1,1)*g(iwig(l,m+1,-m1))
     +                               -r(1,-1)*g(iwig(l,-m-1,m1))
     +                               -r(1,1)*g(iwig(l,-m-1,-m1)))
        enddo
        m = 1
        g(iwig(l+1,-m,-m1-1))
     +  = clmm1(l+1,m,m1+1)*(r(0,-1)*g(iwig(l,-m,m1))
     +                      +r(0,1)*g(iwig(l,-m,-m1)))
     +  + sqrt(0.5d0)*dlmm1(l+1,m,m1+1)*(r(-1,-1)*g(iwig(l,m-1,m1))
     +                              +r(-1,1)*g(iwig(l,m-1,-m1)))
     +  + 0.5d0*dlmm1(l+1,-m,m1+1)*(r(-1,-1)*g(iwig(l,m+1,m1))
     +                             +r(-1,1)*g(iwig(l,m+1,-m1))
     +                             -r(1,-1)*g(iwig(l,-m-1,m1))
     +                             -r(1,1)*g(iwig(l,-m-1,-m1)))
c
c       Eq.(7.9a) of [2].
c
        g(iwig(l+1,0,-m1-1))
     +  = clmm1(l+1,0,m1+1)*(r(0,-1)*g(iwig(l,0,m1))
     +                      +r(0,1)*g(iwig(l,0,-m1)))
     +  - sqrt(0.5d0)*dlmm1(l+1,0,m1+1)*(r(1,-1)*g(iwig(l,1,m1))
     +                                  +r(1,1)*g(iwig(l,1,-m1))
     +                                  +r(-1,-1)*g(iwig(l,-1,m1))
     +                                  +r(-1,1)*g(iwig(l,-1,-m1)))
c
c       Eq.(7.9b) of [2].
c
        m = 1
        g(iwig(l+1,m,-m1-1))
     +  = clmm1(l+1,m,m1+1)*(r(0,-1)*g(iwig(l,m,m1))
     +                      +r(0,1)*g(iwig(l,m,-m1)))
     +  + sqrt(0.5d0)*dlmm1(l+1,m,m1+1)*(r(1,-1)*g(iwig(l,m-1,m1))
     +                                  +r(1,1)*g(iwig(l,m-1,-m1)))
     +  - 0.5d0*dlmm1(l+1,-m,m1+1)*(r(1,-1)*g(iwig(l,m+1,m1))
     +                             +r(1,1)*g(iwig(l,m+1,-m1))
     +                             +r(-1,-1)*g(iwig(l,-m-1,m1))
     +                             +r(-1,1)*g(iwig(l,-m-1,-m1)))
        do m = 2, l+1
          g(iwig(l+1,m,-m1-1))
     +    = clmm1(l+1,m,m1+1)*(r(0,-1)*g(iwig(l,m,m1))
     +                        +r(0,1)*g(iwig(l,m,-m1)))
     +    + 0.5d0*dlmm1(l+1,m,m1+1)*(r(1,-1)*g(iwig(l,m-1,m1))
     +                              +r(1,1)*g(iwig(l,m-1,-m1))
     +                              -r(-1,-1)*g(iwig(l,-m+1,m1))
     +                              -r(-1,1)*g(iwig(l,-m+1,-m1)))
     +    - 0.5d0*dlmm1(l+1,-m,m1+1)*(r(1,-1)*g(iwig(l,m+1,m1))
     +                               +r(1,1)*g(iwig(l,m+1,-m1))
     +                               +r(-1,-1)*g(iwig(l,-m-1,m1))
     +                               +r(-1,1)*g(iwig(l,-m-1,-m1)))
        enddo
c
c       Now do m1 = -l, ..., l
c
        do m1 = -l, l
c
c         Eq.(6.7) of [2]. Note that when m=l blmm1(l+1,-m,m1)=0.
c
          do m = l+1, 2, -1
            g(iwig(l+1,-m,m1)) 
     +      = almm1(l+1,m,m1)*r(0,0)*g(iwig(l,-m,m1))
     +      + 0.5d0*blmm1(l+1,m,m1)*(r(-1,0)*g(iwig(l,m-1,m1))
     +                              +r(1,0)*g(iwig(l,-m+1,m1)))
     +      + 0.5d0*blmm1(l+1,-m,m1)*(r(-1,0)*g(iwig(l,m+1,m1))
     +                               -r(1,0)*g(iwig(l,-m-1,m1)))
          enddo
          m = 1
          g(iwig(l+1,-m,m1)) 
     +    = almm1(l+1,m,m1)*r(0,0)*g(iwig(l,-m,m1))
     +    + sqrt(0.5d0)*blmm1(l+1,m,m1)*r(-1,0)*g(iwig(l,m-1,m1))
     +    + 0.5d0*blmm1(l+1,-m,m1)*(r(-1,0)*g(iwig(l,m+1,m1))
     +                             -r(1,0)*g(iwig(l,-m-1,m1)))
c
c         Eq.(6.3) of [1]
c
          g(iwig(l+1,0,m1))
     +    = almm1(l+1,0,m1)*r(0,0)*g(iwig(l,0,m1))
     +    - blmm1(l+1,0,m1)*(r(1,0)*g(iwig(l,1,m1))
     +                      +r(-1,0)*g(iwig(l,-1,m1)))/sqrt(2.0d0)
c
c         Eq.(6.6) of [1]. Note that when m=l blmm1(l+1,-m,m1)=0.
c
          m = 1
          g(iwig(l+1,m,m1))
     +    = almm1(l+1,m,m1)*r(0,0)*g(iwig(l,m,m1))
     +    + sqrt(0.5d0)*blmm1(l+1,m,m1)*r(1,0)*g(iwig(l,m-1,m1))
     +    - 0.5d0*blmm1(l+1,-m,m1)*(r(1,0)*g(iwig(l,m+1,m1))
     +                             +r(-1,0)*g(iwig(l,-m-1,m1)))
          do m = 2, l+1
            g(iwig(l+1,m,m1))
     +      = almm1(l+1,m,m1)*r(0,0)*g(iwig(l,m,m1))
     +      + 0.5d0*blmm1(l+1,m,m1)*(r(1,0)*g(iwig(l,m-1,m1))
     +                              -r(-1,0)*g(iwig(l,-m+1,m1)))
     +      - 0.5d0*blmm1(l+1,-m,m1)*(r(1,0)*g(iwig(l,m+1,m1))
     +                               +r(-1,0)*g(iwig(l,-m-1,m1)))
          enddo
        enddo
c
c       Finally do m1 = l+1
c
c       Eq.(7.8b) of [1]
c
        m1 = l
        do m = l+1, 2, -1
          g(iwig(l+1,-m,m1+1))
     +    = clmm1(l+1,m,m1+1)*(r(0,1)*g(iwig(l,-m,m1))
     +                        -r(0,-1)*g(iwig(l,-m,-m1)))
     +    + 0.5d0*dlmm1(l+1,m,m1+1)*(r(-1,1)*g(iwig(l,m-1,m1))
     +                              -r(-1,-1)*g(iwig(l,m-1,-m1))
     +                              +r(1,1)*g(iwig(l,-m+1,m1))
     +                              -r(1,-1)*g(iwig(l,-m+1,-m1)))
     +    + 0.5d0*dlmm1(l+1,-m,m1+1)*(r(-1,1)*g(iwig(l,m+1,m1))
     +                               -r(-1,-1)*g(iwig(l,m+1,-m1))
     +                               -r(1,1)*g(iwig(l,-m-1,m1))
     +                               +r(1,-1)*g(iwig(l,-m-1,-m1)))
        enddo
        m = 1
        g(iwig(l+1,-m,m1+1))
     +  = clmm1(l+1,m,m1+1)*(r(0,1)*g(iwig(l,-m,m1))
     +                      -r(0,-1)*g(iwig(l,-m,-m1)))
     +  + sqrt(0.5d0)*dlmm1(l+1,m,m1+1)*(r(-1,1)*g(iwig(l,m-1,m1))
     +                                  -r(-1,-1)*g(iwig(l,m-1,-m1)))
     +  + 0.5d0*dlmm1(l+1,-m,m1+1)*(r(-1,1)*g(iwig(l,m+1,m1))
     +                             -r(-1,-1)*g(iwig(l,m+1,-m1))
     +                             -r(1,1)*g(iwig(l,-m-1,m1))
     +                             +r(1,-1)*g(iwig(l,-m-1,-m1)))
c
c       Eq.(7.5) of [1]
c
        g(iwig(l+1,0,m1+1))
     +  = clmm1(l+1,0,m1+1)*(r(0,1)*g(iwig(l,0,m1))
     +                      -r(0,-1)*g(iwig(l,0,-m1)))
     +  - sqrt(0.5d0)*dlmm1(l+1,0,m1+1)*(r(1,1)*g(iwig(l,1,m1))
     +                                  -r(1,-1)*g(iwig(l,1,-m1))
     +                                  +r(-1,1)*g(iwig(l,-1,m1))
     +                                  -r(-1,-1)*g(iwig(l,-1,-m1)))
c
c       Eq.(7.8a) of [1]
c
        m = 1
        g(iwig(l+1,m,m1+1))
     +  = clmm1(l+1,m,m1+1)*(r(0,1)*g(iwig(l,m,m1))
     +                      -r(0,-1)*g(iwig(l,m,-m1)))
     +  + sqrt(0.5d0)*dlmm1(l+1,m,m1+1)*(r(1,1)*g(iwig(l,m-1,m1))
     +                              -r(1,-1)*g(iwig(l,m-1,-m1)))
     +  - 0.5d0*dlmm1(l+1,-m,m1+1)*(r(1,1)*g(iwig(l,m+1,m1))
     +                             -r(1,-1)*g(iwig(l,m+1,-m1))
     +                             +r(-1,1)*g(iwig(l,-m-1,m1))
     +                             -r(-1,-1)*g(iwig(l,-m-1,-m1)))
        do m = 2, l+1
          g(iwig(l+1,m,m1+1))
     +    = clmm1(l+1,m,m1+1)*(r(0,1)*g(iwig(l,m,m1))
     +                        -r(0,-1)*g(iwig(l,m,-m1)))
     +    + 0.5d0*dlmm1(l+1,m,m1+1)*(r(1,1)*g(iwig(l,m-1,m1))
     +                              -r(1,-1)*g(iwig(l,m-1,-m1))
     +                              -r(-1,1)*g(iwig(l,-m+1,m1))
     +                              +r(-1,-1)*g(iwig(l,-m+1,-m1)))
     +    - 0.5d0*dlmm1(l+1,-m,m1+1)*(r(1,1)*g(iwig(l,m+1,m1))
     +                               -r(1,-1)*g(iwig(l,m+1,-m1))
     +                               +r(-1,1)*g(iwig(l,-m-1,m1))
     +                               -r(-1,-1)*g(iwig(l,-m-1,-m1)))
        enddo
      enddo
      if (irel.eq.2) then
        ! I guess the recursion relations for the transformation
        ! of total momentum functions should be feasible with what
        ! is known for spherical harmonics. However, I have not 
        ! found those recursion relations yet.
        write(iun,*)'Subroutine bld_sh_trans:'
        write(iun,*)'Symmetry operations for total momentum functions'
        write(iun,*)'not implemented yet.'
        call ending
      endif
      end
c
c-----------------------------------------------------------------------
c
      integer function iwig(l,m,m1)
c
c The vector g stores the transformation matrices for all spherical
c harmonics up to l=lmax. The resulting indexing becomes non-trivial
c therefore it seemed reasonable to encapsulate it in a function of
c its own. The expression is essentially stolen from subroutine
c parsop [parsbl.F]. 
c
      implicit none
      integer, intent(in) :: l  ! l quantum number
      integer, intent(in) :: m  ! m quantum number 
      integer, intent(in) :: m1 ! m' quantum number
c
      iwig = 1
      if ((m.gt.l).or.(-m.gt.l).or.(m1.gt.l).or.(-m1.gt.l)) return
      iwig = l*(2*l-1)*(2*l+1)/3+(2*l+1)*(l+m1)+l+m+1
c
      return
      end
c
c-----------------------------------------------------------------------
c
      real*8 function almm1(l,m,m1)
c
c The a^l_{mm'} coefficient of Eq. (6.4) of [1].
c
c [1] J. Ivanic, K. Ruedenberg, "Rotation matrices for real spherical
c     harmonics. Direct determination by recursion", J. Phys. Chem. 100
c     (1996) 6342-6347, doi: 10.1021/jp953350u
c
      implicit none
      integer, intent(in) :: l  ! l quantum number
      integer, intent(in) :: m  ! m quantum number 
      integer, intent(in) :: m1 ! m' quantum number
c
      real*8 e ! enumerator
      real*8 d ! denominator
c
      e = (l+m)*(l-m)
      d = (l+m1)*(l-m1)
      almm1 = sqrt(e/d)
c
      return
      end
c
c-----------------------------------------------------------------------
c
      real*8 function blmm1(l,m,m1)
c
c The b^l_{mm'} coefficient of Eq. (6.5) of [1].
c
c [1] J. Ivanic, K. Ruedenberg, "Rotation matrices for real spherical
c     harmonics. Direct determination by recursion", J. Phys. Chem. 100
c     (1996) 6342-6347, doi: 10.1021/jp953350u
c
      implicit none
      integer, intent(in) :: l  ! l quantum number
      integer, intent(in) :: m  ! m quantum number 
      integer, intent(in) :: m1 ! m' quantum number
c
      real*8 e ! enumerator
      real*8 d ! denominator
c
      e = (l+m)*(l+m-1)
      d = (l+m1)*(l-m1)
      blmm1 = sqrt(e/d)
c
      return
      end
c
c-----------------------------------------------------------------------
c
      real*8 function clmm1(l,m,m1)
c
c The c^l_{mm'} coefficient of Eq. (7.6) of [1].
c
c [1] J. Ivanic, K. Ruedenberg, "Rotation matrices for real spherical
c     harmonics. Direct determination by recursion", J. Phys. Chem. 100
c     (1996) 6342-6347, doi: 10.1021/jp953350u
c
      implicit none
      integer, intent(in) :: l  ! l quantum number
      integer, intent(in) :: m  ! m quantum number 
      integer, intent(in) :: m1 ! m' quantum number
c
      real*8 e ! enumerator
      real*8 d ! denominator
c
      e = (l+m)*(l-m)
      d = (l+m1)*(l+m1-1)
      clmm1 = sqrt(e/d)
c
      return
      end
c
c-----------------------------------------------------------------------
c
      real*8 function dlmm1(l,m,m1)
c
c The d^l_{mm'} coefficient of Eq. (7.7) of [1].
c
c [1] J. Ivanic, K. Ruedenberg, "Rotation matrices for real spherical
c     harmonics. Direct determination by recursion", J. Phys. Chem. 100
c     (1996) 6342-6347, doi: 10.1021/jp953350u
c
      implicit none
      integer, intent(in) :: l  ! l quantum number
      integer, intent(in) :: m  ! m quantum number 
      integer, intent(in) :: m1 ! m' quantum number
c
      real*8 e ! enumerator
      real*8 d ! denominator
c
      e = (l+m)*(l+m-1)
      d = (l+m1)*(l+m1-1)
      dlmm1 = sqrt(e/d)
c
      return
      end
c
c-----------------------------------------------------------------------
