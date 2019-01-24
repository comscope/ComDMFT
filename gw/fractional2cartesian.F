c-----------------------------------------------------------------------
c
!> \brief Convert symmetry operations from fractional to Cartesian
!> coordinates
!>
!>  Transform the rotations and translations of the symmetry
!>  operations from fractional to cartesian coordinates.
!>
!>  By popular demand the symmetry operations in the "ini"-file may
!>  now be listed explicitly in fractional coordinates. The GW code
!>  needs them in Cartesian coordinates (the operators in fractional
!>  coordinates are easier to interpret which is why they are used
!>  in the input file). 
!>
!>  This routine transforms and replaces the rotation matrices
!>  in "rot" and the translation vectors in "trans". The
!>  transformation is specified by the lattice vectors in "rbas".
!>
      subroutine fractional2cartesian(ngroup,rot,trans,rbas)
      use manager_mod
      implicit none
!>
!>    The number of symmetry operations
      integer,          intent(in)    :: ngroup
!>
!>    The rotation matrices
      double precision, intent(inout) :: rot(3,3,ngroup)
!>
!>    The glide vectors (translations)
      double precision, intent(inout) :: trans(3,ngroup)
!>
!>    The lattice vectors
      double precision, intent(in)    :: rbas(3,3)
c
      double precision :: invrbs(3,3) ! Inverse of Transpose of rbas
c
      integer :: ipiv(3) ! Array of pivots
      integer :: info    ! Error code from LAPACK
c
      integer, parameter :: lwork = 18  ! LWORK for LAPACK
      double precision   :: work(lwork) ! WORK array for LAPACK
c
      integer :: ig ! Counter over ngroup symmetry operations
c
c     Setup rbast and invrbs
c
      invrbs = rbas
      call dgetrf(3,3,invrbs,3,ipiv,info)
      if (info.ne.0) then
        if (info.lt.0) then
          write(iun,*)"fractional2cartesian: dgetrf: illegal value of "
     &              //"argument number",-info
        else
          write(iun,*)"fractional2cartesian: dgetrf: singularity "
     &              //"detected",info
        endif
        call ending()
      endif
      call dgetri(3,invrbs,3,ipiv,work,lwork,info)
      if (info.ne.0) then
        if (info.lt.0) then
          write(iun,*)"fractional2cartesian: dgetri: illegal value of "
     &              //"argument number",-info
        else
          write(iun,*)"fractional2cartesian: dgetri: singularity "
     &              //"detected",info
        endif
        call ending()
      endif
c
c     Loop over all symmetry operations and transform them
c
      do ig = 1, ngroup
        call dgemm('n','n',3,3,3,1.0d0,rot(:,:,ig),3,invrbs,3,
     &                           0.0d0,work,3)
        call dgemm('n','n',3,3,3,1.0d0,rbas,3,work,3,
     &                           0.0d0,rot(:,:,ig),3)
        call dgemv('n',3,3,1.0d0,rbas,3,trans(:,ig),1,0.0d0,work,1)
        trans(1:3,ig) = work(1:3)
      enddo
c
      end subroutine fractional2cartesian
c
c-----------------------------------------------------------------------
