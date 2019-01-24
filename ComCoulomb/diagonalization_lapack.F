      subroutine hermitianeigen_cmplxdouble(dimen, w, datamat)
      implicit none
      
      integer, intent(in) :: dimen
      complex*16, intent(inout) :: datamat(dimen, dimen)
      double precision, intent(out) :: w(dimen)
      integer :: lwork, info, j, i, ifail, lwkopt
      complex*16, allocatable :: work(:)
      complex*16 :: worktemp(2*dimen)
      double precision :: rwork(dimen*3), rcondz(dimen), zerrbd(dimen),
     $  eerrbd, eps
      
      call zheev('V','U',dimen,datamat,dimen, w,worktemp,-1,rwork, info)
      lwkopt = worktemp(1)
      allocate(work(lwkopt))
      call zheev('V','U',dimen,datamat,dimen, w,work,lwkopt,rwork, info)
      deallocate(work)
      if (info .ne. 0) then
        write(*,*) 'error in diagonalization'
      end if
      end subroutine hermitianeigen_cmplxdouble


