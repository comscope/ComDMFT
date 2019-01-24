      subroutine pw_from_pb(ind_k,ndim,ndim0,pb,pw,key)
c     key = 1 - for V, W
c     key = 2 - for E^-1 and E
c     key = 3 - for P
c     ------------------------------------------------------

      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: ndim,ndim0,ind_k,key
      complex*16, intent(in) :: pb(ndim0,ndim0)
      complex*16, intent(out) :: pw
      complex*16 :: zdotc
      complex*16, allocatable :: tmp(:)
      allocate(tmp(ndim))
      if(key==1) then
        call zgemv('n',ndim,ndim,(1.d0,0.d0),pb,ndim0,
     &    pw_op(1,ind_k),1,(0.d0,0.d0),tmp,1)
        pw=zdotc(ndim,pw_op(1,ind_k),1,tmp,1)
      else if(key==2) then
        call zgemv('n',ndim,ndim,(1.d0,0.d0),pb,ndim0,
     &    pw_op_tild(1,ind_k),1,(0.d0,0.d0),tmp,1)
        pw=zdotc(ndim,pw_op(1,ind_k),1,tmp,1)
      else if(key==3) then
        call zgemv('n',ndim,ndim,(1.d0,0.d0),pb,ndim0,
     &    pw_op_tild(1,ind_k),1,(0.d0,0.d0),tmp,1)
        pw=zdotc(ndim,pw_op_tild(1,ind_k),1,tmp,1)
      endif
      deallocate(tmp)
      end
