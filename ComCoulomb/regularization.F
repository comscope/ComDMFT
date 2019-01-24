      subroutine regularization(v,ndim,ndim0,c,dev)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      include 'mpif.h'
      character*1,intent(in) :: c
      integer, intent(in) :: ndim,ndim0
      complex*16, intent(inout) :: v(ndim0,ndim0)
      double precision, intent(out) :: dev
      integer :: i,j,k
      double precision, allocatable :: v_diag(:)
      complex*16, allocatable :: a(:,:)
      allocate(a(ndim,ndim),v_diag(ndim))

      call zlacpy('a',ndim,ndim,v,ndim0,a,ndim)
      call timel('*st eig_val_solver regularization***')
C$$$      call eig_val_solver(ndim,ndim,a,v_diag)
      call hermitianeigen_cmplxdouble(ndim,v_diag,a)
      call timel('*dn eig_val_solver regularization***')
c      write(179,*) 'regularization of polarizability'
c      write(179,'(a, f20.12, a, f20.12)')
c     $  'maxval(v_diag)', maxval(v_diag),
c     $  'minval(v_diag)', minval(v_diag)
!      write(179,'(1000f12.6)') (v_diag(i), i=1, ndim)
!      write(100,'(1000f12.6)') (v_diag(i), i=1, ndim)
!      call mpi_barrier(mpi_comm_world,mpierr)

      if(c=='+') then
        dev=v_diag(1)
c        write(179,*) 'less than zero'        
        do i=1,ndim
          if(v_diag(i)<0.0d0) then
c            write(179,*) 'dim', i, v_diag(i)                          
            v_diag(i)=0.d0
          endif
        enddo
      else if(c=='-') then
        dev=v_diag(ndim)
        if(dev>0.d0) then
c          write(179,*) 'larger than zero'
          do i=1,ndim
            if(v_diag(i)>0.0d0) then
c              write(179,*) 'dim', i, v_diag(i)              
              v_diag(i)=0.d0
            endif
          enddo
        endif
      endif
      v=(0.d0,0.d0)
      do j=1,ndim
        do i=1,ndim
          do k=1,ndim
            v(i,j)=v(i,j)+a(i,k)*conjg(a(j,k))*v_diag(k)
          enddo
        enddo
      enddo
      deallocate(a,v_diag)
      end
