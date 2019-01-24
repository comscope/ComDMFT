      subroutine diag_v_coul(v,ndim,vop,qvop,nvop)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer, intent(in) :: ndim
      complex*16, intent(in) :: v(ndim,ndim)
      integer, intent(out) :: nvop
      double precision, intent(out) :: vop(n_pbtot)
      complex*16, intent(out) :: qvop(n_pbtot,n_pbtot)
      integer :: i,i0
      double precision :: am
      double precision, allocatable :: v_diag(:)
      complex*16, allocatable :: a(:,:)
      allocate(a(ndim,ndim),v_diag(ndim))
      a=0.0d0
      vop=0.0d0
      qvop=0.0d0
      call zlacpy('a',ndim,ndim,v,ndim,a,ndim)
      call timel('*start eig_val_solver in diag_v_coul')
C$$$      call eig_val_solver(ndim,ndim,a,v_diag)
      call hermitianeigen_cmplxdouble(ndim,v_diag,a)

      call timel('*done  eig_val_solver in diag_v_coul')
      am=v_diag(ndim)*eps_v_coul
      nvop=0
      do i=1,ndim
        if(v_diag(i)>=am) nvop=nvop+1
      enddo
      i0=ndim-nvop
      do i=1,nvop
        vop(i)=v_diag(i0+i)
        qvop(1:ndim,i)=a(:,i0+i)
c$$$        write(179,*) 'qvop',vop(i)!,qvop(1:dim,i)
      enddo
      deallocate(a,v_diag)
      end
