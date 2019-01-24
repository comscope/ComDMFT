      subroutine rr_mesh_coarse
	use manager_mod
	use parallel_mod
      use sitegen_mod
	use solid_mod
	use units_mod
      implicit none
      integer :: i1,i2,i3,ir,ir0
      real*8 :: pi2,v(3),s
      integer, allocatable :: iref(:)
      real*8, allocatable :: a(:)
      pi2=pi+pi
c -----  Coarse R-vectors inside parallelepiped ------------------
      allocate(rr_coarse(3,nqdiv))
      ir=0
      do i3=0,ndiv(3)-1
        do i2=0,ndiv(2)-1
          do i1=0,ndiv(1)-1
            ir=ir+1
            rr_coarse(:,ir)=i1*rbas(:,1)+i2*rbas(:,2)+i3*rbas(:,3)
          enddo
        enddo
      enddo
c ----------------- List of unique neighbors ------------------------
      allocate(rr_neigh(3,5))
      n_rr_neigh=1
      rr_neigh(:,1)=0.d0
      allocate(a(nqdiv))
      allocate(iref(nqdiv))
      do ir=1,nqdiv
        call zone0_cart(rr_coarse(1,ir),rb0,v)
        a(ir)=dot_product(v,v)
      enddo
      call ordering_all(nqdiv,a,iref)
      s=0.d0
      do ir=2,nqdiv
        ir0=iref(ir)
        if(abs(s-a(ir0))>1.d-10) then
          n_rr_neigh=n_rr_neigh+1
          s=a(ir0)
          call zone0_cart(rr_coarse(1,ir0),rb0,v)
          rr_neigh(:,n_rr_neigh)=v
        endif
        if(n_rr_neigh==5) exit
      enddo
      deallocate(a,iref)    
      end
