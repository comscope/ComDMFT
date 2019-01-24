      subroutine mix_up_down(n,n0,up,upo,dn,dno,t)
	use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
	use units_mod
	implicit none
	integer, intent(in) :: n,n0
	real*8, intent(in) :: upo(n0,n0),dno(n0,n0),t(n,n)
	real*8, intent(inout) :: up(n0,n0),dn(n0,n0)
	real*8, allocatable :: d(:,:,:),s(:,:,:),v2p(:,:),v2po(:,:)
      complex*16, allocatable :: v2(:,:),v2o(:,:),tmp(:,:)
	allocate(d(n,n,2))
	allocate(s(n,n,2))
	allocate(v2(n,n))
	allocate(v2o(n,n))
	allocate(v2p(n,n))
	allocate(v2po(n,n))
	allocate(tmp(n,n))
	call unpack_hermit(v2,dn,n,n0,n,(0.d0,0.d0),(1.d0,0.d0))
	call zgemm('n','c',n,n,n,(1.d0,0.d0),v2,n,t,n,(0.d0,0.d0),tmp,n)
	call zgemm('n','n',n,n,n,(1.d0,0.d0),t,n,tmp,n,(0.d0,0.d0),v2,n)
      call pack_hermit(v2,v2p,n,n,n,0.d0,1.d0)
	call unpack_hermit(v2o,dno,n,n0,n,(0.d0,0.d0),(1.d0,0.d0))
	call zgemm('n','c',n,n,n,(1.d0,0.d0),v2o,n,t,n,(0.d0,0.d0),tmp,n)
	call zgemm('n','n',n,n,n,(1.d0,0.d0),t,n,tmp,n,(0.d0,0.d0),v2o,n)
      call pack_hermit(v2o,v2po,n,n,n,0.d0,1.d0)
	d(:,:,1)=0.5d0*(up(1:n,1:n)-v2p)
	d(:,:,2)=0.5d0*(upo(1:n,1:n)-v2po)
	s(:,:,1)=0.5d0*(up(1:n,1:n)+v2p)
	s(:,:,2)=0.5d0*(upo(1:n,1:n)+v2po)
      d(:,:,1)=adspin*d(:,:,1)+(1.d0-adspin)*d(:,:,2)
      s(:,:,1)=adm_sig*s(:,:,1)+(1.d0-adm_sig)*s(:,:,2)
	up(1:n,1:n)=s(:,:,1)+d(:,:,1)
	v2p=s(:,:,1)-d(:,:,1)
	call unpack_hermit(v2,v2p,n,n,n,(0.d0,0.d0),(1.d0,0.d0))
	call zgemm('n','n',n,n,n,(1.d0,0.d0),v2,n,t,n,(0.d0,0.d0),tmp,n)
	call zgemm('c','n',n,n,n,(1.d0,0.d0),t,n,tmp,n,(0.d0,0.d0),v2,n)
      call pack_hermit(v2,dn,n,n,n0,0.d0,1.d0)
      deallocate(d,s,v2,v2o,v2p,v2po,tmp)
      end
      
      
      
      
      subroutine mixer_e_qp(eps,n,eold,enew)
      use manager_mod
      use parallel_mod
	implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: eps,eold(n)
	real*8, intent(inout) :: enew(n)
	integer :: ndeg_old,ndeg_new,i,i1,i2,j1,j2,j
	real*8 :: ad
	integer, allocatable :: nn_old(:),nn_new(:),n01_old(:,:),
     & 	                      n01_new(:,:),flag(:)
      ad=0.5d0
      allocate(nn_old(n))
      allocate(nn_new(n))
      allocate(n01_old(2,n))
      allocate(n01_new(2,n))
      call equal(1,eold,100.d0*eps,n)
      call equal(1,enew,100.d0*eps,n)
      call degen_list(eps,n,eold,ndeg_old,nn_old,n01_old)
      call degen_list(eps,n,enew,ndeg_new,nn_new,n01_new)
      allocate(flag(ndeg_old))
      if(ndeg_new/=ndeg_old) then
        if(maswrk) write(iun,*)' Error_1 in mixer_e_qp'
        goto 1
      endif
      if(maxval(abs(nn_old(1:ndeg_old)-nn_new(1:ndeg_new)))/=0) then
        if(maswrk) write(iun,*)' Error_2 in mixer_e_qp'
        goto 1
      endif
      flag=0
      do i=1,ndeg_new
        i1=n01_new(1,i)
        i2=n01_new(2,i)
        do j=1,ndeg_old
          if(flag(j)==1) cycle
          j1=n01_old(1,j)
          j2=n01_old(2,j)
          if(nn_old(j)==nn_new(i)) then
            enew(i1:i2)=ad*enew(i1:i2)+(1.d0-ad)*eold(j1:j2)
            flag(j)=1
            goto 2
          else
            cycle
          endif
        enddo
2       continue        
      enddo
1     if(minval(flag)==0) then
        if(maswrk) write(iun,*)' Error in mixer_e_qp: FLAG=0'
      endif
      deallocate(nn_old,nn_new,n01_old,n01_new,flag)
      end