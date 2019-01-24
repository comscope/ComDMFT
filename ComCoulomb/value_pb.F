      subroutine value_pb(rp,k,k0,ubi_r,iat_r,r_transl,val,i_start,ndim,
     &  ndim0,key,s1,nk)
c     key = 0 Product basis     
c     key = 1 Dual Product basis     
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      character*2, intent(out) :: ubi_r
      integer, intent(in) :: k,k0,ndim0,key,nk
      integer, intent(out) :: i_start,ndim,iat_r
      double precision, intent(in) :: rp(3)
      double precision, intent(out) :: r_transl(3)
      complex*16, intent(in) :: s1(nk,nk)
      complex*16, intent(out) :: val(ndim0)
      integer :: iatom,isort,ic,ib,ia,irad,i0,l_i,iei,lfi,ind,m,lmi,i,j,
     &  nn
      double precision :: r00(3),dr0(3),pi2,rmt,f_pb,fas
      complex*16 :: cf
      double precision, allocatable :: ylm(:)
      complex*16, allocatable :: tmp(:)
      pi2=pi+pi
      allocate (ylm(limlpb))
c     -- Figure out if RP belongs to MT or INT --------------
      ubi_r='in'
      nn=3
      lat:do iatom=1,natom
      isort=is(iatom)
      do ic=-nn,nn
        do ib=-nn,nn
          do ia=-nn,nn
            r00=tau(:,iatom)+ic*rbas(:,3)+ib*rbas(:,2)+ia*rbas(:,1)
            dr0=rp-r00
            rmt=par*sqrt(dot_product(dr0,dr0))
            if(rmt<=smt(isort)) then
              ubi_r='mt'
              iat_r=iatom
              r_transl=r00-tau(:,iatom)
              exit lat
            endif
          enddo
        enddo
      enddo
      enddo lat
      if(ubi_r=='mt') then
        i0=0
        ir:do irad=nrad(isort),1,-1
        if(rmt>=r(irad,isort)) then
          i0=max(0,irad-2)
          i0=min(i0,nrad(isort)-3)
          exit ir
        endif
      enddo ir
      if(rmt==0.d0) then
        ylm(1)=1.d0/sqrt(4.d0*pi)
        ylm(2:)=0.d0
      else
        call sphharm(dr0(1)*par/rmt,dr0(2)*par/rmt,dr0(3)*par/rmt,
     &    lmpb(isort),ylm)
      endif
      i_start=iopb(iatom)
      ndim=n_pbmt0(isort)
      do l_i=0,lmpb(isort)
        do iei=1,ntle_pb(l_i,isort)
          lfi=ind_prod0(iei,l_i,isort)
          ind=ind_pb(lfi,isort)
          call fint(fun_pb(ind+i0),r(i0,isort),4,rmt,f_pb)
          do m=-l_i,l_i
            lmi=l_i*(l_i+1)+m+1
            i=ind_prod(iei,lmi,isort)
            val(i)=dcmplx(f_pb*ylm(lmi),0.d0)
          enddo
        enddo
      enddo
      fas=pi2*dot_product(pnt(:,k),r00-tau(:,iatom))
      cf=dcmplx(cos(fas),sin(fas))
      val=cf*val
      else if(ubi_r=='in') then
        i_start=n_pbmt+1
        ndim=nplwgw(k0)
        do i=1,ndim
          j=indpw_gw(i,k)
          r00=pnt(:,k)+gbs(:,j)
          fas=pi2*dot_product(r00,rp)
          val(i)=dcmplx(cos(fas),sin(fas))
        enddo
        if(key==1.and.ndim/=0) then
          allocate(tmp(ndim))
          call zgemv('t',ndim,ndim,(1.d0,0.d0),s1,nk,val,1,
     &      (0.d0,0.d0),tmp,1)
          val(1:ndim)=tmp
          deallocate(tmp)
        endif
      endif
      deallocate(ylm)
      if(ubi_r=='in') iat_r=0
      end
