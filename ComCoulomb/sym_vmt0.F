      subroutine sym_vmt0(iq,iq0,ind_q,v_mm0)
      use atom_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: iq,iq0,ind_q
      integer :: ig,j,i,jatom,jsort,iatom,isort,j0,i0
      double precision :: tt(3),pi2,fas
      complex*16 v_mm0(nlb2,nlb2),cf
      double precision, allocatable :: s0(:,:),sr(:,:),sj(:,:)
      allocate(s0(nlb2,nlb2),sr(nlb2,nlb2),sj(nlb2,nlb2))
      if(iq==iq0) then
        v_mm0=v_mm(:,:,ind_q)
        return
      endif
      ig=k_group(iq)
      sr=0.d0
      sj=0.d0
      do j=1,nlb2
        do i=1,nlb2
          s0(i,j)=dble(v_mm(i,j,ind_q))
        enddo
      enddo
      call sym_pro(ip(1,ig),lmpb,iob2,u(1,ig),sr,s0,nlb2,limlpb)
      do j=1,nlb2
        do i=1,nlb2
          s0(i,j)=ifia(ig)*dimag(v_mm(i,j,ind_q))
        enddo
      enddo
      call sym_pro(ip(1,ig),lmpb,iob2,u(1,ig),sj,s0,nlb2,limlpb)
      do j=1,nlb2
        do i=1,nlb2
          v_mm0(i,j)=dcmplx(sr(i,j),sj(i,j))
        enddo
      enddo
      deallocate(s0,sr,sj)
      pi2=pi+pi
      do jatom=1,natom
        jsort=is(jatom)
        j0=iob2(jatom)
        j=j0+(lmpb(jsort)+1)**2-1
        do iatom=1,natom
          isort=is(iatom)
          i0=iob2(iatom)
          i=i0+(lmpb(isort)+1)**2-1
          tt=tshift(:,iatom,ig)-tshift(:,jatom,ig)
          fas=pi2*dot_product(pnt(:,iq),tt)
          cf=dcmplx(cos(fas),sin(fas))
          v_mm0(i0:i,j0:j)=cf*v_mm0(i0:i,j0:j)
        enddo
      enddo
      end
