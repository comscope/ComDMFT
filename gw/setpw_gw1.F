      subroutine setpw_gw1
	use manager_mod
	use solid_mod
      implicit none
      integer :: iq,k,j,j0,k0,ig,i,gbs_number,j1
      real*8 :: qq(3)
c --------------- K_SYM_0 --------------------------------	
      do k0=1,nqdiv
	  do ig=1,ngroup
	    call rotate(pnt(1,k0),pnt(2,k0),pnt(3,k0),qq(1),qq(2),qq(3),
     &	              u(2,ig),2)  !! direct operation = A*k0
	    call zone1_number(qq,rb0,ndiv,k)
	    k=index_k1(k)
	    k_sym_0(k,ig)=k0
	    qq=pnt(:,k)-qq
          g_sym_0(k,ig)=gbs_number(qq)
	  enddo
	enddo
c --------------- GBS_SYM --------------------------------
      gbs_sym=0	
      do i=1,nplwdif
	  do ig=1,ngroup
          qq=matmul(gbs(:,i),u_opt(:,:,ig)) !!! inv ig-operation
          gbs_sym(i,ig)=gbs_number(qq)
	  enddo
	enddo
c -------- Filling INDPW_GW into another q - points --------------------
	do iq=npnt+1,nqdiv
	  k=i_kref(iq)
	  ig=k_group(iq)
	  do j=1,nplwgw(k)
	    j0=indpw_gw(j,k)
	    j1=g_sym_0(iq,ig)
	    call rotate(gbs(1,j0),gbs(2,j0),gbs(3,j0),qq(1),qq(2),qq(3),
     &		          u(2,ig),2)
          qq=qq-gbs(:,j1)   !! = AG-G_A
	    indpw_gw(j,iq)=gbs_number(qq)
	  enddo  !! over j
	enddo   !! over iq
c -----------------------------------------------------------------------
	iplf_gk=0
	do k=1,nqdiv
	  k0=i_kref(k)
	  do ig=1,nplwgw(k0)
	    i=indpw_gw(ig,k)
	    iplf_gk(i,k)=ig
	  enddo
	enddo
	iplf_bk=0
	do k=1,nqdiv
	  k0=i_kref(k)
	  do ig=1,nbask(k0)/nrel
	    i=indgb(ig,k)
	    iplf_bk(i,k)=ig
	  enddo
	enddo
c ------------------------------------------------------------	
      allocate(nplw_r(npnt))
      if(cutoff<cutoffgw) then
        nplw_rm=nplw_gw
        allocate(indpw_r(nplw_rm,nqdiv))
        nplw_r=nplwgw
        indpw_r=indpw_gw
      else if(cutoff>=cutoffgw) then
        nplw_rm=nbndtm
        allocate(indpw_r(nplw_rm,nqdiv))
        nplw_r=nbask/nrel
        indpw_r=indgb
      endif
      end