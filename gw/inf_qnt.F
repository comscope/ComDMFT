      subroutine inf_qnt(npt,i_krf,k_grp,ind,pt,ipk,k_str,k_str_m)
	use solid_mod
      implicit none
      integer, intent(in) :: npt
      real*8, intent(inout) :: pt(3,nqdiv)
      integer, intent(out) :: i_krf(nqdiv),k_grp(nqdiv),ind(nqdiv),
     &                        ipk(nqdiv,ngroup),k_str(npt),k_str_m
      integer :: i,k,ig,indk,k0,indk0,k1
      real*8 :: qq(3)
	integer, allocatable :: ifla(:)
	allocate(nk_group(nqdiv))
	allocate(kk_group(ngroup,nqdiv))
c -------- Symmetrization k-points into another q - points --------------------
	allocate(ifla(nqdiv))
	ifla=0
	i=npt
	nk_group=0
	do k=1,npt
	  i_krf(k)=k
	  k_grp(k)=1
	  call zone1_number(pt(:,k),rb0,ndiv,indk0)
	  ind(indk0)=k
	  ifla(indk0)=1
	  do ig=1,ngroup
          call rotate(pt(1,k),pt(2,k),pt(3,k),qq(1),qq(2),qq(3),u(2,ig),
     &                2)                     !! direct operation: qq=Rk
	    call zone1_number(qq,rb0,ndiv,indk)
	    if(indk>0) then
	      if(ifla(indk).eq.0) then
	        i=i+1
	        i_krf(i)=k
	        k_grp(i)=ig
	        ind(indk)=i
	        pt(:,i)=qq
	        ifla(indk)=1
	        k1=i
	      else if(ifla(indk)/=0) then
	        k1=ind(indk)
	      endif
	      nk_group(k1)=nk_group(k1)+1
	      kk_group(nk_group(k1),k1)=ig
	    endif
        enddo   !! over igroup
	enddo   !! over k
	deallocate(ifla)
	do ig=1,ngroup
	  do k=1,nqdiv
          call rotate(pt(1,k),pt(2,k),pt(3,k),qq(1),qq(2),qq(3),u(2,ig),
     &                2)
	    call zone1_number(qq,rb0,ndiv,indk)
	    if(indk>0) then
	      do i=1,nqdiv
	        if(ind(indk)==i) then
	          ipk(i,ig)=k
	          exit
	        endif
	      enddo
	    endif
	  enddo
	enddo
	k_str=0
	do k=1,nqdiv
	  k0=i_krf(k)
	  k_str(k0)=k_str(k0)+1
	enddo
	k_str_m=maxval(k_str)
	allocate(k_list(k_str_m,npt))
	k_str=0
	do k=1,nqdiv
	  k0=i_krf(k)
	  k_str(k0)=k_str(k0)+1
	  k_list(k_str(k0),k0)=k
	enddo
	do k=1,npt
	  do i=1,k_str(k)
	    k0=k_list(i,k)
	    ind_k_star_back(k0)=i
	  enddo
	enddo
      end