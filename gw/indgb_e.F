      subroutine indgb_e
	use manager_mod
	use solid_mod
      implicit none
      integer :: iq,k,ig,j,j0,igb1,igb2,igb3
      real*8 :: qq(3),det3,det
c -------- Symmetrization INDGB  into another q - points ----------------------
      det=det3(gbas(1,1),gbas(1,2),gbas(1,3))
	do iq=npnt+1,nqdiv
	  k=i_kref(iq)
	  ig=k_group(iq)
	  do j=1,nbask(k)/nrel
	    j0=indgb(j,k)
	    call rotate(gbs(1,j0),gbs(2,j0),gbs(3,j0),qq(1),qq(2),qq(3),
     &		        u(2,ig),2)    !!  qq=RG
          igb1=idnint(det3(qq,gbas(1,2),gbas(1,3))/det)
          igb2=idnint(det3(gbas(1,1),qq,gbas(1,3))/det)
          igb3=idnint(det3(gbas(1,1),gbas(1,2),qq)/det)
	    indgb(j,iq)=indplw(igb1,igb2,igb3)
	  enddo  !! over j
	enddo   !! over iq
      end