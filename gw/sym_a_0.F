      subroutine sym_a_0(ev_new,k,ev_bnd0,neig,ig,pt,k0)
c     k0 here is always from NPNT set (not npnt)      
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: k,neig,ig,k0
	real*8, intent(in) :: pt(3)
	complex*16, intent(in) :: ev_bnd0(nbasmpw,neig)
	complex*16, intent(out) :: ev_new(nbasmpw,neig)
      integer :: j,ibnd,i1,i2,j0,i,ng,k00
      real*8 :: qq(3),qg(3),pi2,fas,ddot
	complex*16 :: cf
	if(k==k0) then
	  ev_new=ev_bnd0
	  return
	endif
	pi2=pi+pi
c -------- Symmetrization into another q - point ---------------------
      k00=i_kref(k0)
	ng=nbask(k00)/nrel
      i1=0
      i2=0
	do j=1,ng
	  if(irel<2) then
	    do ibnd=1,neig
	      ev_new(j,ibnd)=ev_bnd0(j,ibnd)
	    enddo
	  else if(irel==2) then
	    i1=j
	    i2=ng+j
	    do ibnd=1,neig
	      ev_new(i1,ibnd)=conjg(uj(1,ig))*ev_bnd0(i1,ibnd)
     &                     +conjg(uj(3,ig))*ev_bnd0(i2,ibnd)
	      ev_new(i2,ibnd)=conjg(uj(2,ig))*ev_bnd0(i1,ibnd)
     &                     +conjg(uj(4,ig))*ev_bnd0(i2,ibnd)
 	    enddo
	  endif
	  if(.not.symmorph) then
	    call rotate(shift(1,ig),shift(2,ig),shift(3,ig),qq(1),
     & 	  	          qq(2),qq(3),u(2,ig),1)
	    j0=indgb(j,k0)
          do i=1,3
	      qg(i)=pt(i)+gbs(i,j0)
	    enddo
	    fas=pi2*ddot(3,qg,1,qq,1)
	    cf=dcmplx(cos(fas),-sin(fas))
	    if(irel<2) then
	      do ibnd=1,neig
	        ev_new(j,ibnd)=cf*ev_new(j,ibnd)
	      enddo
	    else if(irel==2) then
	      do ibnd=1,neig
	        ev_new(i1,ibnd)=cf*ev_new(i1,ibnd)
	        ev_new(i2,ibnd)=cf*ev_new(i2,ibnd)
	      enddo
	    endif
	  endif
	enddo  !! over j
      end
