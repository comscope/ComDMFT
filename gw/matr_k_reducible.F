      subroutine matr_k_reducible(ig,ispin,k,k0,ib1,ib2,red,aw1,zw1,aw0,
     &                            zw0,gtild)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: ispin,k,k0,ib1,ib2,ig
      real*8, intent(in) :: gtild(3)
      complex*16, intent(in) :: aw1(nbasmpw,nbndf),zw1(nfun,nbndf),
     &                          aw0(nbasmpw,nbndf),zw0(nfun,nbndf)
      complex*16, intent(out) :: red(ndeg_max,ndeg_max)
	integer :: isort,iatom,l,ie,lm,km,i,lget,je,kmj,in,jn,in1,jn1,mj,
     &	         li,ib,jbas,ibas,j,gbs_number,nb,nbas,ind0,ind1,
     &           j0,i0,jb,k00
	real*8 :: zz,v(3),phase,pi2,v1(3)
	complex*16 :: cc
	complex*16, allocatable :: zw(:,:),aw(:,:),s(:,:),tmp(:)
	k00=i_kref(k)
	nb=ib2-ib1+1
	nbas=nbask(k00)
	pi2=pi+pi
	allocate(zw(nfun,nb))
	allocate(aw(nbasmpw,nb))
c ----- Store transformed coefficients -----------------------------
c ----------------- MT part ---------------------------------
      do ib=1,nb
        call sym_val_left(ig,zw0(1,ib1+ib-1),zw(1,ib),2)
      enddo
      ind1=0.0d0
	do iatom=1,natom
	  isort=is(iatom)
	  ind0=io_lem(iatom)-1   !! zero index for atom iatom
        if(iatom<natom) ind1=io_lem(iatom+1)-1 !! last ind for iatom
	  if(iatom==natom) ind1=nfun !! last index for atom iatom
 	  phase=pi2*dot_product(pnt(:,k),tshift(:,iatom,ig))
	  cc=dcmplx(cos(phase),sin(phase))
	  zw(ind0+1:ind1,:)=cc*zw(ind0+1:ind1,:)
	enddo	
c ----------------- INT part ---------------------------------
      do ib=1,nb
        i=ib1+ib-1
        do j=1,nbas
	    j0=indgb(j,k0)  !! G in GBS-list
          call rotate(gbs(1,j0),gbs(2,j0),gbs(3,j0),v1(1),v1(2),v1(3),
     &                u(2,ig),2)
	    v=pnt(:,k)+v1-gtild
	    phase=-pi2*dot_product(v,shift(:,ig))
	    cc=dcmplx(cos(phase),sin(phase))
	    aw(j,ib)=cc*aw0(j,i)
	  enddo   !! over j
	enddo
c ------------------------------------------------------------------	
      allocate(s(nbas,nbas))
      allocate(tmp(nbas))
c ------ We form the overlap matrix -----------------------------
	do jbas=1,nbas
	  j0=indgb(jbas,k0)  !! G' in GBS-list
        call rotate(gbs(1,j0),gbs(2,j0),gbs(3,j0),v1(1),v1(2),v1(3),
     &              u(2,ig),2)
	  do ibas=1,nbas
	    i0=indgb(ibas,k)  !! G in GBS-list
	    v=v1-gbs(:,i0)-gtild
	    j=gbs_number(v)
	    if(complex_ro) then
		    s(ibas,jbas)=dcmplx(sovr(j),sovi(j))
	    else
		    s(ibas,jbas)=dcmplx(sovr(j),0.d0)
	    endif
	  enddo
	enddo
	red=(0.d0,0.d0)
	do jb=1,nb
        do ib=1,nb
      	  i=ib1+ib-1
c ------- MT contribution ------------------------------------------
	    do iatom=1,natom
	      isort=is(iatom)
	      do lm=1,nrel*(lmb(isort)+1)**2
              if(irel.ne.2) then
		        l=lget(lm)
              else if(irel.eq.2) then
		        call getlimj(lm,l,i,mj,li,0)
	        endif
	        do je=1,ntle(l,isort)
			    in1=1
			    if(augm(je,l,isort)/='LOC') in1=2
			    do jn1=1,in1
	            kmj=io_lem(iatom)-1+indbasa(jn1,je,lm,isort)
	            do ie=1,ntle(l,isort)
			        in=1
			        if(augm(ie,l,isort)/='LOC') in=2
			        do jn=1,in
	                km=io_lem(iatom)-1+indbasa(jn,ie,lm,isort)
		  	          zz=conjg(zw1(km,i))*zw(kmj,jb)
	                red(ib,jb)=red(ib,jb)
     &	                    +zz*ffsmt(jn,jn1,ie,je,l+1,isort,ispin)
	              enddo
	            enddo
	          enddo
	        enddo
	      enddo    !! over lm
	    enddo  !!  over iatom
c ------- Interstitial contribution ------------------------------
	    call zgemm('n','n',nbas,1,nbas,(1.d0,0.d0),s,nbas,aw(1,jb),
     &               nbasmpw,(0.d0,0.d0),tmp,nbas)
          do ibas=1,nbas
            red(ib,jb)=red(ib,jb)+conjg(aw1(ibas,i))*tmp(ibas)
          enddo
	  enddo   !! over ib
	enddo  !! over jb
      deallocate(s,aw,zw,tmp)
      end
