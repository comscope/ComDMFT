      subroutine norma_wf(zb,ev,gb,nbnd,nbas,ac,ispin,indg,pnt_)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: nbnd,nbas,ispin,indg(nbas/nrel)
      real*8, intent(in) :: gb(3,nplwbas),pnt_(3)
      complex*16, intent(in) :: zb(nfun,nbnd),ev(nbasmpw,nbnd)
      real*8, intent(out) :: ac(3,nbnd)
	integer :: isort,iatom,l,ie,lm,km,i,lget,je,kmj,in,jn,in1,jn1,mj,
     &	         li,ib,ibas
	real*8 :: zz
	complex*16, allocatable :: s(:,:),tmp(:,:)
      do ib=1,nbnd	
c ------- MT contribution ------------------------------------------
	  do iatom=1,natom
	    isort=is(iatom)
	    do lm=1,nrel*(lmb(isort)+1)**2
            if(irel.ne.2) then
		      l=lget(lm)
		      li=l+1
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
		  	        zz=conjg(zb(km,ib))*zb(kmj,ib)
	              ac(3,ib)=ac(3,ib)
     &	                    +zz*ffsmt(jn,jn1,ie,je,li,isort,ispin)
	            enddo
	          enddo
	        enddo
	      enddo
	    enddo    !! over lm
	  enddo  !! over iatom
	enddo   !! over ib
c ------- Interstitial contribution ------------------------------
      allocate(s(nbas,nbas))
      allocate(tmp(nbas,nbnd))
c ------ We form the overlap matrix -----------------------------
      call overlap_lapw(gb,pnt_,nbas,s,indg)
	call zgemm('n','n',nbas,nbnd,nbas,(1.d0,0.d0),s,nbas,ev,
     &           nbasmpw,(0.d0,0.d0),tmp,nbas)
      do ib=1,nbnd
        do ibas=1,nbas
          ac(3,ib)=ac(3,ib)+conjg(ev(ibas,ib))*tmp(ibas,ib)
        enddo
      enddo
      deallocate(s,tmp)
      end