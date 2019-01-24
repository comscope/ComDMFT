      subroutine overl_bnd_spin(ov,n,ind_k)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: n,ind_k
      complex*16, intent(out) :: ov(n,n)
	integer :: isort,iatom,l,ie,lm,km,i,lget,je,kmj,in,jn,in1,jn1,mj,
     &	         li,ib,nb,k
	real*8, allocatable :: gb(:,:)
	complex*16, allocatable :: s(:,:),tmp(:,:)
	k=n3_mpi_k(me3_k+1)+ind_k
c ------- MT contribution ------------------------------------------
      allocate(tmp(nfun,n))
      tmp=(0.d0,0.d0)
      do ib=1,n
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
		  	        tmp(kmj,ib)=tmp(kmj,ib)+conjg(z_bnd(km,ib,ind_k,1))
     &		  	               *ffsmt(jn,jn1,ie,je,li,isort,1)
	            enddo
	          enddo
	        enddo
	      enddo
	    enddo    !! over lm
	  enddo  !! over iatom
	enddo   !! over ib
	call zgemm('t','n',n,n,nfun,(1.d0,0.d0),tmp,nfun,
     &	         z_bnd(1,1,ind_k,2),nfun,(0.d0,0.d0),ov,n)
      deallocate(tmp)
c ------- Interstitial contribution ------------------------------
      nb=nbask(k)
      allocate(s(nb,nb))
      allocate(tmp(nb,n))
c ------ We form the overlap matrix -----------------------------
	allocate(gb(3,nb))
      do i=1,nb
	  gb(:,i)=gbs(:,indgb(i,k))
	enddo
      call overlap_lapw(gb,pnt(1,k),nb,s,indgb(1,k))
	call zgemm('n','n',nb,n,nb,(1.d0,0.d0),s,nb,ev_bnd(1,1,ind_k,2),
     &           nbasmpw,(0.d0,0.d0),tmp,nb)
	call zgemm('c','n',n,n,nb,(1.d0,0.d0),ev_bnd(1,1,ind_k,1),nbasmpw,
     &	         tmp,nb,(1.d0,0.d0),ov,n)
      deallocate(s,tmp,gb)
      end
