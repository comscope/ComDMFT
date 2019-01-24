      subroutine sym_sig_rs_rel(f_red,f)
c     Get F for all R-sites from symmetrized (reduced) representation
c     i.e. to get F(R,tEL,t'E'L') from F_RED(EL,E'L',Rtt') 
c     where Rtt' is reduced to compact number of terms
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	complex*16, intent(in) :: f_red(maxel,maxel,nrr_max)
	complex*16, intent(out) :: f(nfun,nfun,ndim3_kk(me3_k+1))
	integer :: j,ig,iatom,isort,jatom,jsort,jnd0,ind0,ir,ir0,iatom0,
     &	         jatom0,l,ie,lm,km,lm0,i,in,jn,li,ist,ii,jj,mj,ind_ir
	complex*16 :: y0j(2*maxb+2),t0j(2*maxb+2)
      complex*16, allocatable :: tmpj(:,:)
      allocate(tmpj(maxel,maxel))
      f=(0.d0,0.d0)
	do jatom=1,natom
	  jsort=is(jatom)
	  jnd0=io_lem(jatom)-1
	  do iatom=1,natom
	    isort=is(iatom)
	    ind0=io_lem(iatom)-1
	    do ind_ir=1,ndim3_kk(me3_k+1)
	      ir=n3_mpi_kk(me3_k+1)+ind_ir
	      ir0=r0_pair(ir,iatom,jatom)
	      iatom0=nrr_red(2,ir0)
	      jatom0=nrr_red(3,ir0)
	      ig=ig_pair(ir,iatom,jatom)
	      do j=1,lfunm(jsort)
              li=0
              do l=0,lmb(isort)
                ist=-1
                if(l.eq.0) ist=1
                do ii=ist,1,2
                  li=li+1
                  jj=l+l+ii
                  do ie=1,ntle(l,isort)
			        in=1
			        if(augm(ie,l,isort)/='LOC') in=2
			        do jn=1,in
	                lm0=0
                      do mj=-jj,jj,2
                        lm0=lm0+1
			            call getlimj(lm,l,ii,mj,li,1)
                        km=indbasa(jn,ie,lm,isort)
	                  t0j(lm0)=f_red(km,j,ir0)
	                enddo
	                call rotate2(t0j,y0j,l,ii,uj(1,ig),1,1)
                      do mj=-jj,jj,2
			            call getlimj(lm,l,ii,mj,li,1)
                        km=indbasa(jn,ie,lm,isort)
	                  lm0=(jj+mj)/2+1
	                  tmpj(km,j)=y0j(lm0)
	                enddo
	              enddo
	            enddo
	          enddo
	        enddo
	      enddo   !! over j
	      do i=1,lfunm(isort)
              li=0
              do l=0,lmb(jsort)
                ist=-1
                if(l.eq.0) ist=1
                do ii=ist,1,2
                  li=li+1
                  jj=l+l+ii
                  do ie=1,ntle(l,jsort)
			        in=1
			        if(augm(ie,l,jsort)/='LOC') in=2
			        do jn=1,in
	                lm0=0
                      do mj=-jj,jj,2
                        lm0=lm0+1
			            call getlimj(lm,l,ii,mj,li,1)
                        km=indbasa(jn,ie,lm,jsort)
	                  t0j(lm0)=tmpj(i,km)
	                enddo
	                call rotate2(t0j,y0j,l,ii,uj(1,ig),1,2)
                      do mj=-jj,jj,2
			            call getlimj(lm,l,ii,mj,li,1)
                        km=indbasa(jn,ie,lm,jsort)
	                  lm0=(jj+mj)/2+1
	                  f(ind0+i,jnd0+km,ind_ir)=
     &				        f(ind0+i,jnd0+km,ind_ir)+y0j(lm0)
	                enddo
	              enddo
	            enddo   !! over ie
	          enddo   !! over ii
	        enddo   !! over l
	      enddo  !! over i
	    enddo
	  enddo
	enddo
	deallocate(tmpj)
      end