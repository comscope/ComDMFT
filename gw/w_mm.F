      subroutine w_mm(ind_tau)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: ind_tau
	integer :: k,k0,ind,iatom,isort,ir0,ind_k,ir,jatom,jsort,ig,
     &	         ind_ir0,j,i,jatom0,iatom0,l1,nl1,l,nl,je,ie,lm1,lm,km1,
     &           km,jat,jat0,iat,iat0
      real*8 :: pi2,v(3),phase
      real*8, allocatable :: ai(:,:,:),aj(:,:,:),a(:,:),b(:,:),c(:,:)
	complex*16, allocatable :: w0(:,:,:,:,:),w1(:,:,:,:,:),
     &	                         cphase(:,:,:,:)
	pi2=pi+pi
	w_rs_mm=0.d0
	do jsort=1,nsort
	  do isort=1,nsort
c ---- Store phase factors ---------------------------------------
          allocate(cphase(nqdiv,nhsort(isort),nhsort(jsort),ngroup))
          do jatom0=1,nhsort(jsort)
            jatom=iat_sort(jatom0,jsort)
            do iatom0=1,nhsort(isort)
              iatom=iat_sort(iatom0,isort)
              do ig=1,ngroup
	          v=tshift(:,iatom,ig)-tshift(:,jatom,ig)
	          do k=1,nqdiv
	            phase=pi2*dot_product(pnt(:,k),v)
	            cphase(k,iatom0,jatom0,ig)=dcmplx(cos(phase),
     & 	                                            sin(phase))
	          enddo
	        enddo
	      enddo
	    enddo
c -----------------------------------------------------------------
          do l1=0,lmpb(jsort)
            nl1=l1+l1+1
c ----- Store transformation matrix -------------------------------
            allocate(aj(nl1,nl1,ngroup))
            ind=l1*(nl1-2)*nl1/3
            do j=1,nl1
              do i=1,nl1
                ind=ind+1
                aj(i,j,1:ngroup)=u(ind,1:ngroup)
              enddo
            enddo
            do ig=1,ngroup
              aj(:,:,ig)=transpose(aj(:,:,ig))
            enddo
c ------------------------------------------------------------------                      
            do l=0,lmpb(isort)
              nl=l+l+1
c ----- Store transformation matrix -------------------------------
              allocate(ai(nl,nl,ngroup))
              ind=l*(nl-2)*nl/3
              do j=1,nl
                do i=1,nl
                  ind=ind+1
                  ai(i,j,1:ngroup)=u(ind,1:ngroup)
                enddo
              enddo
c ------------------------------------------------------------------
	        allocate(w0(nl,nl1,nhsort(isort),nhsort(jsort),npnt))
	        allocate(w1(nqdiv,nl,nl1,nhsort(isort),nhsort(jsort)))
	        do je=1,ntle_pb(l1,jsort)
	          do ie=1,ntle_pb(l,isort)
	            w0=(0.d0,0.d0)
	            do j=1,nl1
	              lm1=l1*l1+j
	              do i=1,nl
	                lm=l*l+i
                      do jatom0=1,nhsort(jsort)
                        jatom=iat_sort(jatom0,jsort)
                        km1=iopb(jatom)+ind_prod(je,lm1,jsort)-1
                        do iatom0=1,nhsort(isort)
                          iatom=iat_sort(iatom0,isort)
                          km=iopb(iatom)+ind_prod(ie,lm,isort)-1
	                    do ind_k=1,ndim3_k(me3_k+1)
	                      k=n3_mpi_k(me3_k+1)+ind_k
                            if(km<km1) then
                              w0(i,j,iatom0,jatom0,k)=0.5d0
     &                  *dcmplx(p_w_solid_q_tau(km,km1,ind_tau,ind_k),
     &                          p_w_solid_q_tau(km1,km,ind_tau,ind_k))
                            else if(km>km1) then
                              w0(i,j,iatom0,jatom0,k)=0.5d0
     &                  *dcmplx(p_w_solid_q_tau(km1,km,ind_tau,ind_k),
     &                         -p_w_solid_q_tau(km,km1,ind_tau,ind_k))
                            else if(km==km1) then
                              w0(i,j,iatom0,jatom0,k)=0.5d0
     &                            *p_w_solid_q_tau(km,km,ind_tau,ind_k)
                            endif
                          enddo
                        enddo
                      enddo
                    enddo   !! over i
                  enddo   !! over j
	            if(nproc_k/=1) call dgop(w0,
     &	            2*npnt*nl*nl1*nhsort(isort)*nhsort(jsort),'  +',
     &                                     comm_pnt)
c ------------------------------------------------------------------
                  w1=(0.d0,0.d0)
                  allocate(a(nl,nl1))
                  allocate(b(nl,nl1))
                  allocate(c(nl,nl1))
	            do ind_k=1,ndim3_kk(me3_k+1)
	              k=n3_mpi_kk(me3_k+1)+ind_k
                    k0=i_kref(k)
	              ig=k_group(k)
	              call zone1_number(pnt(:,k),rb0,ndiv,ind)
                    do jatom0=1,nhsort(jsort)
                      jatom=iat_sort(jatom0,jsort)
                      jat=ip(jatom,ig)
                      jat0=0
                      do j=1,nhsort(jsort)
                        if(jat==iat_sort(j,jsort)) jat0=j
                      enddo
                      do iatom0=1,nhsort(isort)
                        iatom=iat_sort(iatom0,isort)
                        iat=ip(iatom,ig)
                        iat0=0
                        do i=1,nhsort(isort)
                          if(iat==iat_sort(i,isort)) iat0=i
                        enddo
                        a=real(w0(:,:,iat0,jat0,k0))
                        b=matmul(ai(:,:,ig),a)
                        a=matmul(b,aj(:,:,ig))
                        c=imag(w0(:,:,iat0,jat0,k0))
                        b=matmul(ai(:,:,ig),c)
                        c=matmul(b,aj(:,:,ig))
                        w1(ind,:,:,iatom0,jatom0)=
     &                    cphase(k,iatom0,jatom0,ig)*dcmplx(a,c)
                      enddo  !! over iatom0
                    enddo   !! over jatom0
                  enddo   !! over ind_k
                  deallocate(a,b,c)
c -----------------------------------------------------------------
	            if(nproc_k/=1) call dgop(w1,
     &	            2*nqdiv*nl*nl1*nhsort(isort)*nhsort(jsort),'  +',
     &                                     comm_pnt)
                  do jatom0=1,nhsort(jsort)
                    do iatom0=1,nhsort(isort)
	                call fft3(ndiv(1),ndiv(2),ndiv(3),nl*nl1,
     &	                        w1(1,1,1,iatom0,jatom0),1)
	              enddo
	            enddo
	            w1=w1/dfloat(nqdiv)
	            do ind_ir0=1,ndim_nrr(me3_k+1)
	              ir0=n_mpi_nrr(me3_k+1)+ind_ir0
	              ir=nrr_red(1,ir0)
	              iatom=nrr_red(2,ir0)
	              if(is(iatom)/=isort) cycle
	              jatom=nrr_red(3,ir0)
	              if(is(jatom)/=jsort) cycle
                    jatom0=0
                    do j=1,nhsort(jsort)
                      if(jatom==iat_sort(j,jsort)) jatom0=j
                    enddo
                    iatom0=0
                    do i=1,nhsort(isort)
                      if(iatom==iat_sort(i,isort)) iatom0=i
                    enddo
	              do j=1,nl1
	                lm1=l1*l1+j
                      km1=ind_prod(je,lm1,jsort)
	                do i=1,nl
	                  lm=l*l+i
                        km=ind_prod(ie,lm,isort)
	                  w_rs_mm(km,km1,ind_ir0)=w1(ir,i,j,iatom0,jatom0)
	                enddo
	              enddo
	            enddo  !! over ind_ir0
	          enddo   !! over ie
	        enddo   !! over je
	        deallocate(ai,w0,w1)
	      enddo   !! over l
	      deallocate(aj)
	    enddo   !! over l1
	    deallocate(cphase)
	  enddo   !! over isort
	enddo   !! over jsort
      end
