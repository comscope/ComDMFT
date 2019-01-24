      subroutine sigc_mm_k_from_r(ispin,ind_tau,sig_mm,key)
c     key = 0 : (G_tau)*(W_tau) contribution to Sig_C -----------------
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	integer, intent(in) :: ispin,ind_tau,key
	complex*16, intent(in) :: sig_mm(nfun,nfun,ndim3_kk(me3_k+1),2)
	integer :: ind,ind_k,k,nbnd,i_t,j,ind_ir,ir,i
	complex*16, allocatable :: tmp(:,:),tmp1(:,:),tmp2(:,:),p_c(:,:),
     &	                         pcc(:,:,:)
	allocate(p_c(nqdiv,nfun))
	allocate(pcc(nfun,nfun,maxval(ndim3_k)))
	allocate(tmp(nfun,nfun),tmp1(nfun,nbndf))
	allocate(tmp2(nbndf,nbndf))
	do i_t=1,2   !! here is over TAU and BETA-TAU
	  do j=1,nfun
	    p_c=(0.d0,0.d0)
	    do ind_ir=1,ndim3_kk(me3_k+1)
	      ir=n3_mpi_kk(me3_k+1)+ind_ir
	      do i=1,nfun
	        p_c(ir,i)=sig_mm(i,j,ind_ir,i_t)
	      enddo
	    enddo
	    if(nproc_k/=1) call dgop(p_c,2*nqdiv*nfun,'  +',comm_pnt)
	    call fft3(ndiv(1),ndiv(2),ndiv(3),nfun,p_c,-1)
	    do ind_k=1,ndim3_k(me3_k+1)
	      k=n3_mpi_k(me3_k+1)+ind_k
	      call zone1_number(pnt(:,k),rb0,ndiv,ind)
	      do i=1,nfun
	        pcc(i,j,ind_k)=p_c(ind,i)
	      enddo
          enddo  !! over ind_k
	  enddo  !! over j
        do ind_k=1,ndim3_k(me3_k+1)
	    k=n3_mpi_k(me3_k+1)+ind_k
	    nbnd=n_bnd(k,ispin)
	    call zgemm('n','n',nfun,nbnd,nfun,(1.d0,0.d0),pcc(1,1,ind_k),
     &	             nfun,z_bnd(1,1,ind_k,ispin),nfun,(0.d0,0.d0),tmp1,
     &               nfun)
	    call zgemm('c','n',nbnd,nbnd,nfun,(1.d0,0.d0),
     &	             z_bnd(1,1,ind_k,ispin),nfun,tmp1,nfun,(0.d0,0.d0),
     &               tmp2,nbndf)
          if(i_t==1) then
            if(key==0) then
              do j=1,nbnd
                do i=1,j-1
                  sig_c_tau(i,j,1,ind_tau,ind_k,ispin)=
     &              sig_c_tau(i,j,1,ind_tau,ind_k,ispin)+real(tmp2(i,j))
                  sig_c_tau(i,j,2,ind_tau,ind_k,ispin)=
     &              sig_c_tau(i,j,2,ind_tau,ind_k,ispin)+real(tmp2(i,j))
                  sig_c_tau(j,i,1,ind_tau,ind_k,ispin)=
     &              sig_c_tau(j,i,1,ind_tau,ind_k,ispin)+imag(tmp2(i,j))
                  sig_c_tau(j,i,2,ind_tau,ind_k,ispin)=
     &              sig_c_tau(j,i,2,ind_tau,ind_k,ispin)+imag(tmp2(i,j))
                enddo
                sig_c_tau(j,j,1,ind_tau,ind_k,ispin)=
     &            sig_c_tau(j,j,1,ind_tau,ind_k,ispin)+real(tmp2(j,j))
                sig_c_tau(j,j,2,ind_tau,ind_k,ispin)=
     &            sig_c_tau(j,j,2,ind_tau,ind_k,ispin)+real(tmp2(j,j))
              enddo
            endif
          else if(i_t==2) then
            if(key==0) then
              do j=1,nbnd
                do i=1,j-1
                  sig_c_tau(i,j,1,ind_tau,ind_k,ispin)=
     &              sig_c_tau(i,j,1,ind_tau,ind_k,ispin)-real(tmp2(i,j))
                  sig_c_tau(i,j,2,ind_tau,ind_k,ispin)=
     &              sig_c_tau(i,j,2,ind_tau,ind_k,ispin)+real(tmp2(i,j))
                  sig_c_tau(j,i,1,ind_tau,ind_k,ispin)=
     &              sig_c_tau(j,i,1,ind_tau,ind_k,ispin)-imag(tmp2(i,j))
                  sig_c_tau(j,i,2,ind_tau,ind_k,ispin)=
     &              sig_c_tau(j,i,2,ind_tau,ind_k,ispin)+imag(tmp2(i,j))
                enddo
                sig_c_tau(j,j,1,ind_tau,ind_k,ispin)=
     &            sig_c_tau(j,j,1,ind_tau,ind_k,ispin)-real(tmp2(j,j))
                sig_c_tau(j,j,2,ind_tau,ind_k,ispin)=
     &            sig_c_tau(j,j,2,ind_tau,ind_k,ispin)+real(tmp2(j,j))
              enddo
            endif
          endif
	  enddo  !! over ind_k
	enddo  !! over i_t
	deallocate(tmp,tmp1,tmp2,p_c,pcc)
      end
