      subroutine seteny
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
      integer :: ispin,isort,ndimr,li,l,ie,lm,ndim,m,k,ibnd,iatom,
     &	         km,n_bnd0,n_bnd1,ind_k,in,jn,lm0
      real*8 :: tos0,tos1,admixe,enynew,tm,adm
	admixe=0.5d0
	adm=0.5d0
	if(iter>iter_dft) then
	  admixe=admix_e   !! for IDMD=0
	  adm=0.4d0        !! for IDMD=4
	endif
      do ispin=1,nspin
        do isort=1,nsort
          ndimr=nrel*lmb(isort)+1
          lm0=0
	    do li=1,ndimr
            l=0
	      if(irel.le.1) l=li-1
	      if(irel.eq.2) l=li/2
            if(irel.ne.2) ndim=li+li-1
            if(irel.eq.2) ndim=((li+1)/2)*2
	      do ie=1,ntle(l,isort)
			  in=1
			  if(augm(ie,l,isort)/='LOC') in=2
	        tos0=0.d0
	        tos1=0.d0
	        do ind_k=1,ndim3_k(me3_k+1)
	          k=n3_mpi_k(me3_k+1)+ind_k
	          if(augm(ie,l,isort)=='LOC') then
			      n_bnd0=1
			      n_bnd1=n_semic
	          else if(augm(ie,l,isort)/='LOC') then
			      n_bnd0=n_semic+1
			      n_bnd1=n_bnd(k,ispin)
	          endif
	          do ibnd=n_bnd0,n_bnd1
	            tm=0.d0
	            do iatom=1,natom
	              if(is(iatom)==isort) then
	                do jn=1,in
	                  lm=lm0
	                  do m=1,ndim
	                    lm=lm+1
	                    km=io_lem(iatom)+indbasa(jn,ie,lm,isort)-1
	                    tm=tm+conjg(z_bnd(km,ibnd,ind_k,ispin))
     &		 		                 *z_bnd(km,ibnd,ind_k,ispin)
                        enddo
	                enddo
	              endif
	            enddo
	            tm=tm/nhsort(isort)*wgt(k)
 	            tos0=tos0+g_full_00(ibnd,ind_k,ispin)*tm
	            tos1=tos1+e_bnd(ibnd,k,ispin)
     &			           *g_full_00(ibnd,ind_k,ispin)*tm
	          enddo
	        enddo
	        if(nproc_k/=1) then
                call DGOP(tos0,1,'  +',comm_pnt)
                call DGOP(tos1,1,'  +',comm_pnt)
	        endif
              if(idmd(ie,l,isort)==0) then
	          if(abs(tos0).gt.1.d-8) then
                  enynew=tos1/tos0
                  enynew=admixe*enynew
     &			        +(1.d0-admixe)*eny(ie,li,isort,ispin)
                  eny(ie,li,isort,ispin)=enynew
	          endif
	        endif
              if(idmd(ie,l,isort)==4) then
                enynew=chem_pot
                enynew=adm*enynew+(1.d0-adm)*eny(ie,li,isort,ispin)
                eny(ie,li,isort,ispin)=enynew
	        endif
	      enddo   !! over ie
	      lm0=lm0+ndim
	    enddo
	  enddo
	enddo
      end
