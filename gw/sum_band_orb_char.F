      double precision function sum_band_orb_char(ibnd,ind_k,ispin,key)
c     key = 0 For current E_BND band structure (may be DFT, HF, or X)
c     key = 1 For E_QP band structure
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer, intent(in) :: ibnd,ind_k,ispin,key
	logical :: cr1,cr2
	integer :: isort,iatom,l,ie,lm,km,i,lget,je,kmj,in,jn,in1,
     &	         jn1,mj,li,isp
	real*8 :: pdos,zz
      isp=min(ispin,nspin_1)
      pdos=0.d0
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
	      cr2=(correlated(je,l,isort)=='L')
	      if(cr2) then
			  in1=1
			  if(augm(je,l,isort)/='LOC') in1=2
			  do jn1=1,in1
	          kmj=io_lem(iatom)-1+indbasa(jn1,je,lm,isort)
	          do ie=1,ntle(l,isort)
	            cr1=(correlated(ie,l,isort)=='L')
                    zz=0.0d0
	            if(cr1) then
			        in=1
			        if(augm(ie,l,isort)/='LOC') in=2
			        do jn=1,in
	                km=io_lem(iatom)-1+indbasa(jn,ie,lm,isort)
	                if(key==0) then
		  	            zz=conjg(z_bnd(km,ibnd,ind_k,ispin))
     &		  	                *z_bnd(kmj,ibnd,ind_k,ispin)
	                else if(key==1) then
		  	            zz=conjg(z_qp(km,ibnd,ind_k,ispin))
     &		  	                *z_qp(kmj,ibnd,ind_k,ispin)
                      endif
	                pdos=pdos+zz*ffsmt(jn,jn1,ie,je,li,isort,isp)
	              enddo
	            endif
	          enddo
	        enddo
	      endif
	    enddo
	  enddo
	enddo
	sum_band_orb_char=pdos
      end
