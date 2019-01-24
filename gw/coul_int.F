      real*8 function coul_int(i,j,k,l,ispin,jspin)
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: i,j,k,l,ispin,jspin
	integer :: ii,iii,indx,icg1,icg2,icg,lm2,i0,j0,k0,l0,lmi,lmj,lmk,
     &	         lml,ij,kl,ii5,iii5,indx5,icg15,icg25,icg5,lm25,l2,lget
      complex*16 :: cc,cgg
      cgg=(0.0d0,0.0d0)
      i0=i_num(i)
      j0=i_num(j)
      k0=i_num(k)
      l0=i_num(l)
      lmi=lm_num(i)
      lmj=lm_num(j)
      lmk=lm_num(k)
      lml=lm_num(l)
      ij=ind_pair(i0,j0)
      kl=ind_pair(k0,l0)
	cc=(0.d0,0.d0)
	if(irel.le.1) then
	  ii = max0(lmi,lmj)
	  iii = min0(lmi,lmj)
        indx = (ii*(ii-1))/2 + iii
        icg1 = indxcg(indx)
        icg2 = indxcg(indx+1) - 1 
	else
	  icg1=indxcgr(lmi,lmj,1,0)
        icg2=indxcgr(lmi,lmj,2,0)
	endif
      do icg = icg1, icg2
        if(irel.le.1) then
	    lm2 = jcg(icg)
	  else
	    lm2=jcgr(icg)
	  endif
	  if(irel.le.1) then
	    ii5 = max0(lml,lmk)
	    iii5 = min0(lml,lmk)
          indx5 = (ii5*(ii5-1))/2 + iii5
          icg15 = indxcg(indx5)
          icg25 = indxcg(indx5+1) - 1 
	  else
	    icg15=indxcgr(lml,lmk,1,0)
          icg25=indxcgr(lml,lmk,2,0)
	  endif
        do icg5 = icg15, icg25
          if(irel.le.1) then
	      lm25 = jcg(icg5)
	    else
	      lm25=jcgr(icg5)
	    endif
	    if(lm25==lm2) then
	      if(irel.le.1) then
	        cgg=cg(icg)*cg(icg5)
	      else if(irel.eq.2) then
	        cgg=cgr(icg)*cgr(icg5)
	      endif
	      l2=lget(lm2)
	      cc=cc+cgg*dbl_int(ij,kl,l2,ispin,jspin)
	    endif
	  enddo
	enddo
	coul_int=cc
      end
