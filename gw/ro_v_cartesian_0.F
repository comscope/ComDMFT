      subroutine ro_v_cartesian_0(rrr,f_out,ti,iatom,isort,key,ispin)
c     key = 0 Total Density (sum over spins)
c     key = 1 V_Hartree - Z/r
c     key = 2 V_Full = V_Hartree + V_xc - Z/r
	use atom_mod
	use etot_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: key,iatom,isort,ispin
      real*8, intent(in) :: rrr(3),ti(3)
      real*8, intent(out) :: f_out
      integer :: l,m,lm,isym,mt,mt1,iatom1,ig,i0
      real*8 :: r1(3),r2(3),rabs,bb,fr,fr1,sqpi4,c
      real*8, allocatable :: ylm(:)
      f_out=0.d0
      r1=rrr-ti
	rabs=sqrt(dot_product(r1,r1))
      if(rabs>smt(isort)) return
      sqpi4=sqrt(4.d0*pi)
      c=-2.d0*z(isort)*sqpi4
      allocate(ylm(limlpb))
      iatom1=iat_1(isort)
      bb=smt(isort)/(exp(h(isort)*nrad(isort))-1.d0)
      if(iatom/=iatom1) then
        r2=r1
        do ig=1,ngroup
          if(ip(iatom,ig)==iatom1) then
	      call rotate(r2(1),r2(2),r2(3),r1(1),r1(2),r1(3),u(2,ig),1)
	      exit
	    endif
	  enddo
	endif
      call sphharm(r1(1)/rabs,r1(2)/rabs,r1(3)/rabs,lmpb(isort),ylm)
      do l=0,lmpb(isort)
        do m=-l,l
          lm=l*(l+1)+m+1
          if(sym(lm,isort))then
            isym=lmsym(lm,isort)
            i0=log(1.d0+rabs/bb)/h(isort)-2
            i0=max(0,i0)
            i0=min(i0,nrad(isort)-5)
            if(key==0) then
              mt=indmt(isym,isort,1)
              call fint(ro(mt+i0),r(i0,isort),6,rabs,fr)
              if(nspin==2) then
                mt1=indmt(isym,isort,nspin)
                call fint(ro(mt1+i0),r(i0,isort),6,rabs,fr1)
                fr=fr+fr1
              endif
            else if(key==1) then
              mt=indmt(isym,isort,1)
              call fint(v_mt_h(mt+i0),r(i0,isort),6,rabs,fr)
              if(l==0) then
                fr=fr+c/rabs
              endif
            else if(key==2) then
              mt=indmt(isym,isort,1)
              mt1=indmt(isym,isort,ispin)
              call fint(v_mt_h(mt+i0),r(i0,isort),6,rabs,fr)
              call fint(v_mt_xc(mt1+i0),r(i0,isort),6,rabs,fr1)
              fr=fr+fr1
              if(l==0) then
                fr=fr+c/rabs
              endif
            endif
            f_out=f_out+fr*ylm(lm)
          endif
        enddo  !! over m
      enddo   !! over l
      deallocate(ylm)
      end