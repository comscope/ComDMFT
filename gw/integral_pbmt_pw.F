      subroutine integral_pbmt_pw(mt_pw,mt_pw_0,v)
	use atom_mod
	use solid_mod
	use units_mod
      implicit none
      real*8, intent(in) :: v(3)
      complex*16, intent(in) :: mt_pw_0(maxlfpb,nsort)
      complex*16, intent(out) :: mt_pw(n_pbmt)
	integer :: iatom,isort,lm,l,m,ie,km,ind0,ll1,i
      real*8 :: pi4,sqpi4,qq,pi2,phase
      complex*16 :: cc
	real*8, allocatable :: ylm(:)
	allocate (ylm(limlpb))
	pi2=pi+pi
	pi4=4.d0*pi
	sqpi4=sqrt(pi4)
	qq=sqrt(dot_product(v,v))
	if(qq.gt.1.d-9) then
	  call sphharm(v(1)/qq,v(2)/qq,v(3)/qq,maxpb,ylm)
	else
	  ylm=0.d0
	  ylm(1)=1.d0/sqpi4
	endif
	do iatom=1,natom
	  isort=is(iatom)
	  ind0=iopb(iatom)-1
	  phase=pi2*dot_product(v,tau(:,iatom))
	  cc=dcmplx(cos(phase),sin(phase))
        do l=0,lmpb(isort)
	    ll1=l*(l+1)+1
	    do ie=1,ntle_pb(l,isort)
	      i=ind_prod0(ie,l,isort)
            do m=-l,l
	        lm=ll1+m
	        km=ind0+ind_prod(ie,lm,isort)
              mt_pw(km)=cc*ylm(lm)*mt_pw_0(i,isort)
	      enddo
	    enddo  !! over ie
	  enddo !! over l
	enddo   !! over isort
	deallocate (ylm)
      end