      subroutine sckg2
	use atom_mod
	use manager_mod
      implicit none
	integer :: ind,l,j,mj
	real*8 :: u,u1,u2
      allocate(ckg2(2,2,(maxylm+2)*(maxylm+3)))
      ind=0
      do l=0,maxylm+1
        j=l+l+1
        do mj=-j,j,2
          ind=ind+1
          u=dfloat(mj)/dfloat(j)
          u1=sqrt(0.5d0*(1.d0+u))
          u2=sqrt(0.5d0*(1.d0-u))
          ckg2(1,1,ind)=u1
          ckg2(1,2,ind)=-u2
          ckg2(2,1,ind)=u2
          ckg2(2,2,ind)=u1
        enddo
      enddo
      end