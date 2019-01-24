      subroutine integral_m_j(ie,l,isort,qg,res)
	use atom_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: isort,ie,l
      real*8, intent(in) :: qg
      real*8, intent(out) :: res
	integer :: i,j,mtpb,n_period,n,n0,ii,m,m6,i0
      real*8 :: pi2,qq,rv,dqdall,pi2a,h_small,hh,am,h0,ff
	real*8, allocatable :: bess(:),bessd(:),work(:)
	n=nrad(isort)
	if(qg<1.d-10) then
	  if(l==0) then
	    j=ind_prod0(ie,l,isort)
	    mtpb=ind_pb(j,isort)
          allocate(work(0:n))
          do i=0,n
            rv=r(i,isort)
            work(i)=rv*rv*dr(i,isort)*fun_pb(mtpb+i)
          enddo
          res=dqdall(h(isort),work,n)
          deallocate(work)
	  else
	    res=0.d0
	  endif
	else
	  n_period=10
	  allocate (bess(0:l),bessd(0:l))
	  pi2=pi+pi
	  pi2a=pi2/par
	  qq=qg*pi2a
	  h_small=pi2/n_period
c --- We find where the subintervals begin ----------------------
        n0=n
        do i=1,n
          hh=r(i,isort)-r(i-1,isort)
          if(qq*hh>h_small) then
            n0=i-1
            goto 1
          endif
        enddo
c ----- Make n0 multiple of 6 -----------------------------------
1       n0=n0-mod(n0,6)
c ---------------------------------------------------------------
	  j=ind_prod0(ie,l,isort)
	  mtpb=ind_pb(j,isort)
c ---- Part of integral without sub-divisions -------------------
        allocate(work(0:n0))
        do i=0,n0
          rv=r(i,isort)
          call BESSR(qq,rv,bess,bessd,l)
          work(i)=rv*rv*dr(i,isort)*bess(l)*fun_pb(mtpb+i)
        enddo
        res=dqdall(h(isort),work,n0)
        deallocate(work)
c ---- Part of integral with sub-divisions ----------------------
        do ii=n0+1,n
          hh=r(ii,isort)-r(ii-1,isort)
          am=qq*hh/h_small/6
          m=am
          if(dfloat(m)<am) m=m+1
          m6=6*m
          h0=hh/m6
          i0=max(0,ii-3)
          i0=min(i0,n-5)
          allocate(work(0:m6))
          do i=0,m6
            rv=r(ii-1,isort)+h0*i
            call BESSR(qq,rv,bess,bessd,l)
            call fint(fun_pb(mtpb+i0),r(i0,isort),6,rv,ff)
            work(i)=rv*rv*bess(l)*ff
          enddo
          res=res+dqdall(h0,work,m6)
          deallocate(work)
        enddo
	  deallocate (bess,bessd)
	endif
      end