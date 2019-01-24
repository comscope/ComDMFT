      subroutine scg(lmax,c,cindx,js)
c  computes clebsch-gordan coefficients (formerly s104 in asw)
c  but here all is doubleprecision
	use units_mod
      implicit none
	integer, intent(in) :: lmax
	integer, intent(out) :: cindx(*),js(*)
      real*8, intent(out) :: c(*)
      integer s1,s2,s3,t1,t2,t3,nfac,nl,i,ic,lmindx,i1,l1,j1s,j1,m1,n1,
     &	    i2,l2,i31,i32,j2s,k2,j2,m2,n2,m3,mb,nm3,n3,i3,l3
      real*8 :: sr2,fs,q1,t,f102,f100
      real*8, allocatable :: fac(:)
      fs(i)=dfloat(1+4*(i/2)-2*i)
	nfac=100
      allocate(fac(nfac))
      nl=lmax+1
      sr2=dsqrt(2.d0)
      fac(1)=1.d0
      do i=1,nfac-1
        fac(i+1)=dfloat(i)*fac(i)
      enddo
      ic=0
      lmindx=0
      do i1=1,nl
        l1=i1-1
        j1s=2*l1+1
        do j1=1,j1s
          m1=j1-i1
          n1=iabs(m1)
          s1=0
          if(m1.lt.0) s1=1
          t1=0
          if(m1.eq.0) t1=1
          do i2=1,i1
            l2=i2-1
            i31=l1-l2+1
            i32=l1+l2+1
            j2s=2*l2+1
            k2=j1s*j2s
            if(i2.eq.i1) j2s=j1
            do j2=1,j2s
              lmindx=lmindx+1
              cindx(lmindx)=ic+1
              m2=j2-i2
              n2=iabs(m2)
              s2=0
              if(m2.lt.0) s2=1
              t2=0
              if(m2.eq.0) t2=1
              if(m1*m2<0) then
                m3=-n1-n2
                mb=-iabs(n1-n2)
                if(mb==0) nm3=1
                if(mb/=0) nm3=2
              else if(m1*m2==0) then
                m3=m1+m2
                mb=0
                nm3=1
              else if(m1*m2>0) then
                m3=n1+n2
                mb=iabs(n1-n2)
                nm3=2
              endif
5             n3=iabs(m3)
              s3=0
              if(m3.lt.0) s3=1
              t3=0
              if(m3.eq.0) t3=1
              q1=dsqrt(dfloat(k2))*fs(n3+(s1+s2+s3)/2)
     &                            /(2.d0*sr2**(1+t1+t2+t3))
              do i3=i31,i32,2
                l3=i3-1
                if(n3.gt.l3) cycle
                t=0.d0
                if(n1+n2.eq.-n3) t=t+f102(fac,l1,l2,l3)
                if(n1+n2.eq.n3)  t=t+f100(fac,l1,l2,l3,n1,n2,n3)
     &                              *fs(n3+s3)
                if(n1-n2.eq.-n3) t=t+f100(fac,l1,l2,l3,n1,-n2,-n3)
     &                              *fs(n2+s2)
                if(n1-n2.eq.n3)  t=t+f100(fac,l1,l2,l3,-n1,n2,-n3)
     &                              *fs(n1+s1)
                ic=ic+1
                c(ic)=q1*t*f102(fac,l1,l2,l3)
     &                    /(sqrtpi*dsqrt(dfloat(2*l3+1)))
                js(ic)=l3*(l3+1)+m3+1
              enddo
              nm3=nm3-1
              m3=mb
              if(nm3.gt.0) goto 5
            enddo
          enddo
        enddo
      enddo
      cindx(lmindx+1)=ic+1
      deallocate(fac)
      end