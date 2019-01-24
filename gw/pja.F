      real*8 function pja(l,m,k,t,fac)
      integer :: l,m,k,i,j1,j2,i1,i2
      real*8 t,dl,a,p,fac(0:*)
      dl=fac(l+m)/(fac(l)*fac(m))
      i1=l+1
      i2=m+k+l
      p=1.d0-t
      a=1.d0
      do i=l,1,-1
        j1=(i1-i)*(i2+i)
        j2=i*(m+i)
        j2=j2+j2
        a=1.d0-a*p*j1/j2
      enddo
      pja=dl*a
      end
