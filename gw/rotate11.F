	subroutine rotate11(x,y,nlm,u0,key)
	use solid_mod
      implicit none
c	key=1 Y(m)=Sum_m1 D_m_m1 X_m1
c	key=2 Y(m)=Sum_m1 D_m1_m X_m1
	integer :: nlm,key,lmsf,l,nls,m,lm,nls1,m1,l1m1,iwig,lget
      real*8 :: x(nlm),y(nlm),u0(maxwig),sum
	lmsf=lget(nlm)
	do l=0,lmsf
        nls=l*(2*l-1)*(2*l+1)/3+l+1
	  if(key.eq.1) then
          do m=-l,l
            lm=l*(l+1)+m+1
            nls1=nls+m
	      sum=0.d0
            do m1=-l,l
              l1m1=l*(l+1)+m1+1
              iwig=nls1+(2*l+1)*(l+m1)
              sum=sum+u0(iwig)*x(l1m1)
	      enddo
            y(lm)=sum
	    enddo
	  else if(key.eq.2) then
          do m=-l,l
            lm=l*(l+1)+m+1
            nls1=nls+(2*l+1)*(l+m)
	      sum=0.d0
            do m1=-l,l
              l1m1=l*(l+1)+m1+1
              iwig=nls1+m1
              sum=sum+u0(iwig)*x(l1m1)
	      enddo
            y(lm)=sum
	    enddo
	  endif
	enddo
      end