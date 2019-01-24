      subroutine local_orb(work,p0,p1,p2,q0,q1,q2,dr,n,fi,dfi,fidot,
     &	                   dfidot,fidot2,dfidot2,adet,tmp,tmp2,c2,h,
     &                     avv,adv,anorm,fidfid,fid2fid2)
	use manager_mod
	use units_mod
      implicit none
	integer, intent(in) :: n
	real*8, intent(in) :: p0(0:n),q0(0:n),p1(0:n),fidfid,fid2fid2,
     &	                    q1(0:n),p2(0:n),q2(0:n),fi,dfi,fidot,
     &	                    dfidot,fidot2,dfidot2,adet,dr(0:n),c2,h
      real*8, intent(out) :: tmp(0:n),tmp2(0:n),avv,adv,anorm
	integer :: ir
	real*8 :: fidfid2,work(0:n),dqdall
c ------- Normalization, zero value and slope on MT ------
      avv=-(fidot2*dfidot-dfidot2*fidot)/adet
      adv=-(dfidot2*fi-fidot2*dfi)/adet
c -------- Normalization --------------------------------
	do ir=0,n
	  work(ir)=p1(ir)*p2(ir)*dr(ir)
      enddo
	if(irel.ge.1) then
	  do ir=0,n
	    work(ir)=work(ir)+q1(ir)*q2(ir)*dr(ir)/c2
        enddo
	endif
      fidfid2=dqdall(h,work,n)
	anorm=avv*avv+adv*adv*fidfid+fid2fid2+2.d0*adv*fidfid2
	anorm=1.d0/sqrt(anorm)
      do ir=0,n
        tmp(ir)=anorm*(avv*p0(ir)+adv*p1(ir)+p2(ir))
      enddo
      if(irel.ge.1) then
        do ir=0,n
          tmp2(ir)=anorm*(avv*q0(ir)+adv*q1(ir)+q2(ir))
        enddo
      endif
      end