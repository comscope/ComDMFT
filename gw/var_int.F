   	subroutine var_int(ww,p,q,p0,q0,f,g,fg,r,rl,work,nr,h)
	use manager_mod
	use units_mod
	implicit none
	integer :: nr,ir
	real*8 ::  ww(nr),p(nr),q(nr),p0(nr),q0(nr),f(nr),g(nr),fg(nr),
     &	         rl(nr),r(nr),work(nr),c2,h
	c2=clight**2
	do ir=1,nr
	  work(ir)=p(ir)*p0(ir)
	enddo
	if(irel>=1) then
	  do ir=1,nr
	    work(ir)=work(ir)+q(ir)*q0(ir)/c2
	  enddo
	endif
	do ir=1,nr
	  work(ir)=ww(ir)*work(ir)*h
	enddo
	do ir=1,nr
	  f(ir)=work(ir)*rl(ir)
	enddo
	call intalt(nr,f,g)
	fg(1)=0.d0
	do ir=2,nr
	  fg(ir)=g(ir)/rl(ir)/r(ir)
	enddo
	f(1)=0.d0
	do ir=2,nr
	  f(ir)=work(ir)/r(ir)/rl(ir)
	enddo
	call intalt(nr,f,g)
	do ir=1,nr
	  fg(ir)=fg(ir)+(g(nr)-g(ir))*rl(ir)
	enddo
	end
	
	
	
	
	subroutine var_int1(ww,p,f,g,fg,r,rl,work,nr,h)
	implicit none
	integer :: nr,ir
	real*8 ::  ww(nr),p(nr),f(nr),g(nr),fg(nr),rl(nr),r(nr),work(nr),h
	do ir=1,nr
	  work(ir)=ww(ir)*p(ir)*h*r(ir)**2
	enddo
	do ir=1,nr
	  f(ir)=work(ir)*rl(ir)
	enddo
	call intalt(nr,f,g)
	do ir=2,nr
	  fg(ir)=g(ir)/rl(ir)/r(ir)
	enddo
	f(1)=0.d0
	do ir=2,nr
	  f(ir)=work(ir)/r(ir)/rl(ir)
	enddo
	call intalt(nr,f,g)
	do ir=1,nr
	  fg(ir)=fg(ir)+(g(nr)-g(ir))*rl(ir)
	enddo
	end