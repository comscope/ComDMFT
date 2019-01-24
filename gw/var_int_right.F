   	subroutine var_int_right(ww,p,q,f,g,fg,r,rl,work,nr,h)
	use manager_mod
	use units_mod
	implicit none
	integer :: nr,ir
	real*8 ::  ww(nr),p(nr),q(nr),f(nr),g(nr),fg(nr),
     &	         rl(nr),r(nr),work(nr),h
	do ir=1,nr
	  work(ir)=p(ir)*q(ir)
	enddo
	do ir=1,nr
	  work(ir)=ww(ir)*work(ir)*h
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