      subroutine mklegw(n,z,w)
c- quadrature weights for legendre gaussian quadrature
c ----------------------------------------------------------------------
ci inputs
ci   n:   number of mesh points for numerical integration
ci        (on interval -1, 1)
ci   ipr: verbosity
co outputs
co   z,w
cr remarks
cr   integrates 2(n-1) order polynomial exactly on interval (-1,1)
c ----------------------------------------------------------------------
      implicit none
c passed parameters
      integer, intent(in) :: n
      real*8, intent(out) :: z(n),w(n)
      real*8 :: plegn,pnp
c local parameters
      integer :: iroot,in,k
      real*8 :: root,delta,machep,epsmas,pi
      epsmas=1.d-15

      pi = 4 * datan(1d0)
      machep = 10*epsmas
c --- find all the roots of p_n ---
      do iroot = 1, n
        z(iroot) = dcos(pi*(2*iroot-.5d0)/(2*n+1))
        root = z(iroot)
  100   continue
        delta = -plegn(n,root)/pnp(n,root)
        root = root + delta
        if (dabs(delta) .gt. dabs(machep*root)) goto 100
        z(iroot) = root
	enddo
c --- make the weights ---
      do in = 1, n
        w(in) = 2/((1-z(in)**2)*(pnp(n,z(in))**2))
	enddo
c	 ----- Ordering ------------------------------------
	do iroot=1,n
	  delta=z(iroot)
	  k=iroot
	  do in=iroot+1,n
	    if(z(in)<delta) then
	      delta=z(in)
	      k=in
	    endif
	  enddo
	  if(k/=iroot) then
	    z(k)=z(iroot)
	    z(iroot)=delta
	    delta=w(k)
	    w(k)=w(iroot)
	    w(iroot)=delta
	  endif
	enddo
      end