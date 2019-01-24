      subroutine plwovi(g1,g2,g3,pvr,pvi)
	use atom_mod
	use manager_mod
	use solid_mod
	use units_mod
      implicit none
	integer :: isort,iatom
	real*8, intent(in) :: g1,g2,g3
	real*8,intent(out) :: pvr,pvi
	real*8 :: c1,pi2,c2,c4,g,ss,sr,fr,dot,ffi,x,sbes,si
      c1=amegaint/amega
      pi2=2.d0*pi
      c2=4.d0*pi/amega
      c4=pi2/par
      pvr=0.d0
      if(complex_ro) pvi=0.d0
      g=c4*sqrt(g1*g1+g2*g2+g3*g3)
      if(g.lt.1.d-6) then
        pvr=c1
      else
        ss=-c2/g
        sr=0.d0
        if(complex_ro) si=0.d0
        do isort=1,nsort
          fr=0.d0
          if(complex_ro) ffi=0.d0
          do iatom=1,natom
            if(is(iatom).eq.isort) then
              dot=pi2*(g1*tau(1,iatom)+g2*tau(2,iatom)+g3*tau(3,iatom))
              fr=fr+cos(dot)
              if(complex_ro) ffi=ffi+sin(dot)
            endif
          enddo
          x=g*smt(isort)
          sbes=smt(isort)**2*(sin(x)/x-cos(x))/x
          sr=sr+fr*sbes
          if(complex_ro) si=si+ffi*sbes
        enddo
        pvr=ss*sr
        if(complex_ro) pvi=ss*si
      endif
      end