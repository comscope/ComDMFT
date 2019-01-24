      subroutine mxmymz(kin,k,lx,ly,lz)
c- do mirrors in x,y,z if lx,ly,lz=1, respectively
c ----------------------------------------------------------------
ci inputs
ci
co outputs
co
cr remarks
cr
c ----------------------------------------------------------------
      integer :: kin(3),k(3),lx,ly,lz
      k(1) = kin(1)
      k(2) = kin(2)
      k(3) = kin(3)
      if (lx .eq. 1) k(1) = 1-k(1)
      if (ly .eq. 1) k(2) = 1-k(2)
      if (lz .eq. 1) k(3) = 1-k(3)
      end