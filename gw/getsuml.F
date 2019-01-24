      subroutine getsuml(zc,qmt,qmti)
	use atom_mod
	use manager_mod
	use solid_mod
      implicit none
      integer :: iatom,isort,indt,ind0,k1,i1,k,i
      real*8 :: qmt(maxqmt,natom),qmti(maxqmt,natom)
	complex*16 :: zc(nfun),qtmp
      do iatom=1,natom
        isort=is(iatom)
        indt=0
	  ind0=io_lem(iatom)-1 !! zero index of the orbitals in atom IATOM
        do k1=1,lfunm(isort)
          i1=ind0+k1
          do k=1,k1
            i=ind0+k
            indt=indt+1
	      qtmp=conjg(zc(i))*zc(i1)
            qmt(indt,iatom)=dreal(qtmp)
            if(irel.eq.2) qmti(indt,iatom)=dimag(qtmp)
          enddo  !!! over k
        enddo  !!! over k1
      enddo   !!! over iatom
      end