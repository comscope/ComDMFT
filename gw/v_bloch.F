      subroutine v_bloch(a,pt,ispin,ndim,ndim0,nbas,indg,evbn,zbn,p,pj,
     &                   w,wi)
      use atom_mod
	use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
	integer, intent(in) :: ispin,ndim,ndim0,nbas,indg(nbndtm)
	real*8, intent(in) :: pt(3),p(maxel,maxel,natom,nspin),
     &	                    pj(maxel,maxel,natom),w(nplwdif,nspin),
     &                      wi(nplwdif,nspin)
	complex*16, intent(in) :: evbn(nbasmpw,nbndf),zbn(nfun,nbndf)
	real*8, intent(out) :: a(ndim0,ndim0)
c ------- ME from MT -----------------------------------------------
	call mt_me(ispin,a,ndim,ndim0,zbn,p,pj)
c ------- ME from Interstitial --------------------------------------------
	call int_me(a,pt,ispin,ndim,ndim0,nbas,indg,evbn,w,wi)
      end
