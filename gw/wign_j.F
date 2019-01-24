      complex*16 function wign_j(j,m1,m2,al,bet,gam,fac)
      implicit none
      integer :: j,m1,m2
      real*8 :: fac(0:*),al,bet,gam,dm_j,d,a,b
      complex*16 exa,exg
      d=dm_j(j,m1,m2,bet,fac)
      a=m1*al/2.d0
      b=m2*gam/2.d0
      exa=dcmplx(0.d0,a)
      exa=exp(exa)
      exg=dcmplx(0.d0,b)
      exg=exp(exg)
      wign_j=exa*d*exg
      end
