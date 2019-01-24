      subroutine ferm_unpack_tau(f_cmpl,f_dble,n,n0,n1,key)

c     ------- Unpacking from A and B into F(t) and F(b-t)  -------------------
c     G(t)   =   1/2 ( A + B )
c     G(b-t) =   1/2 ( B - A )
c     
c     key = 1  - to get only G(t)
c     key = 2  - to get only G(b-t)
c     key = 3  - to get both of them

      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer, intent(in) :: n,n0,n1,key
      double precision, intent(in) :: f_dble(n1,n1,2)
      complex*16, intent(out) :: f_cmpl(n0,n0,*)
      integer :: i,j
      double precision :: ar,ai,br,bi

      do j=1,n
        do i=1,j-1
          ar=f_dble(i,j,1)
          ai=f_dble(j,i,1)
          br=f_dble(i,j,2)
          bi=f_dble(j,i,2)
          if(key==1) then
            f_cmpl(i,j,1)=0.5d0*dcmplx(ar+br,ai+bi)
          else if(key==2) then
            f_cmpl(i,j,1)=0.5d0*dcmplx(br-ar,bi-ai)
          else if(key==3) then
            f_cmpl(i,j,1)=0.5d0*dcmplx(ar+br,ai+bi)
            f_cmpl(i,j,2)=0.5d0*dcmplx(br-ar,bi-ai)
            f_cmpl(j,i,2)=conjg(f_cmpl(i,j,2))
          endif
          f_cmpl(j,i,1)=conjg(f_cmpl(i,j,1))
        enddo
        ar=f_dble(j,j,1)
        br=f_dble(j,j,2)
        if(key==1) then
          f_cmpl(j,j,1)=0.5d0*(ar+br)
        else if(key==2) then
          f_cmpl(j,j,1)=0.5d0*(br-ar)
        else if(key==3) then
          f_cmpl(j,j,1)=0.5d0*(ar+br)
          f_cmpl(j,j,2)=0.5d0*(br-ar)
        endif
      enddo
      end
