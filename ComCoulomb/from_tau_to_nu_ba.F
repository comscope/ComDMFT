      subroutine from_tau_to_nu_ba(f_nu,n0,f_tau,n1,n)

c     ------- Transform TAU ---> NU  --------------------------
c     Input:
c     F_tau - hermitian, but is stored in square R*8 form = A
c     Output:
c     F_nu is hermitian and stored in A square R*8 form

      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,n0,n1
      double precision, intent(in) :: f_tau(n1,n1,ndim3_tau)
      double precision, intent(out) :: f_nu(n0,n0,ndim3_nu)
      integer :: nn,i,j,ij,ind_tau,nblock,nstep,ind_nu,j0,jm,i_nu,i1,
     $  m0,m
      double precision :: yy1
      double precision, allocatable :: at(:,:),an(:,:)
      nblock=1000
      nstep=max(1,nblock/n)
      nn=n*nstep
      allocate(at(nn,ndim3_tau))
      allocate(an(nn,0:n_nu_d))
      j0=0
      m=interp_nu_d
      m0=m-1      
      do while(j0<n)
        jm=min(nstep,n-j0)
c     -------------------------- COS transformation ------------------
        at=0.d0
        do ind_tau=1,ndim3_tau
          ij=0
          do j=j0+1,j0+jm
            do i=1,n
              ij=ij+1
              at(ij,ind_tau)=f_tau(i,j,ind_tau)
            enddo
          enddo
        enddo
        call from_tau_to_nu_cos_d(an,at,nn)
        ij=0
        do j=j0+1,j0+jm
          do i=1,n
            ij=ij+1
            do ind_nu=1,ndim3_nu
              i_nu=me3_tau*ndim3_nu+ind_nu-1
              if(i_nu<=n_nu_exa) then
                f_nu(i,j,ind_nu)=an(ij,i_nu)
              else if(i_nu<n_nu) then
                i1=ind_nu_int(i_nu)
                call fint(an(ij,i1:i1+m0),w_nu_d(i1),m,w_nu(i_nu),yy1)
                f_nu(i,j,ind_nu)=yy1
              else if(i_nu==n_nu) then
                f_nu(i,j,ind_nu)=an(ij,n_nu_d)
              endif
            enddo
          enddo
        enddo
        j0=j0+jm
      enddo
      deallocate(at,an)
      end
