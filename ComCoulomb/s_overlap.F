      subroutine s_overlap(key,s,n0,n,iref)

c     Calculates interstitial overlap matrix for plane waves ------
c     IREF - reference of the plane waves to GBS array ------------
c     key = 0  Just overlap matrix
c     key = 1  Its inverse

!     if key=0 \int_{unitcell} dr exp{-i(G_1-G_2)r}
!     if key=1 inverse of \int_{unitcell} dr exp{-i(G_1-G_2)r}            

      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod
      implicit none
      include 'mpif.h'
      integer, intent(in) :: key,n,n0,iref(n0)
      complex*16, intent(out) :: s(n0,n0)
      integer :: j,j0,i,i0,ind,info,gbs_number,mpierr
      double precision :: v(3)
      integer, allocatable :: ipiv(:)
      complex*16, allocatable :: work(:)

!      write(100,*)'key',key,'s',s,'n0',n0,'iref',iref,
!     $ 'amega',amega,'complex_ro',complex_ro,'gbs',gbs,
!     $ 'sovi',sovi,'sovr',sovr
!      call mpi_bcast(irel,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)

      s=0.0d0

      do j=1,n
        j0=iref(j)
        do i=1,n
          i0=iref(i)
          v=gbs(:,j0)-gbs(:,i0)
          ind=gbs_number(v)
          if(complex_ro) then
            s(i,j)=dcmplx(sovr(ind),sovi(ind))
          else
            s(i,j)=dcmplx(sovr(ind),0.d0)
          endif
        enddo
      enddo
      s=s*amega
      if(key==1) then
        allocate(ipiv(n),work(3*n))
        call zgetrf(n,n,s,n0,ipiv,info)
        call zgetri(n,s,n0,ipiv,work,3*n,info)
        deallocate(ipiv,work)
      endif

!      write(100,*)'key',key,'s',s,'n0',n0,'iref',iref,
!     $ 'amega',amega,'complex_ro',complex_ro,'gbs',gbs,
!     $ 'sovi',sovi,'sovr',sovr,'gbas',gbas,'indplw',indplw
!      call mpi_bcast(irel,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
      end

