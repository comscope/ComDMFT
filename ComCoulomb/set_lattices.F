      subroutine set_lattices

      use solid_mod

      implicit none
      integer :: nkmx
      double precision, allocatable :: work(:),dlat0(:,:),rlat0(:,:)

      if(allocated(dtau)) deallocate(dtau)
      if(allocated(ipdtau)) deallocate(ipdtau)
      allocate(dtau(3,natom**2),ipdtau(natom,natom))

      call convx(tau,ndtau)
      nkmx=90000
      allocate(dlat0(3,nkmx),rlat0(3,nkmx),work(nkmx))

      call lattc(dlat0,rlat0,nkmx,nkmx,work)
      if(allocated(dlat)) deallocate(dlat)
      if(allocated(rlat)) deallocate(rlat)
      allocate(dlat(3,nkd),rlat(3,nkr))
      dlat=dlat0(:,1:nkd)
      rlat=rlat0(:,1:nkr)
      deallocate(work,dlat0,rlat0)
      end
