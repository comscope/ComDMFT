      subroutine n_overlap(key,ind_k,s,n0pl,npl,n0red,nred,sn)

c     Calculates overlap matrix for reduced PB ------
c     key = 0  Just overlap matrix
c     key = 1  Its inverse

      use atom_mod
      use etot_mod
      use manager_mod
      use solid_mod
      use units_mod

      implicit none
      integer, intent(in) :: key,n0pl,npl,n0red,nred,ind_k
      complex*16, intent(in) :: s(n0pl,n0pl)
      complex*16, intent(out) :: sn(n0red,n0red)
      integer :: i,info
      integer, allocatable :: ipiv(:)
      complex*16, allocatable :: tmp(:,:),tmp1(:,:),work(:)

      allocate(tmp(npl,npl))
      allocate(tmp1(nred,npl))
      tmp=s(1:npl,1:npl)
      do i=1,npl
        tmp(i,i)=tmp(i,i)-(1.d0,0.d0)
      enddo
      call zgemm('c','n',nred,npl,npl,(1.d0,0.d0),
     &  qvopt(n_pbmt+1,1,ind_k),n_pbtot,tmp,npl,(0.d0,0.d0),
     &  tmp1,nred)
      call zgemm('n','n',nred,nred,npl,(1.d0,0.d0),tmp1,nred,
     &  qvopt(n_pbmt+1,1,ind_k),n_pbtot,(0.d0,0.d0),sn,n0red)
      deallocate(tmp,tmp1)
      do i=1,nred
        sn(i,i)=sn(i,i)+(1.d0,0.d0)
      enddo
      if(key==1) then
        allocate(ipiv(nred),work(3*nred))
        call zgetrf(nred,nred,sn,n0red,ipiv,info)
        call zgetri(nred,sn,n0red,ipiv,work,3*nred,info)
        deallocate(ipiv,work)
      endif
      end
