      subroutine pack_hermit(a_compl,a_dble,n,n0,n1,alpha,beta)
      implicit none
      integer, intent(in) :: n,n0,n1
      double precision, intent(in) :: alpha,beta
      complex*16, intent(in) :: a_compl(n0,n0)
      double precision, intent(out) :: a_dble(n1,n1)
      integer :: i,j
      if(alpha==0.d0) a_dble=0.d0
      do j=1,n
        do i=1,j-1
          a_dble(i,j)=alpha*a_dble(i,j)+beta*dble(a_compl(i,j))
          a_dble(j,i)=alpha*a_dble(j,i)+beta*dimag(a_compl(i,j))
        enddo
        a_dble(j,j)=alpha*a_dble(j,j)+beta*dble(a_compl(j,j))
      enddo
      end
      
      
      
      subroutine unpack_hermit(a_compl,a_dble,n,n0,n1,alpha,beta)
!
      implicit none
      integer, intent(in) :: n,n0,n1
      complex*16, intent(in) :: alpha,beta
      double precision, intent(in) :: a_dble(n0,n0)
      complex*16, intent(out) :: a_compl(n1,n1)
      integer :: i,j

      if(alpha==(0.d0,0.d0)) a_compl=(0.d0,0.d0)

      do j=1,n
        do i=1,j-1
          a_compl(i,j)=alpha*a_compl(i,j)
     &      +beta*dcmplx(a_dble(i,j),a_dble(j,i))
          a_compl(j,i)=conjg(a_compl(i,j))
        enddo
        a_compl(j,j)=alpha*a_compl(j,j)+beta*a_dble(j,j)
      enddo
      end

