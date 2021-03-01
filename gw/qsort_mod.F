      MODULE qsort_mod
      IMPLICIT NONE
      PUBLIC :: xorder_qsort,qsortfinal
      PRIVATE :: Partition, qsortc
      
      CONTAINS

      subroutine xorder_qsort(nin,arrayin, nout,barray, iref)
      implicit none
      integer, intent(in) :: nin
      real(8), intent(in) :: arrayin(nin)
      
      integer, intent(out) :: nout
      real(8), intent(out) :: barray(nin)
      integer, intent(out) :: iref(nin)
      
      integer :: ii
      integer :: arrayindex(nin)
      real(8) :: arrayout(nin)
      
      
      call qsortfinal(arrayin, arrayout, arrayindex)
      
      iref(arrayindex(1))=1
      barray(1)=arrayout(1)
      nout=1
      do ii=2, nin
        if (arrayout(ii)-arrayout(ii-1) .lt. 1.0d-8) then
          iref(arrayindex(ii))=nout        
        else
          nout=nout+1
          barray(nout)=arrayout(ii)
          iref(arrayindex(ii))=nout
        endif
      enddo
      end subroutine xorder_qsort
      
      subroutine qsortfinal(arrayin, arrayout, arrayindex)
      implicit none
      real(8), intent(in) :: arrayin(:)
      real(8), intent(out) :: arrayout(Size(arrayin))
      integer, intent(out) :: arrayindex(size(arrayin))
      
      integer :: ii
      
      arrayout=arrayin
      
      do ii=1, size(arrayin)
        arrayindex(ii)=ii
      enddo
      call QsortC(arrayout, arrayindex)
      end subroutine qsortfinal
      
      RECURSIVE SUBROUTINE QsortC(A, index_A)
!     Input/output variables
! A    	: real(precision) array to sort
!     index_A   : integer array indexing elements of A
      implicit none
      REAL(8), INTENT(in out), DIMENSION(:) :: A
      INTEGER, INTENT(in out), DIMENSION(SIZE(A)) :: index_A
      INTEGER :: iq
      
      IF(SIZE(A) > 1) THEN
        CALL Partition(A,index_A, iq)
        CALL QsortC(A(:iq-1),index_A(:iq-1))
        CALL QsortC(A(iq:),index_A(iq:))
      ENDIF
      END SUBROUTINE QsortC
      
      SUBROUTINE Partition(A, index_A, marker)
      implicit none
      REAL(8), INTENT(in out), DIMENSION(:) :: A
      INTEGER, INTENT(in out), DIMENSION(SIZE(A)) :: index_A
      INTEGER, INTENT(out) :: marker
      INTEGER :: i, j
      INTEGER :: index_temp
      REAL(8) :: temp
      REAL(8) :: x              ! pivot point
      x = A(1)
      i= 0
      j= SIZE(A) + 1
      
      DO
        j = j-1
        DO
          IF (A(j) <= x) EXIT
          j = j-1
        END DO
        i = i+1
        DO
          IF (A(i) >= x) EXIT
          i = i+1
        END DO
        IF (i < j) THEN
! exchange A(i) and A(j)
          temp = A(i)
          A(i) = A(j)
          A(j) = temp
          
          index_temp = index_A(i)
          index_A(i) = index_A(j)
          index_A(j) = index_temp
          
        ELSEIF (i == j) THEN
          marker = i+1
          RETURN
        ELSE
          marker = i
          RETURN
        ENDIF
      END DO
      
      END SUBROUTINE Partition
      
      END MODULE 
      
