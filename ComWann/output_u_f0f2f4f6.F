      integer function find_ind(jj, narray,array)

      use comwann_mod
      implicit none
      integer, intent(in) :: jj,narray, array(narray)
      integer :: ii
      
      find_ind=0
      do ii=1, narray
        if (jj .eq. array(ii)) then
          find_ind=ii
        endif
      enddo
      if (find_ind .eq. 0) then
        if (maswrk) write(iun,*) 'problem in find_ind'
c        write(179,*) 'problem in find_ind'
        call ending_qft
      endif
      end




      


