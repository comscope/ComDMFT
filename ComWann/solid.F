      subroutine indexing(ntot,ndivision,divisionarray,flag,n1,n2)
      use comwann_mod
      implicit none
      integer, intent(in) :: ntot,ndivision,
     $  divisionarray(ndivision),flag
      integer, intent(inout) :: n1, n2(ndivision)
      integer :: ii,jj,tempcnt,tmpsize

      tmpsize=1
      do ii=1, ndivision
        tmpsize=tmpsize*divisionarray(ii)
      enddo
      if (tmpsize .ne.  ntot) then
        if (maswrk) then
          write(993,*) 'array_division wrong'
          write(993,'(100(i10,2x))') tmpsize, ntot,
     $    divisionarray
          call flush(993)
        endif
c        write(179,*) 'array_division wrong'        
c        write(179,'(100(i10,2x))') tmpsize, ntot,
c     $    divisionarray
        call ending
      endif
      if (flag .eq. 1) then
        n1=n2(1)
        do ii=2, ndivision
          tempcnt=1
          do jj=1, ii-1
            tempcnt=tempcnt*divisionarray(jj)
          enddo
          n1=n1+(n2(ii)-1)*tempcnt
        enddo
      else
        n2=0
        tempcnt=n1
        do ii=1, ndivision-1
          n2(ii)=tempcnt-((tempcnt-1)/divisionarray(ii))
     $      *divisionarray(ii)
          tempcnt=(tempcnt-n2(ii))/divisionarray(ii)+1
        enddo
        n2(ndivision)=tempcnt
      endif
      end

