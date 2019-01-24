! Here keep fixed format due to primme_f77.h
      module gprimme
      implicit none
      private
      public :: primme

      interface primme
          module procedure zprimme, dprimme
      end interface primme

      contains

      subroutine zprimme(basismax, numemax, maxmatvecs, etol,  
     :                   printlevel,whichevals, method, n, 
     :                   nest, evals, evecs, av)
      implicit none
      integer basismax, numemax, maxmatvecs, printlevel, whichevals, 
     :        method, n, nest
      real(8) etol, evals(numemax)
      complex(8) evecs(n*numemax)
      external :: av

      integer(8) prim
      integer ierr, i
      real(8) rnorms(numemax)
      include 'primme_f77.h'

      ! initialize prim
      call primme_initialize_f77(prim)
      call primme_set_member_f77(prim, primmef77_n, n)
      call primme_set_member_f77(prim, primmef77_numevals, numemax)
      call primme_set_member_f77(prim, primmef77_maxbasissize, 
     :                           basismax)
      call primme_set_member_f77(prim, primmef77_matrixmatvec, av)
! set the method to be used (after n, numevals, and precondition have
! been set. also after basissize is set if desired.)
      call primme_set_method_f77(prim, method, ierr)
      call primme_set_member_f77(prim, primmef77_eps, etol)
      call primme_set_member_f77(prim, primmef77_target, whichevals)
      call primme_set_member_f77(prim, primmef77_printlevel,
     :                           printlevel)
      call primme_set_member_f77(prim, primmef77_maxmatvecs, 
     :                           maxmatvecs)
      call primme_set_member_f77(prim, primmef77_initsize, nest)
! calling the primme solver
      call zprimme_f77(evals, evecs, rnorms, prim, ierr)
      if (ierr /= 0) then
          write(*,'(" zprimme returned with error: ",I0)') ierr
      endif
      call primme_display_stats_f77(prim)
      call primme_free_f77(prim)
      write(*, '(" zprimme:")')
      write(*,'(" dim_H = ",I0)')n
      do i = 1, numemax
          write (*, 9000) i, evals(i),rnorms(i)
      enddo
      return
9000  format (1x,'e(',i1,') = ',f0.16,4x,'residual norm = ', e12.4)
      end subroutine zprimme


      subroutine dprimme(basismax, numemax, maxmatvecs, etol,  
     :                   printlevel,whichevals, method, n, 
     :                   nest, evals, evecs, av)
      implicit none
      integer basismax, numemax, maxmatvecs, printlevel, whichevals, 
     :        method, n, nest
      real(8) etol, evals(numemax)
      real(8) evecs(n*numemax)
      external :: av

      integer(8) prim
      integer ierr, i
      real(8) rnorms(numemax)
      include 'primme_f77.h'

      ! initialize prim
      call primme_initialize_f77(prim)
      call primme_set_member_f77(prim, primmef77_n, n)
      call primme_set_member_f77(prim, primmef77_numevals, numemax)
      call primme_set_member_f77(prim, primmef77_maxbasissize, 
     :                           basismax)
      call primme_set_member_f77(prim, primmef77_matrixmatvec, av)
! set the method to be used (after n, numevals, and precondition have
! been set. also after basissize is set if desired.)
      call primme_set_method_f77(prim, method, ierr)
      call primme_set_member_f77(prim, primmef77_eps, etol)
      call primme_set_member_f77(prim, primmef77_target, whichevals)
      call primme_set_member_f77(prim, primmef77_printlevel,
     :                           printlevel)
      call primme_set_member_f77(prim, primmef77_maxmatvecs, 
     :                           maxmatvecs)
      call primme_set_member_f77(prim, primmef77_initsize, nest)
! calling the primme solver
      call dprimme_f77(evals, evecs, rnorms, prim, ierr)
      if (ierr /= 0) then
          write(*,'(" dprimme returned with error: ",I0)') ierr
      endif
      call primme_display_stats_f77(prim)
      call primme_free_f77(prim)
      write(*, '(" dprimme:")')
      write(*,'(" dim_H = ",I0)')n
      do i = 1, numemax
          write (*, 9000) i, evals(i),rnorms(i)
      enddo
      return
9000  format (1x,'e(',i1,') = ',f0.16,4x,'residual norm = ', e12.4)
      end subroutine dprimme


      end module gprimme
