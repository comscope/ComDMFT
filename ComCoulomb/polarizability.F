      subroutine polarizability(key_u)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use vertex_mod
      use wannier90_mod
      implicit none
      include "mpif.h"
      integer, intent(in) :: key_u
      integer :: ik,ind_tau,ind_k,k,ii,kst,iq,ind
      double precision :: dev0
      double precision, allocatable :: p0(:,:,:)
      complex*16,allocatable :: tmp3(:,:),tmp4(:,:),tmp2(:,:)

!      write(100,*)"iun",iun,"nvopt_max",nvopt_max,"obietto",
!     $ obietto,"p_local_k",p_local_k,"p_w_solid_q_tau",
!     $ p_w_solid_q_tau,"p_wannier",p_wannier,"pol_iter_for_u",
!     $ pol_iter_for_u,"u_mode",u_mode,"ubi",ubi,"we_need_u",
!     $ we_need_u,"weig_rc",weig_rc
!      call mpi_barrier(mpi_comm_world,mpierr)


      if(obietto=='   solid') then
        if (pol_iter_for_u .gt. 0) then
          if (allocated(p_wannier)) deallocate(p_wannier)
          allocate(p_wannier
     $      (nvopt_max,nvopt_max,ndim3_tau,ndim3_k(me3_k+1)))
          p_wannier=0.0d0                      
          if (we_need_u .and.
     $      ((u_mode .eq. 'enpj') 
     $      .or. (u_mode .eq. 'wnds')
     $      .or. (u_mode .eq. 'bnse'))) then
            if (maswrk) write(iun,*)
     $        "polarizability from the correlated bands"
            
            do ii=1, pol_iter_for_u
              call p_r_tau_solid(ii)
              p_wannier=p_wannier+p_w_solid_q_tau
            enddo
          elseif (we_need_u .and.
     $        (u_mode .eq. 'gtop')) then
            call chi0_from_g_correlated(p_wannier)
          endif
C$$$  call p_r_0_cal(1)

          if (p_local_k) then                  
            do ind_k=1, ndim3_k(me3_k+1)
              k=n3_mpi_k(me3_k+1)+ind_k
              p_wannier(:,:,:,ind_k)
     $          =p_wannier(:,:,:,ind_k)*weig_rc(k)
            enddo
          endif          
        endif

c     ---------------------------------------------------------------
c$$$  if (myrevision .eq. 1) then      

c$$$  endif
c$$$  endif

        if (maswrk) write(iun,*)
     $    "polarizability from full bands"
        call p_r_tau_solid(0)
        if(ubi==' gw'.or.ubi=='gwg')
     $    call etot_x_from_pw(2)
      endif
      call timel('**** POLARIZABILITY finished *******')
      end
