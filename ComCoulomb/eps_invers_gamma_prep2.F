      subroutine eps_invers_gamma_prep2
c     Deals with (1-VP)^(-1)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      include "mpif.h"
      integer :: i_len,i_nu,j,i,nd,ii,npbt1,ind_nu,jatom,jsort,jnd0,
     &  iatom,isort,ind0,lm1,jj,lm
      double precision :: a0,an,ai,aa,bb,f,fc
      complex*16 :: cc
      double precision, allocatable :: e0(:),e2(:)
      allocate(e0(0:n_nu),e2(0:n_nu))
      e0=0.d0
      npbt1=n_pbmt+nplwgw(1)
      do ind_nu=1,ndim3_nu
        i_nu=me3_tau*ndim3_nu+ind_nu-1
        call pw_from_pb(1,npbt1,npbt1,eps_invers_gamma(1,1,ind_nu),cc,2)
        e0(i_nu)=dble(cc)
      enddo
c$$$      if(nproc_tau/=1) call mpi_allreduce_(e0,n_nu+1,mpi_sum,comm_k)
      if(nproc_tau/=1) call mpi_allreduce_dble(e0,n_nu+1,mpi_sum,comm_k)      
c     --------------------------------------------------------------------
      do i_nu=1,n_nu
        if(e0(i_nu)<e0(i_nu-1)) then
          e0(i_nu)=e0(i_nu-1)
        endif
      enddo
c     --------------------------------------------------------
      a0=e0(0)
      if(a0<0.d0) then
        an=e0(n_nu)
        do i_nu=0,n_nu
          ai=e0(i_nu)
          e0(i_nu)=an-(ai-an)*an/(a0-an)
        enddo
      endif
c     --------------------------------------------------------
      an=e0(n_nu)
      if(abs(an-1.d0)>1.d-10) then
        a0=e0(0)
        fc=(1.d0-a0)/(an-a0)
        do i_nu=0,n_nu
          ai=e0(i_nu)
          e0(i_nu)=a0+fc*(ai-a0)
        enddo
      endif
c     ----- Plasmon pole smoothing -----------------------------------
      e2=e0(:)
      if(metal) e2(0)=0.d0
      f=e2(0)+0.4d0*(e2(n_nu)-e2(0))
      nd=0
      a0=abs(f-e2(0))
      do i_nu=1,n_nu
        an=abs(f-e2(i_nu))
        if(an<a0) then
          a0=an
          nd=i_nu
        endif
      enddo
      a0=e2(0)/(e2(0)-1.d0)
      an=e2(nd)/(e2(nd)-1.d0)
      aa=(an-a0)/w_nu(nd)**2
      do i_nu=1,n_nu
        bb=w_nu(i_nu)**2*aa+a0
        bb=1.d0/bb
        bb=1.d0-bb
        e2(i_nu)=1.d0/bb
        em1_gg0(1,i_nu)=e2(i_nu)
      enddo
c     -----------------------------------------------------
      i_len=len_trim(allfile)
      if(maswrk) then
        open(3,file=allfile(1:i_len)//'_em1_gg0_fit_nu')
        do i_nu=0,n_nu
          write(3,'(e12.5,2(1x,f18.10))')w_nu(i_nu),e0(i_nu),e2(i_nu)
        enddo
        close(3)
      endif
c     --- Now we correct matrix elements --------------------------------
      do ind_nu=1,ndim3_nu
        i_nu=me3_tau*ndim3_nu+ind_nu-1
        if(i_nu==0.and.metal) then
          a0=0.d0
        else
          a0=e2(i_nu)/e0(i_nu)
        endif
        do jatom=1,natom
          jsort=is(jatom)
          jnd0=iopb(jatom)-1
          do j=1,n_pbmt0(jsort)
            lm1=lm_pbmt(j,jsort)
            if(lm1>1) cycle
            jj=jnd0+j
            do iatom=1,natom
              isort=is(iatom)
              ind0=iopb(iatom)-1
              do i=1,n_pbmt0(isort)
                lm=lm_pbmt(i,isort)
                if(lm>1) cycle
                ii=ind0+i
                eps_invers_gamma(ii,jj,ind_nu)=
     &	          a0*eps_invers_gamma(ii,jj,ind_nu)
              enddo
            enddo
c     -- All elements in interstitial are affected ----------------
            do ii=n_pbmt+1,npbt1
              eps_invers_gamma(ii,jj,ind_nu)=
     &          a0*eps_invers_gamma(ii,jj,ind_nu)
            enddo
          enddo
        enddo
        do jj=n_pbmt+1,npbt1
          do iatom=1,natom
            isort=is(iatom)
            ind0=iopb(iatom)-1
            do i=1,n_pbmt0(isort)
              lm=lm_pbmt(i,isort)
              if(lm>1) cycle
              ii=ind0+i
              eps_invers_gamma(ii,jj,ind_nu)=
     &          a0*eps_invers_gamma(ii,jj,ind_nu)
            enddo
          enddo
c     -- All elements in interstitial are affected ----------------
          do ii=n_pbmt+1,npbt1
            eps_invers_gamma(ii,jj,ind_nu)=
     &        a0*eps_invers_gamma(ii,jj,ind_nu)
          enddo
        enddo
      enddo	
      deallocate(e0,e2)
      end
