      subroutine p_mm_r(p_mm,it,ind_tau,ispin,key,key_g,key_p)
c$$$  calculate eq. 26 if key = 0

c     key = 0 : (G_tau)*(G_beta-tau) contribution to P -----------------
c     key = 1 : (G_tau_UP)*(G_beta-tau_DN) contribution to Xi_T -----
c     key_g = 0 - to use QP Green's function
c     key_g = 1 - to use regular Green's function

      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      include "mpif.h"
      integer, intent(in) :: ispin,it,ind_tau,key,key_g,key_p
      double precision, intent(out) :: p_mm(n_pbmtm,n_pbmtm,nrr_max) 
      integer :: isort,lm3,lm,l,ii,iii,indx,lget,icg1,icg2,icg,km3,lf3,
     &  jatom,jsort,km2,lm2,lf2,l1m1,l1,km1,lm1,lf1,iatom,
     &  km4,lm4,lf4,ie,lf,km,ir0,ind_ir0,isp,isp1
      double precision :: cgg,eps,g1,g2,dum
      complex*16 :: g1c,g2c,cggc,dumc
      double precision, allocatable ::
     $  x(:,:),y(:),y1(:),gr1(:,:,:),gr2(:,:,:)
      complex*16, allocatable :: xc(:,:),yc(:),gc1(:,:,:),gc2(:,:,:)


!      write(100,*)"cg",cg,"cgr",cgr,"fi0",fi0,"ind_prod",ind_prod,
!     $ "ind_prod0",ind_prod0,"indxcg",indxcg,"indxcgr",indxcgr,
!     $ "irel",irel,"is",is,"jcg",jcg,"jcgr",jcgr,"lf_isz",lf_isz,
!     $ "lfunm",lfunm,"lm_isz",lm_isz,"lmpb",lmpb,"maxel",maxel,
!     $ "n_mpi_nrr",n_mpi_nrr,"nrr_max",nrr_max,"nrr_red",nrr_red,
!     $ "nspin",nspin,"nspin_1",nspin_1,"ntle_pb",ntle_pb
!       call mpi_barrier(mpi_comm_world,mpierr)


      if(key==0) then
        isp=min(ispin,nspin_1)
        isp1=isp
      else if(key==1) then
        isp=1
        isp1=nspin
      endif

      if(irel<=1) then
        allocate(gr1(maxel,maxel,ndim_nrr(me3_k+1)))
        allocate(gr2(maxel,maxel,ndim_nrr(me3_k+1)))
        if(key==0) then
!          write(100,*)"before g_rs 1"
          call g_rs_from_ks_mm(gr1,dumc,1,ind_tau,ispin,key_g,key_p)
!          write(100,*)"before g_rs 2"
          call g_rs_from_ks_mm(gr2,dumc,2,ind_tau,ispin,key_g,key_p)
        else if(key==1) then
          call g_rs_from_ks_mm(gr1,dumc,it,ind_tau,1,key_g,key_p)
          call g_rs_from_ks_mm(gr2,dumc,3-it,ind_tau,nspin,key_g,key_p)
        endif		
      else if(irel==2) then
        allocate(gc1(maxel,maxel,ndim_nrr(me3_k+1)))
        allocate(gc2(maxel,maxel,ndim_nrr(me3_k+1)))
        if(key==0) then
          call g_rs_from_ks_mm(dum,gc1,1,ind_tau,ispin,key_g,key_p)
          call g_rs_from_ks_mm(dum,gc2,2,ind_tau,ispin,key_g,key_p)
        endif
      endif

      if(ind_tau==1) 
     &  call timel('** greens function P_R_MM finished *')

c      call log_eachcore('** greens function P_R_MM finished *')


      allocate(y1(n_pbmtm))
      if(irel<=1) allocate(x(maxel,n_pbmtm),y(n_pbmtm))
      if(irel==2) allocate(xc(maxel,n_pbmtm),yc(n_pbmtm))
      eps=1.d-15

      do ind_ir0=1,ndim_nrr(me3_k+1) ! R in eq. 26
        ir0=n_mpi_nrr(me3_k+1)+ind_ir0
        iatom=nrr_red(2,ir0)
        jatom=nrr_red(3,ir0)
        jsort=is(jatom)
        isort=is(iatom)
        do km3=1,lfunm(isort)   ! E3L3 index in eq 26
          lm3=lm_isz(km3,isort)
          lf3=lf_isz(km3,isort)
          if(irel<=1) x=0.d0
          if(irel==2) xc=(0.d0,0.d0) 
          do km4=1,lfunm(jsort) ! E4L4 index in eq 26
            if(irel<=1) g2=gr2(km3,km4,ind_ir0) 
            if(irel==2) then
              g2c=conjg(gc2(km3,km4,ind_ir0))
              g2=abs(g2c)
            endif
            if(abs(g2)>=eps) then
              lm4=lm_isz(km4,jsort)
              lf4=lf_isz(km4,jsort)
              do km2=1,lfunm(jsort) ! E2L2 index in eq 26
                lm2=lm_isz(km2,jsort)
                lf2=lf_isz(km2,jsort)
                if(irel<=1) then
                  ii = max0(lm4,lm2)
                  iii = min0(lm4,lm2)
                  indx = (ii*(ii-1))/2 + iii
                  icg1 = indxcg(indx)
                  icg2 = indxcg(indx+1) - 1
                else
                  icg1=indxcgr(lm2,lm4,1,0)
                  icg2=indxcgr(lm2,lm4,2,0)
                endif
c$$$                if (maswrk) then
c$$$                  if (km3 .eq. 1) then
c$$$                    write(993,'(a, 6i8)') 'icg1, icg2',
c$$$     $                ind_ir0,km4,km2,icg1,icg2, icg2-icg1
c$$$                  endif
c$$$                endif
                do icg = icg1, icg2
                  if(irel<=1) l1m1 = jcg(icg)
                  if(irel==2) l1m1 = jcgr(icg)
                  l1=lget(l1m1)
                  if(l1>lmpb(jsort)) cycle
                  if(irel<=1) then
                    cgg=g2*cg(icg)
                    do ie=1,ntle_pb(l1,jsort)
                      lf1=ind_prod0(ie,l1,jsort)
                      km1=ind_prod(ie,l1m1,jsort)
                      x(km2,km1)=x(km2,km1)+cgg
     &                  *fi0(lf2,lf4,lf1,jsort,isp,isp1) !radial part of <M_{Lk}^t\psi_{E2L2}^t|\psi_{E4L4}^t>
                    enddo
                  else if(irel==2) then
                    cggc=g2c*cgr(icg)
                    do ie=1,ntle_pb(l1,jsort)
                      lf1=ind_prod0(ie,l1,jsort)
                      km1=ind_prod(ie,l1m1,jsort)
                      xc(km2,km1)=xc(km2,km1)+cggc
     &                  *fi0(lf2,lf4,lf1,jsort,isp,isp)
                    enddo
                  endif
                enddo
              enddo             !! over km2
            endif
          enddo                 !! over km4
          do km1=1,lfunm(isort) 
            lm1=lm_isz(km1,isort)
            lf1=lf_isz(km1,isort)
            if(irel<=1) then
              y(1:n_pbmt0(jsort))=0.d0
              do km2=1,lfunm(jsort)
                g1=gr1(km1,km2,ind_ir0) ! G_{tE1L1, t'E2L2}^{R}(\tau),  value at -R is the same as value R
                if(abs(g1)>=eps) then
                  y(1:n_pbmt0(jsort))=y(1:n_pbmt0(jsort))
     &              +g1*x(km2,1:n_pbmt0(jsort))
                endif
              enddo             !! over km2
              ii = max0(lm3,lm1)
              iii = min0(lm3,lm1)
              indx = (ii*(ii-1))/2 + iii
              icg1 = indxcg(indx)
              icg2 = indxcg(indx+1) - 1
            else if(irel==2) then
              yc(1:n_pbmt0(jsort))=(0.d0,0.d0)
              do km2=1,lfunm(jsort)
                g1c=gc1(km1,km2,ind_ir0)
                if(abs(g1c)>=eps) then
                  yc(1:n_pbmt0(jsort))=yc(1:n_pbmt0(jsort))
     &              +g1c*xc(km2,1:n_pbmt0(jsort))
                endif
              enddo             !! over km2
              icg1=indxcgr(lm3,lm1,1,0)
              icg2=indxcgr(lm3,lm1,2,0)
            endif
c$$$            if (maswrk) then
c$$$              write(993,'(a, 6i8)') 'icg1, icg2',
c$$$     $          ind_ir0,km3,km1,icg1,icg2, icg2-icg1
c$$$            endif            
            do icg = icg1, icg2
              if(irel<=1) lm = jcg(icg)
              if(irel==2) lm = jcgr(icg)
              l=lget(lm)
              if(l>lmpb(isort)) cycle
              do ie=1,ntle_pb(l,isort)
                lf=ind_prod0(ie,l,isort)
                km=ind_prod(ie,lm,isort)
                if(irel<=1) then
                  cgg=fi0(lf1,lf3,lf,isort,isp,isp1)*cg(icg) !radial part of <M_{Lk}^t\psi_{E3L3}^t|\psi_{E1L1}^t>*angular part(cg)
                  y1(1:n_pbmt0(jsort))=y(1:n_pbmt0(jsort))*cgg
                else if(irel==2) then
                  cggc=fi0(lf1,lf3,lf,isort,isp,isp)*cgr(icg)
                  y1(1:n_pbmt0(jsort))=yc(1:n_pbmt0(jsort))*cggc
                endif
                do ii=1,n_pbmt0(jsort)
                  p_mm(km,ii,ir0)=p_mm(km,ii,ir0)-y1(ii)
                enddo
              enddo
            enddo               !! over icg
          enddo                 !! over km1
        enddo                   !! over km3
      enddo                     !! over ir0
      deallocate(y1)
      if(irel<=1) then
        deallocate(x,y,gr1,gr2)
      else if(irel==2) then
        deallocate(xc,yc,gc1,gc2)
      endif

      if(ind_tau==1) 
     &  call timel('**  matrix element P_R_MM finished *')
      
c      call log_eachcore('**  matrix element P_R_MM finished *')
      
      
      if(nproc_k/=1) call mpi_allreduce_dble
     $  (p_mm,n_pbmtm**2*nrr_max,mpi_sum,comm_tau)
      end
