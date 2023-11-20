      subroutine check_lda_dmft

      use manager_mod
      use lda_dmft_mod
      use parallel_mod      
      use files_mod
#ifdef USE_HDF5
#ifndef MPI
The code expects to be using MPI if you want to use HDF5
#endif
#endif
#ifdef USE_HDF5
      use hdf5io_mod
#endif

      implicit none
#ifdef MPI
      INCLUDE 'mpif.h'            
      integer :: mpierr
#endif
#ifdef USE_HDF5
      integer(hid_t) :: fid_f =-1
      integer :: h5err
#endif

      logical :: file_exist

      
      if (maswrk) then
        inquire(file="wannier_den_matrix.dat", exist=file_exist)
        if (file_exist) then
          is_lda_dmft=.true.
        else
          is_lda_dmft=.false.
        endif
      endif
#ifdef MPI
      call brdcst(msgint,is_lda_dmft,4,master,MPI_COMM_WORLD)
#endif
      if (is_lda_dmft) then
        nproc_k=nproc
        nproc_tau=1
        iter_dft=2
        iter_hf=0
        iter_qp=0
        iter_gw=0
        restart=.True.
        bandplot=.False.
        dosplot=.False.
#ifdef USE_HDF5
        call fopen(trim(allfile)//'.rst',fid_rst)
        call fwrite(MPI_COMM_WORLD,'dft',"ubi",fid_rst)
        call fwrite(MPI_COMM_WORLD,1,"iteration",fid_rst)
        call fclose(fid_rst)
#else
        if (maswrk) then
          open(3,file=adjustl(trim(checkpoint_dir))//
     $      trim(allfile)//'.rst',form="unformatted")
          write(3)'dft', 1
          close(3)
        endif
#endif

#ifdef USE_HDF5
        call fopen('wannier_den_matrix.dat',fid_f)
        call fread(MPI_COMM_WORLD,dmft_temperature,"temperature",fid_f)
        call fclose(fid_f)
#else
        if (maswrk) then
          open(unit=305, file='wannier_den_matrix.dat',status='old',
     $      form='unformatted')
          read(305) dmft_temperature
          close(305)
        endif
#ifdef MPI
        call mpi_bcast(dmft_temperature,1,mpi_double_precision,
     $    0,MPI_COMM_WORLD,mpierr)
#endif      
#endif  
c$$$        temperature=dmft_temperature
      endif
      end

! initialize_wan_density_matrix reads in wannier_den_matrix.dat
!     wannier_den_matrix.dat can be in either hdf5 or binary format. It contains
!       - dmft_temperature="temperature"
!       - nqdiv_temp = "nqdiv" (ignored) - the number of k-points
!       - n_dmft_bnd="n_dmft_bnd" from nbnd_in_win_co [called 'num_band_in' when saved from comrisb, and called 'n_dmft_bnd' when saved from comlowh] - should be the number of bands
!       - n_dmft_wan="n_dmft_wan" from nbnd_in_win_co - the number of wanniers
!       - dmft_bnd_index="include_band_new"
!       - disentangle_mat="vw_matrix_new"
!       - dmft_den_mat='n_matrix_new' - the density matrix
      subroutine initialize_wan_density_matrix

      use lda_dmft_mod
      use parallel_mod
      use manager_mod
      use atom_mod            
      use solid_mod      
#ifdef USE_HDF5
      use hdf5io_mod
      use hdf5
#endif

      implicit none
#ifdef MPI
      INCLUDE 'mpif.h'
#endif
      integer :: ii,ispin,k,k0,i,mpierr,nqdiv_temp
#ifdef USE_HDF5
      integer :: ng(7)
      integer :: nv(7)
      integer(hid_t) :: fid_f =-1
      integer :: h5err
#endif

#ifdef USE_HDF5

      call fopen('wannier_den_matrix.dat',fid_f)
      call fread(MPI_COMM_WORLD,dmft_temperature,"temperature",fid_f)
!      if (maswrk)  write(iun,*)'dmft_temperature',dmft_temperature
      call fread(MPI_COMM_WORLD,nqdiv_temp,"nqdiv",fid_f)
!      if (maswrk)  write(iun,*)'nqdiv_temp',nqdiv_temp
        allocate(n_dmft_bnd(nqdiv))
        allocate(n_dmft_wan(nqdiv))
      ng(1)=size(n_dmft_bnd,1)
      call fread(MPI_COMM_WORLD,n_dmft_bnd,ng(1),nv(1),
     $                "n_dmft_bnd",fid_f)
!      if (maswrk)  write(iun,*)'n_dmft_bnd',n_dmft_bnd
        n_dmft_bnd_max=maxval(n_dmft_bnd)

      ng(1)=size(n_dmft_wan,1)
      call fread(MPI_COMM_WORLD,n_dmft_wan,ng(1),nv(1),
     $                "n_dmft_wan",fid_f)
!      if (maswrk)  write(iun,*)'n_dmft_bnd',n_dmft_bnd
        n_dmft_wan_max=maxval(n_dmft_wan)

        allocate(dmft_bnd_index(n_dmft_bnd_max,nqdiv))
      ng(1)=size(dmft_bnd_index,1)
      ng(2)=size(dmft_bnd_index,2)
      call fread(MPI_COMM_WORLD,dmft_bnd_index,ng,nv,
     $                "include_band_new",fid_f)
!      if (maswrk)  write(iun,*)'dmft_bnd_index',dmft_bnd_index
        allocate(disentangle_mat(n_dmft_bnd_max,n_dmft_wan_max,nqdiv))
        disentangle_mat=0.0d0
        allocate(dmft_den_mat(n_dmft_wan_max,n_dmft_wan_max,nqdiv))
        dmft_den_mat=0.0d0
      ng(1)=size(disentangle_mat,1)
      ng(2)=size(disentangle_mat,2)
      ng(3)=size(disentangle_mat,3)
      call fread(MPI_COMM_WORLD,disentangle_mat,ng,nv,
     $  "vw_matrix_new",fid_f)
!      if (maswrk) write(iun,*)'disentangle_mat',disentangle_mat
!      if (maswrk) write(iun,*)'size',ng(1),ng(2),ng(3)

      ng(1)=size(dmft_den_mat,1)
      ng(2)=size(dmft_den_mat,2)
      ng(3)=size(dmft_den_mat,3)

      call fread(MPI_COMM_WORLD,dmft_den_mat,ng,nv,
     $ "n_matrix_new",fid_f)
!      if (maswrk) write(iun,*)'dmft_den_mat',dmft_den_mat
!      if (maswrk) write(iun,*)'size',ng(1),ng(2),ng(3)

      call fclose(fid_f)
      if (maswrk) write(iun,*) 'dmft_temperature', temperature
#else

      if (maswrk) then
        open(unit=305, file='wannier_den_matrix.dat',status='old',
     $    form='unformatted')
        read(305) dmft_temperature!
        read(305) nqdiv_temp
        allocate(n_dmft_bnd(nqdiv))
        allocate(n_dmft_wan(nqdiv))        
        read(305) n_dmft_bnd
        n_dmft_bnd_max=maxval(n_dmft_bnd)        !
        read(305) n_dmft_wan
        n_dmft_wan_max=maxval(n_dmft_wan)         !       
        allocate(dmft_bnd_index(n_dmft_bnd_max,nqdiv))
        read(305) dmft_bnd_index                

        allocate(disentangle_mat(n_dmft_bnd_max,n_dmft_wan_max,nqdiv))
        disentangle_mat=0.0d0        
        allocate(dmft_den_mat(n_dmft_wan_max,n_dmft_wan_max,nqdiv))
        dmft_den_mat=0.0d0
        read(305) disentangle_mat
        read(305) dmft_den_mat
        close(305)
c$$$        write(iun,*) 'dmft_temperature', temperature        
      endif

#ifdef MPI      
      call mpi_bcast(dmft_temperature,1,mpi_double_precision,
     $  0,MPI_COMM_WORLD,mpierr)      
      call mpi_bcast(n_dmft_bnd_max,1,mpi_integer,
     $  0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(n_dmft_wan_max,1,mpi_integer,
     $  0,MPI_COMM_WORLD,mpierr)
#endif
      
      if (.not. maswrk) then
        allocate(n_dmft_bnd(nqdiv))
        n_dmft_bnd=0
        allocate(n_dmft_wan(nqdiv))
        n_dmft_wan=0
        allocate(dmft_bnd_index(n_dmft_bnd_max,nqdiv))
        dmft_bnd_index=0
        allocate(dmft_den_mat(n_dmft_wan_max,n_dmft_wan_max,nqdiv))
        dmft_den_mat=0.0d0
        allocate(disentangle_mat(n_dmft_bnd_max,n_dmft_wan_max,nqdiv))
        disentangle_mat=0.0d0        
      endif

#ifdef MPI
      call mpi_bcast(n_dmft_bnd,nqdiv,mpi_integer,
     $  0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(n_dmft_wan,nqdiv,mpi_integer,
     $  0,MPI_COMM_WORLD,mpierr)

      call mpi_bcast(dmft_bnd_index,n_dmft_bnd_max*nqdiv,mpi_integer,
     $  0,MPI_COMM_WORLD,mpierr)
      
      
      call mpi_bcast(dmft_den_mat,n_dmft_wan_max**2*nqdiv,
     $  mpi_double_complex,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(disentangle_mat,
     $  n_dmft_bnd_max*n_dmft_wan_max*nqdiv,
     $  mpi_double_complex,0,MPI_COMM_WORLD,mpierr)
#endif
#endif

c$$$      write(179,*) 'band_index', dmft_bnd_index      
      
c$$$  allocate(is_dmft_bnd_included(nbndf))
c$$$  is_dmft_bnd_included=0
c$$$  do ii=1, n_dmft_bnd
c$$$  is_dmft_bnd_included(dmft_bnd_index(ii))=1
c$$$  enddo

      allocate(ev_bnd_dmft_0
     $  (nbasmpw,nbndf,ndim3_k(me3_k+1),nspin_0))
      ev_bnd_dmft_0=0.0d0
      allocate(z_bnd_dmft_0
     $  (nfun,nbndf,ndim3_k(me3_k+1),nspin_0))
      z_bnd_dmft_0=0.0d0

      ev_bnd_dmft_0=ev_bnd
      z_bnd_dmft_0=z_bnd

      allocate(gfun_dmft_0(maxwf,nspin))
      gfun_dmft_0=0.0d0
      allocate(gfund_dmft_0(maxwf,nspin))      
      gfund_dmft_0=0.0d0

      gfun_dmft_0=gfun
      gfund_dmft_0=gfund      

      allocate(g_full_0(nbndf,nbndf,ndim3_k(me3_k+1),nspin_0))
      g_full_0=0.0d0
          
      end

      
      subroutine dft_basis_rotmat_k(ind_k,k0,ispin,umat)
      use lda_dmft_mod
      use solid_mod
      use atom_mod      
      implicit none

      integer,intent(in) :: ind_k, k0, ispin
      complex*16, intent(out) :: umat(nbndf,nbndf) !<new|0>
      real*8 :: fifi_j(maxlfun,maxlfun,0:2*maxb,nsort)      
      complex*16 :: ff(maxel,maxel,natom)

      umat=0.0d0
      fifi_j=0.0d0
      ff=0.0d0
      
      call fifi_j_prepare_onepoint_2(fifi_j,(/0.0d0,0.0d0,0.0d0/),ispin)
      
      call integral_band_pw_band_mt_wan(fifi_j,ff,(/0.0d0,0.0d0,0.0d0/))

      call integral_band_pw_band_full
     $  (ff,z_bnd(1,1,ind_k,ispin),ev_bnd(1,1,ind_k,ispin),
     $  z_bnd_dmft_0(1,1,ind_k,ispin),
     $  ev_bnd_dmft_0(1,1,ind_k,ispin),
     &  umat,k0,(/0.0d0,0.0d0,0.0d0/),k0,
     $  (/0,0,0/),n_bnd(k0,ispin),n_bnd(k0,ispin))
      end

      subroutine integral_band_pw_band_mt_wan(fifi_j,ff,pn)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      real*8, intent(in) ::
     $  fifi_j(maxlfun,maxlfun,0:2*maxb,nsort),pn(3)
      complex*16, intent(out) :: ff(maxel,maxel,natom)

      integer :: iatom,isort,nl,ind_lm,ndimb,ndimpb,lm1,lm2,l1,l2,ii,
     &  iii,indx,icg1,icg2,i2,mj1,mj2,i1,li1,li2,icg,lmi,lget,
     &  l_i,nm1,nm2,iei,lfi,ie1,km1,lf1,ie2,km2,lf2,ispin,in1,
     &  jn2,in2,jn1
      real*8 :: fas,q(3),qq,sqpi4,pi2
      complex*16 :: c1,c2
      real*8, allocatable :: ylm(:)

      allocate(ylm((2*maxb+1)**2))
      pi2=pi+pi
      sqpi4=sqrt(4.d0*pi)
      q=pn
      qq=sqrt(dot_product(q,q))

      if(qq.gt.1.d-9) then
        call sphharm(q(1)/qq,q(2)/qq,q(3)/qq,2*maxb,ylm)
      else
        ylm=0.d0
        ylm(1)=1.d0/sqpi4
      endif
      ff=(0.d0,0.d0)
      
      do iatom=1,natom
        isort=is(iatom)
        nl=lfunm(isort)
        ndimb=nrel*(lmb(isort)+1)**2
        fas=pi2*dot_product(pn,tau(:,iatom))
        c1=exp(dcmplx(0.d0,-fas))
        do lm1=1,ndimb
          do lm2=1,ndimb
            if(irel<=1) then
              l1=lget(lm1)
              l2=lget(lm2)
              ii = max0(lm1,lm2)
              iii = min0(lm1,lm2)
              indx = (ii*(ii-1))/2 + iii
              icg1 = indxcg(indx)
              icg2 = indxcg(indx+1) - 1 
            else
              call getlimj(lm2,l2,i2,mj2,li2,0)
              call getlimj(lm1,l1,i1,mj1,li1,0)
              icg1=indxcgr(lm1,lm2,1,0)
              icg2=indxcgr(lm1,lm2,2,0)
            endif
            do icg = icg1, icg2
              if(irel<=1) then
                lmi = jcg(icg)
                l_i=lget(lmi)
                c2=c1*dcmplx(0.d0,-1.d0)**l_i*ylm(lmi)*cg(icg)
              else
                lmi=jcgr(icg)
                l_i=lget(lmi)
                c2=c1*dcmplx(0.d0,-1.d0)**l_i*ylm(lmi)*cgr(icg)
              endif
              nm1=ntle(l1,isort)
              nm2=ntle(l2,isort)
              do ie1=1,nm1
                in1=1
                if(augm(ie1,l1,isort)/='LOC') in1=2
                do jn1=1,in1
                  km1=indbasa(jn1,ie1,lm1,isort)
                  lf1=lf_isz(km1,isort)
                  do ie2=1,nm2
                    in2=1
                    if(augm(ie2,l2,isort)/='LOC') in2=2
                    do jn2=1,in2
                      km2=indbasa(jn2,ie2,lm2,isort)
                      lf2=lf_isz(km2,isort)
                      ff(km1,km2,iatom)=
     &                  ff(km1,km2,iatom)+c2
     &                  *fifi_j(lf1,lf2,l_i,isort)
                    enddo
                  enddo         ! over ie2
                enddo
              enddo
            enddo               ! over icg
          enddo                 ! over lm2
        enddo                   ! over lm1
      enddo                     !!! over iatom
      deallocate(ylm)
      end      


      subroutine integral_band_pw_band_full(ff,znew,anew,znew1,anew1,
     $  b_pw_b,k,pn,kq,kq_shift,nbndk,nbndkq)

!     znew and anew at k
!     znew1 and anew1 at kq      
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod

      implicit none
      integer, intent(in) :: k,kq,kq_shift(3),nbndk,nbndkq
      real*8, intent(in) :: pn(3)
      complex*16, intent(in) ::
     $  znew(nfun,nbndf),
     &  anew(nbasmpw,nbndf),
     &  znew1(nfun,nbndf),
     &  anew1(nbasmpw,nbndf),
     &  ff(maxel,maxel,natom)
      complex*16, intent(out) :: b_pw_b(nbndf,nbndf)
      integer :: j,iatom,isort,k0,kq0,ndimb,ind0,ispin,jspin,
     &  igb1,igb2,igb3,ig_tild,jbas,jnd,ia,ib,ic,ibas,ind,ia1,
     &  ib1,ic1,ibas1,ibas0,irl
      real*8 :: v(3),det,det3
      complex*16, allocatable :: s(:,:),tmp(:,:),tmp1(:,:)



      b_pw_b=(0.d0,0.d0)

      k0=i_kref(k)
      kq0=i_kref(kq)

      allocate(s(nbndtm,nbndtm))
      allocate(tmp(nbndtm,nbndf))
      allocate(tmp1(maxel,nbndf))
      s=0.0d0      
c     ----------------- MT contribution ----------------------------------
      do iatom=1,natom
        isort=is(iatom)
        ndimb=lfunm(isort)
        ind0=io_lem(iatom)
        tmp1=0.0d0
        call zgemm('n','n',ndimb,nbndkq,ndimb,
     &    (1.d0,0.d0),ff(1,1,iatom),maxel,
     &    znew1(ind0,1),nfun,(0.d0,0.d0),tmp1,maxel)
        call zgemm('c','n',nbndk,nbndkq,ndimb,
     &    (1.d0,0.d0),znew(ind0,1),nfun,
     &    tmp1,maxel,(1.d0,0.d0),b_pw_b,
     &    nbndf)
      enddo                     !! over iatom
c     ----------------- Interstitial contribution ------------------------
      do jbas=1,nbask(k0)/nrel
        jnd=indgb(jbas,k)
        ia=kq_shift(1)+igbs(1,jnd)
        ib=kq_shift(2)+igbs(2,jnd)
        ic=kq_shift(3)+igbs(3,jnd)
        do ibas=1,nbask(kq0)/nrel
          ind=indgb(ibas,kq)
          ia1=-ia+igbs(1,ind)
          ib1=-ib+igbs(2,ind)
          ic1=-ic+igbs(3,ind)
          j=indplw(ia1,ib1,ic1)
          if(complex_ro) then
            s(jbas,ibas)=dcmplx(sovr(j),sovi(j))
          else
            s(jbas,ibas)=dcmplx(sovr(j),0.d0)
          endif
        enddo
      enddo
      do irl=1,nrel
        ibas1=nbask(kq0)/nrel*(irl-1)
        ibas0=nbask(k0)/nrel*(irl-1)
        tmp=0.0d0
        call zgemm('n','n',nbask(k0)/nrel,nbndkq,nbask(kq0)/nrel,
     &    (1.d0,0.d0),s,nbndtm,anew1(ibas1+1,1),nbasmpw,
     &    (0.d0,0.d0),tmp,nbndtm)
        call zgemm('c','n',nbndk,nbndkq,
     &    nbask(k0)/nrel,(1.d0,0.d0),anew(ibas0+1,1),nbasmpw,
     &    tmp,nbndtm,(1.d0,0.d0),b_pw_b,
     &    nbndf)
      enddo
      deallocate(s,tmp,tmp1)
      end      
      

      subroutine fifi_j_prepare_onepoint_2(fifi_j,pn,ispin)

      use atom_mod
      use manager_mod
      use solid_mod
      use units_mod
      use lda_dmft_mod      
      implicit none
      integer,intent(in) :: ispin
      real*8, intent(in) :: pn(3)
      real*8, intent(out) ::
     $  fifi_j(maxlfun,maxlfun,0:2*maxb,
     &  nsort)
      integer :: isort,ir,l,mt,mt1,iq,lf,lf1
      real*8 :: pi4,qq,rv,dqdall,q(3),c2,sqpi4,pp,pi2a,
     $  bess(0:2*maxb,0:maxnrad),bessd(0:2*maxb),
     $  work(0:maxnrad),work1(0:maxnrad)

      fifi_j=0.0d0

      pi4=4.d0*pi
      pi2a=(pi+pi)/par
      sqpi4=sqrt(pi4)
      c2=clight**2

      q=pn
      qq=sqrt(dot_product(q,q))*pi2a ! bvec
      
      do isort=1,nsort
        do ir=0,nrad(isort)
          rv=r(ir,isort)
          call BESSR(qq,rv,bess(0,ir),bessd,2*lmb(isort))
          rv=rv*rv
          work(ir)=rv*dr(ir,isort)
        enddo
        do l=0,2*lmb(isort)
          do lf1=1,lfun(isort)
            mt1=ind_wf(lf1,isort)
            do lf=1,lfun(isort)
              mt=ind_wf(lf,isort)
              do ir=0,nrad(isort)
                pp=gfun(mt+ir,ispin)*gfun_dmft_0(mt1+ir,ispin)
                if(irel>=1) pp=pp+gfund(mt+ir,ispin)
     &            *gfund_dmft_0(mt1+ir,ispin)/c2
                work1(ir)=work(ir)*bess(l,ir)*pp
              enddo
              fifi_j(lf,lf1,l,isort)=
     &          pi4*dqdall(h(isort),work1,nrad(isort))
            enddo               !! over lf
          enddo                 !! over lf1
        enddo                   !! over l
      enddo                     !! over iatom
      end      



      subroutine check_if_unitary(ndim,ndim0, mat, fail_ind, maxerror)
      implicit none
      integer, intent(in) :: ndim,ndim0
      complex*16, intent(in) :: mat(ndim0,ndim0)
      integer, intent(out) :: fail_ind
      real*8,intent(out) :: maxerror

      integer :: ii,jj
      complex*16 :: mat2(ndim,ndim)
      

      mat2=0.0d0
      call zgemm('n','c',ndim,ndim,ndim,
     $  (1.0d0,0.0d0),mat,ndim0,mat,ndim0,
     $  (0.0d0,0.0d0),mat2,ndim)

      maxerror=0.0d0
      fail_ind=0
      do ii=1, ndim
        do jj=1, ndim
          if (ii .eq. jj) then
            if (cdabs(mat2(ii,jj)-1.0d0) .gt. 1.0d-6) then
              fail_ind=fail_ind+1
              write(179,*) ii,jj,mat2(ii,jj)
            endif
            if (maxerror .lt. cdabs(mat2(ii,jj)-1.0d0)) then
              maxerror=cdabs(mat2(ii,jj)-1.0d0)
            endif
          else
            if (cdabs(mat2(ii,jj)) .gt. 1.0d-6) then
              fail_ind=fail_ind+1
              write(179,*) ii,jj,mat2(ii,jj)                            
            endif
            if (maxerror .lt. cdabs(mat2(ii,jj))) then
              maxerror=cdabs(mat2(ii,jj))
            endif              
          endif
        enddo
      enddo
      end            


      subroutine check_hermitian(mata,matsize1,matsize0,string)

!     assumes that matsize3 < lda and ldb
!     multiplication of dconjg(transpose(a(1:matsize3,1:matsize1))) and b(1:matsize3,1:matsize2)
      
      implicit none
      
      integer,intent(in) :: matsize1,matsize0
      complex*16,intent(in) :: mata(matsize0,matsize0)
      character(len=*), intent(in) :: string
      integer :: ii, jj
      real*8 :: err

      err=1.0d-6
      
      do ii=1, matsize1
        do jj=ii, matsize1        
          if (cdabs(mata(ii,jj)-dconjg(mata(jj,ii))) .gt. err) then
            write(179,'(a, 2i5,2f12.6, 2x, 2f12.6)')
     $        'hermitian fail: '//trim(string),
     $        ii, jj, mata(ii,jj),mata(jj,ii)
            call flush(179)
            write(993,'(a, 2i5,2f12.6, 2x, 2f12.6)')
     $        'hermitian fail: '//trim(string),
     $        ii, jj, mata(ii,jj),mata(jj,ii)
            call flush(993)            
            call ending
          endif
        enddo
      enddo
      end      



      subroutine force_mat_unitary(dim1,dim0,umat)
      implicit none
      integer, intent(in) :: dim1, dim0
      complex*16, intent(inout) :: umat(dim0,dim0)

      integer :: ii,jj,info,lwork
      real*8 :: ss(dim1),rwork(5*dim1)
      complex*16 :: umat_temp(dim1,dim1),uu(dim1,dim1),vt(dim1,dim1),
     $  worktemp(3*dim1)
      complex*16, allocatable :: work(:)

      umat_temp=0.0d0

      do ii=1, dim1
        do jj=1, dim1
          umat_temp(ii,jj)=umat(ii,jj)
        enddo
      enddo

      call zgesvd('A','A',dim1,dim1,umat_temp,dim1,
     $  ss,uu,dim1,vt,dim1,worktemp,-1,rwork,info)
      lwork=max(int(worktemp(1)), 1)
      allocate(work(lwork))
      call zgesvd('A','A',dim1,dim1,umat_temp,dim1,
     $  ss,uu,dim1,vt,dim1,work,lwork,rwork,info)
      deallocate(work)

      umat=0.0d0
      call zgemm('n','n',dim1,dim1,dim1,
     &  (1.d0,0.d0),uu,dim1,vt,dim1,
     &  (0.d0,0.d0),umat,dim0)

      end



c$$$      subroutine print_xsf
c$$$      
c$$$      use parallel_mod
c$$$      use solid_mod
c$$$      use units_mod
c$$$      use atom_mod
c$$$      
c$$$      implicit none
c$$$
c$$$      integer :: ii,jj
c$$$
c$$$      if (maswrk) then
c$$$
c$$$	 open(unit=992, file='coord.xsf')
c$$$	 write(992,'(a)') "CRYSTAL"
c$$$	 write(992,'(a)') "PRIMVEC"
c$$$	 write(992,'(3(f12.7,3x))') (rbas(ii,1)*par*bohr, ii=1,3)
c$$$	 write(992,'(3(f12.7,3x))') (rbas(ii,2)*par*bohr, ii=1,3)
c$$$	 write(992,'(3(f12.7,3x))') (rbas(ii,3)*par*bohr, ii=1,3)
c$$$	 write(992,'(a)') "CONVVEC"
c$$$	 write(992,'(3(f12.7,3x))') (rbas(ii,1)*par*bohr, ii=1,3)
c$$$	 write(992,'(3(f12.7,3x))') (rbas(ii,2)*par*bohr, ii=1,3)
c$$$	 write(992,'(3(f12.7,3x))') (rbas(ii,3)*par*bohr, ii=1,3)
c$$$	 write(992,'(a)') "PRIMCOORD"
c$$$	 write(992,'(2(i5, 2x))') natom, 1
c$$$	 do jj=1, natom
c$$$	   write(992,'(i5, 2x, 3(f12.6,3x))') nint(z(is(jj))),
c$$$     $      (tau(ii,jj)*par*bohr, ii=1,3)
c$$$	 enddo        
c$$$	 close(992)
c$$$	 endif
c$$$      end
c$$$
