      subroutine read_wan_bni_hdf5
      use comwann_mod
      use hdf5io_mod
      use hdf5
      implicit none
      include 'mpif.h'
      integer(hid_t) :: fid_dat = -1
      integer :: itemp =0
      integer :: ng(4)
      integer :: nv(4)
      integer :: jcg_size,jcgr_size,cg_size,
     $ indxcg_size,indxcgr_r1size,indxcgr_r2size,
     $ indxcgr_r3size,indxcgr_r4size,cgr_size,
     $ iatom,i,it,isort,lb,lpb,l,ll,n,ntl(0:10),
     $ iwan,ii,jj,i_len,
     $  cnt1,cnt2,kk
      character*1 :: dummyc
      complex*16 :: ai
      double precision :: dblr,dbli

      call finit()
      i_len=len_trim(rstpath)
      call fopen(rstpath(1:i_len)//'/info.rst',fid_dat)
!      call fread(MPI_COMM_WORLD,itemp,"irel",fid_dat)
!      call fclose(fid_dat)
c      write(993,*)'itemp=',itemp
c      write(179,*)'bk-readini1'
!        read(1)irel,clight
      call fread(MPI_COMM_WORLD,irel,"irel",fid_dat)
      call fread(MPI_COMM_WORLD,clight,"clight",fid_dat)
c        write(993, *) 'irel=', irel
c        write(993, *) 'clight=', clight
!        if(irel==0) irel_core=0
!        irel_core=max(irel,irel_core)
        nrel=1
        if(irel.eq.2) nrel=2
!        nrel_core=1
!        if(irel_core.eq.2) nrel_core=2
!        read(1)allfile
      call fread(MPI_COMM_WORLD,allfile,"allfile",fid_dat)
c        write(993, *) 'allfile=', allfile
!        read(1)par,natom,nsort
      call fread(MPI_COMM_WORLD,par,"par",fid_dat)
      call fread(MPI_COMM_WORLD,natom,"natom",fid_dat)
      call fread(MPI_COMM_WORLD,nsort,"nsort",fid_dat)
c        write(993, *) 'par=', par
c        write(993, *) 'natom=', natom
c        write(993, *) 'nsort=', nsort
c        write(179,*)'bk-readini2'
      allocate(is(natom))
      allocate(tau(3,natom))
      tau=0.0d0
      allocate(txtel(nsort))
      allocate(z(nsort))
      allocate(h(nsort))
      allocate(nrad(nsort))
      allocate(lmb(nsort))
!        read(1)is
      nv = 0
      ng(1)=size(is,1)
      call fread(MPI_COMM_WORLD,is,ng(1),nv(1),"is",fid_dat)
c        write(993, '(a, 21i3)') 'is(1:natom)=',
c     $    (is(iatom),iatom=1,natom)

!        read(1)b_a,c_a
      call fread(MPI_COMM_WORLD,b_a,"b_a",fid_dat)
      call fread(MPI_COMM_WORLD,c_a,"c_a",fid_dat)
c        write(993, *) 'b_a=', b_a
c        write(993, *) 'c_a=', c_a
!        read(1)rbas
      ng(1)=size(rbas,1)
      ng(2)=size(rbas,2)
      call fread(MPI_COMM_WORLD,rbas,ng,nv,"rbas",fid_dat)
c        write(993, '(a, 3f21.16)') 'rbas(1:3, 1)=', (rbas(i, 1), i=1, 3)
c        write(993, '(a, 3f21.16)') 'rbas(1:3, 2)=', (rbas(i, 2), i=1, 3)
c        write(993, '(a, 3f21.16)') 'rbas(1:3, 3)=', (rbas(i, 3), i=1, 3)
!        read(1)tau
      ng(1)=size(tau,1)
      ng(2)=size(tau,2)
      call fread(MPI_COMM_WORLD,tau,ng,nv,"tau",fid_dat)
        do iatom=1,natom
c          write(993, '(a,3f21.16)') 'tau(1:3, iatom)=',
c     $      (tau(i,iatom),i=1,3)
        enddo
!        read(1)nbndf
      call fread(MPI_COMM_WORLD,nbndf,"nbndf",fid_dat)
c        write(993,*) 'nbndf=', nbndf
!        read(1)ndiv
      ng(1)=size(ndiv,1)
      call fread(MPI_COMM_WORLD,ndiv,ng(1),nv(1),"ndiv",fid_dat)
c        write(993, *) 'ndiv(i)=', (ndiv(i), i=1,3)
c     -------- Transform to internal magnetic field --------------
!        read(1)maxb,maxpb,maxntle
      call fread(MPI_COMM_WORLD,maxb,"maxb",fid_dat)
      call fread(MPI_COMM_WORLD,maxpb,"maxpb",fid_dat)
      call fread(MPI_COMM_WORLD,maxntle,"maxntle",fid_dat)
c        write(993,*) 'maxb,pb ntle ',maxb,maxpb,maxntle
      allocate(ntle(0:maxb,nsort))
      allocate(augm(maxntle,0:maxb,nsort))
c      write(179,*)'bk-readini_end6'

!        read(1)txtel
      ng(1)=size(txtel,1)
      call fread(MPI_COMM_WORLD,txtel,ng(1),nv(1),"txtel",fid_dat)

!        read(1)z
      ng(1)=size(z,1)
      call fread(MPI_COMM_WORLD,z,ng(1),nv(1),"z",fid_dat)
!        read(1)h,nrad
      ng(1)=size(h,1)
      call fread(MPI_COMM_WORLD,h,ng(1),nv(1),"h",fid_dat)
      ng(1)=size(nrad,1)
      call fread(MPI_COMM_WORLD,nrad,ng(1),nv(1),"nrad",fid_dat)
!        read(1)lmb
      ng(1)=size(lmb,1)
      call fread(MPI_COMM_WORLD,lmb,ng(1),nv(1),"lmb",fid_dat)
!        read(1)ntle
      ng(1)=size(ntle,1)
      ng(2)=size(ntle,2)
      call fread(MPI_COMM_WORLD,ntle,ng,nv,"ntle",fid_dat)
!        read(1)augm
      ng(1)=size(augm,1)
      ng(2)=size(augm,2)
      ng(3)=size(augm,3)
      call fread(MPI_COMM_WORLD,augm,ng,nv,"augm",fid_dat)

        do isort=1,nsort
c          write(993,*) 'txtel=',txtel(isort) 
c          write(993,*) 'z(isort)=', z(isort)
c          write(993,*) 'h(isort)=', h(isort)
c          write(993,*) 'nrad(isort)=', nrad(isort)
c          write(993,*) 'lmb(isort)=', lmb(isort)
c          write(993, '(a, 10i3)') 'ntle(0:lmb(isort), isort)',
c     $      (ntle(l,isort),l=0,lmb(isort))
c     ------- Reading the type of the functions for L <= L_b ------
          do l=0,lmb(isort)
            do n=1,ntle(l,isort)
c              write(993, *) 'augm(n,l,isort)=', augm(n,l,isort)
            enddo
          enddo
        enddo
c      write(179,*)'bk-readini_end5'
c      write(179,*)'bk-readini_end4'   
      allocate(lfun(nsort))  
      allocate(lfunm(nsort))
!       read(1)lfun
      ng(1)=size(lfun,1)
      call fread(MPI_COMM_WORLD,lfun,ng(1),nv(1),"lfun",fid_dat)
       do isort=1,nsort
c          write(993,*) 'lfun(isort)=', lfun(isort)
       enddo  
!       read(1)lfunm
      ng(1)=size(lfunm,1)
      call fread(MPI_COMM_WORLD,lfunm,ng(1),nv(1),"lfunm",fid_dat)
       do isort=1,nsort
c          write(993,*) 'lfunm(isort)=', lfunm(isort)
       enddo
!       read(1)maxlfun,maxb,maxel,
!     $   maxnrad
      call fread(MPI_COMM_WORLD,maxlfun,"maxlfun",fid_dat)
      call fread(MPI_COMM_WORLD,maxb,"maxb",fid_dat)
      call fread(MPI_COMM_WORLD,maxel,"maxel",fid_dat)
      call fread(MPI_COMM_WORLD,maxnrad,"maxnrad",fid_dat)
c       write(993,*) 'maxlfun,maxb,maxel,maxnrad ',maxlfun,
c     $   maxb,maxel,maxnrad
!       read(1)irel,
!     $   iun,magn,
!     $   nrel,nspin,limlb,nfun,nqdiv
      call fread(MPI_COMM_WORLD,irel,"irel",fid_dat)
      call fread(MPI_COMM_WORLD,iun,"iun",fid_dat)
      iun=993
      call fread(MPI_COMM_WORLD,magn,"magn",fid_dat)
      call fread(MPI_COMM_WORLD,nrel,"nrel",fid_dat)
      call fread(MPI_COMM_WORLD,nspin,"nspin",fid_dat)
      call fread(MPI_COMM_WORLD,limlb,"limlb",fid_dat)
      call fread(MPI_COMM_WORLD,nfun,"nfun",fid_dat)
      call fread(MPI_COMM_WORLD,nqdiv,"nqdiv",fid_dat)
c       write(993,*) 'irel,iun,magn,nrel,nspin,limlb,nfun,nqdiv',
c     $   irel,iun,
c     $   magn,nrel,nspin,limlb,nfun,nqdiv
!       read(1)amega
      call fread(MPI_COMM_WORLD,amega,"amega",fid_dat)
c       write(993,*) 'amega', amega
!       read(1)n_tau
      call fread(MPI_COMM_WORLD,n_tau,"n_tau",fid_dat)
c       write(993,*) 'n_tau',n_tau
!       read(1)npnt
      call fread(MPI_COMM_WORLD,npnt,"npnt",fid_dat)
c       write(993,*) 'npnt',npnt
!       read(1)maxylm
      call fread(MPI_COMM_WORLD,maxylm,"maxylm",fid_dat)
c       write(993,*) 'maxylm',maxylm
!       read(1)cg_size
      call fread(MPI_COMM_WORLD,cg_size,"cg_size",fid_dat)
c       write(993,*) 'cg_size',cg_size
!       read(1)jcg_size
      call fread(MPI_COMM_WORLD,jcg_size,"jcg_size",fid_dat)
c       write(993,*) 'jcg_size',jcg_size
!       read(1)jcgr_size
      call fread(MPI_COMM_WORLD,jcgr_size,"jcgr_size",fid_dat)
c       write(993,*) 'jcgr_size',jcgr_size
!       read(1)indxcg_size
      call fread(MPI_COMM_WORLD,indxcg_size,"indxcg_size",fid_dat)
c       write(993,*) 'indxcg_size',indxcg_size
!       read(1)indxcgr_r1size
      call fread(MPI_COMM_WORLD,indxcgr_r1size,"indxcgr_r1size",
     & fid_dat)
c       write(993,*) 'indxcgr_r1size',indxcgr_r1size
!       read(1)indxcgr_r2size
      call fread(MPI_COMM_WORLD,indxcgr_r2size,"indxcgr_r2size",
     & fid_dat)
c       write(993,*) 'indxcgr_r2size',indxcgr_r2size
!       read(1)indxcgr_r3size
      call fread(MPI_COMM_WORLD,indxcgr_r3size,"indxcgr_r3size",
     & fid_dat)
c       write(993,*) 'indxcgr_r3size',indxcgr_r3size
!       read(1)indxcgr_r4size
      call fread(MPI_COMM_WORLD,indxcgr_r4size,"indxcgr_r4size",
     & fid_dat)
c       write(993,*) 'indxcgr_r4size',indxcgr_r4size       
!       read(1)cgr_size
      call fread(MPI_COMM_WORLD,cgr_size,"cgr_size",fid_dat)
c       write(993,*) 'cgr_size',cgr_size       
!       read(1)complex_ro
      call fread(MPI_COMM_WORLD,complex_ro,"complex_ro",fid_dat)
!       write(993,*) 'complex_ro',complex_ro
!       read(1)maxplw
      call fread(MPI_COMM_WORLD,maxplw,"maxplw",fid_dat)
c       write(993,*) 'maxplw',maxplw
!       read(1)nplwdif
      call fread(MPI_COMM_WORLD,nplwdif,"nplwdif",fid_dat)
c       write(993,*) 'nplwdif',nplwdif
!       read(1)ngroup
      call fread(MPI_COMM_WORLD,ngroup,"ngroup",fid_dat)
c       write(993,*) 'ngroup',ngroup
!       read(1)nbndtm
      call fread(MPI_COMM_WORLD,nbndtm,"nbndtm",fid_dat)
c       write(993,*) 'nbndtm',nbndtm
!       read(1)maxia
      call fread(MPI_COMM_WORLD,maxia,"maxia",fid_dat)
c       write(993,*) 'maxia',maxia
!       read(1)maxib
      call fread(MPI_COMM_WORLD,maxib,"maxib",fid_dat)
c       write(993,*) 'maxib',maxib
!       read(1)maxic
      call fread(MPI_COMM_WORLD,maxic,"maxic",fid_dat)
c       write(993,*) 'maxic',maxic
!       read(1)nplwbas
      call fread(MPI_COMM_WORLD,nplwbas,"nplwbas",fid_dat)
c       write(993,*) 'nplwbas',nplwbas
!       read(1)numplw
      call fread(MPI_COMM_WORLD,numplw,"numplw",fid_dat)
c       write(993,*) 'numplw',numplw
!       read(1)maxwig
      call fread(MPI_COMM_WORLD,maxwig,"maxwig",fid_dat)
c       write(993,*) 'maxwig',maxwig
!       read(1)maxwig_j
      call fread(MPI_COMM_WORLD,maxwig_j,"maxwig_j",fid_dat)
c       write(993,*) 'maxwig_j',maxwig_j
!       read(1)maxwf
      call fread(MPI_COMM_WORLD,maxwf,"maxwf",fid_dat)
c       write(993,*) 'maxwf',maxwf
!       read(1)nspin_0
      call fread(MPI_COMM_WORLD,nspin_0,"nspin_0",fid_dat)
c       write(993,*) 'nspin_0',nspin_0


      
      nbasmpw=nrel*nbndtm
c      write(179,*)'bk-nbasmpw',nbasmpw

c      write(179,*)'bk-readini_end3'
      allocate(ind_wf(maxlfun,nsort))
      allocate(indbas(maxntle,limlb,natom))
      allocate(indbasa(2,maxntle,limlb,nsort))
      allocate(io_lem(natom))
      allocate(lf_isz(maxel,nsort))
      allocate(dr(0:maxnrad,nsort))
      allocate(r(0:maxnrad,nsort))
      allocate(bett((maxylm+1)**2))
      allocate(cg(cg_size))
      allocate(jcg(jcg_size))
      allocate(jcgr(jcgr_size))
      allocate(indxcg(indxcg_size))
      allocate(indxcgr(indxcgr_r1size,indxcgr_r2size,
     $   indxcgr_r3size,0:indxcgr_r4size-1))
      allocate(cgr(cgr_size))
      allocate(gbs(3,maxplw))
      allocate(gbs_sym(nplwdif,ngroup))
      allocate(i_kref(nqdiv))
      allocate(igbs(3,maxplw))
      allocate(indgb(nbndtm,nqdiv))
      allocate(indplw(-maxia:maxia,-maxib:maxib,-maxic:maxic))
      allocate(ip(natom,48))
      allocate(iplf_bk(nplwbas,nqdiv))
      allocate(k_group(nqdiv))
      allocate(nbask(npnt))
      allocate(pnt(3,nqdiv))
      allocate(shift(3,48))
      
      allocate(sovr(numplw))
      if(complex_ro) allocate(sovi(numplw))
 
      allocate(tshift(3,natom,ngroup))
      allocate(u(maxwig,48))
      allocate(uj(maxwig_j,48))
      allocate(li_val(maxlfun,nsort))
 
c      write(179,*)'bk-readini_end2'
!       read(1)ind_wf
      ng(1)=size(ind_wf,1)
      ng(2)=size(ind_wf,2)
      call fread(MPI_COMM_WORLD,ind_wf,ng,nv,"ind_wf",fid_dat)
       do isort=1,nsort
        do ii=1,maxlfun
!         write(993, '(a,i10)')'ind_wf(1:maxlfun,isort)',
!     $      ind_wf(ii,isort)
        enddo
       enddo
!       read(1)indbas  
      ng(1)=size(indbas,1)
      ng(2)=size(indbas,2)
      ng(3)=size(indbas,3)
      call fread(MPI_COMM_WORLD,indbas,ng,nv,"indbas",fid_dat)
       do ii=1,maxntle
        do jj=1,limlb
         do kk=1,natom
!           write(993,'(a,i10)')'indbas',indbas(ii,jj,kk)
         enddo
        enddo
       enddo
!       read(1)io_lem     
      ng(1)=size(io_lem,1)
      call fread(MPI_COMM_WORLD,io_lem,ng(1),nv(1),"io_lem",fid_dat)
       do ii=1,natom
!         write(993,'(a,i10)')'io_lem',io_lem(ii)
       enddo
!       read(1)lf_isz 
      ng(1)=size(lf_isz,1)
      ng(2)=size(lf_isz,2)
      call fread(MPI_COMM_WORLD,lf_isz,ng,nv,"lf_isz",fid_dat)
       do ii=1,maxel
        do jj=1,nsort
!         write(993,'(a,i10)')'lf_isz',lf_isz(ii,jj)
        enddo
       enddo
!       read(1)gbas
      ng(1)=size(gbas,1)
      ng(2)=size(gbas,2)
      call fread(MPI_COMM_WORLD,gbas,ng,nv,"gbas",fid_dat)
       do ii=1,3
        do jj=1,3
!         write(993,'(a,f20.10)')'gbas',gbas(ii,jj)
        enddo
       enddo
!       read(1)dr
      ng(1)=size(dr,1)
      ng(2)=size(dr,2)
      call fread(MPI_COMM_WORLD,dr,ng,nv,"dr",fid_dat)
       do ii=0,maxnrad
        do jj=1,nsort
!         write(993,'(a,f20.10)')'dr',dr(ii,jj)
        enddo
       enddo
!       read(1)indbasa
      ng(1)=size(indbasa,1)
      ng(2)=size(indbasa,2)
      ng(3)=size(indbasa,3)
      ng(4)=size(indbasa,4)
      call fread(MPI_COMM_WORLD,indbasa,ng,nv,"indbasa",fid_dat)
       do ll=1,2
        do ii=1,maxntle
         do jj=1,limlb
          do kk=1,nsort
!           write(993,'(a,i10)')'indbasa',indbasa(ll,ii,jj,kk)
          enddo
         enddo
        enddo
       enddo
!       read(1)r
      ng(1)=size(r,1)
      ng(2)=size(r,2)
      call fread(MPI_COMM_WORLD,r,ng,nv,"r",fid_dat)
       do ii=0,maxnrad
        do jj=1,nsort
!         write(993,'(a,f20.10)')'r',r(ii,jj)
        enddo
       enddo
!       read(1)bett
      ng(1)=size(bett,1)
      call fread(MPI_COMM_WORLD,bett,ng(1),nv(1),"bett",fid_dat)
       do ii=1,(maxylm+1)**2
!         write(993,'(a,f20.10)')'bett',bett(ii)
       enddo
!       read(1)cg
      ng(1)=size(cg,1)
      call fread(MPI_COMM_WORLD,cg,ng(1),nv(1),"cg",fid_dat)
       do ii=1,cg_size
!         write(993,'(a,f20.10)')'cg',cg(ii)
       enddo
!       read(1)jcg
      ng(1)=size(jcg,1)
      call fread(MPI_COMM_WORLD,jcg,ng(1),nv(1),"jcg",fid_dat)
       do ii=1,jcg_size
!         write(993,'(a,f20.10)')'jcg',jcg(ii)
       enddo
!       read(1)jcgr
      ng(1)=size(jcgr,1)
      call fread(MPI_COMM_WORLD,jcgr,ng(1),nv(1),"jcgr",fid_dat)
       do ii=1,jcgr_size
!         write(993,'(a,f20.10)')'jcgr',jcgr(ii)
       enddo
!       read(1)indxcg
      ng(1)=size(indxcg,1)
      call fread(MPI_COMM_WORLD,indxcg,ng(1),nv(1),"indxcg",fid_dat)
       do ii=1,indxcg_size
!         write(993,'(a,f20.10)')'indxcg',indxcg(ii)
       enddo
!       read(1)indxcgr
      ng(1)=size(indxcgr,1)
      ng(2)=size(indxcgr,2)
      ng(3)=size(indxcgr,3)
      ng(4)=size(indxcgr,4)
      call fread(MPI_COMM_WORLD,indxcgr,ng,nv,"indxcgr",fid_dat)
       do ii=1,indxcgr_r1size
        do jj=1,indxcgr_r2size
         do kk=1,indxcgr_r3size
          do ll=1,indxcgr_r4size
!            write(993,'(a,f20.10)')'indxcgr',indxcgr(ii,jj,kk,ll)
          enddo
         enddo
        enddo
       enddo
!       read(1)cgr
      ng(1)=size(cgr,1)
      call fread(MPI_COMM_WORLD,cgr,ng(1),nv(1),"cgr",fid_dat)
       do ii=1,cgr_size
!         write(993,*)'cgr',cgr
       enddo

!       read(1)gbs
      ng(1)=size(gbs,1)
      ng(2)=size(gbs,2)
      call fread(MPI_COMM_WORLD,gbs,ng,nv,"gbs",fid_dat)
       do ii=1,3
        do jj=1,maxplw
!         write(993,'(a,f20.10)')'gbs',gbs(ii,jj)
        enddo
       enddo
!       read(1)gbs_sym
      ng(1)=size(gbs_sym,1)
      ng(2)=size(gbs_sym,2)
      call fread(MPI_COMM_WORLD,gbs_sym,ng,nv,"gbs_sym",fid_dat)
       do ii=1,nplwdif
        do jj=1,ngroup
!         write(993,'(a,i10)')'gbs_sym',gbs_sym(ii,jj)
        enddo
       enddo
!       read(1)i_kref
      ng(1)=size(i_kref,1)
      call fread(MPI_COMM_WORLD,i_kref,ng(1),nv(1),"i_kref",fid_dat)
       do ii=1,nqdiv
!         write(993,'(a,i10)')'i_kref',i_kref(ii)
       enddo
!       read(1)igbs
      ng(1)=size(igbs,1)
      ng(2)=size(igbs,2)
      call fread(MPI_COMM_WORLD,igbs,ng,nv,"igbs",fid_dat)
       do ii=1,3
        do jj=1,maxplw
!         write(993,'(a,i10)')'igbs',igbs(ii,jj)
        enddo
       enddo 
!       read(1)indgb      
      ng(1)=size(indgb,1)
      ng(2)=size(indgb,2)
      call fread(MPI_COMM_WORLD,indgb,ng,nv,"indgb",fid_dat)
       do ii=1,nbndtm
        do jj=1,nqdiv
!         write(993,'(a,i10)')'indgb',indgb(ii,jj)
        enddo
       enddo
!       read(1)indplw
      ng(1)=size(indplw,1)
      ng(2)=size(indplw,2)
      ng(3)=size(indplw,3)
      call fread(MPI_COMM_WORLD,indplw,ng,nv,"indplw",fid_dat)
       do ii=-maxia,maxia
        do jj=-maxib,maxib
         do ll=-maxic,maxic
!           write(993,'(a,i10)')'indplw',indplw(ii,jj,ll)
         enddo
        enddo
       enddo
!       read(1)ip      
      ng(1)=size(ip,1)
      ng(2)=size(ip,2)
      call fread(MPI_COMM_WORLD,ip,ng,nv,"ip",fid_dat) 
       do ii=1,natom
        do jj=1,48
!         write(993,'(a,i10)')'ip',ip(ii,jj)
        enddo
       enddo
!       read(1)iplf_bk
      ng(1)=size(iplf_bk,1)
      ng(2)=size(iplf_bk,2)
      call fread(MPI_COMM_WORLD,iplf_bk,ng,nv,"iplf_bk",fid_dat)
       do ii=1,nplwbas
        do jj=1,nqdiv
!         write(993,'(a,i10)')'iplf_bk',iplf_bk(ii,jj)
        enddo
       enddo
!       read(1)k_group
      ng(1)=size(k_group,1)
      call fread(MPI_COMM_WORLD,k_group,ng(1),nv(1),"k_group",fid_dat)
       do ii=1,nqdiv
!         write(993,'(a,i10)')'k_group',k_group(ii)
       enddo
!       read(1)nbask
      ng(1)=size(nbask,1)
      call fread(MPI_COMM_WORLD,nbask,ng(1),nv(1),"nbask",fid_dat)
       do ii=1,npnt
!         write(993,'(a,i10)')'nbask',nbask(ii)
       enddo
!       read(1)nrdiv
      ng(1)=size(nrdiv,1)
      call fread(MPI_COMM_WORLD,nrdiv,ng(1),nv(1),"nrdiv",fid_dat)
       do ii=1,3
!         write(993,'(a,i10)')'nrdiv',nrdiv(ii)
       enddo
!       read(1)pnt
      ng(1)=size(pnt,1)
      ng(2)=size(pnt,2)
      call fread(MPI_COMM_WORLD,pnt,ng,nv,"pnt",fid_dat)
       do ii=1,3
        do jj=1,nqdiv
!         write(993,'(a,f20.10)')'pnt',pnt(ii,jj)
        enddo
       enddo
!       read(1)shift
      ng(1)=size(shift,1)
      ng(2)=size(shift,2)
      call fread(MPI_COMM_WORLD,shift,ng,nv,"shift",fid_dat)
       do ii=1,3
        do jj=1,48
!         write(993,'(a,f20.10)')'shift',shift(ii,jj)
        enddo
       enddo
!       read(1)sovr
      ng(1)=size(sovr,1)
      call fread(MPI_COMM_WORLD,sovr,ng(1),nv(1),"sovr",fid_dat)
       do ii=1,numplw
!         write(993,'(a,f20.10)')'sovr',sovr(ii)
       enddo
       if(complex_ro) then
!          read(1)sovi
        ng(1)=size(sovi,1)
        call fread(MPI_COMM_WORLD,sovi,ng(1),nv(1),"sovi",fid_dat)
         do ii=1,numplw
!          write(993,'(a,f20.10)')'sovi',sovi(ii)
         enddo
       endif
!       read(1)tshift
      ng(1)=size(tshift,1)
      ng(2)=size(tshift,2)
      ng(3)=size(tshift,3)
      call fread(MPI_COMM_WORLD,tshift,ng,nv,"tshift",fid_dat)
       do ii=1,3
        do jj=1,natom
         do ll=1,ngroup
!           write(993,'(a,f20.10)')'tshift',tshift(ii,jj,ll)
         enddo
        enddo
       enddo
!       read(1)u
      ng(1)=size(u,1)
      ng(2)=size(u,2)
      call fread(MPI_COMM_WORLD,u,ng,nv,"u",fid_dat)
       do ii=1,maxwig
        do jj=1,48
!         write(993,'(a,f20.10)')'u',u(ii,jj)
        enddo
       enddo
!       read(1)uj
      ng(1)=size(uj,1)
      ng(2)=size(uj,2)
      call fread(MPI_COMM_WORLD,uj,ng,nv,"uj",fid_dat)
       do ii=1,maxwig_j
        do jj=1,48
!         write(993,*)'uj',uj(ii,jj)
        enddo
       enddo
!       read(1)li_val
      ng(1)=size(li_val,1)
      ng(2)=size(li_val,2)
      call fread(MPI_COMM_WORLD,li_val,ng,nv,"li_val",fid_dat)
       do ii=1,maxlfun
        do jj=1,nsort
!         write(993,*)'li_val',li_val(ii,jj)
        enddo
       enddo


      call fclose(fid_dat)
      end
