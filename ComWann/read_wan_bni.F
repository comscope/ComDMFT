      subroutine read_wan_bni
      use comwann_mod
      implicit none
      include 'mpif.h'
      integer :: jcg_size,jcgr_size,cg_size,
     $ indxcg_size,indxcgr_r1size,indxcgr_r2size,
     $ indxcgr_r3size,indxcgr_r4size,cgr_size,
     $ iatom,i,isort,l,n,jj,ii,kk,ll,i_len
!     $ iwan,ii,jj,
!     $  cnt1,cnt2,kk
!      character*1 :: dummyc
!      complex*16 :: ai
!      double precision :: dblr,dbli
!      write(179,*)'bk-readini1'
      i_len=len_trim(rstpath)
      if(maswrk) then
        open(1,file=rstpath(1:i_len)//'/info.rst',
     $ status='old',form='unformatted')
        read(1)irel,clight
        write(993, *) 'irel=', irel
        write(993, *) 'clight=', clight
!        if(irel==0) irel_core=0
!        irel_core=max(irel,irel_core)
        nrel=1
        if(irel.eq.2) nrel=2
!        nrel_core=1
!        if(irel_core.eq.2) nrel_core=2
        read(1)allfile
        write(993, *) 'allfile=', allfile
        read(1)par,natom,nsort
        write(993, *) 'par=', par
        write(993, *) 'natom=', natom
        write(993, *) 'nsort=', nsort
!        write(179,*)'bk-readini2'
      endif
      if(goparr) then
        call mpi_bcast(clight,1,mpi_double_precision,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(irel,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(nrel,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
!        call mpi_bcast(nrel_core,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(allfile,72,mpi_character,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(par,1,mpi_double_precision,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(natom,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(nsort,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
c        write(179,*)'bk-readini3'
      endif
      allocate(is(natom))
      allocate(tau(3,natom))
      tau=0.0d0
      allocate(txtel(nsort))
      allocate(z(nsort))
      allocate(h(nsort))
      allocate(nrad(nsort))
      allocate(lmb(nsort))
      if(maswrk) then
        read(1)is
        write(993, '(a, 21i3)') 'is(1:natom)=',
     $    (is(iatom),iatom=1,natom)

        read(1)b_a,c_a
        write(993, *) 'b_a=', b_a
        write(993, *) 'c_a=', c_a
        read(1)rbas
        write(993, '(a, 3f21.16)') 'rbas(1:3, 1)=', (rbas(i, 1), i=1, 3)
        write(993, '(a, 3f21.16)') 'rbas(1:3, 2)=', (rbas(i, 2), i=1, 3)
        write(993, '(a, 3f21.16)') 'rbas(1:3, 3)=', (rbas(i, 3), i=1, 3)
        read(1)tau
        do iatom=1,natom
          write(993, '(a,3f21.16)') 'tau(1:3, iatom)=',
     $      (tau(i,iatom),i=1,3)
        enddo
        read(1)nbndf
        write(993,*) 'nbndf=', nbndf
        read(1)ndiv
        write(993, *) 'ndiv(i)=', (ndiv(i), i=1,3)
        read(1)maxb,maxpb,maxntle
        write(993,*) 'maxb,pb ntle ',maxb,maxpb,maxntle
      endif
      if(goparr) then
        call mpi_bcast(maxb,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(maxpb,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(maxntle,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
      endif
      allocate(ntle(0:maxb,nsort))
      allocate(augm(maxntle,0:maxb,nsort))
c      write(179,*)'bk-readini_end6'
      if(maswrk) then

        read(1)txtel
        read(1)z
        read(1)h,nrad
        read(1)lmb
        read(1)ntle
        read(1)augm
        do isort=1,nsort
          write(993,*) 'txtel=',txtel(isort) 
          write(993,*) 'z(isort)=', z(isort)
          write(993,*) 'h(isort)=', h(isort)
          write(993,*) 'nrad(isort)=', nrad(isort)
          write(993,*) 'lmb(isort)=', lmb(isort)
          write(993, '(a, 10i3)') 'ntle(0:lmb(isort), isort)',
     $      (ntle(l,isort),l=0,lmb(isort))
          do l=0,lmb(isort)
            do n=1,ntle(l,isort)
              write(993, *) 'augm(n,l,isort)=', augm(n,l,isort)
            enddo
          enddo
        enddo
      endif
c      write(179,*)'bk-readini_end5'
      if(goparr) then
        call mpi_bcast(b_a,1,mpi_double_precision,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(c_a,1,mpi_double_precision,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(is,natom,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(rbas,9,mpi_double_precision,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(tau,3*natom,mpi_double_precision,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(nbndf,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(ndiv,3,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(z,nsort,mpi_double_precision,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(h,nsort,mpi_double_precision,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(nrad,nsort,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(augm,3*nsort*(maxb+1)*maxntle,mpi_character,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(lmb,nsort,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(ntle,nsort*(maxb+1),mpi_integer,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(txtel,4*nsort,mpi_character,
     $    0,MPI_COMM_WORLD,mpierr)

      endif
c      write(179,*)'bk-readini_end4'   
      allocate(lfun(nsort))  
      allocate(lfunm(nsort))
      if(maswrk) then
       read(1)lfun
       do isort=1,nsort
          write(993,*) 'lfun(isort)=', lfun(isort)
       enddo  
       read(1)lfunm
       do isort=1,nsort
          write(993,*) 'lfunm(isort)=', lfunm(isort)
       enddo
       read(1)maxlfun,maxb,maxel,
     $   maxnrad
       write(993,*) 'maxlfun,maxb,maxel,maxnrad ',maxlfun,
     $   maxb,maxel,maxnrad
       read(1)irel,
     $   iun,magn,
     $   nrel,nspin,limlb,nfun,nqdiv
       iun=993
       write(993,*) 'irel,iun,magn,nrel,nspin,limlb,nfun,nqdiv',
     $   irel,iun,
     $   magn,nrel,nspin,limlb,nfun,nqdiv
       read(1)amega
       write(993,*) 'amega', amega
       read(1)n_tau
       write(993,*) 'n_tau',n_tau
       read(1)npnt
       write(993,*) 'npnt',npnt
       read(1)maxylm
       write(993,*) 'maxylm',maxylm
       read(1)cg_size
       write(993,*) 'cg_size',cg_size
       read(1)jcg_size
       write(993,*) 'jcg_size',jcg_size
       read(1)jcgr_size
       write(993,*) 'jcgr_size',jcgr_size
       read(1)indxcg_size
       write(993,*) 'indxcg_size',indxcg_size
       read(1)indxcgr_r1size
       write(993,*) 'indxcgr_r1size',indxcgr_r1size
       read(1)indxcgr_r2size
       write(993,*) 'indxcgr_r2size',indxcgr_r2size
       read(1)indxcgr_r3size
       write(993,*) 'indxcgr_r3size',indxcgr_r3size
       read(1)indxcgr_r4size
       write(993,*) 'indxcgr_r4size',indxcgr_r4size       
       read(1)cgr_size
       write(993,*) 'cgr_size',cgr_size       
       read(1)complex_ro
       write(993,*) 'complex_ro',complex_ro
       read(1)maxplw
       write(993,*) 'maxplw',maxplw
       read(1)nplwdif
       write(993,*) 'nplwdif',nplwdif
       read(1)ngroup
       write(993,*) 'ngroup',ngroup
       read(1)nbndtm
       write(993,*) 'nbndtm',nbndtm
       read(1)maxia
       write(993,*) 'maxia',maxia
       read(1)maxib
       write(993,*) 'maxib',maxib
       read(1)maxic
       write(993,*) 'maxic',maxic
       read(1)nplwbas
       write(993,*) 'nplwbas',nplwbas
       read(1)numplw
       write(993,*) 'numplw',numplw
       read(1)maxwig
       write(993,*) 'maxwig',maxwig
       read(1)maxwig_j
       write(993,*) 'maxwig_j',maxwig_j
       read(1)maxwf
       write(993,*) 'maxwf',maxwf
       read(1)nspin_0
       write(993,*) 'nspin_0',nspin_0

      endif

      if(goparr) then
        call mpi_bcast(lfun,nsort,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(lfunm,nsort,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(maxlfun,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(maxel,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(maxnrad,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(iun,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(magn,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(nspin,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(limlb,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(nfun,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(nqdiv,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(amega,1,mpi_double_precision,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(n_tau,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(npnt,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(maxylm,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(cg_size,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(jcg_size,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(jcgr_size,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(indxcg_size,1,mpi_integer,0,MPI_COMM_WORLD,
     $     mpierr)
        call mpi_bcast(cgr_size,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(indxcgr_r1size,1,
     $     mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(indxcgr_r2size,1,
     $     mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(indxcgr_r3size,1,
     $     mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(indxcgr_r4size,1,
     $     mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(maxplw,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(nplwdif,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(ngroup,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(nbndtm,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(complex_ro,1,mpi_logical,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(maxia,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(maxib,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(maxic,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(nplwbas,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(numplw,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(maxwig,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(maxwig_j,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(maxwf,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(nspin_0,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
      endif
      
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
      if(maswrk) then
       read(1)ind_wf
       do isort=1,nsort
        do ii=1,maxlfun
!         write(993, '(a,i10)')'ind_wf(1:maxlfun,isort)',
!     $      ind_wf(ii,isort)
        enddo
       enddo
       read(1)indbas  
       do ii=1,maxntle
        do jj=1,limlb
         do kk=1,natom
!           write(993,'(a,i10)')'indbas',indbas(ii,jj,kk)
         enddo
        enddo
       enddo
       read(1)io_lem     
       do ii=1,natom
!         write(993,'(a,i10)')'io_lem',io_lem(ii)
       enddo
       read(1)lf_isz 
       do ii=1,maxel
        do jj=1,nsort
!         write(993,'(a,i10)')'lf_isz',lf_isz(ii,jj)
        enddo
       enddo
       read(1)gbas
       do ii=1,3
        do jj=1,3
!         write(993,'(a,f20.10)')'gbas',gbas(ii,jj)
        enddo
       enddo
       read(1)dr
       do ii=0,maxnrad
        do jj=1,nsort
!         write(993,'(a,f20.10)')'dr',dr(ii,jj)
        enddo
       enddo
       read(1)indbasa
       do ll=1,2
        do ii=1,maxntle
         do jj=1,limlb
          do kk=1,nsort
!           write(993,'(a,i10)')'indbasa',indbasa(ll,ii,jj,kk)
          enddo
         enddo
        enddo
       enddo
       read(1)r
       do ii=0,maxnrad
        do jj=1,nsort
!         write(993,'(a,f20.10)')'r',r(ii,jj)
        enddo
       enddo
       read(1)bett
       do ii=1,(maxylm+1)**2
!         write(993,'(a,f20.10)')'bett',bett(ii)
       enddo
       read(1)cg
       do ii=1,cg_size
!         write(993,'(a,f20.10)')'cg',cg(ii)
       enddo
       read(1)jcg
       do ii=1,jcg_size
!         write(993,'(a,f20.10)')'jcg',jcg(ii)
       enddo
       read(1)jcgr
       do ii=1,jcgr_size
c$$$         write(993,'(a,f20.10)')'jcgr',jcgr(ii)
       enddo
       read(1)indxcg
       do ii=1,indxcg_size
!         write(993,'(a,f20.10)')'indxcg',indxcg(ii)
       enddo
       read(1)indxcgr
       do ii=1,indxcgr_r1size
        do jj=1,indxcgr_r2size
         do kk=1,indxcgr_r3size
          do ll=1,indxcgr_r4size
!            write(993,'(a,f20.10)')'indxcgr',indxcgr(ii,jj,kk,ll)
          enddo
         enddo
        enddo
       enddo
       read(1)cgr
       do ii=1,cgr_size
!         write(993,*)'cgr',cgr
       enddo

       read(1)gbs
       do ii=1,3
        do jj=1,maxplw
!         write(993,'(a,f20.10)')'gbs',gbs(ii,jj)
        enddo
       enddo
       read(1)gbs_sym
       do ii=1,nplwdif
        do jj=1,ngroup
!         write(993,'(a,i10)')'gbs_sym',gbs_sym(ii,jj)
        enddo
       enddo
       read(1)i_kref
       do ii=1,nqdiv
!         write(993,'(a,i10)')'i_kref',i_kref(ii)
       enddo
       read(1)igbs
       do ii=1,3
        do jj=1,maxplw
!         write(993,'(a,i10)')'igbs',igbs(ii,jj)
        enddo
       enddo 
       read(1)indgb      
       do ii=1,nbndtm
        do jj=1,nqdiv
!         write(993,'(a,i10)')'indgb',indgb(ii,jj)
        enddo
       enddo
       read(1)indplw
       do ii=-maxia,maxia
        do jj=-maxib,maxib
         do ll=-maxic,maxic
!           write(993,'(a,i10)')'indplw',indplw(ii,jj,ll)
         enddo
        enddo
       enddo
       read(1)ip       
       do ii=1,natom
        do jj=1,48
!         write(993,'(a,i10)')'ip',ip(ii,jj)
        enddo
       enddo
       read(1)iplf_bk
       do ii=1,nplwbas
        do jj=1,nqdiv
!         write(993,'(a,i10)')'iplf_bk',iplf_bk(ii,jj)
        enddo
       enddo
       read(1)k_group
       do ii=1,nqdiv
!         write(993,'(a,i10)')'k_group',k_group(ii)
       enddo
       read(1)nbask
       do ii=1,npnt
!         write(993,'(a,i10)')'nbask',nbask(ii)
       enddo
       read(1)nrdiv
       do ii=1,3
!         write(993,'(a,i10)')'nrdiv',nrdiv(ii)
       enddo
       read(1)pnt
       do ii=1,3
        do jj=1,nqdiv
!         write(993,'(a,f20.10)')'pnt',pnt(ii,jj)
        enddo
       enddo
       read(1)shift
       do ii=1,3
        do jj=1,48
!         write(993,'(a,f20.10)')'shift',shift(ii,jj)
        enddo
       enddo
       read(1)sovr
       do ii=1,numplw
!         write(993,'(a,f20.10)')'sovr',sovr(ii)
       enddo
       if(complex_ro) then
          read(1)sovi
         do ii=1,numplw
!          write(993,'(a,f20.10)')'sovi',sovi(ii)
         enddo
       endif
       read(1)tshift
       do ii=1,3
        do jj=1,natom
         do ll=1,ngroup
!           write(993,'(a,f20.10)')'tshift',tshift(ii,jj,ll)
         enddo
        enddo
       enddo
       read(1)u
       do ii=1,maxwig
        do jj=1,48
!         write(993,'(a,f20.10)')'u',u(ii,jj)
        enddo
       enddo
       read(1)uj
       do ii=1,maxwig_j
        do jj=1,48
!         write(993,*)'uj',uj(ii,jj)
        enddo
       enddo
       read(1)li_val
       do ii=1,maxlfun
        do jj=1,nsort
!         write(993,*)'li_val',li_val(ii,jj)
        enddo
       enddo


      endif
      if(goparr) then
        call mpi_bcast(ind_wf,maxlfun*nsort,mpi_integer,
     $    0,MPI_COMM_WORLD,mpierr)


        call mpi_bcast(indbas,maxntle*limlb*natom,mpi_integer,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(io_lem,natom,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(lf_isz,maxel*nsort,mpi_integer,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(gbas,9,mpi_double_precision,0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(dr,(maxnrad+1)*nsort,mpi_double_precision,
     $    0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(r,(maxnrad+1)*nsort,mpi_double_precision,
     $    0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(bett,(maxylm+1)**2,mpi_double_precision,
     $    0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(cg,cg_size,mpi_double_precision,
     $    0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(jcg,jcg_size,mpi_integer,
     $    0,MPI_COMM_WORLD,
     $    mpierr)
c$$$        write(*,*) 'jcgr_size', jcgr_size
        call mpi_bcast(jcgr,jcgr_size,mpi_integer,
     $    0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(indxcg,indxcg_size,mpi_integer,0,
     $    MPI_COMM_WORLD,mpierr)
        call mpi_bcast(indxcgr,
     $     indxcgr_r1size*indxcgr_r2size*indxcgr_r3size*indxcgr_r4size,
     $     mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(cgr,cgr_size,mpi_double_complex,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(gbs,3*maxplw,mpi_double_precision,
     $    0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(gbs_sym,nplwdif*ngroup,mpi_integer,
     $    0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(i_kref,nqdiv,mpi_integer,
     $    0,MPI_COMM_WORLD, mpierr)
        call mpi_bcast(igbs,3*maxplw,mpi_integer,
     $    0,MPI_COMM_WORLD, mpierr)
        call mpi_bcast(indgb,nbndtm*nqdiv,mpi_integer,
     $    0,MPI_COMM_WORLD, mpierr)
        call mpi_bcast(indbasa,2*maxntle*limlb*nsort,mpi_integer,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(indplw,(2*maxia+1)*(2*maxib+1)*(2*maxic+1),
     $     mpi_integer,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(ip,natom*48,mpi_integer,
     $    0,MPI_COMM_WORLD, mpierr)
        call mpi_bcast(iplf_bk,nplwbas*nqdiv,mpi_integer,
     $    0,MPI_COMM_WORLD, mpierr)
        call mpi_bcast(k_group,nqdiv,mpi_integer,
     $    0,MPI_COMM_WORLD, mpierr)
        call mpi_bcast(nbask,npnt,mpi_integer,
     $    0,MPI_COMM_WORLD, mpierr)
        call mpi_bcast(nrdiv,3,mpi_integer,
     $    0,MPI_COMM_WORLD, mpierr)
        call mpi_bcast(pnt,3*nqdiv,mpi_double_precision,
     $    0,MPI_COMM_WORLD,
     $    mpierr)
        call mpi_bcast(shift,3*48,mpi_double_precision,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(sovr,numplw,mpi_double_precision,
     $    0,MPI_COMM_WORLD,mpierr)
        if(complex_ro) then
          call mpi_bcast(sovi,numplw,mpi_double_precision,
     $       0,MPI_COMM_WORLD,mpierr)
        endif
        call mpi_bcast(tshift,3*natom*ngroup,mpi_double_precision,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(u,maxwig*48,mpi_double_precision,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(uj,maxwig_j*48,mpi_double_complex,
     $    0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(li_val,maxlfun*nsort,mpi_integer,
     $    0,MPI_COMM_WORLD, mpierr)

      endif
      close(1)
      end
