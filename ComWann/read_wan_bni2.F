      subroutine read_wan_bni2
      use comwann_mod
      implicit none
      include 'mpif.h'
      integer :: jcg_size,jcgr_size,cg_size,
     $  indxcg_size,indxcgr_r1size,indxcgr_r2size,
     $  indxcgr_r3size,indxcgr_r4size,cgr_size,
     $  iatom,i,isort,l,n,jj,ii,kk,ll,i_len
!     $ iwan,ii,jj,
!     $  cnt1,cnt2,kk
!     character*1 :: dummyc
!     complex*16 :: ai
!     double precision :: dblr,dbli
!     write(179,*)'bk-readini1'
      i_len=len_trim(rstpath)
      if(maswrk) then
        open(1,file=rstpath(1:i_len)//'/info2.rst',
     $    status='old',form='unformatted')
        read(1)maxnsym
        write(993,*) 'maxnsym',maxnsym
        read(1)maxmt
        write(993,*) 'maxmt',maxmt
        read(1)irel_core
        write(993,*) 'irel_core',irel_core
        read(1)nrel_core
        write(993,*) 'nrel_core',nrel_core
        read(1)maxl_cor
        write(993,*) 'maxl_cor',maxl_cor
        read(1)ncormax
        write(993,*) 'ncormax',ncormax
        read(1)max_rl
        write(993,*) 'max_rl',max_rl
        read(1)nspin_1
        write(993,*) 'nspin_1',nspin_1
        read(1)maxmtcor
        write(993,*) 'maxmtcor',maxmtcor
      endif

      if(goparr) then
        call mpi_bcast(maxnsym,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(maxmt,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)        
        call mpi_bcast(irel_core,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(nrel_core,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(maxl_cor,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(ncormax,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(max_rl,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(nspin_1,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
        call mpi_bcast(maxmtcor,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)                
      endif
      
      allocate(smt(nsort))
      smt=0.0d0
      allocate(indmt(maxnsym,nsort,nspin_0))
      indmt=0
      allocate(indfun0(2,maxntle,nrel*maxb+1,nsort))
      indfun0=0
      allocate(iat_1(nsort))
      iat_1=0
      allocate(lc_max(nsort))
      lc_max=0
      allocate(lm_isz(maxel,nsort))
      lm_isz=0
      allocate(ncor(nsort))
      ncor=0
      allocate(indcor(ncormax,nspin_0,nsort))
      indcor=0
      allocate(li_cor(ncormax,nsort))
      li_cor=0
      allocate(rl(0:maxnrad,0:max_rl,nsort))
      rl=0.0d0
      
c     write(179,*)'bk-readini_end2'
      if(maswrk) then
        read(1)smt
        read(1)indmt
        read(1)indfun0
        read(1)iat_1
        read(1)lc_max
        read(1)lm_isz
        read(1)ncor
        read(1)indcor
        read(1)li_cor
        read(1)rl
        close(1)        
      endif
      if(goparr) then
        call mpi_bcast(smt,nsort,
     $    mpi_double_precision,0,MPI_COMM_WORLD, mpierr)
        call mpi_bcast(indmt,maxnsym*nsort*nspin_0,
     $    mpi_integer,0,MPI_COMM_WORLD, mpierr)
        call mpi_bcast(indfun0,2*maxntle*(nrel*maxb+1)*nsort,
     $    mpi_integer,0,MPI_COMM_WORLD, mpierr)
        call mpi_bcast(iat_1,nsort,
     $    mpi_integer,0,MPI_COMM_WORLD, mpierr)
        call mpi_bcast(lc_max,nsort,
     $    mpi_integer,0,MPI_COMM_WORLD, mpierr)
        call mpi_bcast(lm_isz,maxel*nsort,
     $    mpi_integer,0,MPI_COMM_WORLD, mpierr)
        call mpi_bcast(ncor,nsort,
     $    mpi_integer,0,MPI_COMM_WORLD, mpierr)
        call mpi_bcast(indcor,ncormax*nspin_0*nsort,
     $    mpi_integer,0,MPI_COMM_WORLD, mpierr)
        call mpi_bcast(li_cor,ncormax*nsort,
     $    mpi_integer,0,MPI_COMM_WORLD, mpierr)
        call mpi_bcast(rl,(1+maxnrad)*(max_rl+1)*nsort,
     $    mpi_double_precision,0,MPI_COMM_WORLD, mpierr)                                                                                                 
      endif

      end
