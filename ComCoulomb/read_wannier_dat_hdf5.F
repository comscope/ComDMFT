      subroutine read_wannier_dat_hdf5
      use wannier90_mod
      use manager_mod
      use solid_mod
      use parallel_mod
!      use units_mod

      use hdf5io_mod
      use hdf5
      implicit none
      include "mpif.h"
      integer :: iatom,ii,jj,ll,kk,isort,l,n
      integer(hid_t) :: fid_dat = -1
      integer :: h5err
      integer :: ng(4)
      integer :: nv(4)
      integer :: itemp
      integer, allocatable :: tempwindow(:,:)
      complex*16 :: ai
      double precision :: dblr,dbli,trbas(3,3),tgbas(3,3)


      integer :: ik,ibnd,iwan,num_bands_temp,num_wann_temp,
     $  ierr,ndiv_temp(3),irk,j,i,indtmp,m,mp_grid(3)
      double precision :: nelectron_dble, real_lattice_temp(3,3),
     $  recip_lattice_temp(3,3)
      complex*16, allocatable :: v_matrix_temp(:,:,:)

!      call fread(MPI_COMM_WORLD,ubi,"ubi",fid_dat)
      call finit()
!      write(iun,*)'1 fid_dat',fid_dat
!      write(iun,*)'seed_name',seed_name
c      call fopen('wannier.dat',fid_dat)
       call fopen(adjustl(trim(wanpath))//'/wannier.dat',
     $ fid_dat)
!      write(iun,*)'2 fid_dat',fid_dat

      nv=0
 
      ng(1)=size(real_lattice_temp,1)
      ng(2)=size(real_lattice_temp,2)
      call fread(MPI_COMM_WORLD,real_lattice_temp,ng,nv,"rbas",fid_dat)  
c      write(179,*) 'real_lattice',real_lattice_temp
      ng(1)=size(recip_lattice_temp,1)
      ng(2)=size(recip_lattice_temp,2)
      call fread(MPI_COMM_WORLD,recip_lattice_temp,ng,nv,"gbas",fid_dat)
c      write(179,*) 'recip_lattice',recip_lattice_temp
      call fread(MPI_COMM_WORLD,num_bands,"num_bands",fid_dat)
c      write(179,*) 'num_bands',num_bands
      call fread(MPI_COMM_WORLD,num_wann,"num_wann",fid_dat)   
c      write(179,*) 'num_wann',num_wann
      ng(1)=size(mp_grid,1)
      call fread(MPI_COMM_WORLD,mp_grid,ng(1),nv(1),"ndiv",fid_dat)
c      write(179,*) 'mp_grid',mp_grid
        num_kpts=mp_grid(1)*mp_grid(2)*mp_grid(3) !!
        if (.not.allocated(kpt_latt)) then
          allocate(kpt_latt(3,num_kpts),stat=ierr)
          if (ierr/=0)
     $      write(iun,*) 'Error allocating kpt_latt in read_wann_chkpt'
        endif

      ng(1)=size(kpt_latt,1)
      ng(2)=size(kpt_latt,2)
      call fread(MPI_COMM_WORLD,kpt_latt,ng,nv,"kpt_latt",fid_dat)
c      write(179,*) 'kpt_latt',kpt_latt
      allocate(include_bands(num_bands))   
      ng(1)=size(include_bands,1)
      call fread(MPI_COMM_WORLD,include_bands,ng(1),nv(1),
     $                      "include_bands",fid_dat) 
c      write(179,*) 'include_bands',include_bands
      allocate(v_matrix(num_bands,num_wann,nqdiv))
      ng(1)=size(v_matrix,1)
      ng(2)=size(v_matrix,2)
      ng(3)=size(v_matrix,3)
      call fread(MPI_COMM_WORLD,v_matrix,ng,nv,"v_matrix",fid_dat)
c      write(179,*) 'v_matrix',v_matrix
      allocate(eigenvalues(num_bands,nqdiv))
      ng(1)=size(eigenvalues,1)
      ng(2)=size(eigenvalues,2)
      call fread(MPI_COMM_WORLD,eigenvalues,ng,nv,"eigenvalues",
     $ fid_dat)
c      write(179,*) 'eigenvalues',eigenvalues
      call fread(MPI_COMM_WORLD,nrpts,"nrpts",fid_dat)
c      write(179,*) 'nrpts',nrpts
      allocate(irvec(3,nrpts))
      ng(1)=size(irvec,1)
      ng(2)=size(irvec,2)
      call fread(MPI_COMM_WORLD,irvec,ng,nv,"irvec",fid_dat)
c      write(179,*) 'irvec',irvec
      allocate(ndegen(nrpts))
      ng(1)=size(ndegen,1)
      call fread(MPI_COMM_WORLD,ndegen,ng(1),nv(1),
     $                      "ndegen",fid_dat)
c      write(179,*) 'ndegen',ndegen
      if (.not.allocated(u_matrix)) then
          allocate(u_matrix(num_wann,num_wann,num_kpts),stat=ierr)
          if (ierr/=0) write(iun,*)
     $      'Error allocating u_matrix in read_wann_chkpt'
      endif
      ng(1)=size(u_matrix,1)
      ng(2)=size(u_matrix,2)
      ng(3)=size(u_matrix,3)
      call fread(MPI_COMM_WORLD,u_matrix,ng,nv,"u_matrix",fid_dat)
c      write(179,*)'u_matrix',u_matrix      
      call fread(MPI_COMM_WORLD,have_disentangled,
     $"have_disentangled",fid_dat)
c      write(179,*)'have',have_disentangled
      if (have_disentangled) then
          if (.not.allocated(lwindow)) then
            allocate(lwindow(num_bands,num_kpts),stat=ierr)
            if (ierr/=0) write(iun,*)
     $        'Error allocating lwindow in read_wann_chkpt'
          endif
          allocate(tempwindow(num_bands,num_kpts))
          ng(1)=size(tempwindow,1)
          ng(2)=size(tempwindow,2)
          call fread(MPI_COMM_WORLD,tempwindow,ng,nv,"lwindow",fid_dat)
          do ii=1,ng(1)
           do jj=1,ng(2)
             if(tempwindow(ii,jj) .eq. 1) then
              lwindow(ii,jj) = .true.
             else
              lwindow(ii,jj) = .false.
             endif
           enddo
          enddo

c       write(iun,*)'lwindow',lwindow
          if (.not.allocated(u_matrix_opt)) then
            allocate(u_matrix_opt(num_bands,num_wann,num_kpts),
     $        stat=ierr)
            if (ierr/=0) write(iun,*)
     $        'Error allocating u_matrix_opt in read_wann_chkpt'
          endif
          ng(1)=size(u_matrix_opt,1)
          ng(2)=size(u_matrix_opt,2)
          ng(3)=size(u_matrix_opt,3)
          call fread(MPI_COMM_WORLD,u_matrix_opt,ng,nv,
     $           "u_matrix_opt",fid_dat)
c          write(179,*)'u_matrix_opt',u_matrix_opt
      endif


      call fclose(fid_dat)
!      endif
      end

