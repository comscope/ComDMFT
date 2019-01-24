      subroutine read_wannier_dat_hdf5
      use comlowh_mod
!      use manager_mod
!      use atom_mod
!      use solid_mod
!      use parallel_mod
!      use units_mod

      use hdf5io_mod
      use hdf5
      implicit none
      include "mpif.h"
      integer :: ik,ibnd,iwan
      double precision :: nelectron_dble
      integer :: iatom,ii,jj,ll,kk,i,isort,l,n
      integer(hid_t) :: fid_dat = -1
      integer :: h5err
      integer :: ng(4)
      integer :: nv(4)
      integer :: itemp
      complex*16 :: ai
      double precision :: dblr,dbli,trbas(3,3),tgbas(3,3)
!      call fread(MPI_COMM_WORLD,ubi,"ubi",fid_dat)
      call finit()
!      write(iun,*)'1 fid_dat',fid_dat
!      write(iun,*)'seed_name',seed_name
      call fopen('wannier.dat',fid_dat)

!      write(iun,*)'2 fid_dat',fid_dat

      nv=0
 
      ng(1)=size(real_lattice,1)
      ng(2)=size(real_lattice,2)
      call fread(MPI_COMM_WORLD,real_lattice,ng,nv,"rbas",fid_dat)  
!      write(179,*) 'real_lattice',real_lattice
      ng(1)=size(recip_lattice,1)
      ng(2)=size(recip_lattice,2)
      call fread(MPI_COMM_WORLD,recip_lattice,ng,nv,"gbas",fid_dat)
!      write(179,*) 'recip_lattice',recip_lattice
      call fread(MPI_COMM_WORLD,num_bands,"num_bands",fid_dat)
!      write(179,*) 'num_bands',num_bands
      call fread(MPI_COMM_WORLD,num_wann,"num_wann",fid_dat)   
!      write(179,*) 'num_wann',num_wann
      ng(1)=size(ndiv,1)
      call fread(MPI_COMM_WORLD,ndiv,ng(1),nv(1),"ndiv",fid_dat)
!      write(179,*) 'ndiv',ndiv
      nqdiv=ndiv(1)*ndiv(2)*ndiv(3)
      allocate(kpt_latt(3,nqdiv))

      ng(1)=size(kpt_latt,1)
      ng(2)=size(kpt_latt,2)
      call fread(MPI_COMM_WORLD,kpt_latt,ng,nv,"kpt_latt",fid_dat)
!      write(179,*) 'kpt_latt',kpt_latt
      allocate(include_bands(num_bands))   
      ng(1)=size(include_bands,1)
      call fread(MPI_COMM_WORLD,include_bands,ng(1),nv(1),
     $                      "include_bands",fid_dat) 
!      write(179,*) 'include_bands',include_bands
      allocate(v_matrix(num_bands,num_wann,nqdiv))
      ng(1)=size(v_matrix,1)
      ng(2)=size(v_matrix,2)
      ng(3)=size(v_matrix,3)
      call fread(MPI_COMM_WORLD,v_matrix,ng,nv,"v_matrix",fid_dat)
!      write(179,*) 'v_matrix',v_matrix
      allocate(eigenvalues(num_bands,nqdiv))
      ng(1)=size(eigenvalues,1)
      ng(2)=size(eigenvalues,2)
      call fread(MPI_COMM_WORLD,eigenvalues,ng,nv,"eigenvalues",
     $ fid_dat)
!      write(179,*) 'eigenvalues',eigenvalues
      call fread(MPI_COMM_WORLD,nrpts,"nrpts",fid_dat)
!      write(179,*) 'nrpts',nrpts
      allocate(irvec(3,nrpts))
      ng(1)=size(irvec,1)
      ng(2)=size(irvec,2)
      call fread(MPI_COMM_WORLD,irvec,ng,nv,"irvec",fid_dat)
!      write(179,*) 'irvec',irvec
      allocate(ndegen(nrpts))
      ng(1)=size(ndegen,1)
      call fread(MPI_COMM_WORLD,ndegen,ng(1),nv(1),
     $                      "ndegen",fid_dat)
!      write(179,*) 'ndegen',ndegen
      call fclose(fid_dat)
!      endif


      nelectron_dble=0.0
      do ik=1, nqdiv
        do ibnd=1, num_bands
          do iwan=1, num_wann
            if (eigenvalues(ibnd,ik) .lt. 0.0d0) then
              nelectron_dble=nelectron_dble+1.0d0/dble(nqdiv)*
     $          cdabs(v_matrix(ibnd,iwan,ik))**2
            endif
          enddo
        enddo
      enddo

c      nelectron=dble(nint(nelectron_dble*num_cell))/num_cell
c     $  +delta_charge
c      write(iun,'(a, 3f12.6)')
c     $  'nelectron in one cell', nelectron_dble,
c     $  delta_charge, nelectron

      if (is_spinorbit .eq. 1) then
        nelectron=dble(nint(nelectron_dble*num_cell))/num_cell
     $    +delta_charge
      else
        nelectron=dble(nint(nelectron_dble*num_cell*2))/2.0d0/num_cell
     $    +delta_charge
      endif
      write(iun,'(a, 3f12.6)')
     $  'nelectron in one cell', nelectron_dble,
     $  delta_charge, nelectron

      end

