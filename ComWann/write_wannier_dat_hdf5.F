      subroutine write_wannier_dat_hdf5
      use comwann_mod
!      use manager_mod
!      use atom_mod
!      use solid_mod
!      use parallel_mod
!      use units_mod

      use hdf5io_mod
      use hdf5
      implicit none
      include "mpif.h"
      integer :: iatom,ii,jj,ll,kk,i,isort,l,n
      integer(hid_t) :: fid_dat = -1
      integer :: h5err
      integer :: ng(4)
      integer :: nv(4)
      integer :: itemp
      integer, allocatable :: tempwindow(:,:)
      complex*16 :: ai
      double precision :: dblr,dbli,trbas(3,3),tgbas(3,3)
      logical have_disentangled
!      call fwrite(MPI_COMM_WORLD,ubi,"ubi",fid_dat)
      call finit()
c      write(iun,*)'fid_dat',fid_dat
      call fopen(adjustl(trim(seed_name))//'.dat',fid_dat)
c      write(iun,*)'fid_dat',fid_dat

      trbas=rbas*par*bohr
      tgbas=gbas*2.0d0*pi/par/bohr
      nv=0
 
      ng(1)=size(trbas,1)
      ng(2)=size(trbas,2)
      call fwrite(MPI_COMM_WORLD,trbas,ng,nv,"rbas",fid_dat)  
!      write(180,*)'trbas',trbas

      ng(1)=size(tgbas,1)
      ng(2)=size(tgbas,2)
      call fwrite(MPI_COMM_WORLD,tgbas,ng,nv,"gbas",fid_dat)
!      write(180,*)'tgbas',tgbas
      call fwrite(MPI_COMM_WORLD,num_bands,"num_bands",fid_dat)
!      write(180,*)'num_bands',num_bands
      call fwrite(MPI_COMM_WORLD,num_wann,"num_wann",fid_dat)   
!      write(180,*)'num_wann',num_wann
      ng(1)=size(ndiv,1)
      call fwrite(MPI_COMM_WORLD,ndiv,ng(1),nv(1),"ndiv",fid_dat)
!      write(180,*)'ndiv',ndiv
      ng(1)=size(kpt_latt,1)
      ng(2)=size(kpt_latt,2)
      call fwrite(MPI_COMM_WORLD,kpt_latt,ng,nv,"kpt_latt",fid_dat)
!      write(180,*)'kpt_latt',kpt_latt   
      ng(1)=size(include_bands,1)
      call fwrite(MPI_COMM_WORLD,include_bands,ng(1),nv(1),
     $                      "include_bands",fid_dat) 
!      write(180,*)'include_bands',include_bands
      ng(1)=size(v_matrix,1)
      ng(2)=size(v_matrix,2)
      ng(3)=size(v_matrix,3)
      call fwrite(MPI_COMM_WORLD,v_matrix,ng,nv,"v_matrix",fid_dat)
!      write(180,*)'v_matrix',v_matrix
      ng(1)=size(eigenvalues,1)
      ng(2)=size(eigenvalues,2)
      call fwrite(MPI_COMM_WORLD,eigenvalues,ng,nv,"eigenvalues",
     $ fid_dat)
!      write(180,*)'eigenvalues',eigenvalues
      call fwrite(MPI_COMM_WORLD,nrpts,"nrpts",fid_dat)
!      write(180,*)'nrpts',nrpts
      ng(1)=size(irvec,1)
      ng(2)=size(irvec,2)
      call fwrite(MPI_COMM_WORLD,irvec,ng,nv,"irvec",fid_dat)
!      write(180,*)'irvec',irvec
      ng(1)=size(ndegen,1)
      call fwrite(MPI_COMM_WORLD,ndegen,ng(1),nv(1),
     $                      "ndegen",fid_dat)
!      write(180,*)'ndegen',ndegen
      have_disentangled=.true.
      call fwrite(MPI_COMM_WORLD,have_disentangled,
     $"have_disentangled",fid_dat)

      ng(1)=size(u_matrix,1)
      ng(2)=size(u_matrix,2)
      ng(3)=size(u_matrix,3)
      call fwrite(MPI_COMM_WORLD,u_matrix,ng,nv,"u_matrix",fid_dat)

      if (have_disentangled) then

      ng(1)=size(u_matrix_opt,1)
      ng(2)=size(u_matrix_opt,2)
      ng(3)=size(u_matrix_opt,3)
      call fwrite(MPI_COMM_WORLD,u_matrix_opt,ng,nv,
     $"u_matrix_opt",fid_dat) 

      ng(1)=size(lwindow,1)
      ng(2)=size(lwindow,2)
      allocate(tempwindow(ng(1),ng(2)))
      do ii=1,ng(1)
       do jj=1,ng(2)
         if(lwindow(ii,jj) .eq. .true.) then
          tempwindow(ii,jj) = 1
         else
          tempwindow(ii,jj) = 0
         endif
       enddo
      enddo

      call fwrite(MPI_COMM_WORLD,tempwindow,ng,nv,"lwindow",fid_dat)
      deallocate(tempwindow)
      endif

      call fclose(fid_dat)
!      endif
      end

