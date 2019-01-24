      subroutine output_eigen_states(ch0,key,ispin)
c	key = 0  from DFT
c	key = 1  from HF
c	key = 2  from X
c	key = 3  from QP
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      character*4, intent(in) :: ch0
	integer, intent(in) :: key,ispin
	integer :: i_len,k,ibnd,n,i,ind_k
	real*8 :: e,sum_band_orb_char
	real*8, allocatable :: pdos(:,:)
	allocate(pdos(nbndf,npnt))
	pdos=0.d0
      do ind_k=1,ndim3_k(me3_k+1)
	  k=n3_mpi_k(me3_k+1)+ind_k
	  n=n_bnd(k,ispin)
	  do ibnd=1,n
	    if(key<=2) then
	      pdos(ibnd,k)=sum_band_orb_char(ibnd,ind_k,ispin,0)
	    else if(key==3) then
	      pdos(ibnd,k)=sum_band_orb_char(ibnd,ind_k,ispin,1)
	    endif
	  enddo
	enddo
	if(nproc_k/=1) call dgop(pdos,nbndf*npnt,'  +',comm_pnt)
c ---------------------------------------------------------	
	if(maswrk) then
	  i_len=len_trim(allfile)
	  if(key==0) open(3,file=allfile(1:i_len)//'_dft.eig')
	  if(key==1) open(3,file=allfile(1:i_len)//'_hf.eig')
	  if(key==2) open(3,file=allfile(1:i_len)//ch0//'_x.eig')
	  if(key==3) open(3,file=allfile(1:i_len)//ch0//'_qp.eig')
	  if(ispin==1) write(3,*)
     &	     ' Eigen states relative to Chemical Potential (eV)'
	  if(ispin==2) call goto_end(3)
c         write(3,*)'npnt',npnt,'pnt',pnt,'sym_label',sym_label,
c     &    'n_bnd',n_bnd,'e_bnd',e_bnd,'pdos',pdos,'nbndf',nbndf,
c     &    'comm_pnt',comm_pnt
        do k=1,npnt
	    n=n_bnd(k,ispin)
	    write(3,'(a9,i1,a5,i4,4x,a1,3f9.5)')' Ispin = ',ispin,
     &	    ' k = ',k,sym_label(k),(pnt(i,k),i=1,3)
	    do ibnd=1,n
	      if(key<=2) then
	        e=e_bnd(ibnd,k,ispin)
	      else if(key==3) then
	        e=e_qp(ibnd,k,ispin)
	      endif
	      if(key<=2) then
	        write(3,10)ibnd,(e-chem_pot)*evolt/2,pdos(ibnd,k)
	      else if(key==3) then
	        write(3,10)ibnd,(e-chem_pot_qp)*evolt/2,pdos(ibnd,k)
	      endif
	    enddo
	  enddo
	  close(3)
	endif
	deallocate(pdos)
10	format(i5,1x,e15.8,1x,f6.3)
11	format(i5,1x,e15.8,1x,a1,1x,f6.3)
      end
