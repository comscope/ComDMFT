      subroutine sym_sig_r_from_nonsym(f_in,f_out)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	complex*16, intent(in) :: f_in(nfun,nfun,nqdiv)
	complex*16, intent(out) :: f_out(nfun,nfun,nqdiv)
	integer :: j,ig,iatom,isort,jatom,jsort,i,ip0,jp0,iatom0,jatom0,
     &	         ip00,jp00,ir,ir0
      real*8 :: v(3),v1(3)
      complex*16, allocatable :: tmp(:,:),tmp1(:)
      allocate(tmp(maxel,maxel))
      allocate(tmp1(maxel))
	f_out=(0.d0,0.d0)
	do jatom=1,natom
	  jsort=is(jatom)
	  jp0=io_lem(jatom)-1
	  do iatom=1,natom
	    isort=is(iatom)
	    ip0=io_lem(iatom)-1
	    do ir=1,nqdiv
	      do ig=1,ngroup
	        v=rr_coarse(:,ir)+tshift(:,iatom,ig)-tshift(:,jatom,ig)
              call rotate(v(1),v(2),v(3),v1(1),v1(2),v1(3),u(2,ig),1)
	        call zone1_number(v1,gbas,ndiv,ir0)
	        iatom0=ip(iatom,ig)
	        jatom0=ip(jatom,ig)
	        jp00=io_lem(jatom0)-1
	        ip00=io_lem(iatom0)-1
	        do i=1,lfunm(isort)
	          call sym_left_tt(ig,
     &                           f_in(ip00+i,jp00+1:jp00+lfunm(jsort),
     &                                ir0),tmp(i,1:lfunm(jsort)),2,
     &                           lfunm(jsort),jsort)
              enddo
	        do j=1,lfunm(jsort)
	          call sym_left_tt(ig,tmp(1:lfunm(isort),j),
     &                           tmp1(1:lfunm(isort)),2,lfunm(isort),
     &                           isort)
                f_out(ip0+1:ip0+lfunm(isort),jp0+j,ir)=
     &            f_out(ip0+1:ip0+lfunm(isort),jp0+j,ir)
     &           +tmp1(1:lfunm(isort))
              enddo
	      enddo    !! over ig
	    enddo   !! over ir
	  enddo  !! over iatom
	enddo   !! over jatom
	deallocate(tmp,tmp1)
	f_out=f_out/ngroup
      end