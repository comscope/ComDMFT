      subroutine sym_mm_rr_boson(f_in,f_out)
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	real*8, intent(in) :: f_in(n_pbmt,n_pbmt,nqdiv)
	real*8, intent(out) :: f_out(n_pbmt,n_pbmt,nqdiv)
	integer :: j,ig,iatom,isort,jatom,jsort,i,ip0,jp0,iatom0,jatom0,
     &	         ip00,jp00,ir,ir0
      real*8 :: v(3),v1(3)
      real*8, allocatable :: tmp(:,:),tmp1(:)
      allocate(tmp(n_pbmtm,n_pbmtm))
      allocate(tmp1(n_pbmtm))
	f_out=0.d0
	do jatom=1,natom
	  jsort=is(jatom)
	  jp0=iopb(jatom)-1
	  do iatom=1,natom
	    isort=is(iatom)
	    ip0=iopb(iatom)-1
	    do ir=1,nqdiv
	      do ig=1,ngroup
	        v=rr_coarse(:,ir)+tshift(:,iatom,ig)-tshift(:,jatom,ig)
              call rotate(v(1),v(2),v(3),v1(1),v1(2),v1(3),u(2,ig),1)
	        call zone1_number(v1,gbas,ndiv,ir0)
	        iatom0=ip(iatom,ig)
	        jatom0=ip(jatom,ig)
	        jp00=iopb(jatom0)-1
	        ip00=iopb(iatom0)-1
	        do i=1,n_pbmt0(isort)
	          call sym_bos_tt_r(ig,
     &                          f_in(ip00+i,jp00+1:jp00+n_pbmt0(jsort),
     &                               ir0),tmp(i,1:n_pbmt0(jsort)),
     &                          n_pbmt0(jsort),jsort)
              enddo
	        do j=1,n_pbmt0(jsort)
	          call sym_bos_tt_r(ig,tmp(1:n_pbmt0(isort),j),
     &                          tmp1(1:n_pbmt0(isort)),n_pbmt0(isort),
     &                          isort)
                f_out(ip0+1:ip0+n_pbmt0(isort),jp0+j,ir)=
     &            f_out(ip0+1:ip0+n_pbmt0(isort),jp0+j,ir)
     &           +tmp1(1:n_pbmt0(isort))
              enddo
	      enddo    !! over ig
	    enddo   !! over ir
	  enddo  !! over iatom
	enddo   !! over jatom
	deallocate(tmp,tmp1)
	f_out=f_out/ngroup
      end