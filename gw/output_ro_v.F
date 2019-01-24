      subroutine output_ro_v
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
	integer :: isort,i_len,nn,isym,mt,mt1,ir,i
	real*8 :: sqpi4
	real*8, allocatable :: tmp(:,:)
	if(.not.maswrk) return
	sqpi4=sqrt(4.d0*pi)
	allocate(tmp(0:maxnrad,5))
	i_len=len_trim(allfile)
	open(3,file=allfile(1:i_len)//'_Electronic_density')
	do isort=1,nsort
	  write(3,*)' Isort ',isort
	  nn=min(5,nsym(isort))
	  do isym=1,nn
	    if(nspin==1) then
            mt=indmt(isym,isort,1)
            do ir=0,nrad(isort)
              tmp(ir,isym)=ro(mt+ir)
            enddo
	    else if(nspin==2) then
            mt=indmt(isym,isort,1)
            mt1=indmt(isym,isort,2)
            do ir=0,nrad(isort)
              tmp(ir,isym)=ro(mt+ir)+ro(mt1+ir)
            enddo
          endif
        enddo
        do ir=0,nrad(isort)
          write(3,'(6(1x,e12.5))')r(ir,isort),(tmp(ir,i),i=1,nn)
        enddo
	enddo
	close(3)
	open(3,file=allfile(1:i_len)//'_Hartree_potential')
	do isort=1,nsort
	  write(3,*)' Isort ',isort
	  nn=min(5,nsym(isort))
	  do isym=1,nn
          mt=indmt(isym,isort,1)
          if(isym/=1) then
            do ir=0,nrad(isort)
              tmp(ir,isym)=v_mt_h(mt+ir)
            enddo
          else if(isym==1) then
            tmp(0,isym)=v_mt_h(mt)
            do ir=1,nrad(isort)
              tmp(ir,isym)=v_mt_h(mt+ir)
     &                    -2.d0*z(isort)/r(ir,isort)*sqpi4
            enddo
          endif
        enddo
        do ir=0,nrad(isort)
          write(3,'(6(1x,e12.5))')r(ir,isort),(tmp(ir,i),i=1,nn)
        enddo
	enddo
	close(3)
	open(3,file=allfile(1:i_len)//'_Exchange_Correlation_Potential')
	do isort=1,nsort
	  write(3,*)' Isort ',isort
	  nn=min(5,nsym(isort))
	  do isym=1,nn
          mt=indmt(isym,isort,1)
          do ir=0,nrad(isort)
            tmp(ir,isym)=v_mt_xc(mt+ir)
          enddo
        enddo
        do ir=0,nrad(isort)
          write(3,'(6(1x,e12.5))')r(ir,isort),(tmp(ir,i),i=1,nn)
        enddo
	enddo
	close(3)
	deallocate(tmp)
      end