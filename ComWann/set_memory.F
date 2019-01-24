      subroutine set_memory
      use comwann_mod

      implicit none
      integer :: iatom,isort,
!     $  tmparr(3),
!     $  k0,k1,k2,k3,
!     $  in_out(2,2,2,ndiv(1),ndiv(2),ndiv(3)),
     $  j,j1,j0
      character*20 :: ch_atom
!      double precision ::
!     $  tempvec3(3),weig_rc_tmp(ndiv(1),ndiv(2),ndiv(3)),dist_rc_tmp,
!     $  maxdev,gtild(3),v(3)
!      character*300 :: ch,tempchar
!      double precision, allocatable :: tmppnt(:,:),dist_r(:)
      integer, external :: gbs_number

      allocate(ffsmt(2,2,maxntle,maxntle,nrel*maxb+1,nsort,nspin))
      allocate(gfun(maxwf,nspin))
      gfun=0.0d0
      allocate(gfund(maxwf,nspin))
      gfund=0.0d0
      allocate(e_bnd(nbndf,npnt,nspin_0))
      
      allocate(n_bnd(npnt,nspin_0))

      allocate(e_qp(nbndf,npnt,nspin_0))
      e_qp=0.0d0
      allocate(ev_bnd(nbasmpw,nbndf,ndim_k(me+1),nspin_0))
      ev_bnd=0.0d0
      allocate(z_bnd(nfun,nbndf,ndim_k(me+1),nspin_0))
      z_bnd=0.0d0
      allocate(ev_qp(nbasmpw,nbndf,ndim_k(me+1),nspin_0))
      ev_qp=0.0d0
      allocate(z_qp(nfun,nbndf,ndim_k(me+1),nspin_0))
      z_qp=0.0d0

!     for wannier input
      allocate(atom_symbols(natom))
      do iatom=1, natom
        isort=is(iatom)
        if (iatom .lt. 10) then
          write(ch_atom,'(i1)') iatom
          atom_symbols(iatom)=trim(txtel(isort))//'_'//
     $      trim(ch_atom)//repeat(" ",16)
          
        elseif (iatom .lt. 100) then
          write(ch_atom,'(i2)') iatom
          atom_symbols(iatom)=trim(txtel(isort))//'_'//
     $      trim(ch_atom)//repeat(" ",15)
        else
          write(ch_atom,'(i3)') iatom
          atom_symbols(iatom)=trim(txtel(isort))//'_'//
     $      trim(ch_atom)//repeat(" ",14)
        endif
      enddo

      allocate(kpt_latt(3,nqdiv))
      
      kpt_latt=0.0d0
!      indtemp=0
      kpt_latt=matmul(transpose(rbas),pnt)



      end
