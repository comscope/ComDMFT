!******************************************************************************
! Copyright c 2013, The Ames Laboratory, Iowa State University, and Rutgers
! University*.  All rights reserved.
!
! This software was authored by Yongxin Yao, Nicola Lanata*, Gabriel Kotliar*,
! Cai-Zhuang Wang, and Kai-Ming Ho, at The Ames Laboratory and
! Rutgers University and was supported by the U.S.
! Department of Energy (DOE), Office of Science,
! Basic Energy Sciences, Materials Science and Engineering Division.
! The Ames Laboratory is operated by Iowa State University for DOE
! under U.S. Government contract DE-AC02-07CH11358.
! The U.S. Government has the rights to use, reproduce, and
! distribute this software.
! NEITHER THE GOVERNMENT, THE AMES LABORATORY, IOWA STATE UNIVERSITY,
! NOR RUTGERS UNIVERSITY MAKES ANY WARRANTY,
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
! If software is modified to produce derivative works,
! such modified software should be clearly marked,
! so as to not confuse it with the version available from
! The Ames Laboratory and Rutgers University.
!
! Additionally, redistribution and use in source and binary forms,
! with or without modification,
! are permitted provided that the following conditions are met:
!
!     Redistribution of source code must retain the above copyright notice,
!     this list of conditions, and the following disclaimer.
!
!     Redistribution in binary form must reproduce the above copyright notice,
!     this list of conditions, and the following disclaimer
!     in the documentation and/or other materials provided with distribution.
!
!     Neither the name of The Ames Laboratory, Iowa State University,
!     Rutgers University, the U.S. Government, nor the names of
!     its contributors may be used to endorse or promote products derived
!     from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE AMES LABORATORY, IOWA STATE UNIVERSITY,
! RUTGERS UNIVERSITY, AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
! THE IMPLIED WARRANTIES OF MERCHANTABILITY
! AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
! IN NO EVENT SHALL THE GOVERNMENT, THE AMES LABORATORY,
! IOWA STATE UNIVERSITY, RUTGERS UNIVERSITY, OR CONTRIBUTORS BE LIABLE
! FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
! HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
! WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
! OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
! EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!******************************************************************************

module bandstru
    use gprec
    use gmpi
    use gutil
    use gconstant, only: d0,d1,z0,z1,zi,iu_kgen
    use ghdf5_base
    use warehouse
    use gtime
    implicit none

    type band_stru
        !< iso=2 => soc
        integer ispin_in,nspin_in
        integer n_frozen
        real(q) ehybrd,eband(2),ets2(2)
        !< ne(3). 1: total number of bands, 2-3: correlated bands interval
        integer,pointer :: ne(:,:,:)
        real(q),pointer :: ek(:,:,:) ! kohn-sham / gutz eigen values
        real(q),allocatable :: normso(:,:,:) ! list of spin-up/dn weights.
        complex(q),pointer :: r(:,:,:),la1(:,:,:),nc_phy(:,:,:),nks(:,:,:), &
                &d0(:,:,:),nrl(:,:,:) ! for onsite local 1pdm
        complex(q),pointer :: psi0_b(:,:,:) ! <bare psi0 | basis orbital>
        complex(q),pointer :: vk (:,:,:,:) ! eigen-vector
        complex(q),pointer :: hk0(:,:,:,:) ! lda bare dispersion
        integer nmax ! nmax: miximal number of bands
        integer nmaxin !< maximal number of bands inside the enrgy window.
        real(q),pointer :: ferwe(:,:,:),ferwer(:,:,:) 
        !< fermi-weight, convention including kpt%wt and spin degeneracy
        real(q) :: eflda,nelet,nelel,nelec,mup_dn=0._q
        real(q) :: ef=0._q,nelet_frozen=0._q
        real(q) ebmin,ebmax,ebwidth
        complex(q),pointer::cwt(:,:,:)   ! occupation matrix for \rho
    end type band_stru
      
    type k_points
        integer::dim=-1,diml=-1,ismear,icor,ensemble=0
        real(q),pointer :: wt(:)=>null()
        real(q) twt,delta,cordeg
        character::file_name*512
    end type k_points
      
    type sym_info
        integer::nop,mode=0,ie
    end type sym_info
      
    type (band_stru) :: bnd
    type (sym_info)  :: sym
    type (k_points)  :: kpt
    

    contains


    subroutine set_bnd_info(io)
    integer,intent(in)::io

    integer iso,isp
    logical lexist

    bnd%nelet_frozen=sum(wh%fz(:)%nelect)
    bnd%n_frozen=sum(wh%fz(:)%nsorb)

    call gh5_open_r('GPARAMBANDS.h5',f_id)
    call gh5_read(iso,'/iso',f_id)
    call gh5_read(bnd%ispin_in,'/ispin',f_id)
    call gh5_read(kpt%dim,'/kptdim',f_id)
    call gh5_read(bnd%nmax,'/nbmax',f_id)
    wh%ispin_in=bnd%ispin_in
    if(io>0)then
        write(io,'(" ispin_in=",i2," iso=",i2, &
                &" ispin=",i2," ispo=",i2)')bnd%ispin_in, &
                &wh%iso,wh%ispin,wh%ispo
    endif
    if(iso==2.and.wh%ispin==2)then
        allocate(bnd%normso(bnd%nmax,2,kpt%dim))
        call gh5_read(bnd%normso,bnd%nmax,2,kpt%dim,'/normso',f_id)
    endif
    ! modify ispin_in for iso=2
    bnd%nspin_in=max(1, bnd%ispin_in/iso)
    allocate(bnd%ne(3,kpt%dim,bnd%nspin_in))
    call h5lexists_f(f_id,'/NE_LIST_SPIN1', &
            &lexist,gh5_err)
    if(lexist)then
        ! upper-case name: fortran convention
        do isp=1,bnd%nspin_in
            call gh5_read(bnd%ne(:,:,isp),3,kpt%dim,'/NE_LIST_SPIN'//&
                    &trim(int_to_str(isp)),f_id)
        enddo
    else
        bnd%ne(1,:,:)=bnd%nmax
        bnd%ne(2,:,:)=1
        bnd%ne(3,:,:)=bnd%nmax
    endif

    call h5lexists_f(f_id,'/ensemble',lexist,gh5_err)
    if(lexist)then
        call gh5_read(kpt%ensemble,'/ensemble',f_id)
    endif

    allocate(kpt%wt(kpt%dim))
    call gh5_read(kpt%wt,kpt%dim,'/kptwt',f_id)
    call gh5_read(kpt%ismear,'/ismear',f_id)
    call gh5_read(kpt%delta,'/delta',f_id)
    call gh5_read(bnd%nelet,'/nelectron',f_id)
    call h5lexists_f(f_id,'/kptfname',lexist,gh5_err)
    if(lexist)then
        call gh5_read(kpt%file_name,512,'/kptfname',f_id)
    endif
    call gh5_read(sym%nop,'/symnop',f_id)
    call gh5_read(sym%ie,'/symie',f_id)
    call gh5_read_wh_h1e(bnd%nspin_in)
    call gh5_close(f_id)
    call rotate_h1e_list()
    call calc_herm_matrices_pp(wh%h1e,'h1e',wh%hm,.false.,io,-1)
    kpt%twt=sum(kpt%wt)
    kpt%wt=kpt%wt/kpt%twt

    ! Consistence check
    if(wh%iso /= iso)then
        stop ' error in set_bnd_info: wh%iso /= bnd%iso!'
    endif

    bnd%nmaxin=maxval(bnd%ne(3,:,:)-bnd%ne(2,:,:)+1)
    if(io>0)then
        write(io,'(" min/max(bnd%ne(1,:,:))=",2i8)')minval(bnd%ne(1,:,:)), &
                &maxval(bnd%ne(1,:,:))
        write(io,'(" min/max(bnd%ne(2,:,:))=",2i8)')minval(bnd%ne(2,:,:)), &
                &maxval(bnd%ne(2,:,:))
        write(io,'(" min/max(bnd%ne(3,:,:))=",2i8)')minval(bnd%ne(3,:,:)), &
                &maxval(bnd%ne(3,:,:))
    endif

    bnd%nelec=0
    do isp=1,bnd%nspin_in
        bnd%nelec=bnd%nelec+sum((bnd%ne(2,:,isp)-1)*kpt%wt)*wh%rspo
    enddo
    bnd%nelel=bnd%nelet-bnd%nelec
    if(io>0)then
        write(io,'(" bnd%nmax=",i4)')bnd%nmax
        write(io,'(" valence electrons: total=",f8.1)')bnd%nelet
        write(io,'("         correlated block=",f8.1)')bnd%nelel
    endif
    call set_kpt_icor(io)
    return
      
    end subroutine set_bnd_info
    

    subroutine set_kpt_diml(nkpl)
    integer,intent(in)::nkpl

    kpt%diml=nkpl
    return
    
    end subroutine set_kpt_diml


    subroutine alloc_bnd()
    integer i
    
    if(kpt%diml<0)then
        write(0,'(" error in alloc_bnd: kpt%diml not set yet!")')
        stop
    endif
    allocate(bnd%ferwe (bnd%nmax,kpt%dim,wh%nspin))
    allocate(bnd%ferwer(bnd%nmax,kpt%dim,wh%nspin))
    allocate(bnd%ek(bnd%nmax,kpt%dim,wh%nspin))
    bnd%ek=0
    allocate(bnd%hk0(bnd%nmaxin**2,sym%nop,kpt%diml,bnd%nspin_in))
    allocate(bnd%vk (bnd%nmaxin**2,sym%nop,kpt%diml, &
            &wh%nspin))
    allocate(bnd%psi0_b(bnd%nmaxin**2,kpt%diml,bnd%nspin_in))
    bnd%psi0_b=0
    allocate(bnd%r(wh%nasotot,wh%nasotot,wh%nspin))
    bnd%r=0
    forall(i=1:wh%nasotot) bnd%r(i,i,:)=1._q
    allocate(bnd%d0 (wh%nasotot,wh%nasotot,wh%nspin))
    allocate(bnd%la1(wh%nasotot,wh%nasotot,wh%nspin))
    bnd%la1=0
    allocate(bnd%nc_phy(wh%nasotot,wh%nasotot,wh%nspin))
    allocate(bnd%nks(wh%nasotot,wh%nasotot,wh%nspin))
    allocate(bnd%nrl(wh%nasotot,wh%nasotot,wh%nspin))
    return
      
    end subroutine alloc_bnd
    

    subroutine read_bare_hamiltonian()
    integer ivec,ikp,ikpl,iks,nbtot,isym,nbands,isp
    complex(q),pointer::p_2d(:,:)
      
    ikpl=0
    do ivec=1,gp%nvec
        call gh5_open_r(file_name(gp%kvec(ivec,1),'BAREHAM'), f_id)
        do iks=1,gp%kvec(ivec,2)
            ikp=gp%kvec(ivec,3)+iks
            if(gp%ipar==1)then ! openmp
                ikpl=ikp
            else
                ikpl=ikpl+1
            endif
            do isp=1,bnd%nspin_in
                nbtot=bnd%ne(1,ikp,isp)
                call gh5_read(bnd%ek(1:nbtot,ikp,isp),nbtot,'/IKP_'// &
                        &trim(int_to_str(ikp))//'/ek0_spin'//&
                        &trim(int_to_str(isp)),f_id)
            enddo
            if(wh%iso==1.and.bnd%nspin_in==1.and.wh%nspin==2)then
                bnd%ek(1:nbtot,ikp,2)=bnd%ek(1:nbtot,ikp,1)
            endif
            do isp=1,bnd%nspin_in
                nbands=bnd%ne(3,ikp,isp)-bnd%ne(2,ikp,isp)+1
                p_2d(1:nbands,1:nbands)=>bnd%psi0_b(1:nbands**2,ikpl,isp)
                call gh5_read(p_2d,nbands,nbands, &
                        &'/IKP_'//trim(int_to_str(ikp))// &
                        &'/T_PSIK0_TO_HK0_BASIS_SPIN'// &
                        &trim(int_to_str(isp)),f_id)
                do isym=1,sym%nop
                    p_2d(1:nbands,1:nbands)=>bnd%hk0(1:nbands**2,isym, &
                            &ikpl,isp)
                    call gh5_read(p_2d,nbands,nbands, &
                            &'/IKP_'//trim(int_to_str(ikp))// &
                            &'/ISYM_'//trim(int_to_str(isym))//&
                            &'/HK0_SPIN'//trim(int_to_str(isp)),f_id)
                enddo
            enddo
        enddo
        call gh5_close(f_id)
    enddo
    nullify(p_2d)
#ifdef mpi_mode
    call dsum_all_mpi(bnd%ek,bnd%nmax*kpt%dim*wh%nspin)
#endif
    return
      
    end subroutine read_bare_hamiltonian
    

    ! rotate bare hamiltonian to symmetry-adapted basis.
    subroutine rotate_bare_hamiltonian()
    integer ivec,ikp,ikpl,iks,isym,nbands,isp
    complex(q),pointer::p_2d(:,:)

    ! check if necessary
    if(.not.associated(wh%db2sab))return

    ikpl=0
    do ivec=1,gp%nvec
        do iks=1,gp%kvec(ivec,2)
            ikp=gp%kvec(ivec,3)+iks
            if(gp%ipar==1)then
                ikpl=ikp
            else
                ikpl=ikpl+1
            endif
            do isp=1,bnd%nspin_in
                nbands=bnd%ne(3,ikp,isp)-bnd%ne(2,ikp,isp)+1
                p_2d(1:nbands,1:nbands)=>bnd%psi0_b(1:nbands**2,ikpl,isp)
                call site_wise_rtrans(p_2d,nbands,1)
                do isym=1,sym%nop
                    p_2d(1:nbands,1:nbands)=>bnd%hk0(1:nbands**2, &
                            &isym,ikpl,isp)
                    call site_wise_rtrans(p_2d,nbands,1)
                    p_2d=transpose(p_2d)
                    call site_wise_rtrans(p_2d,nbands,-1)
                    p_2d=transpose(p_2d)
                enddo
            enddo
        enddo
    enddo
    nullify(p_2d)
    return
      
    end subroutine rotate_bare_hamiltonian


    subroutine site_wise_rtrans(a,n,mode)
    integer,intent(in)::n,mode
    complex(q),intent(inout)::a(n,n)

    integer nbase,i,naso
    complex(q),pointer::p_trans(:,:)
    complex(q),target::utrans(wh%nasomax**2)

    nbase=0
    do i=1,wh%num_imp
        naso=wh%co(i)%dimso
        p_trans(1:naso,1:naso) => utrans(1:naso**2)
        if(wh%iso==1)then
            p_trans=wh%co(i)%db2sab(1:naso,1::2)
        else
            p_trans = wh%co(i)%db2sab
        endif
        if(mode<0)then
            p_trans=conjg(p_trans)
        endif
        call anmxbmm("n",a(:,nbase+1:nbase+naso), &
                &p_trans,n,naso)
        nbase=nbase+naso
    enddo
    nullify(p_trans)
    return

    end subroutine site_wise_rtrans


    subroutine rotate_h1e_list()
    integer i,na2,na
    complex(q),pointer::p_h1e(:,:)
    complex(q),target::h1e(wh%na2max**2)

    ! check if necessary
    if(.not.associated(wh%db2sab))return

    do i=1,wh%num_imp
        na2=wh%co(i)%dim2; na=wh%co(i)%dim
        if(wh%iso==1.and.(.not.associated(wh%db2sab)))then
            ! {{orbs}_up, {orbs}_dn} -> {{orb_up, orb_dn}}
            p_h1e(1:na2,1:na2) => h1e(1:na2**2)
            p_h1e=wh%co(i)%h1e
            wh%co(i)%h1e(1::2,1::2)=p_h1e(1:na,1:na)
            wh%co(i)%h1e(2::2,2::2)=p_h1e(1+na:,1+na:)
            wh%co(i)%h1e(1::2,2::2)=0
            wh%co(i)%h1e(2::2,1::2)=0
        else
            call uhau(wh%co(i)%h1e,wh%co(i)%db2sab,na2,na2)
        endif
    enddo
    nullify(p_h1e)
    return

    end subroutine rotate_h1e_list


    subroutine rm_h1e_from_bare_hamiltonian()
    integer i,ivec,isym,ikp,iks,ikpl,isp
    integer naso,nstep,nbase,nasot,nbands
    complex(q) h1e(wh%nasotot,wh%nasotot,bnd%nspin_in)
    complex(q),pointer::p_hk0(:,:)

    h1e=0
    nbase=1
    do i=1,wh%num_imp
        naso=wh%co(i)%dimso
        nstep=wh%co(i)%dim2/naso
        do isp=1,bnd%nspin_in
            h1e(nbase:nbase+naso-1,nbase:nbase+naso-1,isp)= &
                    &wh%co(i)%h1e(isp::nstep,isp::nstep)
        enddo
        nbase=nbase+naso
    enddo
    nasot=wh%nasotot

    ikpl=0
    do ivec=1,gp%nvec
        do iks=1,gp%kvec(ivec,2)
            ikp=gp%kvec(ivec,3)+iks
            if(gp%ipar==1)then
                ikpl=ikp
            else
                ikpl=ikpl+1
            endif
            do isp=1,bnd%nspin_in
                nbands=bnd%ne(3,ikp,isp)-bnd%ne(2,ikp,isp)+1
                do isym=1,sym%nop
                    p_hk0(1:nbands,1:nbands)=>bnd%hk0(1:nbands**2, &
                            &isym,ikpl,isp)
                    p_hk0(1:nasot,1:nasot)=p_hk0(1:nasot,1:nasot)-h1e(:,:,isp)
                enddo
            enddo
        enddo
    enddo
    nullify(p_hk0)
    return

    end subroutine rm_h1e_from_bare_hamiltonian

    !*************************************************************************  
    !< calc band energy
    !*************************************************************************
    subroutine calc_ebnd()
    integer ikp,ikpl,iks,ivec,isp,ispp,nemin,nemax,nbands,i,j
    complex(q),pointer::p_a(:,:),p_v(:,:),p_psi0b(:,:)
    complex(q),target::a(bnd%nmaxin**2)
     
    if(wh%ispin/=2.or.wh%iso/=2)then
        bnd%eband=bnd%ets2
        do ikp=1,kpt%dim
            do isp=1,wh%nspin
                ispp=min(isp,bnd%nspin_in)
                bnd%eband(isp)=bnd%eband(isp)+ &
                        &sum(bnd%ek(:,ikp,isp)*bnd%ferwe(:,ikp,isp))
            enddo
        enddo
    else
        bnd%eband=0
        ikpl=0
        do ivec=1,gp%nvec
            do iks=1,gp%kvec(ivec,2)
                ikp=gp%kvec(ivec,3)+iks
                if(gp%ipar==1)then
                    ikpl=ikp
                else
                    ikpl=ikpl+1
                endif
                do isp=1,2
                    nemin=bnd%ne(2,ikp,1)
                    nemax=bnd%ne(3,ikp,1)
                    if(nemin>1)then
                        bnd%eband(isp)=bnd%eband(isp)+sum(&
                                &bnd%ek(:nemin-1,ikp,1)* &
                                &bnd%ferwe(:nemin-1,ikp,1)* &
                                &bnd%normso(:nemin-1,isp,ikp))
                    endif
                    nbands=nemax-nemin+1
                    p_a(1:nbands,1:nbands)=>a(1:nbands**2)
                    p_psi0b(1:nbands,1:nbands)=>bnd%psi0_b(1:nbands**2,ikpl,1)
                    p_v(1:nbands,1:nbands)=>bnd%vk(1:nbands**2,sym%ie,ikpl,1)
                    call zgemm('n','n',nbands,nbands,nbands,z1,p_psi0b, &
                            &nbands,p_v,nbands,z0,p_a,nbands)
                    do i=1,nbands; do j=1,nbands
                        bnd%eband(isp)=bnd%eband(isp)+ &
                                &bnd%ek(nemin-1+i,ikp,1)*p_a(i,j)* &
                                &conjg(p_a(i,j))*bnd%ferwe(nemin-1+i,ikp,1)* &
                                &bnd%normso(nemin-1+j,isp,ikp)
                    enddo; enddo
                enddo
            enddo
        enddo
#ifdef mpi_mode
        call sum_all_mpi(bnd%eband,2)
#endif
        bnd%eband=bnd%eband+bnd%ets2/2
        nullify(p_a,p_v,p_psi0b)
    endif
    return
      
    end subroutine calc_ebnd


    subroutine calc_e_hybrd()

    bnd%ehybrd=real(sum(wh%r*(wh%d0)),q)
    return

    end subroutine calc_e_hybrd


    !*************************************************************************
    !< check magnetic moment for the case of spin-polarized calculations.
    !*************************************************************************
    subroutine chk_mag_moment(io)
    integer,intent(in)::io

    integer isp,ikp
    real(q) nel(2)

    if(io<0)return
    if(wh%nspin==1.or.wh%iso==2)return
    nel=0
    do isp=1,wh%nspin
        do ikp=1,kpt%dim
            nel(isp) = nel(isp) + sum(bnd%ferwe(:,ikp,isp))*kpt%wt(ikp)
        enddo
    enddo
    write(io,'(" total (quasi-particle) magnetic moment = ", f12.4)') &
            &nel(2)- nel(1)
    return

    end subroutine chk_mag_moment


    !*************************************************************************
    !> f_n <psi_n|a> r_{a,A} r^+_{B,b} <b|psi_n>
    !! =r^+_{B,b} <b|psi_n> f_n <psi_n|a> r_{a,A}
    !*************************************************************************
    subroutine calc_nabr_1k(nabr,vk,ferwe,nbands,wtk,wtk0,isp)
    integer,intent(in)       :: nbands,isp
    real(q),intent(in)       :: wtk,wtk0,ferwe(nbands)
    complex(q),intent(in)    :: vk(nbands,nbands)
    complex(q),intent(inout) :: nabr(nbands,nbands)

    integer ia,nasot
    complex(q) r(wh%nasotot,wh%nasotot)

    !< r^+_{B,b}<b|psi>f<psi|a>
    call calc_nabr1_1k(nabr,vk,ferwe,nbands,wtk,isp)
    nasot=wh%nasotot
    r=bnd%r(:,:,isp)
    !< r^+_{B,b}<b|psi>f<psi|a>r(a,A)
    call anmxbmm('n',nabr(:,1:nasot),r,nbands,nasot)
    nabr=transpose(nabr) ! to {A,B}
    nabr(1:nasot,1:nasot)=nabr(1:nasot,1:nasot)+ &
            &(-bnd%nrl(:,:,isp)+bnd%nc_phy(:,:,isp))*wtk0
    return

    end subroutine calc_nabr_1k


    !*************************************************************************
    !> r^+ applied to right side only, for d
    !! r^+ <b|psi>f<psi|a>. 
    !*************************************************************************
    subroutine calc_nabr1_1k(nabr,vk,ferwe,nbands,wtk,isp)
    integer,intent(in)       :: nbands,isp
    real(q),intent(in)       :: wtk,ferwe(nbands)
    complex(q),intent(in)    :: vk(nbands,nbands)
    complex(q),intent(out) :: nabr(nbands,nbands)

    integer nasot
    complex(q) zr(wh%nasotot,wh%nasotot),nabr_sub(wh%nasotot,nbands)

    call calc_nab_1k(nabr,vk,ferwe,nbands,nbands,wtk)
    nasot=wh%nasotot
    zr=bnd%r(:,:,isp)
    nabr_sub=nabr(1:nasot,:)
    ! r^+ <b|psi>f<psi|a>
    call annxb('c',zr,nabr_sub,nasot,nbands,nabr(1:nasot,:))
    return

    end subroutine calc_nabr1_1k
     

    !*************************************************************************
    !D0^{t}_{A,a} 
    != f_n <psi_n|a> h_{A,B} r^+_{B,b} <b|psi_n>
    != h_{A,B} r^+_{B,b} <b|psi_n> f_n <psi_n|a>
    !*************************************************************************
    subroutine add_da_1k(hk0,vk,ferwe,nbands,wtk,isp)
    integer,intent(in)    :: nbands,isp
    real(q),intent(in)    :: wtk,ferwe(nbands)
    complex(q),intent(in) :: hk0(nbands,nbands),vk(nbands,nbands)

    complex(q) nabr(nbands,nbands),zd(wh%nasotot,wh%nasotot), &
            &hk0_(wh%nasotot,nbands)

    call calc_nabr1_1k(nabr,vk,ferwe,nbands,wtk,isp) ! r^+ <b|psi>f<psi|a>
    hk0_=hk0(1:wh%nasotot,:)
    call zgemm('n','n',wh%nasotot,wh%nasotot,nbands,z1, &
            &hk0_,wh%nasotot,nabr(:,1:wh%nasotot),nbands,z0,zd,wh%nasotot)
    ! D0 can be real or complex type .
    bnd%d0(:,:,isp)=bnd%d0(:,:,isp)+zd
    return

    end subroutine add_da_1k


    !< calculate n_{ab} at a single k-point. 
    subroutine calc_nab_1k(nab,vk,ferwe,nbasis,nbands,wtk)
    integer,intent(in)     :: nbasis,nbands
    real(q),intent(in)     :: wtk,ferwe(nbands)
    complex(q),intent(in)  :: vk(nbasis,nbands)
    complex(q),intent(out) :: nab(nbasis,nbasis)

    integer ib
    complex(q) vf(nbasis,nbands)

    do ib=1,nbands
        vf(:,ib)=vk(:,ib)*ferwe(ib)*wtk
    enddo
    call zgemm('n','c',nbasis,nbasis,nbands,z1,vf,nbasis,vk,nbasis,z0, &
            &nab,nbasis) ! <b|psi>f<psi|a>
    return

    end subroutine calc_nab_1k


    subroutine calc_da0()
    integer isym,ivec,iks,ikp,ikpl,isp,ispp,nbands,nemin,nemax
    real(q) wtk
    complex(q),pointer :: p_vk(:,:),p_hk0(:,:)
    real(q),pointer :: ferwe(:)

    bnd%d0=0
    wtk=1._q/sym%nop/wh%rspo
    ikpl=0
    do ivec=1,gp%nvec
        do iks=1,gp%kvec(ivec,2)
            ikp=gp%kvec(ivec,3)+iks
            if(gp%ipar==1)then
                ikpl=ikp
            else
                ikpl=ikpl+1
            endif
            do isp=1,wh%nspin
                ispp=min(isp,bnd%nspin_in)
                nemin=bnd%ne(2,ikp,ispp)
                nemax=bnd%ne(3,ikp,ispp)
                nbands=nemax-nemin+1
                ferwe=>bnd%ferwe(nemin:nemax,ikp,isp)
                do isym=1,sym%nop
                    p_hk0(1:nbands,1:nbands)=>bnd%hk0(1:nbands**2, &
                            &isym,ikpl,ispp)
                    p_vk(1:nbands,1:nbands)=>bnd%vk(1:nbands**2,isym,ikpl,isp)
                    call add_da_1k(p_hk0,p_vk,ferwe,nbands,wtk,isp)
                enddo                    
            enddo
        enddo
    enddo
#ifdef mpi_mode
    call zsum_all_mpi(bnd%d0,wh%nasotot*wh%nasotot*wh%nspin)
#endif
    nullify(p_vk,p_hk0,ferwe)

    ! back to D0_{a,A}
    do isp=1,wh%nspin
        bnd%d0(:,:,isp)=transpose(bnd%d0(:,:,isp))
    enddo
    call map_wh_bnd_matrix(wh%d0,bnd%d0,.true.)
    return

    end subroutine calc_da0


    subroutine calc_nks()
    integer ikpl,ivec,iks,ikp,isym,isp,ispp,nbands,nbase,&
            &n1,n2,i,naso,nemin,nemax
    real(q) wtk
    real(q),pointer :: p_ferwe(:)
    complex(q),pointer :: p_vk(:,:),p_vk_sub(:,:),p_nks(:,:)
    complex(q),target::nks(wh%nasomax*wh%nasomax)
    complex(q),target::vk_sub(wh%nasomax*bnd%nmaxin)
   
    complex(q) tmp(bnd%nmaxin,bnd%nmaxin)

    wh%nks=0
    wtk=1._q/sym%nop/wh%rspo
    ikpl=0
    do ivec=1,gp%nvec
        do iks=1,gp%kvec(ivec,2)
            ikp=gp%kvec(ivec,3)+iks
            if(gp%ipar==1)then
                ikpl=ikp
            else
                ikpl=ikpl+1
            endif
            do isp=1,wh%nspin
                ispp=min(isp,bnd%nspin_in)
                nemin=bnd%ne(2,ikp,ispp)
                nemax=bnd%ne(3,ikp,ispp)
                nbands=nemax-nemin+1
                p_ferwe=>bnd%ferwe(nemin:nemax,ikp,isp)
                do isym=1,sym%nop
                    p_vk(1:nbands,1:nbands)=>bnd%vk(1:nbands**2,isym,ikpl,isp)
                    nbase=0
                    do i=1,wh%num_imp
                        naso=wh%co(i)%dimso
                        n1=1+(isp-1)*naso; n2=naso*isp
                        p_vk_sub(1:naso,1:nbands)=>vk_sub(1:naso*nbands)
                        p_vk_sub=p_vk(nbase+1:nbase+naso,1:nbands)
                        p_nks(1:naso,1:naso)=>nks(1:naso*naso)
                        call calc_nab_1k(p_nks,p_vk_sub,p_ferwe,naso, &
                                &nbands,wtk)
                        wh%co(i)%nks(n1:n2,n1:n2)=wh%co(i)%nks(n1:n2,n1:n2) +&
                                &transpose(p_nks)
                        nbase=nbase+naso
                    enddo
                enddo
            enddo
        enddo
    enddo
    nullify(p_ferwe,p_vk,p_vk_sub,p_nks)
#ifdef mpi_mode
    call sum_all_mpi(wh%nks,wh%na2112)
#endif
    do i=1,wh%num_imp
        call co_nks_patch_order(wh%co(i),wh%ispo,wh%iso)
    enddo
    return

    end subroutine calc_nks


    !*************************************************************************
    !> f_n <psi_n|a> r_{a,A} r^+_{B,b} <b|psi_n>
    !! =r^+_{B,b} <b|psi_n> f_n <psi_n|a> r_{a,A}
    !*************************************************************************
    subroutine calc_rnrl()
    integer isp,n
      
    n=wh%nasotot
    do isp=1,wh%nspin
        !< <b|psi_n> f_n <psi_n|a> r_{a,A}
        call zgemm('t','n',n,n,n,z1,bnd%nks(:,:,isp),n,bnd%r(:,:,isp),n, &
                &z0,bnd%nrl(:,:,isp),n)
        !< r^+_{B,b}<b|psi_n> f_n <psi_n|a> r_{a,A}
        call annxb('c',bnd%r(:,:,isp),bnd%nrl(:,:,isp),n,n)
        !< bring back to {A,B}
        bnd%nrl(:,:,isp)=transpose(bnd%nrl(:,:,isp))
    enddo
    return
      
    end subroutine calc_rnrl
      

    subroutine map_wh_bnd_matrix(x,y,lback)
    complex(q) x(wh%na2112),y(wh%nasotot,wh%nasotot,wh%nspin)
    logical,intent(in)::lback

    integer i,isp
    integer na2,naso,nbase,ibase
    complex(q),target::buf(wh%na2max**2)
    complex(q),pointer::p_buf(:,:)

    nbase=0
    ibase=0
    do i=1,wh%num_imp
        na2=wh%co(i)%dim2
        naso=wh%co(i)%dimso
        p_buf(1:na2,1:na2)=>buf(1:na2**2)
        if(.not.lback)then
            buf(1:na2**2)=x(ibase+1:ibase+na2*na2)
            if(wh%iso==1)call orbital_spin_trans(p_buf,na2,.false.)
            do isp=1,wh%nspin
                y(nbase+1:nbase+naso,nbase+1:nbase+naso,isp)= &
                        &p_buf((isp-1)*naso+1:isp*naso,(isp-1)*naso+1:isp*naso)
            enddo
        else
            buf=0
            do isp=1,wh%nspin
                p_buf((isp-1)*naso+1:isp*naso,(isp-1)*naso+1:isp*naso)= &
                        &y(nbase+1:nbase+naso,nbase+1:nbase+naso,isp)
            enddo
            if(wh%ispo==1)p_buf(1+naso:na2,1+naso:na2)=p_buf(1:naso,1:naso)
            if(wh%iso==1)call orbital_spin_trans(p_buf,na2,.true.)
            x(ibase+1:ibase+na2*na2)=buf(1:na2**2)
        endif
        nbase=nbase+naso
        ibase=ibase+na2*na2
    enddo
    return

    end subroutine map_wh_bnd_matrix


    !*************************************************************************
    subroutine calc_band_all(io)
    integer,intent(in)::io

    integer isym,ivec,iks,ikp,ikpl,nbands,isp,ispp
    complex(q),pointer :: hk(:,:)
    real(q),pointer :: ek_list(:,:,:)
    complex(q) :: r(wh%nasotot,wh%nasotot),la1(wh%nasotot,wh%nasotot)
     
    call set_time_point(1,2)

#ifdef mpi_mode
    allocate(ek_list(bnd%nmax,kpt%dim,wh%nspin))
    ek_list=0
#else
    ek_list=>bnd%ek
#endif

    ikpl=0
    !$omp parallel do &
    !$omp &firstprivate(ikpl) &
    !$omp &private(ivec,isym,isp,ispp,iks,ikp,nbands,hk,r,la1) &
    !$omp &schedule(static,1)
    do ivec=1,gp%nvec
        do iks=1,gp%kvec(ivec,2)
            ikp=gp%kvec(ivec,3)+iks
            if(gp%ipar==1)then
                ikpl=ikp
            else
                ikpl=ikpl+1
            endif
            do isym=1,sym%nop
                do isp=1,wh%nspin
                    ispp=min(isp,bnd%nspin_in)
                    nbands=bnd%ne(3,ikp,ispp)-bnd%ne(2,ikp,ispp)+1
                    hk(1:nbands,1:nbands)=>bnd%vk(1:nbands**2,isym,ikpl,isp)
                    r=bnd%r(:,:,isp)
                    la1=bnd%la1(:,:,isp)
                    bnd%vk(1:nbands**2,isym,ikpl,isp)=&
                            &bnd%hk0(1:nbands**2,isym,ikpl,ispp)
                    if(isym==1)then
                        ! save the energies outside of the energy window.
#ifdef mpi_mode
                        ek_list(:,ikp,isp)=bnd%ek(:,ikp,isp)
#endif
                        call calc_band_1k(hk,r,la1,nbands,wh%nasotot, &
                                &ek_list(bnd%ne(2,ikp,ispp): &
                                &bnd%ne(3,ikp,ispp),ikp,isp))
                    else
                        call calc_band_1k(hk,r,la1,nbands,wh%nasotot)   
                    endif
                enddo
            enddo !isym
        enddo ! iks
    enddo ! ivec
    !$omp end parallel do
    nullify(hk)
     
#ifdef mpi_mode
    call dsum_all_mpi(ek_list,bnd%nmax*kpt%dim*wh%nspin)
    bnd%ek=ek_list
    deallocate(ek_list)
#endif
    nullify(ek_list)

    call set_time_point(2,2)
    call print_time_usage('calc_band_all',2,io)
    return 

    end subroutine calc_band_all


    !*************************************************************************
    subroutine calc_band_1k(hk,r,la1,nbands,nasot,ek)
    integer,intent(in)::nbands,nasot
    complex(q),intent(in)::r(nasot,nasot),la1(nasot,nasot)
    complex(q),intent(inout)::hk(nbands,nbands)
    real(q),optional,intent(out)::ek(nbands)

    real(q) w(nbands)

    call calc_hamil_1k(hk,r,la1,nbands,nasot)
    call hermev('v','u',hk,w,nbands)
    if(present(ek))then
        ek=w
    endif
    return

    end subroutine calc_band_1k


    ! Here bnd%hk0 will be replaced by the renormalized occupation matrix
    ! of the bare Kohn-Sham bands.
    subroutine calc_kswt(io)
    integer,intent(in)::io

    integer ivec,iks,ikp,ikpl,nbands,isp,i,ispp,nsp,naso,nbase
    real(q) sumwt(2),maxoffdiag,fac_irspo
    complex(q),pointer :: p_vk(:,:),p_kswt(:,:),p_uk(:,:),p_utrans(:,:)
    complex(q),target::utrans(wh%nasomax**2)
    real(q),pointer :: p_ferwe(:)
     
    call set_time_point(1,2)
    sumwt=0
    maxoffdiag=0
    ikpl=0
    fac_irspo=1._q/wh%rspo

    !$omp parallel do &
    !$omp &firstprivate(ikpl) &
    !$omp &private(ivec,isp,ispp,iks,ikp,i,nbands,nbase,nsp, &
    !$omp         &p_kswt,p_ferwe,p_vk,p_uk,utrans,p_utrans, &
    !$omp         &naso,f_id) &
    !$omp &schedule(static,1) &
    !$omp &reduction(+:sumwt) &
    !$omp &reduction(max:maxoffdiag)
    do ivec=1,gp%nvec
        do iks=1,gp%kvec(ivec,2)
            ikp=gp%kvec(ivec,3)+iks
            if(gp%ipar==1)then
                ikpl=ikp
            else
                ikpl=ikpl+1
            endif
            do ispp=1,bnd%nspin_in
                if(bnd%nspin_in==1)then
                    nsp=wh%nspin
                else
                    nsp=ispp
                endif
                nbands=bnd%ne(3,ikp,ispp)-bnd%ne(2,ikp,ispp)+1
                ! do not need hk0 anymore, store kswt instead.
                p_kswt(1:nbands,1:nbands)=>bnd%hk0(1:nbands**2,sym%ie, &
                        &ikpl,ispp)
                p_kswt=0
                p_uk(1:nbands,1:nbands)=>bnd%psi0_b(1:nbands**2,ikpl,ispp)
                do isp=ispp,nsp
                    p_ferwe=>bnd%ferwe(bnd%ne(2,ikp,ispp): &
                            &bnd%ne(3,ikp,ispp),ikp,isp)
                    p_vk(1:nbands,1:nbands)=>bnd%vk(1:nbands**2,sym%ie,&
                            &ikpl,isp)
                    call calc_kswt_1k(p_kswt,p_vk,p_uk,p_ferwe,nbands, &
                            &fac_irspo,kpt%wt(ikp),isp)
                    sumwt(1)=sumwt(1)+sum(p_ferwe)
                enddo ! isp
                p_kswt=p_kswt*wh%rspo

                do i=1,nbands
                    sumwt(2)=sumwt(2)+real(p_kswt(i,i),q)
                    if(i==nbands)exit
                    maxoffdiag=max(maxoffdiag, &
                            &maxval(abs(p_kswt(i+1:,i)))/kpt%wt(ikp))
                enddo
            enddo
        enddo ! iks
    enddo ! ivec
    !$omp end parallel do
    nullify(p_vk,p_uk,p_ferwe,p_kswt,p_utrans)
      
#ifdef mpi_mode
    call gmpi_barrier()
    call sum_all_mpi(sumwt,2)
    call max_all_mpi(maxoffdiag)
#endif

    if(io>0)then
        write(io,'(" correlated subspace, sum_ferwt = ",f0.7, &
                &" sum_kswt = ",f0.7)')sumwt
        write(io,'("     max off diagonal of psik occ. mat. = ", &
                &f0.6)')maxoffdiag
        sumwt(1) = sumwt(1)-sumwt(2)
        if(abs(sumwt(1))>1.e-4_q)then
            write(0,'(" warning: too large sum_ferwt-sum_kswt = ", &
                    &f0.5,"!")')sumwt(1)
            write(io,'(" warning: too large sum_ferwt-sum_kswt = ", &
                    &f0.5,"!")')sumwt(1)
        endif
    endif

    call set_time_point(2,2)
    call print_time_usage('calc_kswt',2,io)
    return 

    end subroutine calc_kswt


    !*************************************************************************
    subroutine calc_kswt_1k(kswt,vk,uk,ferwe,nbands,wtk,wtk0,isp)
    integer,intent(in)::nbands,isp
    real(q),intent(in)::ferwe(nbands),wtk,wtk0
    complex(q),intent(in)::vk(nbands,nbands),uk(nbands,nbands)
    complex(q),intent(inout)::kswt(nbands,nbands)

    complex(q) nabr(nbands,nbands)

    call calc_nabr_1k(nabr,vk,ferwe,nbands,wtk,wtk0,isp) ! ~ f<psi|a><b|psi>
    nabr=transpose(nabr) ! ~ form <a|psi>f<psi|b>
    call uhau(nabr,uk,nbands,nbands,trul='n',trur='c')
    kswt=kswt+nabr ! density matrix <a|psi>f_psi<psi|b>
    return

    end subroutine calc_kswt_1k


    !*************************************************************************
    subroutine calc_hamil_1k(hk,r,la1,nbands,nasot)
    integer,intent(in)::nbands,nasot
    complex(q),intent(in)::r(nasot,nasot),la1(nasot,nasot)
    complex(q),intent(inout)::hk(nbands,nbands)

    integer ia
    complex(q) hk_sub(nasot,nbands)

    hk_sub = hk(1:nasot,:)
    call annxb('n',r,hk_sub,nasot,nbands,hk(1:nasot,:))
    call anmxbmm('c',hk(:,1:nasot),r,nbands,nasot) ! rhr^+
    hk(1:nasot,1:nasot)=hk(1:nasot,1:nasot)+la1
    return

    end subroutine calc_hamil_1k


    !*************************************************************************
    !< calc total magnetic moment
    !*************************************************************************
    subroutine calc_mup_dn(io)
    integer,intent(in)::io

    integer ik
      
    if(wh%nspin==1)return
    bnd%mup_dn=0
    do ik=1,kpt%dim
        bnd%mup_dn=bnd%mup_dn+sum(bnd%ferwe(:,ik,1)-bnd%ferwe(:,ik,2))
    enddo
    if(io>0)then
        write(io,'(" total magnetic moment: ",f8.3)')bnd%mup_dn
    endif
    return
      
    end subroutine calc_mup_dn
    

    subroutine calc_corr_ebwidth(io)
    integer,intent(in)::io
      
    integer ik,isp
    real(q) emin,emax
      
    if(io<0)return
    emin=100._q; emax=-100._q
    do ik=1,kpt%dim
        do isp=1,bnd%nspin_in
            emin=min(emin,bnd%ek(bnd%ne(2,ik,isp),ik,isp))
            emax=max(emax,bnd%ek(bnd%ne(3,ik,isp),ik,isp))
        enddo
    enddo
    write(io,'(" correlated block: emin/emix=",2f10.4)')emin,emax
    bnd%ebmin=emin; bnd%ebmax=emax
    bnd%ebwidth=emax-emin
    return
      
    end subroutine calc_corr_ebwidth
   

    subroutine set_kpt_icor(io)
    integer,intent(in)::io
      
    if(kpt%ismear/=-5)return
    kpt%icor=1
    if(abs(kpt%delta)>100._q)kpt%icor=0
    kpt%cordeg=-1.d-6
    if(kpt%delta>100._q)then
        kpt%cordeg=-kpt%delta+100._q
    elseif(kpt%delta>0._q)then
        kpt%cordeg=-kpt%delta
    elseif(kpt%delta<-100._q)then
        kpt%cordeg= kpt%delta-100._q
    elseif(kpt%delta<0._q)then
        kpt%cordeg= kpt%delta
    endif
    if(abs(kpt%cordeg)>.01_q)kpt%cordeg=-1.d-6
    if(io>0)write(io,'(" bz integration with tetra method, icor=",i2)')&
            &kpt%icor
    if(kpt%cordeg<0._q.and.io>0)write(io,'(" equal occupancy of degenerate &
        &states, tol=",e8.1)')-kpt%cordeg
    return
      
    end subroutine set_kpt_icor
      
    !*************************************************************************
    subroutine gutz_fermi(io)
    integer,intent(in)::io
      
    if(kpt%ismear==-5)then
        call gutz_fermi_tetra_w2k()
    elseif(kpt%ismear==0.or.kpt%ismear==-1)then
        call get_ef_fun()
    else
        write(0,'(" ismear=",i2)')kpt%ismear
        stop ' error: unsupported ismear!'
    endif
    call adjust_efermi_bandgap(io)
    if(io>0)write(io,'(" gutz fermi level=",f16.8)')bnd%ef
    return
      
    end subroutine gutz_fermi

    subroutine adjust_efermi_bandgap(io)
    integer,intent(in)::io

    integer ihomo,i
    real(q) ecmin,evmax

    ! check band gap
    ihomo=0
    do i=1,bnd%ne(3,1,1)
        if(bnd%ek(i,1,1)<bnd%ef)then
            ihomo=i
        else
            exit
        endif
    enddo
    ! conduction band minimum
    ecmin=minval(bnd%ek(ihomo+1,:,:))
    ! valence band maximum
    evmax=maxval(bnd%ek(ihomo,:,:))
    if(evmax<bnd%ef.and.ecmin>bnd%ef)then
        if(io>0)write(io,'(" band gap identified!")')
        bnd%ef=(evmax+ecmin)/2
    endif
    return

    end subroutine adjust_efermi_bandgap
      
    !*************************************************************************
    subroutine gutz_fermi_tetra_w2k()
    integer,parameter::nw=250000
    integer :: iw(nw)
    real(q),allocatable::eb(:,:,:),e_(:,:),weight_(:)
    integer,allocatable::nehelp(:,:)
    real(q)elecn,ef
    integer nkpt,nemax,jspin,iso,nbmax
    integer ik,isp,ispp,n1,n2
      
    nkpt=kpt%dim; nemax=bnd%nmax; iso=wh%iso
    allocate(nehelp(nkpt,2),eb(nemax,nkpt,2))
    eb=3._q
    do isp=1,2
        ispp=min(isp,bnd%nspin_in)
        nehelp(:,isp)=bnd%ne(1,:,ispp)
        ispp=min(isp,wh%nspin)
        do ik=1,nkpt
            eb(1:nehelp(ik,isp),ik,isp)=bnd%ek(1:nehelp(ik,isp),ik,ispp)
        enddo
    enddo
    elecn=bnd%nelet-1.d-10-bnd%nelet_frozen
    jspin=wh%ispin
    allocate(e_(nemax*jspin,nkpt)); e_=0
    e_(1:nemax,:) = eb(:,:,1)
    if(jspin==2)then
        e_(1+nemax:,:) = eb(:,:,2)
    endif

    allocate(weight_(nemax*jspin*nkpt)); weight_=0
    open(iu_kgen,file=kpt%file_name,status='old')
    ! consistent check
    read(iu_kgen,*)n1
    rewind(iu_kgen)
    if(n1/=nkpt)then
        bnd%ef=0
        return
    endif
    call dos(nemax*jspin,nkpt,e_,weight_,elecn/2.d0*iso*jspin,ef,iw,nw)
    call eweigh(ef,weight_,nemax,nkpt,jspin,nehelp,eb,nbmax)
    close(iu_kgen)
    do isp=1,wh%nspin
        do ik=1,nkpt
            n1=(ik-1)*nemax*jspin+(isp-1)*nemax+1
            n2=(ik-1)*nemax*jspin+ isp   *nemax
            bnd%ferwe(:,ik,isp)=weight_(n1:n2)*wh%rspo  ! convention
        enddo
    enddo
    bnd%ets2=0; bnd%ef=ef
    deallocate(eb,e_,weight_,nehelp)
    return
      
    end subroutine gutz_fermi_tetra_w2k
      
    !*************************************************************************
    subroutine bnd_modify_frozen()
    integer ik,i,isp,ispp,ntop
      
    if(bnd%n_frozen==0)return
    do ik=1,kpt%dim
        do isp=1,wh%nspin
            ispp=min(isp,bnd%nspin_in)
            ntop=bnd%ne(3,ik,ispp)
            do i=ntop-bnd%n_frozen*wh%iso/2+1,ntop
                if(abs(bnd%ek(i,ik,isp)-30._q)>1.e-6_q)then
                    write(0,'(" fetal error in bnd_modify_frozen (!=30): &
                            &fixed band level = ",2f12.4)')bnd%ek(i,ik,isp)
                    stop
                endif
                bnd%ek(i,ik,isp)=bnd%ef
                bnd%ferwe(i,ik,isp)=real(bnd%nelet_frozen,q)/bnd%n_frozen* &
                        &wh%rspo*kpt%wt(ik)
            enddo
        enddo
    enddo
    return
      
    end subroutine bnd_modify_frozen
    

    subroutine set_fermi_weight(mu)
    implicit none
    real(q) mu
      
    integer isp,ispp,ikp,ib
    real(q) dt
      
    bnd%ferwe=0
    do isp=1,wh%nspin
        ispp=min(isp,bnd%nspin_in)
        do ikp=1,kpt%dim
            do ib=1,bnd%ne(1,ikp,ispp)
                dt=(bnd%ek(ib,ikp,isp)-mu)/kpt%delta
                select case(kpt%ismear)
                case(-1)
                    bnd%ferwe(ib,ikp,isp)=fermi_fun(dt)
                case(0)
                    bnd%ferwe(ib,ikp,isp)=gauss_fun(dt)
                case default
                    stop ' error in nelect_fun: kpt%ismear not supported!'
                end select
                bnd%ferwe(ib,ikp,isp)=bnd%ferwe(ib,ikp,isp)*wh%rspo*kpt%wt(ikp)
            enddo
        enddo
    enddo
    return
      
    end subroutine set_fermi_weight
    

      
end module bandstru
  
!*****************************************************************************
subroutine get_ef_fun()
    use gprec
    use bandstru
    implicit none
    real(q) emin,emax
    real(q),parameter::eps=1.e-10_q
    real(q),external::rtbis,dif_nele_fun
      
    if(kpt%ensemble==0)then ! Canonical ensemble
        bnd%ferwe=0
        if(abs(bnd%nelet-bnd%nelet_frozen)<1.e-6_q)then
            bnd%ef=0; return
        endif
        emin=minval(bnd%ek); emax=maxval(bnd%ek)
        bnd%ef=rtbis(dif_nele_fun,emin,emax,eps)
    else
        call set_fermi_weight(bnd%ef)
    endif
      
    select case(kpt%ismear)
    case(-1)
        call calc_entropy_fermi()
    case(0)
        call calc_correction_gauss()
    end select
    return
      
end subroutine get_ef_fun
  
!*****************************************************************************
subroutine calc_correction_gauss()
    use gprec
    use gconstant
    use bandstru
    implicit none
    integer ikp,isp,ispp,ib
    real(q) eta,de
      
    bnd%ets2=0
    do isp=1,wh%nspin
        eta=0
        ispp=min(isp,bnd%nspin_in)
        do ikp=1,kpt%dim
            do ib=1,bnd%ne(1,ikp,ispp)
                de=(bnd%ek(ib,ikp,isp)-bnd%ef)/kpt%delta
                de=de*de
                if(de<15._q)eta=eta+0.5*kpt%delta*exp(-de)*kpt%wt(ikp)*wh%rspo
            enddo
        enddo
        eta=-eta*2._q/sqrt(pi)
        bnd%ets2(isp)=eta/2._q
    enddo
    return
      
end subroutine calc_correction_gauss
  
!*****************************************************************************
subroutine calc_entropy_fermi()
    use gprec
    use gconstant
    use bandstru
    implicit none
    integer ikp,isp,ispp,ib
    real(q) entr,focc,f1,f2,eint
      
    bnd%ets2=0
    do isp=1,wh%nspin
        entr=0
        ispp=min(isp,bnd%ispin_in)
        do ikp=1,kpt%dim
            do ib=1,bnd%ne(1,ikp,ispp)
                focc=bnd%ferwe(ib,ikp,isp)/wh%rspo/kpt%wt(ikp)
                f1=focc; f2=1._q-focc
                if(f1>0._q.and.f2>0._q)then
                    eint=f1*log(f1)+f2*log(f2)
                    entr=entr+eint*kpt%wt(ikp)*wh%rspo
                endif
            enddo
        enddo
        bnd%ets2(isp)=kpt%delta*entr
    enddo
    return
      
end subroutine calc_entropy_fermi
  
!*****************************************************************************
function dif_nele_fun(mu)
    use gprec
    use bandstru
    use gutil
    implicit none
    real(q) mu,dif_nele_fun
      
    real(q) nelet
      
    call set_fermi_weight(mu)
    nelet=sum(bnd%ferwe)
    dif_nele_fun=nelet-(bnd%nelet-bnd%nelet_frozen)
    return
      
end function dif_nele_fun
  
!*****************************************************************************
function rtbis(func,x1,x2,tol)
    use gprec
    implicit none
    real(q) x1,x2,tol,dx,xmid,fmid,f,rtbis
    integer j
    integer,parameter::jmax=500
    real(q),external::func
      
    fmid=func(x2)
    f=func(x1)
    if(f*fmid.ge.0._q) goto 900
    if(f.lt.0.)then
        rtbis=x1
        dx=x2-x1
    else
        rtbis=x2
        dx=x1-x2
    endif
    do j=1,jmax
        dx=dx*.5_q
        xmid=rtbis+dx
        fmid=func(xmid)
        if(fmid.le.0._q)rtbis=xmid
        if(abs(fmid)<tol)then
            rtbis=xmid
            return
        endif
    enddo
    goto 910
900 if(abs(fmid).lt.tol)then
        rtbis=x2
        return
    else if(abs(f).lt.tol)then
        rtbis=x1
        return
    endif
    write(0,'(" fmid and f=",2f16.8)')fmid,f
    stop ' error in rtbis: root must be bracketed for bisection.'
910 stop ' error in rtbis: too many bisections in rtbis'
      
end function rtbis
