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

module warehouse
    use gprec
    use gutil
    use gconstant, only: zi,small
    use corrorb
    use hdf5
    use ghdf5_base
    implicit none

    type frozen
        integer :: nsorb,nelect
        integer,pointer :: idx_orb(:)
    end type frozen
      
    type matrix_basis
        integer :: dimhst=0,dimhsmax=0,dimhs1123
        !< dimension of the matrix basis for each atom
        integer,pointer :: dim_hs(:)
        complex(q),pointer :: hs(:)=>null() ! matrix basis
        integer,pointer :: m_struct(:)=>null() ! self-energy structure
    end type matrix_basis

    type ware_house
        integer :: num_imp,ntyp
        integer,allocatable :: ityp_imp(:), & !< type of the impurity.
                &imap_imp(:), & !< equiv. impurity map index. imap_imp(i)<=i.
                &na2_imp(:), & !< impurity spin-orbital dimension list
                !< Index to whether calculate it on local core.
                &local_imp(:),naso_imp(:), &
                &nval_bot_ityp(:),nval_top_ityp(:)
        real(q),allocatable :: et1(:),eu2(:)
        integer :: na2112,na2max,nasomax,nasotot
        integer :: iso,ispin,ispin_in,ispo,rspo,nspin,ivext=1
        type(corr_orb),allocatable::co(:)
        type(frozen),pointer::fz(:)

        real(q) :: ef,edcla1,symerr=0._q
        integer :: r_factor=2
        complex(q),pointer :: r(:),r0(:),d(:),d0(:),la1(:),la2(:),z(:), &
                &copy(:),nks(:),nc_var(:),nc_phy(:),h1e(:),isimix(:), &
                &vext(:)=>null(),&
                &sx(:)=>null(),sy(:)=>null(),sz(:)=>null(), &
                &lx(:)=>null(),ly(:)=>null(),lz(:)=>null(),r_coef(:), &
                ! additional rotations for bare hamiltonian.
                ! complex 
                &db2sab(:)=>null() ! complex spherical 
        real(q),pointer :: nks_coef(:),la1_coef(:),ncv_coef(:)
        !> hermitian matrix basis, (hermitian) maxtrix basis with
        !! selective mott localization.
        type (matrix_basis) :: hm, hm_l, hm_r
    end type ware_house

    type(ware_house) :: wh
      
    contains


    subroutine init_warehouse(io)
    integer,intent(in)::io

    logical lexist
   
    call gh5_open_r('GPARAM.h5',f_id)
    !< Get num_imp, ityp_imp, imap_imp and na2_imp
    call set_impurities_index()
    call write_impurities_index(io)
    call gh5_read_wh_hs('/dim_hs_imp','/HS','/SIGMA_STRUCT',wh%hm)
    call set_wh_db2sab(io)
    call set_wh_sl_vec()
    allocate(wh%co(wh%num_imp))
    call set_v2e_list(io)
    call gh5_close(f_id)

    inquire(file='GVEXT.h5',exist=lexist)
    if(lexist)then
        call gh5_open_r('GVEXT.h5',f_id)
        call set_wh_vext(io)
        call gh5_read(wh%ivext, '/givext', f_id)
        call gh5_close(f_id)
    endif

    inquire(file='GMOTT.h5',exist=lexist)
    if(lexist)call gh5_open_r('GMOTT.h5',f_id)
    call set_wh_hs('/dim_hs_r_imp',"/HS_R","/SIGMA_STRUCT_R",wh%hm_r,lexist)
    call set_wh_hs('/dim_hs_l_imp',"/HS_L","/SIGMA_STRUCT_L",wh%hm_l,lexist)
    call set_mott_localized_orbitals(lexist)
    if(lexist)call gh5_close(f_id)
    call output_mott_localized_orbitals(io)
    call alloc_warehouse()
    call link_co_warehouse()
    return

    end subroutine init_warehouse


    !< Impurities should be ordered according to types.
    subroutine set_v2e_list(io)
    integer,intent(in)::io

    integer i,imap,na2

    do i=1,wh%num_imp
        imap=wh%imap_imp(i)
        if(imap/=i)then
            wh%co(i)%v2e=>wh%co(imap)%v2e
            wh%co(i)%v_j2e=>wh%co(imap)%v_j2e
        else
            na2=wh%na2_imp(i)
            allocate(wh%co(i)%v2e(na2,na2,na2,na2))
            call gh5_read(wh%co(i)%v2e,na2,na2,na2,na2,'/IMPURITY_'// &
                    &trim(int_to_str(i))//'/V2E',f_id)
            if(io>0)then
                write(io,'(" imp = ",I2," v2e(1,1,1,1) = ",f0.2)')i,&
                    &real(wh%co(i)%v2e(1,1,1,1))
            endif
            allocate(wh%co(i)%v_j2e(na2,na2,na2,na2))
            call get_coul_exchange(na2,wh%co(i)%v2e,wh%co(i)%v_j2e)
        endif
    enddo
    return

    end subroutine set_v2e_list


    subroutine calc_la1_hf_list()
    integer i

    do i=1,wh%num_imp
        call calc_co_la1_hf(wh%co(i))
    enddo
    return

    end subroutine calc_la1_hf_list


    subroutine link_co_warehouse()
    integer i,na2,na22,naso,nasoo,dim_hs,ibase,jbase

    ibase=0; jbase=0
    do i=1,wh%num_imp
        na2=wh%na2_imp(i)
        na22=na2*na2
        wh%co(i)%dim2=na2; wh%co(i)%dim=na2/2
        naso=na2*wh%iso/2
        wh%co(i)%dimso=naso
        nasoo=wh%co(i)%dimso**2
        wh%co(i)%dim_hs_l = wh%hm_l%dim_hs(i)
        wh%co(i)%nks   (1:na2,1:na2) => wh%nks   (ibase+1:ibase+na22)
        wh%co(i)%isimix(1:na2,1:na2) => wh%isimix(ibase+1:ibase+na22)
        wh%co(i)%nc_var(1:na2,1:na2) => wh%nc_var(ibase+1:ibase+na22)
        wh%co(i)%nc_phy(1:na2,1:na2) => wh%nc_phy(ibase+1:ibase+na22)
        wh%co(i)%r     (1:na2,1:na2) => wh%r     (ibase+1:ibase+na22)
        wh%co(i)%r0    (1:na2,1:na2) => wh%r0    (ibase+1:ibase+na22)
        wh%co(i)%z     (1:na2,1:na2) => wh%z     (ibase+1:ibase+na22)
        wh%co(i)%d0    (1:na2,1:na2) => wh%d0    (ibase+1:ibase+na22)
        wh%co(i)%d     (1:na2,1:na2) => wh%d     (ibase+1:ibase+na22)
        wh%co(i)%la1   (1:na2,1:na2) => wh%la1   (ibase+1:ibase+na22)
        wh%co(i)%la2   (1:na2,1:na2) => wh%la2   (ibase+1:ibase+na22)
        wh%co(i)%h1e   (1:na2,1:na2) => wh%h1e   (ibase+1:ibase+na22)

        if(associated(wh%vext))then
            wh%co(i)%vext(1:na2,1:na2) => wh%vext(ibase+1:ibase+na22)
        endif

        if(associated(wh%sx))then
            wh%co(i)%sx(1:na2,1:na2) => wh%sx(ibase+1:ibase+na22)
            wh%co(i)%sy(1:na2,1:na2) => wh%sy(ibase+1:ibase+na22)
        endif
        if(associated(wh%sz))then
            wh%co(i)%sz(1:na2,1:na2) => wh%sz(ibase+1:ibase+na22)
        endif
        if(associated(wh%lx))then
            wh%co(i)%lx(1:na2,1:na2) => wh%lx(ibase+1:ibase+na22)
            wh%co(i)%ly(1:na2,1:na2) => wh%ly(ibase+1:ibase+na22)
        endif
        if(associated(wh%lz))then
            wh%co(i)%lz(1:na2,1:na2) => wh%lz(ibase+1:ibase+na22)
        endif
        if(associated(wh%db2sab))then
            wh%co(i)%db2sab(1:na2,1:na2) => &
                &wh%db2sab(ibase+1:ibase+na22)
        endif

        wh%co(i)%m_struct(1:na2,1:na2) => wh%hm%m_struct(ibase+1:ibase+na22)

        dim_hs=wh%hm_l%dim_hs(i)
        if(dim_hs>0)then
            wh%co(i)%hs_l(1:na2,1:na2,1:dim_hs) => wh%hm_l%hs(jbase+1: &
                    &jbase+na22*dim_hs)
        endif
        jbase=jbase+na22*dim_hs
        ibase=ibase+na22
    enddo
    return

    end subroutine link_co_warehouse


    ! \sum_ij \lambda_ij <c^\dagger_i c_j>
    subroutine calc_edcla1()
    
    wh%edcla1=real(sum(wh%la1*wh%nks),q)
    return

    end subroutine calc_edcla1


    subroutine set_wh_hs_dimtot(mb)
    type(matrix_basis),intent(inout)::mb

    integer i

    mb%dimhsmax=maxval(mb%dim_hs)
    mb%dimhst=0
    mb%dimhs1123=0
    do i=1,wh%num_imp
        mb%dimhs1123=mb%dimhs1123+wh%na2_imp(i)**2*mb%dim_hs(i)
        if(wh%imap_imp(i)==i)then
            mb%dimhst=mb%dimhst+mb%dim_hs(i)
        endif
    enddo
    return

    end subroutine set_wh_hs_dimtot


    subroutine set_impurities_index()
    logical lexist

    call gh5_read(wh%iso,'/iso',f_id)
    call gh5_read(wh%ispin,'/ispin',f_id)
    call gh5_read(wh%num_imp,'/num_imp',f_id)
    allocate(wh%ityp_imp(wh%num_imp), wh%imap_imp(wh%num_imp), &
            &wh%na2_imp(wh%num_imp),wh%naso_imp(wh%num_imp))
    call gh5_read(wh%ityp_imp,wh%num_imp,'/ITYP_IMP',f_id)
    call gh5_read(wh%imap_imp,wh%num_imp,'/IMAP_IMP',f_id)
    call gh5_read(wh%na2_imp,wh%num_imp,'/na2_imp',f_id)

    wh%ispo=max(wh%iso,wh%ispin)
    wh%rspo=3-wh%ispo
    wh%nspin=max(1,wh%ispin/wh%iso)
    wh%na2max=maxval(wh%na2_imp)
    wh%nasomax=wh%na2max*wh%iso/2
    wh%na2112=sum(wh%na2_imp**2)
    wh%nasotot=sum(wh%na2_imp)*wh%iso/2
    wh%naso_imp=wh%na2_imp*wh%iso/2
    wh%ntyp=maxval(wh%ityp_imp)

    call h5lexists_f(f_id,'/nval_bot_ityp',lexist,gh5_err)
    if(lexist)then
        allocate(wh%nval_bot_ityp(wh%ntyp),wh%nval_top_ityp(wh%ntyp))
        call gh5_read(wh%nval_bot_ityp,wh%ntyp,'/nval_bot_ityp',f_id)
        call gh5_read(wh%nval_top_ityp,wh%ntyp,'/nval_top_ityp',f_id)
    endif
    return

    end subroutine set_impurities_index


    subroutine wh_set_local_imp(local_imp)
    integer,intent(in)::local_imp(wh%num_imp)

    allocate(wh%local_imp(wh%num_imp))
    wh%local_imp=local_imp
    return

    end subroutine wh_set_local_imp


    subroutine write_impurities_index(io)
    integer,intent(in)::io

    if(io<0)return
    write(io,'(" total number of impurities = ",I4)')wh%num_imp
    write(io,'(" impurity type indices:")')
    write(io,'(4x,10I4)')wh%ityp_imp
    write(io,'(" impurity imap indices:")')
    write(io,'(4x,10I4)')wh%imap_imp
    write(io,'(" impurity num_spin_orbitals:")')
    write(io,'(4x,10I4)')wh%na2_imp
    return

    end subroutine write_impurities_index


    subroutine alloc_warehouse()
      
    allocate(wh%r(wh%na2112), wh%r0(wh%na2112), &
            &wh%z(wh%na2112), wh%d(wh%na2112), &
            &wh%d0(wh%na2112), wh%la1(wh%na2112), &
            &wh%la2(wh%na2112), wh%nks(wh%na2112), &
            &wh%isimix(wh%na2112), wh%nc_var(wh%na2112), &
            &wh%nc_phy(wh%na2112), wh%copy(wh%na2112), &
            &wh%h1e(wh%na2112), &
            &wh%nks_coef(wh%hm_l%dimhst),wh%la1_coef(wh%hm_l%dimhst), &
            &wh%ncv_coef(wh%hm_l%dimhst),wh%r_coef(wh%hm_r%dimhst), &
            &wh%et1(wh%num_imp), wh%eu2(wh%num_imp))
    wh%r=0; wh%r0=0; wh%z=0; wh%d=0; wh%d0=0; wh%la1=0; wh%la2=0
    wh%nks=0; wh%isimix=0; wh%nc_var=0; wh%nc_phy=0
    return
      
    end subroutine alloc_warehouse
     

    subroutine init_wh_x(io,r_default)
    integer,intent(in)::io
    real(q),intent(in)::r_default

    logical lexist

    inquire(file='WH_RL_INP.h5',exist=lexist)
    if(lexist)then
        call gh5_read_wh_rl()
    else
        wh%la1=wh%h1e
        call set_diagonal_r(r_default)
        if(wh%ivext==0.and.associated(wh%vext))then
            wh%la1=wh%la1+wh%vext
        endif
    endif
    call modify_r_la1_frozen(0, 30._q)
    call calc_r_pp(io,'r-inp')
    call calc_la1_pp(io,'la1-inp')
    return

    end subroutine init_wh_x


    subroutine set_diagonal_r(r_default)
    real(q),intent(in)::r_default

    integer i,j

    do i=1,wh%num_imp
        do j=1,wh%na2_imp(i)
            wh%co(i)%r(j,j)=r_default
        enddo
    enddo
    return

    end subroutine set_diagonal_r


    !*************************************************************************
    subroutine gh5_wrt_wh_matrix_list(apath, a)
    character(*),intent(in)::apath
    complex(q),target,intent(in)::a(wh%na2112)
      
    integer i,na2,na22,ibase
    complex(q),pointer::p_a(:,:)
      
    ibase=0
    do i = 1,wh%num_imp
        na2=wh%na2_imp(i)
        na22=na2**2
        p_a(1:na2,1:na2)=>a(ibase+1:ibase+na22)
        call gh5_write(p_a,na2,na2,'/IMPURITY_'//trim(int_to_str(i))// &
                &apath,f_id)
        ibase=ibase+na22
    enddo
    nullify(p_a)
    return
      
    end subroutine gh5_wrt_wh_matrix_list


    subroutine calc_isimix(io)
    integer,intent(in)::io

    integer i,imap

    do i=1,wh%num_imp
        imap=wh%imap_imp(i)
        if(imap==i)then
            call calc_co_isimix(wh%co(i))
        else
            wh%co(i)%isimix=wh%co(imap)%isimix
        endif
    enddo
    if(io>0)then
        call output_matrices('isimix-var',wh%isimix,wh%na2112,wh%num_imp, &
                &wh%na2_imp,io,0)
    endif
    return

    end subroutine calc_isimix


    subroutine eval_sl_vec_all(mode,io)
    integer,intent(in)::mode,io

    integer i
    character name_*12

    if(.not.associated(wh%sx))return
    if(mode==1)then
        name_='var-hf-only'
    else
        name_='physical'
    endif
    do i=1,wh%num_imp
        call eval_co_sl_vec(wh%co(i),mode)
        if(io>0)then
            write(io,'(" imp =",i4," s_xyz(",a12,") =",3f12.5)')i,name_, &
                &wh%co(i)%s_val(:,mode)
            write(io,'(10x," l_xyz(",a12,") =",3f12.5)')name_, &
                &wh%co(i)%l_val(:,mode)
        endif
    enddo
    return

    end subroutine eval_sl_vec_all


    subroutine gh5_read_wh_h1e(nspin_in)
    integer,intent(in)::nspin_in

    integer i,naso,na22,ibase
    complex(q),allocatable::h1e(:,:)

    wh%h1e=0
    do i=1,wh%num_imp
        naso=wh%naso_imp(i)
        if(i>1)then
            if(wh%naso_imp(i-1)/=naso)then
                deallocate(h1e)
                allocate(h1e(naso,naso))
            endif
        else
            allocate(h1e(naso,naso))
        endif
        call gh5_read(h1e,naso,naso,'/IMPURITY_'//trim(int_to_str(i))// &
                &"/H1E_SPIN1",f_id)
        wh%co(i)%h1e(1:naso,1:naso)=h1e
        if(wh%iso==1)then
            if(nspin_in==2)then
                call gh5_read(h1e,naso,naso,'/IMPURITY_' &
                        &//trim(int_to_str(i))//"/H1E_SPIN2",f_id)
            endif
            wh%co(i)%h1e(1+naso:,1+naso:)=h1e
        endif
    enddo
    deallocate(h1e)
    return
      
    end subroutine gh5_read_wh_h1e


    subroutine gh5_read_wh_matrix_list(apath, dim_imp, a)
    character(*),intent(in)::apath
    integer,intent(in)::dim_imp(*)
    complex(q),target,intent(out)::a(*)
      
    integer i,na2,na22,ibase
    complex(q),pointer::p_a(:,:)
      
    ibase=0
    do i = 1,wh%num_imp
        na2=dim_imp(i)
        na22=na2**2
        p_a(1:na2,1:na2)=>a(ibase+1:ibase+na22)
        call gh5_read(p_a,na2,na2,'/IMPURITY_'//trim(int_to_str(i))// &
                &apath,f_id)
        ibase=ibase+na22
    enddo
    nullify(p_a)
    return
      
    end subroutine gh5_read_wh_matrix_list
    

    !*************************************************************************
    subroutine gh5_wrt_wh_rl(fname)
    character(*) fname

    call gh5_open_w(fname, f_id)
    call gh5_create_impurity_groups(f_id)
    call gh5_wrt_wh_matrix_list('/R',wh%r)
    call gh5_wrt_wh_matrix_list('/LAMBDA',wh%la1)
    call gh5_close(f_id)
    return
      
    end subroutine gh5_wrt_wh_rl
      
      
    !*************************************************************************
    subroutine gh5_read_wh_rl()
    logical lexist

    call gh5_open_r('WH_RL_INP.h5', f_id)
    call gh5_read_wh_matrix_list('/LAMBDA',wh%na2_imp,wh%la1)
    call h5lexists_f(f_id,'/IMPURITY_1/R',lexist,gh5_err)
    if(lexist)then
        call gh5_read_wh_matrix_list('/R',wh%na2_imp,wh%r)
    else
        call set_diagonal_r(0.999_q)
    endif
    call gh5_close(f_id)
    return
      
    end subroutine gh5_read_wh_rl
    

    !*************************************************************************
    ! read [s_x, s_y, s_z] and [l_x, l_y, l_z]
    !*************************************************************************
    subroutine set_wh_sl_vec()
    integer i
    logical lexist
   
    call h5lexists_f(f_id,'/IMPURITY_1/SZ',lexist,gh5_err)
    if(lexist)then
        allocate(wh%sz(wh%na2112))
        call gh5_read_wh_matrix_list('/SZ',wh%na2_imp,wh%sz)
    endif
    call h5lexists_f(f_id,'/IMPURITY_1/LZ',lexist,gh5_err)
    if(lexist)then
        allocate(wh%lz(wh%na2112))
        call gh5_read_wh_matrix_list('/LZ',wh%na2_imp,wh%lz)
    endif
    call h5lexists_f(f_id,'/IMPURITY_1/SX',lexist,gh5_err)
    if(lexist)then
        allocate(wh%sx(wh%na2112),wh%sy(wh%na2112), &
                &wh%lx(wh%na2112),wh%ly(wh%na2112))
        call gh5_read_wh_matrix_list('/SX',wh%na2_imp,wh%sx)
        call gh5_read_wh_matrix_list('/SY',wh%na2_imp,wh%sy)
        call gh5_read_wh_matrix_list('/LX',wh%na2_imp,wh%lx)
        call gh5_read_wh_matrix_list('/LY',wh%na2_imp,wh%ly)
    endif
    return
      
    end subroutine set_wh_sl_vec
    

    subroutine gh5_read_wh_hs(dpath,hpath,spath,mb)
    character(*),intent(in)::dpath,hpath,spath
    type(matrix_basis),intent(inout)::mb

    integer na2,dim_hs,na223,i,ibase,jbase
    complex(q),pointer::mb_hs(:,:,:)
    integer,pointer::mb_m_struct(:,:)
    
    allocate(mb%dim_hs(wh%num_imp))
    call gh5_read(mb%dim_hs,wh%num_imp,dpath,f_id)
    call set_wh_hs_dimtot(mb)
    allocate(mb%hs(mb%dimhs1123),mb%m_struct(wh%na2112))
    mb%m_struct=0

    ibase=0; jbase=0
    do i=1,wh%num_imp
        na2=wh%na2_imp(i)
        dim_hs=mb%dim_hs(i)
        na223=na2**2*dim_hs
        if(na223>0)then
            mb_hs(1:na2,1:na2,1:dim_hs)=>mb%hs(ibase+1:ibase+na223)
            call gh5_read(mb_hs,na2,na2,dim_hs,'/IMPURITY_'// &
                    &trim(int_to_str(i))//hpath,f_id)
            ibase=ibase+na223
            mb_m_struct(1:na2,1:na2)=>mb%m_struct(jbase+1:jbase+na2**2)
            call gh5_read(mb_m_struct,na2,na2,'/IMPURITY_'// &
                    &trim(int_to_str(i))//spath,f_id)
        endif
        jbase=jbase+na2**2
    enddo
    nullify(mb_hs,mb_m_struct)
    return
      
    end subroutine gh5_read_wh_hs
     

    !*************************************************************************
    subroutine set_wh_hs(dpath,hpath,spath,mb,lopen)
    character(*),intent(in) :: dpath,hpath,spath
    type(matrix_basis),intent(inout)::mb
    logical,intent(in)::lopen

    logical lexist
     
    if(lopen)then
        call h5lexists_f(f_id,dpath,lexist,gh5_err)
    else
        lexist=.false.
    endif
    if(lexist)then
        call gh5_read_wh_hs(dpath,hpath,spath,mb)
    else
        mb%dim_hs=>wh%hm%dim_hs
        mb%hs=>wh%hm%hs
        mb%dimhsmax=wh%hm%dimhsmax
        mb%dimhst=wh%hm%dimhst
        mb%dimhs1123=wh%hm%dimhs1123
    endif
    return
      
    end subroutine set_wh_hs


    ! Convention: 
    ! iso=1 -> {{complext spherical harmonics (CSH)}_up}
    ! to symmetry adapted basis (SAB);
    ! iso=2 -> relativistic harmonmics to SAB;
    subroutine set_wh_db2sab(io)
    integer,intent(in)::io

    logical lexist
   
    call h5lexists_f(f_id,'/IMPURITY_1/DB_TO_SAB',lexist,gh5_err)
    if(lexist)then
        allocate(wh%db2sab(wh%na2112))
        call gh5_read_wh_matrix_list('/DB_TO_SAB',wh%na2_imp,wh%db2sab)
        if(io>0)then
            write(io,'(" transformation db_to_sab read in.")')
        endif
    endif
    return
      
    end subroutine set_wh_db2sab
     

    ! read in a list of local external potentials.
    subroutine set_wh_vext(io)
    integer,intent(in)::io

    allocate(wh%vext(wh%na2112))
    call gh5_read_wh_matrix_list('/VEXT',wh%na2_imp,wh%vext)
    if(io>0)then
        write(io,'(" list of local vext read in.")')
        call output_matrices('v_ext',wh%vext,wh%na2112,wh%num_imp, &
                &wh%na2_imp,io,0)
    endif
    return
 
    end subroutine set_wh_vext


    !*************************************************************************
    ! Mott localized spin-orbital information
    !*************************************************************************
    subroutine set_mott_localized_orbitals(lopen)
    logical,intent(in)::lopen

    integer i,i_
    logical lexist
      
    allocate(wh%fz(wh%num_imp))
    if(lopen)then
        call h5lexists_f(f_id,'/IMPURITY_1/num_mott_orbitals',lexist,gh5_err)
    else
        lexist=.false.
    endif
    if(.not.lexist)then
        wh%fz(:)%nsorb=0; wh%fz(:)%nelect=0
    else
        do i=1,wh%num_imp
            call gh5_read(wh%fz(i)%nsorb,'/IMPURITY_'//trim(int_to_str(i))// &
                    &'/num_mott_orbitals', f_id)
            call gh5_read(wh%fz(i)%nelect,'/IMPURITY_'//trim(int_to_str(i))// &
                    &'/num_mott_electrons', f_id)
            if(wh%fz(i)%nsorb>0)then
                allocate(wh%fz(i)%idx_orb(wh%fz(i)%nsorb))
                call gh5_read(wh%fz(i)%idx_orb,wh%fz(i)%nsorb, &
                        &'/IMPURITY_'//trim(int_to_str(i))//&
                        &'/mott_orbital_indices', f_id)
            endif
        enddo
    endif
    return
      
    end subroutine set_mott_localized_orbitals
    

    !*************************************************************************
    subroutine output_mott_localized_orbitals(io)
    integer,intent(in)::io
      
    integer i
      
    if(io<0)return
    write(io,'(" mott localized orbital info:")')
    do i=1,wh%num_imp
        write(io,'(" i=",i3," nsorb=",i2," nelect=",i2)')i, &
                &wh%fz(i)%nsorb,wh%fz(i)%nelect
        if(wh%fz(i)%nsorb>0)then
            write(io,'("    idx_orb=",14i3)')wh%fz(i)%idx_orb
        endif
    enddo
    return
      
    end subroutine output_mott_localized_orbitals
    

    !*************************************************************************
    subroutine modify_r_la1_frozen(mode,val)
    integer,intent(in)::mode
    real(q),intent(in)::val

    integer i,j,j_
      
    do i=1,wh%num_imp
        do j=1,wh%fz(i)%nsorb
            j_=wh%fz(i)%idx_orb(j)
            wh%co(i)%r(j_,:)=0
            wh%co(i)%la1(j_,:)=0; wh%co(i)%la1(:,j_)=0
            if(mode>0)then
                wh%co(i)%la1(j_,j_)=val
            endif
        enddo
    enddo
    return
      
    end subroutine modify_r_la1_frozen


    !*************************************************************************
    !< post processing: print and symmetizations.
    !*************************************************************************
    subroutine calc_nks_pp(io,lchk)
    integer,intent(in)::io
    logical,optional,intent(in)::lchk

    real(q) maxerr

    if(io>0)then
        call output_matrices('nks-unsym',wh%nks,wh%na2112,wh%num_imp, &
                &wh%na2_imp,io,1)
    endif
    wh%copy=wh%nks
    call symm_dm_across_atoms(wh%nks)
    ! Get expansion coefficients only
    call hm_expand_all_herm(wh%nks,wh%nks_coef,wh%hm_l,1,.true.)
    ! Symmetrization (note non-zero diagonals)
    call hm_expand_all_sym(wh%nks,wh%hm,.true.,.true.)
    maxerr=maxval(abs(wh%nks-wh%copy))
    wh%symerr=max(wh%symerr,maxerr)
    if(io>0)then
        call output_matrices('nks-sym',wh%nks,wh%na2112,wh%num_imp, &
                &wh%na2_imp,io,1)
        write(io,'(" max error due to symmetrization = ",f12.6)') maxerr
    endif
    call chk_eigens_matrix_list('nks',wh%nks,wh%na2112,wh%num_imp, &
            &wh%na2_imp,io,0)
    call calc_nks_tot(io)
    if(.not.present(lchk))then
        call chk_nks_tiny()
    endif
    return

    end subroutine calc_nks_pp


    subroutine symm_dm_across_atoms(dm)
    complex(q),intent(inout)::dm(wh%na2112)

    integer i,j,isum,nbase,na22
    complex(q) dms(wh%na2max*wh%na2max)

    dms=0; isum=0; nbase=0
    do i=1,wh%num_imp
        na22=wh%na2_imp(i)**2
        nbase=nbase+na22
        if(isum==0)then
            dms(1:na22)=dm(nbase-na22+1:nbase)
        else
            dms(1:na22)=dms(1:na22)+dm(nbase-na22+1:nbase)
        endif
        isum=isum+1
        if(i<wh%num_imp)then
            if(wh%imap_imp(i)==wh%imap_imp(i+1))cycle
        endif
        if(isum>1)then
            dms=dms/isum
            do j=1,isum
                dm(nbase-j*na22+1:nbase-(j-1)*na22)=dms(1:na22)
            enddo
        endif
        isum=0; dms=0
    enddo

    end subroutine symm_dm_across_atoms


    subroutine calc_nks_tot(io)
    integer,intent(in)::io

    integer i
    real(q) res

    do i=1,wh%num_imp
        if(wh%ispin_in==1)then
            ! case of lda
            res = trace_a(wh%co(i)%nks,wh%co(i)%dim2)
            wh%co(i)%net=res/2
        else
            ! case of lsda
            call calc_co_net(wh%co(i))
        endif
    enddo
    if(io>0)then
        write(io,'(" nele_loc total:")')
        write(io,'(4x,2(2f10.5,3x))')(wh%co(i)%net,i=1,wh%num_imp)
    endif
    return

    end subroutine calc_nks_tot


    subroutine calc_ncvar_pp(io)
    integer,intent(in)::io

    real(q) maxerr

    if(io>0)then
        call output_matrices('ncv-unsym',wh%nc_var,wh%na2112,wh%num_imp, &
                &wh%na2_imp,io,1)
    endif
    wh%copy=wh%nc_var
    ! Get expansion coefficients only
    call hm_expand_all_herm(wh%nc_var,wh%ncv_coef,wh%hm_l,1,.true.)
    ! Symmetrization (note non-zero diagonals)
    call hm_expand_all_sym(wh%nc_var,wh%hm,.true.,.true.)
    maxerr=maxval(abs(wh%nc_var-wh%copy))
    wh%symerr=max(wh%symerr,maxerr)
    if(io>0)then
        call output_matrices('ncv-sym',wh%nc_var,wh%na2112,wh%num_imp, &
                &wh%na2_imp,io,1)
        write(io,'(" max error due to symmetrization = ",f12.6)')maxerr
    endif
    return

    end subroutine calc_ncvar_pp


    subroutine chk_nks_tiny()
    integer i,i1
    real(q) res

    do i=1,wh%num_imp
        do i1=1,wh%co(i)%dim2
            res=real(wh%co(i)%nks(i1,i1),q)
            if(res<small)then
                write(0,'(" warning in chk_nks_tiny: too small diagonal &
                        &elements",2f16.5)')wh%co(i)%nks(i1,i1)
            endif
            wh%co(i)%nks(i1,i1)=max(real(wh%co(i)%nks(i1,i1),q),small)
        enddo
    enddo
    return

    end subroutine chk_nks_tiny


    subroutine calc_ncphy_pp(io)
    integer,intent(in)::io

    integer i
    real(q) maxerr

    if(io>0)then
        call output_matrices('ncp-unsym',wh%nc_phy,wh%na2112,wh%num_imp, &
                &wh%na2_imp,io,1)
    endif
    wh%copy=wh%nc_phy
    call hm_expand_all_sym(wh%nc_phy,wh%hm,.true.,.true.)
    maxerr=maxval(abs(wh%nc_phy-wh%copy))
    wh%symerr=max(wh%symerr,maxerr)
    if(io>0)then
        call output_matrices('ncp-sym',wh%nc_phy,wh%na2112,wh%num_imp, &
                &wh%na2_imp,io,1)
        write(io,'(" max error due to symmetrization = ",f12.6)')maxerr
    endif
    call renorm_ncphy(io)
    if(io>0)then
        call output_matrices('ncp-renorm',wh%nc_phy,wh%na2112,wh%num_imp, &
                &wh%na2_imp,io,1)
    endif
    return

    end subroutine calc_ncphy_pp


    ! for wh%h1e, etc.
    subroutine calc_herm_matrices_pp(matrices,sname,mb,ltrans,io,mode)
    complex(q),intent(inout)::matrices(wh%na2112)
    character,intent(in)::sname*3
    type(matrix_basis),intent(in)::mb
    logical,intent(in)::ltrans
    integer,intent(in)::io,mode

    real(q) maxerr

    if(io>0)then
        call output_matrices(sname//'-unsym',matrices,wh%na2112,wh%num_imp, &
                &wh%na2_imp,io,mode)
    endif
    wh%copy=matrices
    call symm_dm_across_atoms(matrices)
    call hm_expand_all_sym(matrices,mb,ltrans,.true.)
    maxerr=maxval(abs(matrices-wh%copy))
    wh%symerr=max(wh%symerr,maxerr)
    if(io>0)then
        call output_matrices(sname//'-sym',matrices,wh%na2112,wh%num_imp, &
                &wh%na2_imp,io,mode)
        write(io,'(" max error due to symmetrization = ",f12.6)')maxerr
    endif
    return

    end subroutine calc_herm_matrices_pp


    subroutine calc_r01_pp(io)
    integer,intent(in)::io

    integer i

    ! r0 -> r
    forall(i=1:wh%num_imp)wh%co(i)%r=matmul(transpose(wh%co(i)%isimix), &
            &wh%co(i)%r0)

    if(io>0)then
        call output_matrices('r0-out',wh%r0,wh%na2112,wh%num_imp, &
                &wh%na2_imp,io,0)
    endif

    call  calc_r_pp(io,'r-out')

    ! r-> z=r^\dagger r
    do i=1,wh%num_imp
        call annxb('c',wh%co(i)%r,wh%co(i)%r,wh%co(i)%z, &
                &wh%na2_imp(i),wh%na2_imp(i))
    enddo

    if(io>0)then
        call output_matrices('z-out-sym',wh%z,wh%na2112,wh%num_imp, &
                &wh%na2_imp,io,0)
    endif

    call chk_eigens_matrix_list('z',wh%z,wh%na2112,wh%num_imp, &
            &wh%na2_imp,io,0)
    return

    end subroutine calc_r01_pp


    subroutine calc_r_pp(io,sname)
    integer,intent(in)::io
    character,intent(in)::sname*5

    integer i
    real(q) maxerr

    if(io>0)then
        call output_matrices(sname,wh%r,wh%na2112,wh%num_imp, &
                &wh%na2_imp,io,0)
    endif
    wh%copy=wh%r
    call hm_expand_all_general(wh%r,wh%r_coef,wh%hm_r,0,.false.)
    maxerr=maxval(abs(wh%r-wh%copy))
    wh%symerr=max(wh%symerr,maxerr)
    if(io>0)then
        call output_matrices(sname//'-sym',wh%r,wh%na2112,wh%num_imp, &
                &wh%na2_imp,io,0)
        write(io,'(" max error due to symmetrization = ",f12.6)')maxerr
    endif
    return

    end subroutine calc_r_pp


    subroutine calc_la1_pp(io,sname)
    integer,intent(in)::io
    character,intent(in)::sname*7

    integer i
    real(q) maxerr

    if(io>0)then
        call output_matrices(sname,wh%la1,wh%na2112,wh%num_imp, &
                &wh%na2_imp,io,0)
    endif
    wh%copy=wh%la1
    call hm_expand_all_herm(wh%la1,wh%la1_coef,wh%hm_l,0,.false.)
    maxerr=maxval(abs(wh%la1-wh%copy))
    wh%symerr=max(wh%symerr,maxerr) 
    if(io>0)then
        call output_matrices(sname//'-sym',wh%la1,wh%na2112,wh%num_imp, &
                &wh%na2_imp,io,0)
        write(io,'(" max error due to symmetrization = ",f12.6)')maxerr
    endif
    return

    end subroutine calc_la1_pp


    subroutine calc_lambdac(io)
    integer,intent(in) :: io

    integer i

    do i=1,wh%num_imp
        call calc_co_lambdac(wh%co(i))
    enddo

    if(io>0)then
        call output_matrices('la2',wh%la2,wh%na2112,wh%num_imp, &
                &wh%na2_imp,io,0)
    endif
    return

    end subroutine calc_lambdac


    !*************************************************************************
    subroutine calc_da0_pp(io)
    integer,intent(in)::io

    integer i
    real(q) maxerr

    if(io>0)then
        call output_matrices('d0-unsym',wh%d0,wh%na2112,wh%num_imp, &
                &wh%na2_imp,io,0)
    endif
    wh%copy=wh%d0
    call symm_dm_across_atoms(wh%d0)
    call hm_expand_all_sym(wh%d0,wh%hm_r,.false.,.false.)
    maxerr=maxval(abs(wh%d0-wh%copy))
    wh%symerr=max(wh%symerr,maxerr)
    if(io>0)then
        call output_matrices('d0-sym',wh%d0,wh%na2112,wh%num_imp, &
                &wh%na2_imp,io,0)
        write(io,'(" max error due to symmetrization = ",f12.6)')maxerr
    endif
    return

    end subroutine calc_da0_pp


    !*************************************************************************
    subroutine calc_da(io)
    integer,intent(in)::io

    call calc_da0_pp(io)
    call d0_to_d()
    if(io>0)then
       call output_matrices('d-sym',wh%d,wh%na2112,wh%num_imp, &
            &wh%na2_imp,io,0)
    endif
    return

    end subroutine calc_da


    !*************************************************************************
    !< mode>0: get c; mode<0: get a; mode=0: symmetrization (get c then a)
    !*************************************************************************
    subroutine hm_expand_all_herm(a,c,mb,mode,ltrans)
    type(matrix_basis),intent(in)::mb
    complex(q),target,intent(inout)::a(wh%na2112)
    real(q),intent(inout)::c(mb%dimhst)
    integer,intent(in)::mode
    logical,intent(in)::ltrans

    integer i,ibase,imap,abase,hsbase,dim_hs,na2,na22,ibase_list(wh%num_imp)
    complex(q),pointer::p_a(:,:),p_hs(:,:,:)
    complex(q) c1(mb%dimhsmax)

    ibase=0
    abase=0
    hsbase=0
    do i=1,wh%num_imp
        na2=wh%na2_imp(i)
        na22=na2**2
        dim_hs=mb%dim_hs(i)
        imap=wh%imap_imp(i)
        if(imap==i)then
            ibase_list(i)=ibase
        else
            ibase_list(i)=ibase_list(imap)
        endif
        p_a(1:na2,1:na2)=>a(abase+1:abase+na22)
        if(dim_hs>0)then
            p_hs(1:na2,1:na2,1:dim_hs)=>mb%hs(hsbase+1:hsbase+na22*dim_hs)
            if(mode<=0)then
                c1(1:dim_hs)=c(ibase_list(i)+1:ibase_list(i)+dim_hs)
            endif
            call get_hm_expand(p_a,p_hs,na2,dim_hs,c1(1:dim_hs), &
                    &mode,ltrans,.true.)
            if(imap==i)then
                if(mode>=0)then
                    c(ibase+1:ibase+dim_hs)=c1(1:dim_hs)
                endif
                ibase=ibase+dim_hs
            endif
        endif
        abase=abase+na22
        hsbase=hsbase+na22*dim_hs
    enddo
    return

    end subroutine hm_expand_all_herm


    subroutine hm_expand_all_general(a,c,mb,mode,ltrans)
    type(matrix_basis),intent(in)::mb
    complex(q),target,intent(inout)::a(wh%na2112)
    complex(q),intent(inout)::c(mb%dimhst)
    integer,intent(in)::mode
    logical,intent(in)::ltrans

    integer i,imap,ibase,abase,hsbase,nbase,dim_hs,na2,na22
    complex(q),pointer::p_a(:,:),p_hs(:,:,:)

    ibase=0
    abase=0
    hsbase=0
    do i=1,wh%num_imp
        na2=wh%na2_imp(i)
        na22=na2**2
        dim_hs=mb%dim_hs(i)
        if(dim_hs>0)then
            p_a(1:na2,1:na2)=>a(abase+1:abase+na22)
            p_hs(1:na2,1:na2,1:dim_hs)=>mb%hs(hsbase+1:hsbase+na22*dim_hs)
            imap=wh%imap_imp(i)
            if(imap==i)then
                call get_hm_expand(p_a,p_hs,na2,dim_hs, &
                        &c(ibase+1:ibase+dim_hs),mode,ltrans,.false.)
                ibase=ibase+dim_hs
            elseif(mode<=0)then
                nbase=sum(wh%na2_imp(1:imap-1)**2)
                a(abase+1:abase+na22)=a(nbase+1:nbase+na22)
            endif
        elseif(mode<0)then
            a(abase+1:abase+na22)=0
        endif
        abase=abase+na22
        hsbase=hsbase+na22*dim_hs
    enddo
    return

    end subroutine hm_expand_all_general


    subroutine hm_expand_all_sym(a,mb,ltrans,lherm)
    complex(q),target,intent(inout)::a(wh%na2112)
    type(matrix_basis),intent(in)::mb
    logical,intent(in)::ltrans,lherm

    integer i,abase,hsbase,dim_hs,na2,na22
    complex(q),pointer::p_a(:,:),p_hs(:,:,:)

    abase=0
    hsbase=0
    do i=1,wh%num_imp
        na2=wh%na2_imp(i)
        na22=na2**2
        dim_hs=mb%dim_hs(i)
        if(dim_hs>0)then
            p_a(1:na2,1:na2)=>a(abase+1:abase+na22)
            p_hs(1:na2,1:na2,1:dim_hs)=>mb%hs(hsbase+1:hsbase+na22*dim_hs)
            call get_hm_expand(p_a,p_hs,na2,dim_hs,ltrans,lherm)
        endif
        abase=abase+na22
        hsbase=hsbase+na22*dim_hs
    enddo
    return

    end subroutine hm_expand_all_sym


    !*************************************************************************
    ! renormalize nc_phy
    !*************************************************************************
    subroutine renorm_ncphy(io)
    integer,intent(in)::io

    integer i
    real(q) sum_ks,sum_phy,sum_var

    do i=1,wh%num_imp
        sum_ks =trace_a(wh%co(i)%nks   , wh%co(i)%dim2)
        sum_phy=trace_a(wh%co(i)%nc_phy, wh%co(i)%dim2)
        sum_var=trace_a(wh%co(i)%nc_var, wh%co(i)%dim2)
        if(io>0)then
            write(io,'(" imp=",i3," sum_phy-sum_var=",f16.8)') &
                    &i,sum_phy-sum_var
            write(io,'(7x       ," sum_phy-sum_ks =",f16.8, &
                    &" would be renormalized!")')sum_phy-sum_ks
        endif
        if(sum_phy<1.e-16_q)then
            wh%co(i)%nc_phy = 0._q
        else
            wh%co(i)%nc_phy=wh%co(i)%nc_phy/sum_phy*sum_ks ! renormalized
        endif
    enddo
    return

    end subroutine renorm_ncphy


    !*************************************************************************
    subroutine d0_to_d()
    integer i

    do i=1,wh%num_imp
        call co_d0_to_d(wh%co(i))
    enddo
    return

    end subroutine d0_to_d


    subroutine gh5_create_impurity_groups(fid)
    integer(hid_t),intent(in)::fid

    integer i

    if(fid<0)return
    do i = 1, wh%num_imp
        call gh5_create_group("/IMPURITY_"//trim(int_to_str(i)),f_id)
    enddo
    return

    end subroutine gh5_create_impurity_groups


    subroutine calc_et1_list()
    integer i

    do i=1,wh%num_imp
        wh%et1(i)=sum(wh%co(i)%h1e*wh%co(i)%nc_phy)
    enddo
    return

    end subroutine calc_et1_list


    subroutine calc_eu2_hf_list()
    integer i

    do i=1,wh%num_imp
        call calc_co_eu2_hf(wh%co(i), wh%eu2(i))
    enddo
    return

    end subroutine calc_eu2_hf_list


    subroutine output_energies_wh(io)
    integer,intent(in)::io

    integer i

    write(io,'(" impurity-wise interaction energy:")')
    write(io,'(4x,5f12.5)')(wh%eu2(i), i=1,wh%num_imp)
    write(io,'(" total impurity interaction energy = ",f0.7)') &
            &sum(wh%eu2)
    write(io,'(" impurity-wise one-body part energy:")')
    write(io,'(4x,5f12.5)')(wh%et1(i), i=1,wh%num_imp)
    write(io,'(" total impurity one-body part energy = ",f0.7)') &
            &sum(wh%et1)
    write(io,'(" total lambda-dc energy = ",f0.7)') wh%edcla1
    return

    end subroutine output_energies_wh



end module warehouse
