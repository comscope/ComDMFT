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

module gspci
    ! sparse-matrix based CI, with impurity-bath basis.
    use gprec
    use gutil
    use gprimme
    use sparse
    use ghdf5_base
    implicit none

    type gspci_mem
        integer::norb,norb2,nstates=0,nval_bot,nval_top
        type(dcoo_matrix),pointer::mcoo(:,:),nvcoo(:,:),npcoo(:,:)
        integer,pointer::bs(:),ibs(:),idx(:)
    end type gspci_mem

    integer,allocatable::mem_binom(:,:)
    type(gspci_mem)::dmem
    type(zcsr_matrix)::ucsr

    complex(q),pointer::p_d(:,:),p_la2(:,:)

    private
    public::dmem,gspci_gh5exe,av1_gspci

    contains


    subroutine gspci_gh5exe()
    integer i,ierr,norb,norb2,nval_bot,nval_top,nv
    integer,allocatable::m_struct(:,:)
    character*32 str
    real(q) etot,de
    logical lexist
    complex(q),allocatable::h1e(:,:),daalpha(:,:),lambdac(:,:),v2e(:,:,:,:), &
            &dm(:,:)
    complex(q),pointer::v(:)=>null()

    call get_command_argument(1, value=str, status=ierr)
    if(ierr/=0)then
        write(0,'(" Error in gspci_gh5start: missing inline argument!")')
        stop 1 
    endif
    read(str,*)i
    call gh5_init()
    call gh5_open_r('EMBED_HAMIL_'//trim(int_to_str(i))//'.h5',f_id)
    call gh5_read(norb,'/na2',f_id)
    call gh5_read(nval_bot,'/nval_bot',f_id)
    call gh5_read(nval_top,'/nval_top',f_id)
    allocate(m_struct(norb,norb))
    call gh5_read(m_struct,norb,norb,'/sigma_struct',f_id)
    norb2=norb*2
    allocate(h1e(norb,norb),daalpha(norb,norb),lambdac(norb,norb), &
            &v2e(norb,norb,norb,norb),dm(norb2,norb2))
    call gh5_read(h1e,norb,norb,'/H1E',f_id)
    call gh5_read(daalpha,norb,norb,'/D',f_id)
    call gh5_read(lambdac,norb,norb,'/LAMBDA',f_id)
    call gh5_read(v2e,norb,norb,norb,norb,'/V2E',f_id)
    call gh5_close(f_id)

    call init_gspci_mem(norb,nval_bot,nval_top,m_struct)

    ! Check previous solution vector
    inquire(file='EMBED_HAMIL_RES_'//trim(int_to_str(i))//'.h5',exist=lexist)
    if(lexist)then
        allocate(v(dmem%nstates))
        call gh5_open_r('EMBED_HAMIL_RES_'//trim(int_to_str(i))//'.h5',f_id)
        call gh5_read(nv,'/dimv',f_id)
        if(nv==dmem%nstates)then
            allocate(v(nv))
            call gh5_read(v,nv,'/evec',f_id)
        endif
        call gh5_close(f_id)
    endif

    call solve_hembed_spci_drive(h1e,daalpha,lambdac,v2e,v,dm,etot)
    
    de=trace_a(lambdac,norb) 
    etot=etot-de

    call gh5_open_w('EMBED_HAMIL_RES_'//trim(int_to_str(i))//'.h5',f_id)
    call gh5_write(etot,'/emol',f_id)
    call gh5_write(dm,norb2,norb2,'/DM',f_id)
    call gh5_write(dmem%nstates,'/dimv',f_id)
    call gh5_write(v,dmem%nstates,'/evec',f_id)
    call gh5_close(f_id)
    call gh5_end()

    deallocate(m_struct,h1e,daalpha,lambdac,v2e,dm,v)
    nullify(v)
    return

    end subroutine gspci_gh5exe


    subroutine gspci_init(n)
    integer,intent(in)::n

    allocate(mem_binom(n,0:n))
    call set_binoms(n,mem_binom)
    return

    end subroutine gspci_init


    subroutine init_gspci_mem(norb,nval_bot,nval_top,m_struct)
    integer,intent(in)::norb,nval_bot,nval_top,m_struct(norb,norb)
    
    call gspci_init(norb)
    call set_full_fock_states_spci(dmem,norb,nval_bot,nval_top)
    call set_mncoo_spci(dmem,norb,nval_bot,nval_top,m_struct)
    return

    end subroutine init_gspci_mem


    subroutine set_full_fock_states_spci(dmem,norb,nval_bot,nval_top)
    type(gspci_mem),intent(inout)::dmem
    integer,intent(in)::norb,nval_bot,nval_top

    integer i

    dmem%norb=norb; dmem%norb2=norb*2
    dmem%nval_bot=nval_bot
    dmem%nval_top=nval_top
    call set_fock_state_indices(norb,mem_binom,dmem%idx,dmem%bs, &
            &dmem%ibs)

    dmem%nstates=0
    do i=nval_bot,nval_top
        dmem%nstates=dmem%nstates+ &
                &(dmem%idx(i+1)-dmem%idx(i))**2
    enddo
    return

    end subroutine set_full_fock_states_spci


    subroutine solve_hembed_spci_drive(h1e,daalpha,lab,v2e,v, &
            &dm,etot)
    complex(q),target,intent(in)::h1e(dmem%norb,dmem%norb), &
            &daalpha(dmem%norb,dmem%norb), &
            &lab(dmem%norb,dmem%norb), &
            &v2e(dmem%norb,dmem%norb,dmem%norb,dmem%norb)
    complex(q),pointer::v(:)
    complex(q),intent(out)::dm(dmem%norb2,dmem%norb2)
    real(q),intent(out)::etot

    integer method
    real(q) w(dmem%norb2)

    if(ucsr%nrow==0)then
        call calc_ucsr_spci(dmem,ucsr,h1e,v2e,dmem%norb)
    endif

    p_d=>daalpha
    p_la2=>lab

    if(dmem%norb2<=8)then
        method=0
    else
        method=1
    endif

    if(.not.associated(v))then
        allocate(v(dmem%nstates))
        w(1)=0
    else
        w(1)=1
    endif

    call solve_hembed_spci(v,w(1:2),method=method)
    etot=w(1)
    call calc_dm_spci(v,dm)
    nullify(p_d,p_la2)
    return

    end subroutine solve_hembed_spci_drive


    subroutine solve_hembed_spci(v,w,method)
    complex(q),intent(inout)::v(dmem%nstates)
    integer,intent(in),optional::method
    real(q),intent(inout)::w(2)
    integer::i_method=1
    external::av_gspci

    if(present(method))then
        i_method=method
    endif

    if(i_method==0)then
        call lapack_diag_spci(v,dmem%nstates,w)
    else
#ifdef debug_mode
        call primme_diag_chk(dmem%nstates,10,av_gspci)
#endif
        call primme_diag(v,dmem%nstates,w,av_gspci)
    endif
    return

    end subroutine solve_hembed_spci


    !*************************************************************************
    subroutine lapack_diag_spci(v,n,w)
    integer,intent(in) :: n
    complex(q),intent(out) :: v(n)

    real(q),intent(out) :: w(2)

    integer i
    real(q) w_(n)
    complex(q) z(n,n),v1(n)

    z=0
    call setup_hdns_spci(z,n)
    call hermev('v','l',z,w_,n)
    v=z(:,1)
    w=w_(:2)
    return

    end subroutine lapack_diag_spci


    subroutine set_mncoo_spci(dmem,norb,nbot,ntop,m_struct)
    type(gspci_mem),intent(inout)::dmem
    integer,intent(in)::norb,nbot,ntop,m_struct(norb,norb)

    integer i,j

    allocate(dmem%mcoo(norb,norb),dmem%nvcoo(norb,norb),dmem%npcoo(norb,norb))
    do i=1,norb
        do j=1,norb
            if(m_struct(i,j)<=0)cycle

            ! c_j^\dagger f_i, note act to left
            dmem%mcoo(j,i)%nrow=dmem%nstates
            dmem%mcoo(j,i)%ncol=dmem%nstates
            call set_dcoo_ij_spci(dmem%mcoo(j,i),j,i+norb,.false., &
                    &.true.,norb,nbot,ntop,dmem%idx,dmem%bs,dmem%ibs)

            if(i<j)cycle ! hermitian
            ! f_j f_i^\dagger, note act to left
            dmem%nvcoo(j,i)%nrow=dmem%nstates
            dmem%nvcoo(j,i)%ncol=dmem%nstates
            call set_dcoo_ij_spci(dmem%nvcoo(j,i),j+norb,i+norb,.true., &
                    &.false.,norb,nbot,ntop,dmem%idx,dmem%bs,dmem%ibs)

            ! c_i^\dagger c_j  
            dmem%npcoo(i,j)%nrow=dmem%nstates
            dmem%npcoo(i,j)%ncol=dmem%nstates
            call set_dcoo_ij_spci(dmem%npcoo(i,j),i,j,.false., &
                    &.true.,norb,nbot,ntop,dmem%idx,dmem%bs,dmem%ibs)
        enddo
    enddo
    return

    end subroutine set_mncoo_spci


    subroutine set_dcoo_ij_spci(dcoo,icf,jcf,iact,jact,norb,nbot,ntop, &
            &idx,bs,ibs)
    type(dcoo_matrix),intent(inout)::dcoo
    logical,intent(in)::iact,jact
    integer,intent(in)::icf,jcf,norb,nbot,ntop,idx(0:),bs(:),ibs(0:)

    integer irun,ival,istate,i1,i2,jstate,nnz,nbase,dval, &
            &itmp,isgn,mask1,ibs1,ibs2,nfs

    dval=0
    if(icf<=norb)then
        if(iact)then
            dval=dval+1
        else
            dval=dval-1
        endif
    endif
    if(jcf<=norb)then
        if(jact)then
            dval=dval+1
        else
            dval=dval-1
        endif
    endif
    mask1=ishft(1,norb)-1
    do irun=1,2
        if(irun==2)then
            dcoo%nnz=nnz
            allocate(dcoo%i(nnz),dcoo%j(nnz),dcoo%a(nnz))
        endif
        nnz=0; nbase=0; istate=0
        do ival=nbot,ntop
            if(dval==1)then
                if(ival==ntop)then
                    exit
                elseif(ival==nbot)then
                    nbase=(idx(nbot+1)-idx(nbot))**2
                endif
            elseif(dval==-1.and.ival==nbot)then
                istate=istate+(idx(nbot+1)-idx(nbot))**2
                cycle
            elseif(abs(dval)>1)then
                write(0,'(" Error in set_dcoo_ij_spci!")')
                stop 2
            endif
            nfs=idx(ival+dval+1)-idx(ival+dval)
            do i1=idx(ival),idx(ival+1)-1
                do i2=idx(norb-ival),idx(norb-ival+1)-1
                    istate=istate+1
                    itmp=ior(bs(i1),ishft(bs(i2),norb))
                    isgn=1
                    call act_state(itmp,icf-1,iact,isgn)
                    if(isgn==0)cycle
                    call act_state(itmp,jcf-1,jact,isgn)
                    if(isgn==0)cycle
                    nnz=nnz+1
                    if(irun==1)cycle
                    ibs1=ibs(iand(itmp,mask1))
                    ibs2=ibs(ishft(itmp,-norb))
                    jstate=nbase+(ibs1-idx(ival+dval))*nfs+ &
                            &ibs2-idx(norb-ival-dval)+1
                    dcoo%i(nnz)=istate
                    dcoo%j(nnz)=jstate
                    dcoo%a(nnz)=isgn
                enddo
            enddo
            nbase=nbase+nfs**2
        enddo
    enddo
    return

    end subroutine set_dcoo_ij_spci


    subroutine calc_ucsr_spci(dmem,zcsr,h1e,v2e,norb)
    type(gspci_mem),intent(in)::dmem
    type(zcsr_matrix),intent(inout)::zcsr
    integer,intent(in)::norb
    complex(q),intent(in)::h1e(norb,norb),v2e(norb,norb,norb,norb)

    complex(q) z_row(maxval(dmem%idx(dmem%nval_bot+1:dmem%nval_top+1) &
            &-dmem%idx(dmem%nval_bot:dmem%nval_top)))
    integer nstates,nnz,irun,ival,istate,i,jstate,j,p,q,r,s, &
            &itmp(4),isgn(4)

    nstates=dmem%idx(dmem%nval_top+1)-dmem%idx(dmem%nval_bot)
    zcsr%nrow=nstates
    zcsr%ncol=nstates
    allocate(zcsr%i(nstates+1)); zcsr%i(1)=1

    do irun=1,2
        if(irun==2)then
            allocate(zcsr%j(nnz),zcsr%a(nnz))
        endif
        nnz=0; istate=0
        do ival=dmem%nval_bot,dmem%nval_top
            do i=dmem%idx(ival),dmem%idx(ival+1)-1
                istate=istate+1
                z_row=0
                do p=1,norb
                    ! <v|p^\dagger
                    isgn(1)=1
                    itmp(1)=dmem%bs(i)
                    call act_state(itmp(1),p-1,.false.,isgn(1))
                    if(isgn(1)==0)cycle

                    ! one-body part
                    do q=1,norb
                        if(abs(h1e(p,q))<1.d-10)cycle
                        ! <v|p^\dagger q
                        isgn(2)=isgn(1)
                        itmp(2)=itmp(1)
                        call act_state(itmp(2),q-1,.true.,isgn(2))
                        if(isgn(2)==0)cycle
                        jstate=dmem%ibs(itmp(2))
                        if(jstate>i)cycle
                        jstate=jstate-dmem%idx(ival)+1
                        z_row(jstate)=z_row(jstate)+isgn(2)*h1e(p,q)
                    enddo

                    ! two-body
                    do q=1,p ! take care of factor 1/2
                        ! <v|p^\dagger q^\dagger
                        isgn(2)=isgn(1)
                        itmp(2)=itmp(1)
                        call act_state(itmp(2),q-1,.false.,isgn(2))
                        if(isgn(2)==0)cycle
                        do r=1,norb
                            ! <v|p^\dagger q^\dagger r 
                            isgn(3)=isgn(2)
                            itmp(3)=itmp(2)
                            call act_state(itmp(3),r-1,.true.,isgn(3))
                            if(isgn(3)==0)cycle
                            do s=1,norb
                                if(abs(v2e(p,s,q,r))<1.d-10)cycle
                                ! <v|p^\dagger q^\dagger r s
                                isgn(4)=isgn(3)
                                itmp(4)=itmp(3)
                                call act_state(itmp(4),s-1,.true.,isgn(4))
                                if(isgn(4)==0)cycle
                                jstate=dmem%ibs(itmp(4))
                                if(jstate>i)cycle
                                jstate=jstate-dmem%idx(ival)+1
                                z_row(jstate)=z_row(jstate)+isgn(4)* &
                                        &v2e(p,s,q,r)
                            enddo
                        enddo
                    enddo
                enddo
                if(irun==1)then
                    nnz=nnz+count(abs(z_row(1:i-dmem%idx(ival)+1))>1.d-10)
                    zcsr%i(istate+1)=nnz+1
                    cycle
                else
                    do j=dmem%idx(ival),i
                        jstate=j-dmem%idx(ival)+1
                        if(abs(z_row(jstate))<=1.d-10)cycle
                        nnz=nnz+1
                        zcsr%j(nnz)=j-dmem%idx(dmem%nval_bot)+1
                        zcsr%a(nnz)=z_row(jstate)
                    enddo
                endif
            enddo
        enddo
    enddo
    return

    end subroutine calc_ucsr_spci


    subroutine setup_hdns_spci(a,n)
    integer,intent(in)::n
    complex(q),intent(out)::a(n,n)

    integer i,j,inz,r,s,nfk,nbase,i1,j1,i1_,i2,ival

    a=0

    ! D_ij
    do i=1,dmem%norb
        do j=1,dmem%norb
            if(abs(p_d(i,j))<1.d-10)cycle
            do inz=1,dmem%mcoo(j,i)%nnz
                r=dmem%mcoo(j,i)%i(inz)
                s=dmem%mcoo(j,i)%j(inz)
                a(r,s)=a(r,s)+dmem%mcoo(j,i)%a(inz)*p_d(i,j)
                a(s,r)=a(s,r)+dmem%mcoo(j,i)%a(inz)*conjg(p_d(i,j))
            enddo
        enddo
    enddo

    ! lambda_c (la2)
    do i=1,dmem%norb
        do j=1,i
            if(abs(p_la2(i,j))<1.d-10)cycle
            do inz=1,dmem%nvcoo(j,i)%nnz
                r=dmem%nvcoo(j,i)%i(inz)
                s=dmem%nvcoo(j,i)%j(inz)
                a(r,s)=a(r,s)+dmem%nvcoo(j,i)%a(inz)*p_la2(i,j)
                if(i==j)cycle
                a(s,r)=a(s,r)+dmem%nvcoo(j,i)%a(inz)*conjg(p_la2(i,j))
            enddo
        enddo
    enddo

    ! h_loc (ucsr)
    nbase=0
    do ival=dmem%nval_bot,dmem%nval_top
        nfk=dmem%idx(ival+1)-dmem%idx(ival)
        do i1=dmem%idx(ival),dmem%idx(ival+1)-1
            i1_=i1-dmem%idx(dmem%nval_bot)+1
            do inz=ucsr%i(i1_),ucsr%i(i1_+1)-1
                j1=ucsr%j(inz)+dmem%idx(dmem%nval_bot)-1
                do i2=1,nfk
                    r=nbase+(i1-dmem%idx(ival))*nfk+i2
                    s=nbase+(j1-dmem%idx(ival))*nfk+i2
                    a(r,s)=a(r,s)+ucsr%a(inz)
                    if(r==s)cycle
                    a(s,r)=a(s,r)+conjg(ucsr%a(inz))
                enddo
            enddo
        enddo
        nbase=nbase+nfk**2
    enddo 
    return

    end subroutine setup_hdns_spci


    subroutine act_ucsr_spci(v1,v2)
    complex(q),intent(in)::v1(*)
    complex(q),intent(inout)::v2(*)

    integer nbase,ival,nfk,i1,i1_,i2,j1,inz,r,s

    nbase=0
    do ival=dmem%nval_bot,dmem%nval_top
        nfk=dmem%idx(ival+1)-dmem%idx(ival)
        do i1=dmem%idx(ival),dmem%idx(ival+1)-1
            i1_=i1-dmem%idx(dmem%nval_bot)+1
            do inz=ucsr%i(i1_),ucsr%i(i1_+1)-1
                j1=ucsr%j(inz)+dmem%idx(dmem%nval_bot)-1
                do i2=1,nfk
                    r=nbase+(i1-dmem%idx(ival))*nfk+i2
                    s=nbase+(j1-dmem%idx(ival))*nfk+i2
                    v2(r)=v2(r)+ucsr%a(inz)*v1(s)
                    if(r==s)cycle
                    v2(s)=v2(s)+conjg(ucsr%a(inz))*v1(r)
                enddo
            enddo
        enddo
        nbase=nbase+nfk**2
    enddo
    return

    end subroutine act_ucsr_spci


    subroutine av1_gspci(v1,v2)
    complex(q),intent(in)::v1(*)
    complex(q),intent(out)::v2(*)

    integer i,j

    ! D_ij
    do i=1,dmem%norb
        do j=1,dmem%norb
            if(abs(p_d(i,j))<1.d-10)cycle
            call zs_dcoomux (p_d(i,j),dmem%mcoo(j,i),v1,v2)
            call zs_dcoohmux(conjg(p_d(i,j)),dmem%mcoo(j,i),v1,v2)
        enddo
    enddo

    ! lambda_c (la2)
    do i=1,dmem%norb
        do j=1,i
            if(abs(p_la2(i,j))<1.d-10)cycle
            call zs_dcoomux (p_la2(i,j),dmem%nvcoo(j,i),v1,v2)
            if(i==j)cycle
            call zs_dcoohmux(conjg(p_la2(i,j)),dmem%nvcoo(j,i),v1,v2)
        enddo
    enddo

    ! h_loc (ucsr)
    call act_ucsr_spci(v1,v2)
    return

    end subroutine av1_gspci   


    subroutine calc_dm_spci(v,dm)
    complex(q),intent(in)::v(dmem%nstates) 
    complex(q),intent(out)::dm(dmem%norb2,dmem%norb2)

    integer i,j,i_,j_

    ! c_i^\dagger c_j
    do i=1,dmem%norb
        do j=1,i
            call zvh_dcoo_zv(dmem%npcoo(i,j),v,dm(i,j))
            if(i==j)cycle
            dm(j,i)=conjg(dm(i,j))
        enddo
    enddo

    ! f_j f_i^\dagger = \delta_i,j - f_i^\dagger f_j
    do i=1,dmem%norb
        i_=i+dmem%norb
        do j=1,i
            j_=j+dmem%norb
            call zvh_dcoo_zv(dmem%nvcoo(j,i),v,dm(i_,j_))
            if(i_==j_)then
                dm(i_,j_)=1-dm(i_,j_)
            else
                dm(i_,j_)=-dm(i_,j_)
                dm(j_,i_)=conjg(dm(i_,j_))
            endif
        enddo
    enddo

    ! c_j^\dagger f_i
    do i=1,dmem%norb
        i_=i+dmem%norb
        do j=1,dmem%norb
            call zvh_dcoo_zv(dmem%mcoo(j,i),v,dm(j,i_))
            dm(i_,j)=conjg(dm(j,i_))
        enddo
    enddo 
    return

    end subroutine calc_dm_spci


end module gspci


subroutine av_gspci(v1,v2,k,primme)
    use gprec
    use gspci
    use gconstant
    implicit none
    integer,intent(in)::k
    integer(8),intent(in)::primme
    complex(q),intent(in)::v1(*)
    complex(q),intent(out)::v2(*)

    integer n,i

    n=dmem%nstates
    do i=1,k
        v2(n*(i-1)+1:n*i)=0
        call av1_gspci(v1(n*(i-1)+1),v2(n*(i-1)+1))
    enddo
    return

end subroutine av_gspci
