#ifdef ncf_csr
subroutine setup_ncab_op()
    use gprec
    use gspci
    use sparse
    use gutil
    implicit none
    integer ia1,ia2,isum,istate,nbase,ival,itmp1,isgn,nfs,nfs_l,i1,i2,ibs
    integer jdx(dmem%nstates)
    real(q) aval(dmem%nstates)

    allocate(dmem%ncabop(dmem%norb,dmem%norb))
    do ia1=1,dmem%norb; do ia2=1,dmem%norb
        dmem%ncabop(ia1,ia2)%nrow=dmem%nstates
        dmem%ncabop(ia1,ia2)%ncol=dmem%nstates
        allocate(dmem%ncabop(ia1,ia2)%i(dmem%nstates+1))
        dmem%ncabop(ia1,ia2)%i(1)=1
        nbase=0; istate=0; isum=1
        ! < state | c_ia1^\dagger c_ia2
        do ival=dmem%nval_bot,dmem%nval_top
            nfs=dmem%idx(ival+1)-dmem%idx(ival)
            nfs_l=dmem%idx_l(dmem%norb-ival+1)-dmem%idx_l(dmem%norb-ival)
            if(nfs_l<=0)cycle
            do i1=dmem%idx(ival),dmem%idx(ival+1)-1
                itmp1=dmem%bs(i1)
                isgn=1 
                call act_state(itmp1,ia1-1,.false.,isgn)
                if(isgn/=0)then
                    call act_state(itmp1,ia2-1,.true.,isgn)
                    if(isgn/=0)then
                        ibs=dmem%ibs(itmp1)
                        if(ibs<=0)isgn=0
                    endif
                endif
                do i2=dmem%idx_l(dmem%norb-ival), &
                        &dmem%idx_l(dmem%norb-ival+1)-1
                    istate=istate+1
                    if(isgn/=0)then
                        jdx(isum)=nbase+(ibs-dmem%idx(ival))*nfs_l+ &
                                &i2-dmem%idx_l(dmem%norb-ival)+1
                        aval(isum)=real(isgn,q)
                        isum=isum+1
                    endif
                    dmem%ncabop(ia1,ia2)%i(istate+1)=isum
                enddo
            enddo
            nbase=nbase+nfs*nfs_l
        enddo
        isum=isum-1
        allocate(dmem%ncabop(ia1,ia2)%j(isum),dmem%ncabop(ia1,ia2)%a(isum))
        dmem%ncabop(ia1,ia2)%j=jdx(1:isum)
        dmem%ncabop(ia1,ia2)%a=aval(1:isum)
    enddo; enddo
    return

end subroutine setup_ncab_op


subroutine setup_nfab_op()
    use gprec
    use gspci
    use sparse
    use gutil
    implicit none
    integer ia1,ia2,isum,istate,nbase,ival,itmp2,isgn,nfs,nfs_l,i1,i2,ibs
    integer jdx(dmem%nstates)
    real(q) aval(dmem%nstates)

    allocate(dmem%nfabop(dmem%norb,dmem%norb))
    do ia1=1,dmem%norb; do ia2=1,dmem%norb
        dmem%nfabop(ia1,ia2)%nrow=dmem%nstates
        dmem%nfabop(ia1,ia2)%ncol=dmem%nstates
        allocate(dmem%nfabop(ia1,ia2)%i(dmem%nstates+1))
        dmem%nfabop(ia1,ia2)%i(1)=1
        nbase=0; istate=0; isum=1
        ! < state | c_ia1^\dagger c_ia2
        do ival=dmem%nval_bot,dmem%nval_top
            nfs=dmem%idx(ival+1)-dmem%idx(ival)
            nfs_l=dmem%idx_l(dmem%norb-ival+1)-dmem%idx_l(dmem%norb-ival)
            if(nfs_l<=0)cycle
            do i1=dmem%idx(ival),dmem%idx(ival+1)-1
                do i2=dmem%idx_l(dmem%norb-ival), &
                        &dmem%idx_l(dmem%norb-ival+1)-1
                    istate=istate+1
                    itmp2=dmem%bs_l(i2)
                    isgn=1 
                    call act_state(itmp2,ia1-1,.false.,isgn)
                    if(isgn/=0)then
                        call act_state(itmp2,ia2-1,.true.,isgn)
                        if(isgn/=0)then
                            ibs=dmem%ibs_l(itmp2)
                            if(ibs<=0)isgn=0
                        endif
                    endif
                    if(isgn/=0)then
                        jdx(isum)=nbase+(i1-dmem%idx(ival))*nfs_l+ &
                                &ibs-dmem%idx_l(dmem%norb-ival)+1
                        aval(isum)=real(isgn,q)
                        isum=isum+1
                    endif
                    dmem%nfabop(ia1,ia2)%i(istate+1)=isum
                enddo
            enddo
            nbase=nbase+nfs*nfs_l
        enddo
        isum=isum-1
        allocate(dmem%nfabop(ia1,ia2)%j(isum),dmem%nfabop(ia1,ia2)%a(isum))
        dmem%nfabop(ia1,ia2)%j=jdx(1:isum)
        dmem%nfabop(ia1,ia2)%a=aval(1:isum)
    enddo; enddo
    return

end subroutine setup_nfab_op


subroutine chk_dm_spci()
    use gprec
    use gspci
    use sparse
    implicit none
    integer i,j,i_,j_
    complex(q) zes

    ! c_i^\dagger c_j
    do i=1,dmem%norb
        do j=1,dmem%norb
            call zvh_dcsr_zv(dmem%ncabop(i,j),dmem%v,zes)
            if(abs(zes-dmem%dm(i,j))>1.d-5)then
                write(0,'(" Error in nc(",i0,",",i0,"):(",f0.6,",",f0.6, &
                        &") vs (",f0.6,",",f0.6,")")')&
                        &i,j,dmem%dm(i,j),zes
            endif
        enddo
    enddo

    ! f_j f_i^\dagger = \delta_i,j - f_i^\dagger f_j
    do i=1,dmem%norb
        i_=i+dmem%norb
        do j=1,dmem%norb
            j_=j+dmem%norb
            call zvh_dcsr_zv(dmem%nfabop(i,j),dmem%v,zes)
            if(abs(zes-dmem%dm(i_,j_))>1.d-5)then
                write(0,'(" Error in nf(",i0,",",i0,"):(",f0.6,",",f0.6, &
                        &") vs (",f0.6,",",f0.6,")")')&
                        &i,j,dmem%dm(i_,j_),zes
            endif
        enddo
    enddo
    return

end subroutine chk_dm_spci
#endif ncf_csr


#ifdef constraint_s2
subroutine set_s2_spci()
    use gprec
    use gspci
    use gutil
    implicit none
    integer irun,ival,istate,i1,i2,jstate,nnz,nbase,nhalf, &
            &itmp,mask1,ibs1,ibs2,nfs,nfs_l,ia1,ia2,inz0, &
            &itmp1,itmp2
    real(q) sz1,sz2

    nhalf=dmem%norb/2
    dmem%s2op%nrow=dmem%nstates
    dmem%s2op%ncol=dmem%nstates
    allocate(dmem%s2op%i(dmem%nstates+1)); dmem%s2op%i(1)=1
    mask1=ishft(1,dmem%norb)-1
    do irun=1,2
        if(irun==2)then
            allocate(dmem%s2op%j(nnz),dmem%s2op%a(nnz))
        endif
        nnz=0; nbase=0; istate=0
        do ival=dmem%nval_bot,dmem%nval_top
            nfs=dmem%idx(ival+1)-dmem%idx(ival)
            nfs_l=dmem%idx_l(dmem%norb-ival+1)-dmem%idx_l(dmem%norb-ival)
            do i1=dmem%idx(ival),dmem%idx(ival+1)-1
                if(irun==2)then
                    call calc_state_sz(dmem%bs(i1),dmem%norb,sz1,0)
                endif
                do i2=dmem%idx_l(dmem%norb-ival), &
                        &dmem%idx_l(dmem%norb-ival+1)-1
                    if(irun==2)then
                        call calc_state_sz(dmem%bs_l(i2),dmem%norb,sz2,0)
                    endif
                    istate=istate+1
                    nnz=nnz+1 ! sz^2+sz
                    if(irun==2)then
                        inz0=nnz
                        dmem%s2op%j(nnz)=istate
                        dmem%s2op%a(nnz)=(sz1+sz2)**2+sz1+sz2
                    endif
                    itmp=ior(dmem%bs(i1),ishft(dmem%bs_l(i2),dmem%norb))
                    do ia1=1,dmem%norb
                        itmp1=itmp
                        ! s_ia1^+
                        if(btest(itmp1,2*ia1-1))then
                            itmp1=ibclr(itmp1,2*ia1-1)
                        else
                            cycle
                        endif
                        if(btest(itmp1,2*ia1-2))then
                            cycle
                        else
                            itmp1=ibset(itmp1,2*ia1-2)
                        endif
                        do ia2=1,dmem%norb
                            ! s_ia2^- s_ia1^+
                            itmp2=itmp1
                            if(ia2==ia1)then
                                if(irun==2)then
                                    dmem%s2op%a(inz0)=dmem%s2op%a(inz0)+1._q
                                endif
                                cycle
                            ! Only keep upper trianglar part
                            elseif(ia1<=nhalf.and.ia2>ia1.and.ia2<=nhalf)then
                                cycle
                            elseif(ia1>nhalf.and.ia2<=nhalf)then
                                cycle
                            elseif(ia1>nhalf.and.ia2>ia1.and.ia2>nhalf)then
                                cycle
                            else
                                if(btest(itmp2,2*ia2-2))then
                                    itmp2=ibclr(itmp2,2*ia2-2)
                                else
                                    cycle
                                endif
                                if(btest(itmp2,2*ia2-1))then
                                    cycle
                                else
                                    itmp2=ibset(itmp2,2*ia2-1)
                                endif
                            endif
                            nnz=nnz+1
                            if(irun==1)cycle
                            ibs1=dmem%ibs(iand(itmp2,mask1))
                            ibs2=dmem%ibs_l(ishft(itmp2,-dmem%norb))
                            jstate=nbase+(ibs1-dmem%idx(ival))*nfs_l+ &
                                    &ibs2-dmem%idx_l(dmem%norb-ival)+1
                            dmem%s2op%j(nnz)=jstate
                            dmem%s2op%a(nnz)=1._q
                        enddo
                    enddo
                    dmem%s2op%i(istate+1)=nnz+1
                enddo
            enddo
            nbase=nbase+nfs*nfs_l
        enddo
    enddo
    return

end subroutine set_s2_spci


subroutine chk_eval_s2(lstop)
    use gprec
    use gspci
    use sparse
    implicit none
    logical,intent(in)::lstop

    complex(q) eval

    eval=0
    call zvh_sydcsr_zv(dmem%s2op,dmem%v,eval)
    write(0,*) "<S^2> = ", eval
    if(lstop)then
        if(abs(eval)>1.d-8)then
            write(0,'(" Warning: <S^2> is not zero!")')
        endif
    endif
    return

end subroutine chk_eval_s2


subroutine av1_gspci_s2(v1,v2)
    use gprec
    use gspci
    use sparse
    implicit none
    complex(q),intent(in)::v1(*)
    complex(q),intent(inout)::v2(*)

    ! s2op
    if(abs(dmem%lambda_j2)>1.d-12)then
        call dcsr_sylamuzx_sk(dmem%lambda_j2,dmem%s2op,v1,v2)
    endif
    return

end subroutine av1_gspci_s2


subroutine add_hdns_spci_s2(a,n)
    use gprec
    use gspci
    implicit none
    integer,intent(in)::n
    complex(q),intent(inout)::a(n,n)

    integer i,inz,j

    ! s2, here we search for s=0 solution.
    if(abs(dmem%lambda_j2)>1.d-12)then
        do i=1,dmem%s2op%nrow
            do inz=dmem%s2op%i(i),dmem%s2op%i(i+1)-1
                j=dmem%s2op%j(inz)
                a(i,j)=a(i,j)+dmem%s2op%a(inz)*dmem%lambda_j2
                if(i==j)cycle
                a(j,i)=a(j,i)+dmem%s2op%a(inz)*dmem%lambda_j2
            enddo
        enddo
    endif
    return

end subroutine add_hdns_spci_s2
#endif constraint_s2


#ifdef constraint_dsym
subroutine chk_eval_dsym(lstop)
    use gprec
    use gutil
    use gspci
    use ghdf5_base
    use ghdf5
    implicit none
    logical lstop

    logical lexist
    integer ival,nfs,nfs_l,nbase,dimu,j
    complex(q) zes
    complex(q),allocatable::phi(:,:),phj(:,:),u(:,:)
    complex(q),external::zdotc

    ! dsym projector
    inquire(file='GLROT.h5', exist=lexist)
    if(.not.lexist)return
    call gh5_open_r('GLROT.h5',f_id)
    call gh5_read(dimu,'/IMPURITY_'//trim(int_to_str(dmem%imp))// &
            &'/dim_rot',f_id)
    if(dimu>0)then
        nbase=0
        zes=0
        do ival=dmem%nval_bot,dmem%nval_top
            nfs=dmem%idx(ival+1)-dmem%idx(ival)
            nfs_l=dmem%idx_l(dmem%norb-ival+1)-dmem%idx_l(dmem%norb-ival)
            if(nfs>1)then ! skip trivial case
                allocate(phi(nfs,nfs),phj(nfs,nfs),u(nfs,nfs))
                phj=0
                call phi_vec_to_mat(dmem%v(nbase+1:nbase+nfs*nfs_l), &
                        &phi,nfs,dmem%bs_l(dmem%idx_l(dmem%norb-ival): &
                        &dmem%idx_l(dmem%norb-ival+1)-1), &
                        &nfs_l,dmem%norb,dmem%ibs,dmem%idx(ival))
                do j=1,dimu
                    ! read jth rotation matrix
                    call gh5_read(u,nfs,nfs,'/IMPURITY_'// &
                            &trim(int_to_str(dmem%imp))//&
                            &'/valence_block_'//trim(int_to_str(ival))// &
                            &'/ROT_'//trim(int_to_str(j)), f_id)
                    call uhau(phi,u,nfs,nfs,uhau=phj)
                enddo
                zes=zes-zdotc(nfs*nfs,phj,1,phi,1)
                deallocate(phi,phj,u)
            elseif(nfs*nfs_l==1)then
                zes=zes-dmem%v(nbase+1)*conjg(dmem%v(nbase+1))*dimu
            endif
            nbase=nbase+nfs*nfs_l
        enddo
        zes=zes/dimu+1
        write(0,*) "<1-P_dsym> = ", zes
        if(lstop)then
            if(abs(zes)>1.d-8)then
                write(0,'(" Warning: <1-P_dsym> is not zero!")')
            endif
        endif
    endif
    call gh5_close(f_id)
    return


end subroutine chk_eval_dsym


! v2 = v2 + lambda_dsym/dim_{Pr}*\sum_{Pr}{(1-Pr)*v1}
subroutine av1_gspci_dsym(v1,v2)
    use gprec
    use gutil
    use gspci
    use ghdf5_base
    use ghdf5
    implicit none
    complex(q),intent(in)::v1(*)
    complex(q),intent(inout)::v2(*)

    integer ival,nfs,nfs_l,nbase,dimu,j
    complex(q),allocatable::phi(:,:),phj(:,:),u(:,:)

    ! dsym projector
    if(abs(dmem%lambda_dsym)<1.d-12)return
    call gh5_open_r('GLROT.h5',f_id)
    call gh5_read(dimu,'/IMPURITY_'//trim(int_to_str(dmem%imp))// &
            &'/dim_rot',f_id)
    if(dimu>0)then
        nbase=0
        do ival=dmem%nval_bot,dmem%nval_top
            nfs=dmem%idx(ival+1)-dmem%idx(ival)
            nfs_l=dmem%idx_l(dmem%norb-ival+1)-dmem%idx_l(dmem%norb-ival)
            if(nfs>1)then ! skip trivial case
                allocate(phi(nfs,nfs),phj(nfs,nfs),u(nfs,nfs))
                phj=0
                call phi_vec_to_mat(v1(nbase+1:nbase+nfs*nfs_l), &
                        &phi,nfs,dmem%bs_l(dmem%idx_l(dmem%norb-ival): &
                        &dmem%idx_l(dmem%norb-ival+1)-1), &
                        &nfs_l,dmem%norb,dmem%ibs,dmem%idx(ival))
                do j=1,dimu
                    ! read jth rotation matrix
                    call gh5_read(u,nfs,nfs,'/IMPURITY_'// &
                            &trim(int_to_str(dmem%imp))//&
                            &'/valence_block_'//trim(int_to_str(ival))// &
                            &'/ROT_'//trim(int_to_str(j)), f_id)
                    call uhau(phi,u,nfs,nfs,uhau=phj)
                enddo
                phj=-dmem%lambda_dsym/dimu*phj+dmem%lambda_dsym*phi
                call phi_mat_to_vec(v2(nbase+1:nbase+nfs*nfs_l), & 
                        &phj,nfs,dmem%bs_l(dmem%idx_l(dmem%norb-ival): &
                        &dmem%idx_l(dmem%norb-ival+1)-1), &
                        &nfs_l,dmem%norb,dmem%ibs,dmem%idx(ival))
                deallocate(phi,phj,u)
            endif
            nbase=nbase+nfs*nfs_l
        enddo
    endif
    call gh5_close(f_id)
    return


end subroutine av1_gspci_dsym
#endif constraint_dsym


#ifdef nc_csr
! In physical subspace.
subroutine setup_npcoo_op1()
    use gprec
    use gspci
    use sparse
    use gutil
    implicit none
    integer nstates,ia1,ia2,isum,nbase,ival,itmp1,isgn,i1,i2,ibs

    allocate(dmem%ncabop(dmem%norb,dmem%norb))

    nstates=dmem%idx(dmem%nval_top)-dmem%idx(dmem%nval_bot)
    do ia1=1,dmem%norb; do ia2=1,ia1
        call alloc_dcoo(dmem%npcoo(ia1,ia2),nstates,nstates,nstates)
        nbase=dmem%idx(dmem%nval_bot)-1
        isum=0
        ! < state | c_ia1^\dagger c_ia2
        do ival=dmem%nval_bot,dmem%nval_top
            do i1=dmem%idx(ival),dmem%idx(ival+1)-1
                itmp1=dmem%bs(i1)
                isgn=1 
                call act_state(itmp1,ia1-1,.false.,isgn)
                if(isgn==0)cycle
                call act_state(itmp1,ia2-1,.true.,isgn)
                if(isgn==0)cycle
                ibs=dmem%ibs(itmp1)
                isum=isum+1
                dmem%npcoo(ia1,ia2)%i(isum)=i1-nbase
                dmem%npcoo(ia1,ia2)%j(isum)=ibs-nbase
                dmem%npcoo(ia1,ia2)%a(isum)=real(isgn,q)
            enddo
        enddo
        dmem%npcoo(ia1,ia2)%nnz=isum
    enddo; enddo
    return

end subroutine setup_ncab_op
#endif nc_csr


subroutine chk_eval_j2_sumway(cm_vec,label)
    use gprec
    use gspci
    use sparse
    implicit none
    character(*),intent(in)::label
    complex(q),intent(in)::cm_vec(dmem%norb,dmem%norb,3)

    integer i
    complex(q) eval

    eval=0
    do i=1,3
        call calc_eval_gspci_j2op(cm_vec(:,:,i),dmem%v,eval)
    enddo
    write(0,*) label,"(sum_way) = ", eval
    return

end subroutine chk_eval_j2_sumway


subroutine av1_gspci_j2_sumway(lambda_j2,cm_vec,v1,v2)
    use gprec
    use gspci
    use sparse
    implicit none
    real(q),intent(in)::lambda_j2
    complex(q),intent(in)::cm_vec(dmem%norb,dmem%norb,3)
    complex(q),intent(in)::v1(*)
    complex(q),intent(inout)::v2(*)

    complex(q) cm(dmem%norb,dmem%norb)
    integer i

    if(abs(lambda_j2)>1.d-12)then
        do i=1,3
            cm=cm_vec(:,:,i)*sqrt(lambda_j2)
            call av1_gspci_j2op(cm,v1,v2)
        enddo
    endif
    return

end subroutine av1_gspci_j2_sumway


! Using equation J^2 = J-J+ + Jz^2 + Jz
subroutine av1_gspci_j2_npzway(lambda_j2,cm_n,cm_z,v1,v2)
    use gprec
    use gspci
    use sparse
    implicit none
    real(q),intent(in)::lambda_j2
    complex(q),intent(in)::cm_n(dmem%norb,dmem%norb),cm_z(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(*)
    complex(q),intent(inout)::v2(*)

    if(abs(lambda_j2)>1.d-12)then
        call av1_gspci_jz2_pjz(lambda_j2,cm_z,v1,v2)
        call av1_gspci_jnjp(lambda_j2,cm_n,v1,v2)
    endif
    return

end subroutine av1_gspci_j2_npzway


subroutine modify_m_struct(svec,lvec)
    use gprec
    use gspci
    implicit none
    complex(q),intent(in)::svec(dmem%norb,dmem%norb,3), &
            &lvec(dmem%norb,dmem%norb,3)

    integer i,j

    do i=1,dmem%norb; do j=1,dmem%norb
        if(dmem%m_struct(i,j)>0)cycle
        if(sum(abs(svec(i,j,:)))>1.d-12.or. &
                &sum(abs(svec(j,i,:)))>1.d-12.or. &
                &sum(abs(lvec(i,j,:)))>1.d-12.or. &
                &sum(abs(lvec(j,i,:)))>1.d-12)then
            dmem%m_struct(i,j)=101
        endif
    enddo; enddo
    return

end subroutine modify_m_struct


subroutine set_full_fock_states_l_spci()
    use gprec
    use gspci
    implicit none
    integer i

    if(dmem%norb_mott>0)then
        call set_dmem_mott()
    else
        dmem%bs_l=>dmem%bs; dmem%ibs_l=>dmem%ibs
        dmem%idx_l=>dmem%idx
    endif

    dmem%nstates=0
    do i=dmem%nval_bot,dmem%nval_top
        dmem%nstates=dmem%nstates+ &
                &(dmem%idx(i+1)-dmem%idx(i))* &
                &(dmem%idx_l(dmem%norb-i+1)-dmem%idx_l(dmem%norb-i))
    enddo
    return

end subroutine set_full_fock_states_l_spci


subroutine set_dmem_mott()
    use gprec
    use gspci
    use gutil
    implicit none
    integer i,n,nbs,isum,sume

    nbs=ishft(1,dmem%norb)
    allocate(dmem%idx_l(0:dmem%norb+1),dmem%bs_l(nbs),dmem%ibs_l(0:nbs-1))
    dmem%ibs_l=0
    dmem%idx_l(0)=1; isum=1
    do n=0,dmem%norb
        if(n<=dmem%norb-dmem%nelect_mott)then
            do i=dmem%idx(n),dmem%idx(n+1)-1
                call sum_empty_slots_fs(dmem%bs(i),dmem%norb_mott, &
                        &dmem%iorb_mott,sume)
                if(sume==dmem%nelect_mott)then
                    dmem%bs_l(isum)=dmem%bs(i)
                    dmem%ibs_l(dmem%bs(i))=isum
                    isum=isum+1
                endif
            enddo
        endif
        dmem%idx_l(n+1)=isum
    enddo
    return

end subroutine set_dmem_mott


subroutine solve_hembed_spci_drive()
    use gprec
    use gspci
    implicit none
    integer method
    real(q) w(dmem%norb2)

    if(dmem%norb2<=4)then
        method=0
    else
        method=1
    endif

    if(.not.associated(dmem%v))then
        allocate(dmem%v(dmem%nstates))
        w(1)=0
    else
        w(1)=1
    endif
    call solve_hembed_spci(w(1:2),method)
    dmem%etot=w(1)
    return

end subroutine solve_hembed_spci_drive


subroutine solve_hembed_spci(w,method)
    use gprec
    use gspci
    use gutil
    implicit none
    real(q),intent(inout)::w(2)
    integer method
    external::av_gspci

    if(method==0)then
        call lapack_diag_spci(dmem%v,dmem%nstates,w)
    else
        call primme_diag(dmem%v,dmem%nstates,w,av_gspci)
    endif
    return

end subroutine solve_hembed_spci


subroutine solve_hembed_spci_chk()
    use gspci
    use gutil
    implicit none

    integer nchk
    external::av_gspci

    nchk=min(10,dmem%nstates-1)
    call zprimme_diag_chk(dmem%nstates,nchk,av_gspci)

end subroutine solve_hembed_spci_chk

subroutine lapack_diag_spci(v,n,w)
    use gprec
    use gspci
    use gutil
    implicit none
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


subroutine set_mncoo_spci()
    use gprec
    use gspci
    implicit none
    integer i,j

    allocate(dmem%mcoo(dmem%norb,dmem%norb), &
            &dmem%nvcoo(dmem%norb,dmem%norb), &
            &dmem%npcoo(dmem%norb,dmem%norb))
    do i=1,dmem%norb
        do j=1,dmem%norb

            if(dmem%m_struct(i,j)<=0)cycle
            ! c_j^\dagger f_i, note act to left
            dmem%mcoo(j,i)%nrow=dmem%nstates
            dmem%mcoo(j,i)%ncol=dmem%nstates
            call set_dcoo_ij_spci(dmem%mcoo(j,i),j,i+dmem%norb,.false.,.true.)

            if(i<j)cycle ! hermitian

            ! f_j f_i^\dagger, note act to left
            dmem%nvcoo(j,i)%nrow=dmem%nstates
            dmem%nvcoo(j,i)%ncol=dmem%nstates
            call set_dcoo_ij_spci(dmem%nvcoo(j,i),j+dmem%norb,i+dmem%norb, &
                    &.true.,.false.)

            ! c_i^\dagger c_j  
            dmem%npcoo(i,j)%nrow=dmem%nstates
            dmem%npcoo(i,j)%ncol=dmem%nstates
            call set_dcoo_ij_spci(dmem%npcoo(i,j),i,j,.false.,.true.)
        enddo
    enddo
    return

end subroutine set_mncoo_spci


subroutine set_dcoo_ij_spci(dcoo,icf,jcf,iact,jact)
    use gprec
    use gspci
    use sparse
    use gutil
    implicit none
    type(dcoo_matrix),intent(inout)::dcoo
    logical,intent(in)::iact,jact
    integer,intent(in)::icf,jcf

    integer irun,ival,istate,i1,i2,jstate,nnz,nbase,dval, &
            &itmp,isgn,mask1,ibs1,ibs2,nfs,nfs_l

    dval=0
    if(icf<=dmem%norb)then
        if(iact)then
            dval=dval+1
        else
            dval=dval-1
        endif
    endif
    if(jcf<=dmem%norb)then
        if(jact)then
            dval=dval+1
        else
            dval=dval-1
        endif
    endif
    mask1=ishft(1,dmem%norb)-1
    do irun=1,2
        if(irun==2)then
            dcoo%nnz=nnz
            allocate(dcoo%i(nnz),dcoo%j(nnz),dcoo%a(nnz))
        endif
        nnz=0; nbase=0; istate=0
        do ival=dmem%nval_bot,dmem%nval_top
            if(dval==1)then
                if(ival==dmem%nval_top)then
                    exit
                elseif(ival==dmem%nval_bot)then
                    nbase=(dmem%idx(dmem%nval_bot+1)-&
                            & dmem%idx(dmem%nval_bot))&
                            &*(dmem%idx_l(dmem%norb-dmem%nval_bot+1)- &
                            & dmem%idx_l(dmem%norb-dmem%nval_bot))
                endif
            elseif(dval==-1.and.ival==dmem%nval_bot)then
                istate=istate+(dmem%idx(dmem%nval_bot+1)-&
                        & dmem%idx(dmem%nval_bot)) &
                        &*(dmem%idx_l(dmem%norb-dmem%nval_bot+1)-&
                        & dmem%idx_l(dmem%norb-dmem%nval_bot))
                cycle
            elseif(abs(dval)>1)then
                write(0,'(" Error in set_dcoo_ij_spci!")')
                stop 2
            endif
            nfs=dmem%idx(ival+dval+1)-dmem%idx(ival+dval)
            nfs_l=dmem%idx_l(dmem%norb-(ival+dval)+1)- &
                    &dmem%idx_l(dmem%norb-(ival+dval))
            if(nfs_l<=0)then
                istate=istate+(dmem%idx(ival+1)-dmem%idx(ival))* &
                        &(dmem%idx_l(dmem%norb-ival+1)- &
                        &dmem%idx_l(dmem%norb-ival))
                cycle
            endif
            do i1=dmem%idx(ival),dmem%idx(ival+1)-1
                do i2=dmem%idx_l(dmem%norb-ival), &
                        &dmem%idx_l(dmem%norb-ival+1)-1
                    istate=istate+1
#ifdef DEBUG
                    if(istate<=0.or.istate>dmem%nstates)then
                        write(0,'(" Error in set_dcoo_ij_spci with &
                                &istate = ", i0)')istate
                        stop 3
                    endif
#endif
                    itmp=ior(dmem%bs(i1),ishft(dmem%bs_l(i2),dmem%norb))
                    isgn=1
                    call act_state(itmp,icf-1,iact,isgn)
                    if(isgn==0)cycle
                    call act_state(itmp,jcf-1,jact,isgn)
                    if(isgn==0)cycle
                    ibs2=dmem%ibs_l(ishft(itmp,-dmem%norb))
                    if(ibs2<=0)cycle
                    nnz=nnz+1
                    if(irun==1)cycle
                    ibs1=dmem%ibs(iand(itmp,mask1))
                    jstate=nbase+(ibs1-dmem%idx(ival+dval))*nfs_l+ &
                            &ibs2-dmem%idx_l(dmem%norb-ival-dval)+1
#ifdef DEBUG
                    if(jstate<=0.or.jstate>dmem%nstates)then
                        write(0,'(" Error in set_dcoo_ij_spci with &
                                &jstate = ", i0)')jstate
                        stop 4
                    endif
#endif
                    dcoo%i(nnz)=istate
                    dcoo%j(nnz)=jstate
                    dcoo%a(nnz)=isgn
                enddo
            enddo
            nbase=nbase+nfs*nfs_l
        enddo
    enddo
    return

end subroutine set_dcoo_ij_spci


subroutine calc_ucsr_spci()
    use gprec, only:dp=>q
    use gspci
    use gutil
    implicit none
    complex(dp) z_row(maxval(dmem%idx(dmem%nval_bot+1:dmem%nval_top+1) &
            &-dmem%idx(dmem%nval_bot:dmem%nval_top)))
    integer nstates,nnz,irun,ival,istate,i,jstate,j,p,q,r,s, &
            &itmp(4),isgn(4)

    nstates=dmem%idx(dmem%nval_top+1)-dmem%idx(dmem%nval_bot)
    dmem%ucsr%nrow=nstates
    dmem%ucsr%ncol=nstates
    allocate(dmem%ucsr%i(nstates+1)); dmem%ucsr%i(1)=1

    do irun=1,2
        if(irun==2)then
            allocate(dmem%ucsr%j(nnz),dmem%ucsr%a(nnz))
        endif
        nnz=0; istate=0
        do ival=dmem%nval_bot,dmem%nval_top
            do i=dmem%idx(ival),dmem%idx(ival+1)-1
                istate=istate+1
                z_row=0
                do p=1,dmem%norb
                    ! <v|p^\dagger
                    isgn(1)=1
                    itmp(1)=dmem%bs(i)
                    call act_state(itmp(1),p-1,.false.,isgn(1))
                    if(isgn(1)==0)cycle

                    ! one-body part
                    do q=1,dmem%norb
                        if(abs(dmem%h1e(p,q))<1.d-10)cycle
                        ! <v|p^\dagger q
                        isgn(2)=isgn(1)
                        itmp(2)=itmp(1)
                        call act_state(itmp(2),q-1,.true.,isgn(2))
                        if(isgn(2)==0)cycle
                        jstate=dmem%ibs(itmp(2))
                        if(jstate>i)cycle
                        jstate=jstate-dmem%idx(ival)+1
                        z_row(jstate)=z_row(jstate)+isgn(2)*dmem%h1e(p,q)
                    enddo

                    ! two-body
                    do q=1,p ! take care of factor 1/2
                        ! <v|p^\dagger q^\dagger
                        isgn(2)=isgn(1)
                        itmp(2)=itmp(1)
                        call act_state(itmp(2),q-1,.false.,isgn(2))
                        if(isgn(2)==0)cycle
                        do r=1,dmem%norb
                            ! <v|p^\dagger q^\dagger r 
                            isgn(3)=isgn(2)
                            itmp(3)=itmp(2)
                            call act_state(itmp(3),r-1,.true.,isgn(3))
                            if(isgn(3)==0)cycle
                            do s=1,dmem%norb
                                if(abs(dmem%v2e(p,s,q,r))<1.d-10)cycle
                                ! <v|p^\dagger q^\dagger r s
                                isgn(4)=isgn(3)
                                itmp(4)=itmp(3)
                                call act_state(itmp(4),s-1,.true.,isgn(4))
                                if(isgn(4)==0)cycle
                                jstate=dmem%ibs(itmp(4))
                                if(jstate>i)cycle
                                jstate=jstate-dmem%idx(ival)+1
                                z_row(jstate)=z_row(jstate)+isgn(4)* &
                                        &dmem%v2e(p,s,q,r)
                            enddo
                        enddo
                    enddo
                enddo
                if(irun==1)then
                    nnz=nnz+count(abs(z_row(1:i-dmem%idx(ival)+1))>1.d-10)
                    dmem%ucsr%i(istate+1)=nnz+1
                    cycle
                else
                    do j=dmem%idx(ival),i
                        jstate=j-dmem%idx(ival)+1
                        if(abs(z_row(jstate))<=1.d-10)cycle
                        nnz=nnz+1
                        dmem%ucsr%j(nnz)=j-dmem%idx(dmem%nval_bot)+1
                        dmem%ucsr%a(nnz)=z_row(jstate)
                    enddo
                endif
            enddo
        enddo
    enddo
    return

end subroutine calc_ucsr_spci


subroutine setup_hdns_spci_dlh(a,n)
    use gprec
    use gspci
    implicit none
    integer,intent(in)::n
    complex(q),intent(out)::a(n,n)

    integer i,j,inz,r,s,nfs,nfs_l,nbase,i1,j1,i1_,i2,ival

    a=0

    ! D_ij
    do i=1,dmem%norb
        do j=1,dmem%norb
            if(abs(dmem%daalpha(i,j))<1.d-10)cycle
            do inz=1,dmem%mcoo(j,i)%nnz
                r=dmem%mcoo(j,i)%i(inz)
                s=dmem%mcoo(j,i)%j(inz)
                a(r,s)=a(r,s)+dmem%mcoo(j,i)%a(inz)*dmem%daalpha(i,j)
                a(s,r)=a(s,r)+dmem%mcoo(j,i)%a(inz)*conjg(dmem%daalpha(i,j))
            enddo
        enddo
    enddo

    ! lambda_c (la2)
    do i=1,dmem%norb
        do j=1,i
            if(abs(dmem%lambdac(i,j))<1.d-10)cycle
            do inz=1,dmem%nvcoo(j,i)%nnz
                r=dmem%nvcoo(j,i)%i(inz)
                s=dmem%nvcoo(j,i)%j(inz)
                a(r,s)=a(r,s)+dmem%nvcoo(j,i)%a(inz)*dmem%lambdac(i,j)
                if(i==j)cycle
                a(s,r)=a(s,r)+dmem%nvcoo(j,i)%a(inz)*conjg(dmem%lambdac(i,j))
            enddo
        enddo
    enddo

    ! h_loc (ucsr)
    nbase=0
    do ival=dmem%nval_bot,dmem%nval_top
        nfs=dmem%idx(ival+1)-dmem%idx(ival)
        nfs_l=dmem%idx_l(dmem%norb-ival+1)-dmem%idx_l(dmem%norb-ival)
        if(nfs_l<=0)cycle
        do i1=dmem%idx(ival),dmem%idx(ival+1)-1
            i1_=i1-dmem%idx(dmem%nval_bot)+1
            do inz=dmem%ucsr%i(i1_),dmem%ucsr%i(i1_+1)-1
                j1=dmem%ucsr%j(inz)+dmem%idx(dmem%nval_bot)-1
                do i2=1,nfs_l
                    r=nbase+(i1-dmem%idx(ival))*nfs_l+i2
                    s=nbase+(j1-dmem%idx(ival))*nfs_l+i2
                    a(r,s)=a(r,s)+dmem%ucsr%a(inz)
                    if(r==s)cycle
                    a(s,r)=a(s,r)+conjg(dmem%ucsr%a(inz))
                enddo
            enddo
        enddo
        nbase=nbase+nfs*nfs_l
    enddo
    return

end subroutine setup_hdns_spci_dlh


subroutine av1_gspci_dlh(v1,v2)
    use gprec
    use gspci
    use sparse
    implicit none
    complex(q),intent(in)::v1(*)
    complex(q),intent(inout)::v2(*)

    integer i,j

    ! D_ij
    do i=1,dmem%norb
        do j=1,dmem%norb
            if(abs(dmem%daalpha(i,j))<1.d-10)cycle
            call zs_dcoomux (dmem%daalpha(i,j),dmem%mcoo(j,i),v1,v2)
            call zs_dcoohmux(conjg(dmem%daalpha(i,j)),dmem%mcoo(j,i),v1,v2)
        enddo
    enddo

    ! lambda_c (la2)
    do i=1,dmem%norb
        do j=1,i
            if(abs(dmem%lambdac(i,j))<1.d-10)cycle
            call zs_dcoomux (dmem%lambdac(i,j),dmem%nvcoo(j,i),v1,v2)
            if(i==j)cycle
            call zs_dcoohmux(dmem%lambdac(j,i),dmem%nvcoo(j,i),v1,v2)
        enddo
    enddo

    ! h_loc (ucsr)
    call act_ucsr_spci(v1,v2)
    return

end subroutine av1_gspci_dlh


subroutine act_ucsr_spci(v1,v2)
    use gprec
    use gspci
    implicit none
    complex(q),intent(in)::v1(*)
    complex(q),intent(inout)::v2(*)

    integer nbase,ival,nfs,nfs_l,i1,i1_,i2,j1,inz,r,s

    nbase=0
    do ival=dmem%nval_bot,dmem%nval_top
        nfs=dmem%idx(ival+1)-dmem%idx(ival)
        nfs_l=dmem%idx_l(dmem%norb-ival+1)-dmem%idx_l(dmem%norb-ival)
        if(nfs_l<=0)cycle
        do i1=dmem%idx(ival),dmem%idx(ival+1)-1
            i1_=i1-dmem%idx(dmem%nval_bot)+1
            do inz=dmem%ucsr%i(i1_),dmem%ucsr%i(i1_+1)-1
                j1=dmem%ucsr%j(inz)+dmem%idx(dmem%nval_bot)-1
                do i2=1,nfs_l
                    r=nbase+(i1-dmem%idx(ival))*nfs_l+i2
                    s=nbase+(j1-dmem%idx(ival))*nfs_l+i2
                    v2(r)=v2(r)+dmem%ucsr%a(inz)*v1(s)
                    if(r==s)cycle
                    v2(s)=v2(s)+conjg(dmem%ucsr%a(inz))*v1(r)
                enddo
            enddo
        enddo
        nbase=nbase+nfs*nfs_l
    enddo
    return

end subroutine act_ucsr_spci


subroutine calc_dm_spci()
    use gprec
    use gspci
    use sparse
    implicit none
    integer i,j,i_,j_

    dmem%dm=0

    ! c_i^\dagger c_j
    do i=1,dmem%norb
        do j=1,i
            call zvh_dcoo_zv(dmem%npcoo(i,j),dmem%v,dmem%dm(i,j))
            if(i==j)cycle
            dmem%dm(j,i)=conjg(dmem%dm(i,j))
        enddo
    enddo

    ! f_j f_i^\dagger = \delta_i,j - f_i^\dagger f_j
    do i=1,dmem%norb
        i_=i+dmem%norb
        do j=1,i
            j_=j+dmem%norb
            call zvh_dcoo_zv(dmem%nvcoo(j,i),dmem%v,dmem%dm(i_,j_))
            if(i_==j_)then
                dmem%dm(i_,j_)=1-dmem%dm(i_,j_)
            else
                dmem%dm(i_,j_)=-dmem%dm(i_,j_)
                dmem%dm(j_,i_)=conjg(dmem%dm(i_,j_))
            endif
        enddo
    enddo

    ! c_j^\dagger f_i
    do i=1,dmem%norb
        i_=i+dmem%norb
        do j=1,dmem%norb
            call zvh_dcoo_zv(dmem%mcoo(j,i),dmem%v,dmem%dm(j,i_))
            dmem%dm(i_,j)=conjg(dmem%dm(j,i_))
        enddo
    enddo
    return

end subroutine calc_dm_spci


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


! A component of angular momentum vector acting on |v1>
! |v2> += \sum_{a,b}{cm_{a,b} (nc_{a,b} + nf_{a,b})} |v1>
subroutine av1_gspci_jop(cm,v1,v2)
    use gprec
    use gspci
    use sparse
    implicit none
    complex(q),intent(in)::cm(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(dmem%nstates)
    complex(q),intent(inout)::v2(dmem%nstates)

    integer i,j

    do i=1,dmem%norb; do j=1,i
        if(abs(cm(i,j))>1.d-12)then
            call zs_dcoomux( cm(i,j),dmem%npcoo(i,j),v1,v2)
            call zs_dcoomux(-cm(i,j),dmem%nvcoo(j,i),v1,v2)
            if(i==j)then
                v2=v2+v1*cm(i,j)
            endif
        endif
        if(i==j)cycle
        if(abs(cm(j,i))>1.d-12)then
            call zs_dcoohmux( cm(j,i),dmem%npcoo(i,j),v1,v2)
            call zs_dcoohmux(-cm(j,i),dmem%nvcoo(j,i),v1,v2)
        endif
    enddo; enddo
    return

end subroutine av1_gspci_jop


! Expectation value.
! <v2| \sum_{a,b}{cm_{a,b} (nc_{a,b} + nf_{a,b})} |v1>
subroutine zv2h_gspci_jop_v1(cm,v1,v2,zes)
    use gprec
    use gspci
    use sparse
    implicit none
    complex(q),intent(in)::cm(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(dmem%nstates),v2(dmem%nstates)
    complex(q),intent(inout)::zes
    complex(q),external::zdotc

    integer i,j

    do i=1,dmem%norb; do j=1,i
        if(abs(cm(i,j))>1.d-12)then
            call zv2h_sdcoo_zv1( cm(i,j),dmem%npcoo(i,j),v1,v2,zes)
            call zv2h_sdcoo_zv1(-cm(i,j),dmem%nvcoo(j,i),v1,v2,zes)
            if(i==j)then
                zes=zes+zdotc(dmem%nstates,v2,1,v1,1)*cm(i,j)
            endif
        endif
        if(i==j)cycle
        if(abs(cm(j,i))>1.d-12)then
            call zv2h_sdcooh_zv1( cm(j,i),dmem%npcoo(i,j),v1,v2,zes)
            call zv2h_sdcooh_zv1(-cm(j,i),dmem%nvcoo(j,i),v1,v2,zes)
        endif
    enddo; enddo
    return

end subroutine zv2h_gspci_jop_v1


! The square of one component of angular momentum operator scting on |v1>
! |v2> += (s/l/j)_{x/y/z}^2 |v1>
subroutine av1_gspci_j2op(cm,v1,v2)
    use gprec
    use gspci
    use sparse
    implicit none
    complex(q),intent(in)::cm(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(dmem%nstates)
    complex(q),intent(inout)::v2(dmem%nstates)

    complex(q) v1p(dmem%nstates)

    v1p=0
    call av1_gspci_jop(cm,v1,v1p)
    call av1_gspci_jop(cm,v1p,v2)
    return

end subroutine av1_gspci_j2op


! |v2> += lambda_j2 * (J_z^2 + Jz) |v1>
subroutine av1_gspci_jz2_pjz(lambda_j2,cm_z,v1,v2)
    use gprec
    use gspci
    use sparse
    implicit none
    real(q),intent(in)::lambda_j2
    complex(q),intent(in)::cm_z(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(dmem%nstates)
    complex(q),intent(inout)::v2(dmem%nstates)

    complex(q) v1p(dmem%nstates)
    complex(q) cm(dmem%norb,dmem%norb)

    v1p=0
    cm=cm_z*sqrt(lambda_j2)
    call av1_gspci_jop(cm,v1,v1p)
    v2=v2+sqrt(lambda_j2)*v1p
    call av1_gspci_jop(cm,v1p,v2)
    return

end subroutine av1_gspci_jz2_pjz


! |v2> += lambda_j2 * (J-J+) |v1>
subroutine av1_gspci_jnjp(lambda_j2,cm_n,v1,v2)
    use gprec
    use gspci
    use sparse
    implicit none
    real(q),intent(in)::lambda_j2
    complex(q),intent(in)::cm_n(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(dmem%nstates)
    complex(q),intent(inout)::v2(dmem%nstates)

    complex(q) v1p(dmem%nstates)
    complex(q) cm(dmem%norb,dmem%norb)

    v1p=0
    cm=conjg(transpose(cm_n))*sqrt(lambda_j2)
    call av1_gspci_jop(cm,v1,v1p)
    cm=cm_n*sqrt(lambda_j2)
    call av1_gspci_jop(cm,v1p,v2)
    return

end subroutine av1_gspci_jnjp


! The expectation value of square of one component of angular momentum 
! operator with respect to |v1>, i.e.,
! <v1| (s/l/j)_{x/y/z}^2 |v1>
subroutine calc_eval_gspci_j2op(cm,v1,zes)
    use gprec
    use gspci
    use sparse
    implicit none
    complex(q),intent(in)::cm(dmem%norb,dmem%norb),v1(dmem%nstates)
    complex(q),intent(inout)::zes

    complex(q) v1p(dmem%nstates)

    v1p=0
    call av1_gspci_jop(cm,v1,v1p)
    call zv2h_gspci_jop_v1(cm,v1p,v1,zes)
    return

end subroutine calc_eval_gspci_j2op


subroutine calc_save_rho_cp_blks(imp)
    use gprec
    use sparse
    use gutil
    use ghdf5_base
    use ghdf5
    use gspci
    implicit none
    integer,intent(in)::imp

    integer ival,nstates,i,j
    complex(q),allocatable::rho(:,:)
    type(dcoo_matrix)::npcoo(dmem%norb,dmem%norb)

    call gh5_open_w('EMBED_HAMIL_ANALYSIS_'//trim(int_to_str(imp))//'.h5',f_id)
    do ival=dmem%nval_bot,dmem%nval_top
        nstates=dmem%idx(ival+1)-dmem%idx(ival)
        allocate(rho(nstates,nstates))
        call calc_reduced_rho_nblk(rho,nstates,ival)
        if(maxval(abs(rho))>1.d-16)then
            call gh5_create_group('/valence_block_'// &
                    &trim(int_to_str(ival)),f_id)
            call gh5_write(rho,nstates,nstates,'/valence_block_'// &
                    &trim(int_to_str(ival))//'/RHO',f_id)
            call calc_npcoo_nblk(npcoo,ival)
            do i=1,dmem%norb; do j=1,i
                if(npcoo(i,j)%nnz==0)cycle
                call gh5_write_compound(npcoo(i,j),'/valence_block_'// &
                        &trim(int_to_str(ival))//'/NP_'//trim(int_to_str(i))// &
                        &'_'//trim(int_to_str(j)),f_id)
                call dealloc_dcoo(npcoo(i,j))
            enddo; enddo
        endif
        deallocate(rho)
    enddo
    call gh5_close(f_id)


end subroutine calc_save_rho_cp_blks


subroutine calc_npcoo_nblk(npcoo,ival)
    use gprec
    use gspci
    use sparse
    use gutil
    implicit none
    integer,intent(in)::ival
    type(dcoo_matrix),intent(out)::npcoo(dmem%norb,dmem%norb)

    integer nstates,ia1,ia2,i1,isum,itmp1,isgn,ibs,nbase

    nstates=dmem%idx(ival+1)-dmem%idx(ival)
    nbase=dmem%idx(ival)-1
    do ia1=1,dmem%norb; do ia2=1,ia1
        call alloc_dcoo(npcoo(ia1,ia2),nstates,nstates,nstates)
        isum=0
        ! < state | c_ia1^\dagger c_ia2
        do i1=dmem%idx(ival),dmem%idx(ival+1)-1
            itmp1=dmem%bs(i1)
            isgn=1
            call act_state(itmp1,ia1-1,.false.,isgn)
            if(isgn==0)cycle
            call act_state(itmp1,ia2-1,.true.,isgn)
            if(isgn==0)cycle
            ibs=dmem%ibs(itmp1)
            isum=isum+1
            npcoo(ia1,ia2)%i(isum)=i1-nbase
            npcoo(ia1,ia2)%j(isum)=ibs-nbase
            npcoo(ia1,ia2)%a(isum)=real(isgn,q)
        enddo
        npcoo(ia1,ia2)%nnz=isum
    enddo; enddo
    return

end subroutine calc_npcoo_nblk


! Generate reduced manybody denity matrix in valence block n.
subroutine calc_reduced_rho_nblk(rho,n,ival)
    use gprec
    use gspci
    implicit none
    integer,intent(in)::n,ival
    complex(q),intent(out)::rho(n,n)
    
    integer i,i_,j,j_,nbase,nfs_l
    complex(q),external::zdotc

    rho=0
    nfs_l=dmem%idx_l(dmem%norb-ival+1)-dmem%idx_l(dmem%norb-ival)
    if(nfs_l<=0)return

    nbase=0
    do i=dmem%nval_bot,ival-1
        nbase=nbase+(dmem%idx(i+1)-dmem%idx(i))* &
                &(dmem%idx_l(dmem%norb-i+1)-dmem%idx_l(dmem%norb-i))
    enddo
    do i=dmem%idx(ival),dmem%idx(ival+1)-1
        i_=i-dmem%idx(ival)+1
        do j=dmem%idx(ival),i
            j_=j-dmem%idx(ival)+1
            rho(i_,j_)=zdotc(nfs_l,dmem%v(nbase+(j-dmem%idx(ival))*nfs_l+1),1,&
                    &dmem%v(nbase+(i-dmem%idx(ival))*nfs_l+1),1)
            if(i==j)cycle
            rho(j_,i_)=conjg(rho(i_,j_))
        enddo
    enddo
    return


end subroutine calc_reduced_rho_nblk


subroutine calc_save_phi_matrix_blks(imp)
    use gprec
    use sparse
    use gutil
    use ghdf5_base
    use ghdf5
    use gspci
    implicit none
    integer,intent(in)::imp

    integer ival,nfs,nfs_l,nbase
    complex(q),allocatable::phi(:,:)

    call gh5_open_w('EMBED_HAMIL_PHIMAT_'//trim(int_to_str(imp))//'.h5',f_id)
    nbase=0
    do ival=dmem%nval_bot,dmem%nval_top
        nfs=dmem%idx(ival+1)-dmem%idx(ival)
        nfs_l=dmem%idx_l(dmem%norb-ival+1)-dmem%idx_l(dmem%norb-ival)
        allocate(phi(nfs,nfs))
        call phi_vec_to_mat(dmem%v(nbase+1:nbase+nfs*nfs_l), &
                &phi,nfs,dmem%bs_l(dmem%idx_l(dmem%norb-ival): &
                &dmem%idx_l(dmem%norb-ival+1)-1), &
                &nfs_l,dmem%norb,dmem%ibs,dmem%idx(ival))
        if(maxval(abs(phi))>1.d-16)then
            call gh5_create_group('/valence_block_'// &
                    &trim(int_to_str(ival)),f_id)
            call gh5_write(phi,nfs,nfs,'/valence_block_'// &
                    &trim(int_to_str(ival))//'/PHI',f_id)
        endif
        deallocate(phi)
        nbase=nbase+nfs*nfs_l
    enddo
    call gh5_close(f_id)


end subroutine calc_save_phi_matrix_blks

