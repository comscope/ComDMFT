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
    real(q) w(dmem%norb2)

    if(.not.associated(dmem%v))then
        allocate(dmem%v(dmem%nstates))
        w(1)=0
    else
        w(1)=1
    endif
    call solve_hembed_spci(w(1:2))
    dmem%etot=w(1)
    return

end subroutine solve_hembed_spci_drive


subroutine solve_hembed_spci(w)
    use gprec
    use gspci
    use gutil
    implicit none
    real(q),intent(inout)::w(2)
    external::av_gspci

    call primme_diag(dmem%v,dmem%nstates,w,av_gspci)
    return

end subroutine solve_hembed_spci


! <v1|c_i^\dagger c_j)|v1>
subroutine expval_gspci_npij(cid,cj,v1,npij)
    use gprec
    use gspci
    use gutil
    implicit none
    integer,intent(in)::cid,cj
    complex(q),intent(inout)::npij
    complex(q),intent(in)::v1(*)

    integer ival,i1,nbase,itmp,isgn1,ibs1,nfs,nfs_l
    complex(q),external::zdotc

    nbase=0
    do ival=dmem%nval_bot,dmem%nval_top
        nfs=dmem%idx(ival+1)-dmem%idx(ival)
        nfs_l=dmem%idx_l(dmem%norb-ival+1)- &
                &dmem%idx_l(dmem%norb-ival)
        if(nfs_l<=0)cycle
        do i1=dmem%idx(ival),dmem%idx(ival+1)-1
            itmp=dmem%bs(i1)
            isgn1=1
            call act_state(itmp,cj-1,.false.,isgn1)
            if(isgn1==0)cycle
            call act_state(itmp,cid-1,.true.,isgn1)
            if(isgn1==0)cycle
            ibs1=dmem%ibs(itmp)
            npij=npij+zdotc(nfs_l,v1(nbase+(ibs1-dmem%idx(ival))*nfs_l+1), &
                    &1,v1(nbase+(i1-dmem%idx(ival))*nfs_l+1),1)*isgn1
        enddo
        nbase=nbase+nfs*nfs_l
    enddo
    return

end subroutine expval_gspci_npij


! <v2|cij c_i^\dagger c_j + cji c_j^\dagger c_i)|v1>
subroutine expval_gspci_npij2(cij,cji,cid,cj,v1,v2,npij)
    use gprec
    use gspci
    use gutil
    implicit none
    integer,intent(in)::cid,cj
    complex(q),intent(in)::cij,cji
    complex(q),intent(inout)::npij
    complex(q),intent(in)::v1(*),v2(*)

    integer ival,i1,nbase,itmp,isgn1,ibs1,nfs,nfs_l
    complex(q),external::zdotc

    nbase=0
    do ival=dmem%nval_bot,dmem%nval_top
        nfs=dmem%idx(ival+1)-dmem%idx(ival)
        nfs_l=dmem%idx_l(dmem%norb-ival+1)- &
                &dmem%idx_l(dmem%norb-ival)
        if(nfs_l<=0)cycle
        do i1=dmem%idx(ival),dmem%idx(ival+1)-1
            itmp=dmem%bs(i1)
            isgn1=1
            call act_state(itmp,cj-1,.false.,isgn1)
            if(isgn1==0)cycle
            call act_state(itmp,cid-1,.true.,isgn1)
            if(isgn1==0)cycle
            ibs1=dmem%ibs(itmp)
            if(abs(cij)>1.d-12)then
                npij=npij+zdotc(nfs_l,v2(nbase+(ibs1-dmem%idx(ival))*nfs_l+1), &
                        &1,v1(nbase+(i1-dmem%idx(ival))*nfs_l+1),1)*isgn1*cij
            endif
            if(cid==cj)cycle
            if(abs(cji)>1.d-12)then
                npij=npij+zdotc(nfs_l,v2(nbase+(i1-dmem%idx(ival))*nfs_l+1), &
                        &1,v1(nbase+(ibs1-dmem%idx(ival))*nfs_l+1),1)*isgn1*cji
            endif
        enddo
        nbase=nbase+nfs*nfs_l
    enddo
    return

end subroutine expval_gspci_npij2


! |v2> += (cij c_i^\dagger c_j + cji c_j^\dagger c_i)|v1>
subroutine av1_gspci_npij(cij,cji,cid,cj,v1,v2)
    use gprec
    use gspci
    use gutil
    implicit none
    integer,intent(in)::cid,cj
    complex(q),intent(in)::cij,cji
    complex(q),intent(in)::v1(*)
    complex(q),intent(inout)::v2(*)

    integer ival,i1,nbase,itmp,isgn1,ibs1,nfs,nfs_l

    nbase=0
    do ival=dmem%nval_bot,dmem%nval_top
        nfs=dmem%idx(ival+1)-dmem%idx(ival)
        nfs_l=dmem%idx_l(dmem%norb-ival+1)- &
                &dmem%idx_l(dmem%norb-ival)
        if(nfs_l<=0)cycle
        do i1=dmem%idx(ival),dmem%idx(ival+1)-1
            itmp=dmem%bs(i1)
            isgn1=1
            call act_state(itmp,cj-1,.false.,isgn1)
            if(isgn1==0)cycle
            call act_state(itmp,cid-1,.true.,isgn1)
            if(isgn1==0)cycle
            ibs1=dmem%ibs(itmp)
            call zaxpy(nfs_l,cij*isgn1,v1(nbase+(i1-dmem%idx(ival))*nfs_l+1),&
                    &1,v2(nbase+(ibs1-dmem%idx(ival))*nfs_l+1),1) 
            
        enddo
        nbase=nbase+nfs*nfs_l
    enddo
    return

end subroutine av1_gspci_npij


! <v1|f_i^\dagger c_j)|v1>
subroutine expval_gspci_mij(fid,cj,v1,mij)
    use gprec
    use gspci
    use gutil
    implicit none
    integer,intent(in)::fid,cj
    complex(q),intent(inout)::mij
    complex(q),intent(in)::v1(*)

    integer ival,istate,i1,i2,jstate,nbase, &
            &itmp,isgn1,isgn2,ibs1,ibs2,nfs1,nfs2,nfs_l1,nfs_l2

    nbase=0
    istate=(dmem%idx(dmem%nval_bot+1)-&
            & dmem%idx(dmem%nval_bot)) &
            &*(dmem%idx_l(dmem%norb-dmem%nval_bot+1)-&
            & dmem%idx_l(dmem%norb-dmem%nval_bot))
    do ival=dmem%nval_bot+1,dmem%nval_top
        nfs1=dmem%idx(ival)-dmem%idx(ival-1)
        nfs_l1=dmem%idx_l(dmem%norb-ival+2)- &
                &dmem%idx_l(dmem%norb-ival+1)
        nfs2=dmem%idx(ival+1)-dmem%idx(ival)
        nfs_l2=dmem%idx_l(dmem%norb-ival+1)- &
                &dmem%idx_l(dmem%norb-ival)
        if(nfs_l1<=0)then
            istate=istate+nfs2*nfs_l2
            cycle
        endif
        if(nfs_l2<=0)then
            nbase=nbase+nfs1*nfs_l1
            cycle
        endif
        do i1=dmem%idx(ival),dmem%idx(ival+1)-1
            itmp=dmem%bs(i1)
            isgn1=1
            call act_state(itmp,cj-1,.false.,isgn1)
            if(isgn1==0)then
                istate=istate+nfs_l2
                cycle
            endif
            if(mod(ival-1,2)==1)then
                isgn1=-isgn1  ! additional sign for f_i^\dagger
            endif
            ibs1=dmem%ibs(itmp)
            do i2=dmem%idx_l(dmem%norb-ival), &
                    &dmem%idx_l(dmem%norb-ival+1)-1
                istate=istate+1
                if(abs(v1(istate))<1.d-16)cycle
                itmp=dmem%bs_l(i2)
                isgn2=isgn1
                call act_state(itmp,fid-1,.true.,isgn2)
                if(isgn2==0)cycle
                ibs2=dmem%ibs_l(itmp)
                if(ibs2<=0)cycle
                jstate=nbase+(ibs1-dmem%idx(ival-1))*nfs_l1+ &
                        &ibs2-dmem%idx_l(dmem%norb-ival+1)+1
                mij=mij+conjg(v1(jstate))*v1(istate)*isgn2
            enddo
        enddo
        nbase=nbase+nfs1*nfs_l1
    enddo
    return

end subroutine expval_gspci_mij


! |v2> += (dij^* f_i^\dagger c_j + dij c_j^\dagger f_i)|v1>
subroutine av1_gspci_mij(dij,fid,cj,v1,v2)
    use gprec
    use gspci
    use gutil
    implicit none
    integer,intent(in)::fid,cj
    complex(q),intent(in)::dij
    complex(q),intent(in)::v1(*)
    complex(q),intent(inout)::v2(*)

    integer ival,istate,i1,i2,jstate,nbase, &
            &itmp,isgn1,isgn2,ibs1,ibs2,nfs1,nfs2,nfs_l1,nfs_l2

    nbase=0
    istate=(dmem%idx(dmem%nval_bot+1)-&
            & dmem%idx(dmem%nval_bot)) &
            &*(dmem%idx_l(dmem%norb-dmem%nval_bot+1)-&
            & dmem%idx_l(dmem%norb-dmem%nval_bot))
    do ival=dmem%nval_bot+1,dmem%nval_top
        nfs1=dmem%idx(ival)-dmem%idx(ival-1)
        nfs_l1=dmem%idx_l(dmem%norb-ival+2)- &
                &dmem%idx_l(dmem%norb-ival+1)
        nfs2=dmem%idx(ival+1)-dmem%idx(ival)
        nfs_l2=dmem%idx_l(dmem%norb-ival+1)- &
                &dmem%idx_l(dmem%norb-ival)
        if(nfs_l1<=0)then
            istate=istate+nfs2*nfs_l2
            cycle
        endif
        if(nfs_l2<=0)then
            nbase=nbase+nfs1*nfs_l1
            cycle
        endif
        do i1=dmem%idx(ival),dmem%idx(ival+1)-1
            itmp=dmem%bs(i1)
            isgn1=1
            call act_state(itmp,cj-1,.false.,isgn1)
            if(isgn1==0)then
                istate=istate+nfs_l2
                cycle
            endif
            if(mod(ival-1,2)==1)then
                isgn1=-isgn1  ! additional sign for f_i^\dagger
            endif
            ibs1=dmem%ibs(itmp)
            do i2=dmem%idx_l(dmem%norb-ival), &
                    &dmem%idx_l(dmem%norb-ival+1)-1
                istate=istate+1
                if(abs(v1(istate))<1.d-16)cycle
                itmp=dmem%bs_l(i2)
                isgn2=isgn1
                call act_state(itmp,fid-1,.true.,isgn2)
                if(isgn2==0)cycle
                ibs2=dmem%ibs_l(itmp)
                if(ibs2<=0)cycle
                jstate=nbase+(ibs1-dmem%idx(ival-1))*nfs_l1+ &
                        &ibs2-dmem%idx_l(dmem%norb-ival+1)+1
                v2(jstate)=v2(jstate)+conjg(dij)*v1(istate)*isgn2
                v2(istate)=v2(istate)+dij*v1(jstate)*isgn2
            enddo
        enddo
        nbase=nbase+nfs1*nfs_l1
    enddo
    return

end subroutine av1_gspci_mij


! <v1|f_j f_i^\dagger|v1>
subroutine expval_gspci_nvij(fid,fj,v1,nvij)
    use gprec
    use gspci
    use gutil
    implicit none
    integer,intent(in)::fid,fj
    complex(q),intent(in)::v1(*)
    complex(q),intent(inout)::nvij

    integer ival,istate,i1,i2,jstate,nbase, &
            &itmp,isgn,ibs1,ibs2,nfs,nfs_l

    nbase=0
    do ival=dmem%nval_bot,dmem%nval_top
        nfs=dmem%idx(ival+1)-dmem%idx(ival)
        nfs_l=dmem%idx_l(dmem%norb-ival+1)- &
                &dmem%idx_l(dmem%norb-ival)
        if(nfs_l<=0)then
            cycle
        endif
        do i2=dmem%idx_l(dmem%norb-ival), &
                &dmem%idx_l(dmem%norb-ival+1)-1
            itmp=dmem%bs_l(i2)
            isgn=1
            call act_state(itmp,fid-1,.true.,isgn)
            if(isgn==0)cycle
            call act_state(itmp,fj-1,.false.,isgn)
            if(isgn==0)cycle
            ibs2=dmem%ibs_l(itmp)
            if(ibs2<=0)cycle
            do i1=dmem%idx(ival),dmem%idx(ival+1)-1
                istate=nbase+(i1-dmem%idx(ival))*nfs_l+ &
                        &i2-dmem%idx_l(dmem%norb-ival)+1
                if(abs(v1(istate))<1.d-16)cycle
                jstate=nbase+(i1-dmem%idx(ival))*nfs_l+ &
                        &ibs2-dmem%idx_l(dmem%norb-ival)+1
                nvij=nvij+conjg(v1(jstate))*v1(istate)*isgn
            enddo
        enddo
        nbase=nbase+nfs*nfs_l
    enddo
    return

end subroutine expval_gspci_nvij


! <v1|f_j f_i^\dagger|v1>
subroutine expval_gspci_nvij2(cij,cji,fid,fj,v1,v2,nvij)
    use gprec
    use gspci
    use gutil
    implicit none
    integer,intent(in)::fid,fj
    complex(q),intent(in)::cij,cji
    complex(q),intent(in)::v1(*),v2(*)
    complex(q),intent(inout)::nvij

    integer ival,istate,i1,i2,jstate,nbase,itmp,isgn,ibs1,ibs2,nfs,nfs_l

    nbase=0
    do ival=dmem%nval_bot,dmem%nval_top
        nfs=dmem%idx(ival+1)-dmem%idx(ival)
        nfs_l=dmem%idx_l(dmem%norb-ival+1)- &
                &dmem%idx_l(dmem%norb-ival)
        if(nfs_l<=0)then
            cycle
        endif
        do i2=dmem%idx_l(dmem%norb-ival), &
                &dmem%idx_l(dmem%norb-ival+1)-1
            itmp=dmem%bs_l(i2)
            isgn=1
            call act_state(itmp,fid-1,.true.,isgn)
            if(isgn==0)cycle
            call act_state(itmp,fj-1,.false.,isgn)
            if(isgn==0)cycle
            ibs2=dmem%ibs_l(itmp)
            if(ibs2<=0)cycle
            do i1=dmem%idx(ival),dmem%idx(ival+1)-1
                istate=nbase+(i1-dmem%idx(ival))*nfs_l+ &
                        &i2-dmem%idx_l(dmem%norb-ival)+1
                if(abs(v1(istate))<1.d-16)cycle
                jstate=nbase+(i1-dmem%idx(ival))*nfs_l+ &
                        &ibs2-dmem%idx_l(dmem%norb-ival)+1
                nvij=nvij+conjg(v2(jstate))*v1(istate)*isgn*cij
                if(fid==fj)cycle
                nvij=nvij+conjg(v2(istate))*v1(jstate)*isgn*cji
            enddo
        enddo
        nbase=nbase+nfs*nfs_l
    enddo
    return

end subroutine expval_gspci_nvij2


! |v2> += (cij f_j f_i^\dagger|v1>
subroutine av1_gspci_nvij(cij,cji,fid,fj,v1,v2)
    use gprec
    use gspci
    use gutil
    implicit none
    integer,intent(in)::fid,fj
    complex(q),intent(in)::cij,cji
    complex(q),intent(in)::v1(*)
    complex(q),intent(inout)::v2(*)

    integer ival,istate,i1,i2,jstate,nbase, &
            &itmp,isgn,ibs1,ibs2,nfs,nfs_l

    nbase=0
    do ival=dmem%nval_bot,dmem%nval_top
        nfs=dmem%idx(ival+1)-dmem%idx(ival)
        nfs_l=dmem%idx_l(dmem%norb-ival+1)- &
                &dmem%idx_l(dmem%norb-ival)
        if(nfs_l<=0)then
            cycle
        endif
        do i2=dmem%idx_l(dmem%norb-ival), &
                &dmem%idx_l(dmem%norb-ival+1)-1
            itmp=dmem%bs_l(i2)
            isgn=1
            call act_state(itmp,fid-1,.true.,isgn)
            if(isgn==0)cycle
            call act_state(itmp,fj-1,.false.,isgn)
            if(isgn==0)cycle
            ibs2=dmem%ibs_l(itmp)
            if(ibs2<=0)cycle
            do i1=dmem%idx(ival),dmem%idx(ival+1)-1
                istate=nbase+(i1-dmem%idx(ival))*nfs_l+ &
                        &i2-dmem%idx_l(dmem%norb-ival)+1
                if(abs(v1(istate))<1.d-16)cycle
                jstate=nbase+(i1-dmem%idx(ival))*nfs_l+ &
                        &ibs2-dmem%idx_l(dmem%norb-ival)+1
                v2(jstate)=v2(jstate)+cij*v1(istate)*isgn
                if(fid==fj)cycle
                v2(istate)=v2(istate)+cji*v1(jstate)*isgn
            enddo
        enddo
        nbase=nbase+nfs*nfs_l
    enddo
    return

end subroutine av1_gspci_nvij


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


subroutine av1_gspci_dlh(v1,v2)
    use gprec
    use gspci
    use sparse
    implicit none
    complex(q),intent(in)::v1(*)
    complex(q),intent(out)::v2(*)

    integer i,j

    ! D_ij
    do i=1,dmem%norb
        do j=1,dmem%norb
            if(abs(dmem%daalpha(i,j))<1.d-10)cycle
            call av1_gspci_mij(dmem%daalpha(i,j),i,j,v1,v2)
        enddo
    enddo

    ! lambda_c (la2)
    do i=1,dmem%norb
        do j=1,i
            if(abs(dmem%lambdac(i,j))<1.d-10)cycle
            call av1_gspci_nvij(dmem%lambdac(i,j),dmem%lambdac(j,i),i,j,v1,v2)
        enddo
    enddo

    ! h_loc (ucsr)
    call av1_gspci_ucsr(v1,v2)
    return

end subroutine av1_gspci_dlh


subroutine av1_gspci_ucsr(v1,v2)
    use gprec
    use gspci
    implicit none
    complex(q),intent(in)::v1(*)
    complex(q),intent(inout)::v2(*)

    integer nbase,ival,nfs,nfs_l,i1,i1_,j1,inz

    nbase=0
    do ival=dmem%nval_bot,dmem%nval_top
        nfs=dmem%idx(ival+1)-dmem%idx(ival)
        nfs_l=dmem%idx_l(dmem%norb-ival+1)-dmem%idx_l(dmem%norb-ival)
        if(nfs_l<=0)cycle
        do i1=dmem%idx(ival),dmem%idx(ival+1)-1
            i1_=i1-dmem%idx(dmem%nval_bot)+1
            do inz=dmem%ucsr%i(i1_),dmem%ucsr%i(i1_+1)-1
                j1=dmem%ucsr%j(inz)+dmem%idx(dmem%nval_bot)-1
                call zaxpy(nfs_l,dmem%ucsr%a(inz), &
                        &v1(nbase+(j1-dmem%idx(ival))*nfs_l+1),1,&
                        &v2(nbase+(i1-dmem%idx(ival))*nfs_l+1),1)
                if(i1==j1)cycle
                call zaxpy(nfs_l,conjg(dmem%ucsr%a(inz)), &
                        &v1(nbase+(i1-dmem%idx(ival))*nfs_l+1),1,&
                        &v2(nbase+(j1-dmem%idx(ival))*nfs_l+1),1)
            enddo
        enddo
        nbase=nbase+nfs*nfs_l
    enddo
    return

end subroutine av1_gspci_ucsr


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
            call expval_gspci_npij(i,j,dmem%v,dmem%dm(i,j))
            if(i==j)cycle
            dmem%dm(j,i)=conjg(dmem%dm(i,j))
        enddo
    enddo

    ! f_j f_i^\dagger = \delta_i,j - f_i^\dagger f_j
    do i=1,dmem%norb
        i_=i+dmem%norb
        do j=1,i
            j_=j+dmem%norb
            call expval_gspci_nvij(i,j,dmem%v,dmem%dm(i_,j_))
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
            call expval_gspci_mij(i,j,dmem%v,dmem%dm(j,i_))
            dmem%dm(i_,j)=dmem%dm(j,i_)
            dmem%dm(j,i_)=conjg(dmem%dm(j,i_))
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
            call av1_gspci_npij(cm(i,j),cm(j,i),i,j,v1,v2)
            call av1_gspci_nvij(-cm(i,j),-cm(j,i),i,j,v1,v2)
            if(i==j)then
                v2=v2+v1*cm(i,j)
            endif
        endif
    enddo; enddo
    return

end subroutine av1_gspci_jop


! Expectation value.
! <v2| \sum_{a,b}{cm_{a,b} (nc_{a,b} + nf_{a,b})} |v1>
! Possible a bug here since <J^2> not real. fixed in onfly2
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
        if(abs(cm(i,j))<1.d-12.and.abs(cm(j,i))<1.d-12)cycle
        call expval_gspci_npij2( cm(i,j), cm(j,i),i,j,v1,v2,zes)
        call expval_gspci_nvij2(-cm(i,j),-cm(j,i),i,j,v1,v2,zes)
        if(i==j)then
            zes=zes+zdotc(dmem%nstates,v2,1,v1,1)*cm(i,j)
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
