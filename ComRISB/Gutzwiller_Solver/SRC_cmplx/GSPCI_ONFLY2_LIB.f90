! Using equation J^2 = J-J+ + Jz^2 + Jz
subroutine av1_gspci_j2_npzway(lambda_j2,cm_n,cm_z,v1,v2)
    use gprec
    use gspci
    use sparse
    implicit none
    real(q),intent(in)::lambda_j2
#ifdef EmbReal
    real(q),intent(in)::cm_n(dmem%norb,dmem%norb),cm_z(dmem%norb,dmem%norb)
    real(q),intent(in)::v1(*)
    real(q),intent(inout)::v2(*)
#else
    complex(q),intent(in)::cm_n(dmem%norb,dmem%norb),cm_z(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(*)
    complex(q),intent(inout)::v2(*)
#endif

    if(abs(lambda_j2)>1.d-10)then
        call av1_gspci_jz2_pjz(lambda_j2,cm_z,v1,v2)
        call av1_gspci_jnjp(lambda_j2,cm_n,v1,v2)
    endif
    return

end subroutine av1_gspci_j2_npzway


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
    write(0,'(" dim_phi = ", i0)')dmem%nstates
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
subroutine expval_gspci_npij(v1,npij)
    use gprec
    use gspci
    use gutil
    implicit none
#ifdef EmbReal
    real(q),intent(inout)::npij(dmem%norb,dmem%norb)
    real(q),intent(in)::v1(*)
#else
    complex(q),intent(inout)::npij(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(*)
#endif

    integer ival,i1,i,j,nbase,itmp,isgn1,ibs1,nfs,nfs_l
#ifdef EmbReal
    real(q),external::ddot
#else
    complex(q),external::zdotc
#endif

    nbase=0
    do ival=dmem%nval_bot,dmem%nval_top
        nfs=dmem%idx(ival+1)-dmem%idx(ival)
        nfs_l=dmem%idx_l(dmem%norb-ival+1)- &
                &dmem%idx_l(dmem%norb-ival)
        if(nfs_l<=0)cycle
        do i1=dmem%idx(ival),dmem%idx(ival+1)-1
            do i=1,dmem%norb; do j=1,i
                itmp=dmem%bs(i1)
                isgn1=1
                call act_state(itmp,j-1,.false.,isgn1)
                if(isgn1==0)cycle
                call act_state(itmp,i-1,.true.,isgn1)
                if(isgn1==0)cycle
                ibs1=dmem%ibs(itmp)
                npij(i,j)=npij(i,j)+ &
#ifdef EmbReal
                        &ddot( &
#else                            
                        &zdotc( &
#endif                            
                        &nfs_l,v1(nbase+(ibs1-dmem%idx(ival))*nfs_l+1), &
                        &1,v1(nbase+(i1-dmem%idx(ival))*nfs_l+1),1) &
                        &*isgn1
            enddo; enddo
        enddo
        nbase=nbase+nfs*nfs_l
    enddo

    ! Complex conjugate part
    do i=1,dmem%norb; do j=1+i,dmem%norb
#ifdef EmbReal
        npij(i,j)=npij(j,i)
#else
        npij(i,j)=conjg(npij(j,i))
#endif
    enddo; enddo
    return

end subroutine expval_gspci_npij


! <v2| \sum_ij cij c_i^\dagger c_j + cji c_j^\dagger c_i)|v1>
subroutine expval_gspci_npij2(cm,v1,v2,zes)
    use gprec
    use gspci
    use gutil
    implicit none
#ifdef EmbReal
    real(q),intent(in)::cm(dmem%norb,dmem%norb)
    real(q),intent(inout)::zes
    real(q),intent(in)::v1(*),v2(*)
#else
    complex(q),intent(in)::cm(dmem%norb,dmem%norb)
    complex(q),intent(inout)::zes
    complex(q),intent(in)::v1(*),v2(*)
#endif

    integer ival,i1,i,j,nbase,itmp,isgn1,ibs1,nfs,nfs_l
#ifdef EmbReal
    real(q),external::ddot
#else
    complex(q),external::zdotc
#endif

    nbase=0
    do ival=dmem%nval_bot,dmem%nval_top
        nfs=dmem%idx(ival+1)-dmem%idx(ival)
        nfs_l=dmem%idx_l(dmem%norb-ival+1)- &
                &dmem%idx_l(dmem%norb-ival)
        if(nfs_l<=0)cycle
        do i1=dmem%idx(ival),dmem%idx(ival+1)-1
            do i=1,dmem%norb; do j=1,i
                if(abs(cm(i,j))<1.d-10.and.abs(cm(j,i))<1.d-10)cycle
                itmp=dmem%bs(i1)
                isgn1=1
                call act_state(itmp,j-1,.false.,isgn1)
                if(isgn1==0)cycle
                call act_state(itmp,i-1,.true.,isgn1)
                if(isgn1==0)cycle
                ibs1=dmem%ibs(itmp)
                if(abs(cm(i,j))>1.d-10)then
                    zes=zes+ &
#ifdef EmbReal
                            &ddot( &
#else
                            &zdotc( &
#endif
                            &nfs_l,v2(nbase+(ibs1-dmem%idx(ival))*nfs_l+1), &
                            &1,v1(nbase+(i1-dmem%idx(ival))*nfs_l+1),1) &
                            &*isgn1*cm(i,j)
                endif
                if(i==j)cycle
                if(abs(cm(j,i))>1.d-10)then
                    zes=zes+ &
#ifdef EmbReal
                            &ddot( &
#else
                            &zdotc( &
#endif
                            &nfs_l,v2(nbase+(i1-dmem%idx(ival))*nfs_l+1), &
                            &1,v1(nbase+(ibs1-dmem%idx(ival))*nfs_l+1),1) &
                            &*isgn1*cm(j,i)
                endif
            enddo; enddo
        enddo
        nbase=nbase+nfs*nfs_l
    enddo
    return

end subroutine expval_gspci_npij2


! |v2> += (\sum_ij cij c_i^\dagger c_j|v1>
subroutine av1_gspci_npij(cm,v1,v2)
    use gprec
    use gspci
    use gutil
    implicit none
#ifdef EmbReal
    real(q),intent(in)::cm(dmem%norb,dmem%norb)
    real(q),intent(in)::v1(*)
    real(q),intent(inout)::v2(*)
#else
    complex(q),intent(in)::cm(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(*)
    complex(q),intent(inout)::v2(*)
#endif

    integer ival,i1,i,j,nbase,itmp,isgn1,ibs1,nfs,nfs_l

    nbase=0
    do ival=dmem%nval_bot,dmem%nval_top
        nfs=dmem%idx(ival+1)-dmem%idx(ival)
        nfs_l=dmem%idx_l(dmem%norb-ival+1)- &
                &dmem%idx_l(dmem%norb-ival)
        if(nfs_l<=0)cycle
        do i1=dmem%idx(ival),dmem%idx(ival+1)-1
            do i=1,dmem%norb; do j=1,i
                if(abs(cm(i,j))<1.d-10.and.abs(cm(j,i))<1.d-10)cycle
                itmp=dmem%bs(i1)
                isgn1=1
                call act_state(itmp,j-1,.false.,isgn1)
                if(isgn1==0)cycle
                call act_state(itmp,i-1,.true.,isgn1)
                if(isgn1==0)cycle
                ibs1=dmem%ibs(itmp)
                if(abs(cm(i,j))>1.d-10)then
#ifdef EmbReal
                    call daxpy(nfs_l,cm(i,j)*isgn1, &
#else
                    call zaxpy(nfs_l,cm(i,j)*isgn1, &
#endif
                            &v1(nbase+(i1-dmem%idx(ival))*nfs_l+1),1,&
                            &v2(nbase+(ibs1-dmem%idx(ival))*nfs_l+1),1) 
                endif
                if(i==j)cycle
                if(abs(cm(j,i))>1.d-10)then
#ifdef EmbReal
                    call daxpy(nfs_l,cm(j,i)*isgn1, &
#else
                    call zaxpy(nfs_l,cm(j,i)*isgn1, &
#endif
                            &v1(nbase+(ibs1-dmem%idx(ival))*nfs_l+1),1,&
                            &v2(nbase+(i1-dmem%idx(ival))*nfs_l+1),1) 
                endif
            enddo; enddo
        enddo
        nbase=nbase+nfs*nfs_l
    enddo
    return

end subroutine av1_gspci_npij


! <v1|f_i^\dagger c_j)|v1>
subroutine expval_gspci_mij(v1,mij)
    use gprec
    use gspci
    use gutil
    implicit none
#ifdef EmbReal
    real(q),intent(inout)::mij(dmem%norb,dmem%norb)
    real(q),intent(in)::v1(*)
#else
    complex(q),intent(inout)::mij(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(*)
#endif

    integer ival,istate,i1,i2,i,j,jstate,nbase,mbase, &
            &itmp,isgn1,isgn2,ibs1,ibs2,nfs1,nfs2,nfs_l1,nfs_l2

    nbase=0
    mbase=(dmem%idx(dmem%nval_bot+1)-&
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
            mbase=mbase+nfs2*nfs_l2
            cycle
        endif
        if(nfs_l2<=0)then
            nbase=nbase+nfs1*nfs_l1
            cycle
        endif
        do i1=dmem%idx(ival),dmem%idx(ival+1)-1
            do i=1,dmem%norb; do j=1,dmem%norb
                itmp=dmem%bs(i1)
                isgn1=1
                call act_state(itmp,j-1,.false.,isgn1)
                if(isgn1==0)cycle
                if(mod(ival-1,2)==1)then
                    isgn1=-isgn1  ! additional sign for f_i^\dagger
                endif
                ibs1=dmem%ibs(itmp)
                do i2=dmem%idx_l(dmem%norb-ival), &
                        &dmem%idx_l(dmem%norb-ival+1)-1
                    istate=mbase+(i1-dmem%idx(ival))*nfs_l2+i2- &
                            &dmem%idx_l(dmem%norb-ival)+1
                    if(abs(v1(istate))<1.d-16)cycle
                    itmp=dmem%bs_l(i2)
                    isgn2=isgn1
                    call act_state(itmp,i-1,.true.,isgn2)
                    if(isgn2==0)cycle
                    ibs2=dmem%ibs_l(itmp)
                    if(ibs2<=0)cycle
                    jstate=nbase+(ibs1-dmem%idx(ival-1))*nfs_l1+ &
                            &ibs2-dmem%idx_l(dmem%norb-ival+1)+1
#ifdef EmbReal
                    mij(i,j)=mij(i,j)+v1(jstate)*v1(istate)*isgn2
#else
                    mij(i,j)=mij(i,j)+conjg(v1(jstate))*v1(istate)*isgn2
#endif
                enddo
            enddo; enddo
        enddo
        mbase=mbase+nfs2*nfs_l2
        nbase=nbase+nfs1*nfs_l1
    enddo
    return

end subroutine expval_gspci_mij


! |v2> += \sum_ij (dij^* f_i^\dagger c_j + dij c_j^\dagger f_i)|v1>
subroutine av1_gspci_mij(v1,v2)
    use gprec
    use gspci
    use gutil
    implicit none
#ifdef EmbReal
    real(q),intent(in)::v1(*)
    real(q),intent(inout)::v2(*)
#else
    complex(q),intent(in)::v1(*)
    complex(q),intent(inout)::v2(*)
#endif

    integer ival,istate,i1,i2,jstate,nbase,nnz,inz,nnz1, &
            &itmp,isgn2,ibs2,nfs1,nfs2,nfs_l1,nfs_l2
    integer,allocatable::isgn1(:),is(:),js(:),ibs1(:),is1(:)
#ifdef EmbReal
    real(q),allocatable::dij(:)
#else
    complex(q),allocatable::dij(:)
#endif

    nnz=count(abs(dmem%daalpha)<1.d-10)
    allocate(isgn1(nnz),is(nnz),is1(nnz),js(nnz),ibs1(nnz),dij(nnz))
    nnz=0
    do i1=1,dmem%norb; do i2=1,dmem%norb
        if(abs(dmem%daalpha(i1,i2))<1.d-10)cycle
        nnz=nnz+1
        is(nnz)=i1; js(nnz)=i2
    enddo; enddo

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
            nnz1=1
            do inz=1,nnz
                itmp=dmem%bs(i1)
                isgn1(nnz1)=1
                call act_state(itmp,js(inz)-1,.false.,isgn1(nnz1))
                if(isgn1(nnz1)==0)cycle
                ibs1(nnz1)=dmem%ibs(itmp)
                is1(nnz1)=is(inz)
                dij(nnz1)=dmem%daalpha(is(inz),js(inz))
                nnz1=nnz1+1
            enddo
            nnz1=nnz1-1
            if(nnz1<=0)then
                istate=istate+nfs_l2
                cycle
            endif
            if(mod(ival-1,2)==1)then
                isgn1=-isgn1  ! additional sign for f_i^\dagger
            endif
            do i2=dmem%idx_l(dmem%norb-ival), &
                    &dmem%idx_l(dmem%norb-ival+1)-1
                istate=istate+1
                if(abs(v1(istate))<1.d-16)cycle
                do inz=1,nnz1
                    itmp=dmem%bs_l(i2)
                    isgn2=isgn1(inz)
                    call act_state(itmp,is1(inz)-1,.true.,isgn2)
                    if(isgn2==0)cycle
                    ibs2=dmem%ibs_l(itmp)
                    if(ibs2<=0)cycle
                    jstate=nbase+(ibs1(inz)-dmem%idx(ival-1))*nfs_l1+ &
                            &ibs2-dmem%idx_l(dmem%norb-ival+1)+1
#ifdef EmbReal
                    v2(jstate)=v2(jstate)+dij(inz)*v1(istate)*isgn2
#else
                    v2(jstate)=v2(jstate)+conjg(dij(inz))*v1(istate)*isgn2
#endif
                    v2(istate)=v2(istate)+dij(inz)*v1(jstate)*isgn2
                enddo
            enddo
        enddo
        nbase=nbase+nfs1*nfs_l1
    enddo
    deallocate(isgn1,is,is1,js,ibs1,dij)
    return

end subroutine av1_gspci_mij


! |v2> += P |v1>
!      += \sum_k v_k P c_i^\dagger ... f_j^\dagger |0>
!      += \sum_k v_k \sum_n Pn c_i^\dagger Pn^\dagger ... 
!                           Pn f_j^\dagger Pn^\dagger |0>
!      += \sum_k v_k \sum_n \sum_i' u_ii' c_i'^\dagger ...
!                           \sum_j' u_jj' f_j'^\dagger |0>
subroutine av1_gspci_pdsym(lambda_p,nsym,plist,v1,v2)
    use gprec
    use gspci
    use gutil
    implicit none
    integer,intent(in)::nsym
    real(q),intent(in)::lambda_p
#ifdef EmbReal
    real(q),intent(in)::v1(*),plist(dmem%norb,dmem%norb,nsym)
    real(q),intent(inout)::v2(*)
#else
    complex(q),intent(in)::v1(*),plist(dmem%norb,dmem%norb,nsym)
    complex(q),intent(inout)::v2(*)
#endif

    integer ival,istate,i1,i2,jstate,nbase,itmp,nfs,nfs_l,ni_skip,nj_skip, &
            &isym,j1,j2,k,ivalp,jbase
    integer ia(dmem%norb),ib(dmem%norb),ib_skip(dmem%norb), &
            &ja(dmem%norb),jb(dmem%norb),jb_skip(dmem%norb)
    real(q) dtmp(dmem%norb,dmem%norb)
#ifdef EmbReal
    real(q) coef1,coef2
#else
    complex(q) coef1,coef2
#endif
    real(q) lambda_

    lambda_=lambda_p/nsym
    nbase=0
    do ival=dmem%nval_bot,dmem%nval_top
        ivalp=dmem%norb-ival
        nfs=dmem%idx(ival+1)-dmem%idx(ival)
        nfs_l=dmem%idx_l(ivalp+1)-dmem%idx_l(ivalp)
        if(nfs_l<=0)cycle
        if(ival==0)then
            v2(nbase+1)=v2(nbase+1)+lambda_*v1(nbase+1)*nsym
            nbase=nbase+1
            cycle
        endif
        do i1=dmem%idx(ival),dmem%idx(ival+1)-1
            itmp=dmem%bs(i1)
            call set_occupied_orbitals(itmp,ia,dmem%norb)
            do isym=1,nsym
                dtmp=abs(plist(:,:,isym))
                call set_skip_orbitals(dtmp,dmem%norb, &
                        &ia(1:ival),ival,ib_skip,ni_skip)
                lj1: do j1=dmem%idx(ival),dmem%idx(ival+1)-1
                    itmp=dmem%bs(j1)
                    do k=1,ni_skip
                        if(btest(itmp,ib_skip(k)-1))cycle lj1
                    enddo
                    coef1=1._q
                    call set_occupied_orbitals(itmp,ib,dmem%norb)
                    call calc_fock_state_coef(plist(:,:,isym),dmem%norb, &
                            &ia,ib,ival,coef1)
                    if(abs(coef1)<1.d-12)cycle lj1
                    istate=nbase+(i1-dmem%idx(ival))*nfs_l
                    jbase=nbase+(j1-dmem%idx(ival))*nfs_l
                    do i2=dmem%idx_l(ivalp),dmem%idx_l(ivalp+1)-1
                        istate=istate+1
                        if(abs(v1(istate))<1.d-16)cycle
                        itmp=dmem%bs_l(i2)
                        call set_occupied_orbitals(itmp,ja,dmem%norb)
                        call set_skip_orbitals(dtmp,dmem%norb, &
                                &ja(1:ivalp),ivalp,jb_skip,nj_skip)
                        lj2: do j2=dmem%idx_l(ivalp),dmem%idx_l(ivalp+1)-1
                            itmp=dmem%bs_l(j2)
                            do k=1,nj_skip
                                if(btest(itmp,jb_skip(k)-1))cycle lj2
                            enddo
                            coef2=coef1
                            call set_occupied_orbitals(itmp,jb,dmem%norb)
                            call calc_fock_state_coef(plist(:,:,isym), &
                                    &dmem%norb,ja,jb,ivalp,coef2)
                            if(abs(coef2)<1.d-16)cycle lj2
                            jstate=jbase+j2-dmem%idx_l(ivalp)+1
                            v2(jstate)=v2(jstate)+coef2*v1(istate)*lambda_
                        enddo lj2
                    enddo
                enddo lj1
            enddo
        enddo
        nbase=nbase+nfs*nfs_l
    enddo
    return

end subroutine av1_gspci_pdsym


subroutine expval_gspci_pdsym(plist,nsym,v1,zes)
    use gprec
    use gspci
    use gutil
    implicit none
    integer,intent(in)::nsym
#ifdef EmbReal
    real(q),intent(in)::v1(*),plist(dmem%norb,dmem%norb,nsym)
    real(q),intent(out)::zes
#else
    complex(q),intent(in)::v1(*),plist(dmem%norb,dmem%norb,nsym)
    complex(q),intent(out)::zes
#endif

    integer ival,istate,i1,i2,jstate,nbase,itmp,nfs,nfs_l,ni_skip,nj_skip, &
            &isym,j1,j2,k,ivalp,jbase
    integer ia(dmem%norb),ib(dmem%norb),ib_skip(dmem%norb), &
            &ja(dmem%norb),jb(dmem%norb),jb_skip(dmem%norb)
    real(q) dtmp(dmem%norb,dmem%norb)
#ifdef EmbReal
    real(q) coef1,coef2
#else
    complex(q) coef1,coef2
#endif

    zes=0
    nbase=0
    do ival=dmem%nval_bot,dmem%nval_top
        ivalp=dmem%norb-ival
        nfs=dmem%idx(ival+1)-dmem%idx(ival)
        nfs_l=dmem%idx_l(ivalp+1)-dmem%idx_l(ivalp)
        if(nfs_l<=0)cycle
        if(ival==0)then
#ifdef EmbReal
            zes=zes+v1(nbase+1)*v1(nbase+1)*nsym
#else
            zes=zes+v1(nbase+1)*conjg(v1(nbase+1))*nsym
#endif
            nbase=nbase+1
            cycle
        endif
        do i1=dmem%idx(ival),dmem%idx(ival+1)-1
            itmp=dmem%bs(i1)
            call set_occupied_orbitals(itmp,ia,dmem%norb)
            do isym=1,nsym
                dtmp=abs(plist(:,:,isym))
                call set_skip_orbitals(dtmp,dmem%norb, &
                        &ia(1:ival),ival,ib_skip,ni_skip)
                lj1: do j1=dmem%idx(ival),dmem%idx(ival+1)-1
                    itmp=dmem%bs(j1)
                    do k=1,ni_skip
                        if(btest(itmp,ib_skip(k)-1))cycle lj1
                    enddo
                    coef1=1._q
                    call set_occupied_orbitals(itmp,ib,dmem%norb)
                    call calc_fock_state_coef(plist(:,:,isym),dmem%norb, &
                            &ia,ib,ival,coef1)
                    if(abs(coef1)<1.d-12)cycle lj1
                    istate=nbase+(i1-dmem%idx(ival))*nfs_l
                    jbase=nbase+(j1-dmem%idx(ival))*nfs_l
                    do i2=dmem%idx_l(ivalp),dmem%idx_l(ivalp+1)-1
                        istate=istate+1
                        if(abs(v1(istate))<1.d-16)cycle
                        itmp=dmem%bs_l(i2)
                        call set_occupied_orbitals(itmp,ja,dmem%norb)
                        call set_skip_orbitals(dtmp,dmem%norb, &
                                &ja(1:ivalp),ivalp,jb_skip,nj_skip)
                        lj2: do j2=dmem%idx_l(ivalp),dmem%idx_l(ivalp+1)-1
                            itmp=dmem%bs_l(j2)
                            do k=1,nj_skip
                                if(btest(itmp,jb_skip(k)-1))cycle lj2
                            enddo
                            coef2=coef1
                            call set_occupied_orbitals(itmp,jb,dmem%norb)
                            call calc_fock_state_coef(plist(:,:,isym), &
                                    &dmem%norb,ja,jb,ivalp,coef2)
                            if(abs(coef2)<1.d-16)cycle lj2
                            jstate=jbase+j2-dmem%idx_l(ivalp)+1
#ifdef EmbReal
                            zes=zes+coef2*v1(istate)*v1(jstate)
#else
                            zes=zes+coef2*v1(istate)*conjg(v1(jstate))
#endif
                        enddo lj2
                    enddo
                enddo lj1
            enddo
        enddo
        nbase=nbase+nfs*nfs_l
    enddo
    zes=zes/nsym
    return

end subroutine expval_gspci_pdsym


subroutine chk_expval_gspci_pdsym(plist,nsym,v1)
    use gprec
    use gspci
    implicit none
    integer,intent(in)::nsym
#ifdef EmbReal
    real(q),intent(in)::v1(*),plist(dmem%norb,dmem%norb,nsym)
#else
    complex(q),intent(in)::v1(*),plist(dmem%norb,dmem%norb,nsym)
#endif

    integer i
#ifdef EmbReal
    real(q) ::zes
#else
    complex(q) ::zes
#endif

    do i=1,nsym
        call expval_gspci_pdsym(plist(:,:,i),1,v1,zes)
        write(0,'(" isym = ",i2," <v|pojector_sym|v> =",2f18.10)')i,zes
    enddo
    return

end subroutine chk_expval_gspci_pdsym


! <v1|f_i^\dagger f_j|v1>
subroutine expval_gspci_nvij(v1,nvij)
    use gprec
    use gspci
    use gutil
    implicit none
#ifdef EmbReal
    real(q),intent(in)::v1(*)
    real(q),intent(inout)::nvij(dmem%norb,dmem%norb)
#else
    complex(q),intent(in)::v1(*)
    complex(q),intent(inout)::nvij(dmem%norb,dmem%norb)
#endif

    integer ival,i1,i2,i,j,nbase,itmp,isgn,ibs2,nfs,nfs_l,ifsr,ifsl
#ifdef EmbReal
    real(q),allocatable::v2(:),v2c(:)
    real(q),external::ddot
#else
    complex(q),allocatable::v2(:),v2c(:)
    complex(q),external::zdotc
#endif

    itmp=maxval(dmem%idx(dmem%nval_bot+1:)- &
            &dmem%idx(dmem%nval_bot:dmem%nval_top-1))
    allocate(v2(itmp),v2c(itmp))

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
            do i=1,dmem%norb; do j=1,i
                itmp=dmem%bs_l(i2)
                isgn=1
                call act_state(itmp,j-1,.false.,isgn)
                if(isgn==0)cycle
                call act_state(itmp,i-1,.true.,isgn)
                if(isgn==0)cycle
                ibs2=dmem%ibs_l(itmp)
                if(ibs2<=0)cycle
                ifsr=nbase+i2  -dmem%idx_l(dmem%norb-ival)+1
                ifsl=nbase+ibs2-dmem%idx_l(dmem%norb-ival)+1
                v2(1:nfs)=v1(ifsr:ifsr+(nfs-1)*nfs_l:nfs_l)
                if(ibs2==i2)then
                    nvij(i,j)=nvij(i,j)+ &
#ifdef EmbReal
                            &ddot(nfs,v2(1),1,v2(1),1) &
#else
                            &zdotc(nfs,v2(1),1,v2(1),1) &
#endif
                            &*isgn
                else
                    v2c(1:nfs)=v1(ifsl:ifsl+(nfs-1)*nfs_l:nfs_l)
                    nvij(i,j)=nvij(i,j)+ &
#ifdef EmbReal
                            &ddot(nfs,v2c(1),1,v2(1),1) &
#else
                            &zdotc(nfs,v2c(1),1,v2(1),1) &
#endif                            
                            &*isgn
                endif
            enddo; enddo
        enddo
        nbase=nbase+nfs*nfs_l
    enddo
    deallocate(v2,v2c)

    ! Complex conjugate part
    do i=1,dmem%norb; do j=1+i,dmem%norb
#ifdef EmbReal
        nvij(i,j)=nvij(j,i)
#else
        nvij(i,j)=conjg(nvij(j,i))
#endif
    enddo; enddo
    return

end subroutine expval_gspci_nvij


! <v2|\sum_ij c_ij f_i^\dagger f_j|v1>
subroutine expval_gspci_nvij2(cm,v1,v2,zes)
    use gprec
    use gspci
    use gutil
    implicit none
#ifdef EmbReal
    real(q),intent(in)::cm(dmem%norb,dmem%norb)
    real(q),intent(in)::v1(*),v2(*)
    real(q),intent(inout)::zes
#else
    complex(q),intent(in)::cm(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(*),v2(*)
    complex(q),intent(inout)::zes
#endif

    integer ival,i1,i2,i,j,nbase,itmp,isgn,ibs2,nfs,nfs_l,ifsl,ifsr
#ifdef EmbReal
    real(q),allocatable::v1p(:),v2p(:)
    real(q),external::ddot
#else
    complex(q),allocatable::v1p(:),v2p(:)
    complex(q),external::zdotc
#endif

    itmp=maxval(dmem%idx(dmem%nval_bot+1:)- &
            &dmem%idx(dmem%nval_bot:dmem%nval_top-1))
    allocate(v1p(itmp),v2p(itmp))

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
            do i=1,dmem%norb; do j=1,i
                if(abs(cm(i,j))<1.d-10.and.abs(cm(j,i))<1.d-10)cycle
                itmp=dmem%bs_l(i2)
                isgn=1
                call act_state(itmp,j-1,.false.,isgn)
                if(isgn==0)cycle
                call act_state(itmp,i-1,.true.,isgn)
                if(isgn==0)cycle
                ibs2=dmem%ibs_l(itmp)
                if(ibs2<=0)cycle
                ifsr=nbase+i2  -dmem%idx_l(dmem%norb-ival)+1
                ifsl=nbase+ibs2-dmem%idx_l(dmem%norb-ival)+1
                if(abs(cm(i,j))>1.d-10)then
                    v1p(1:nfs)=v1(ifsr:ifsr+(nfs-1)*nfs_l:nfs_l)
                    v2p(1:nfs)=v2(ifsl:ifsl+(nfs-1)*nfs_l:nfs_l)
                    zes=zes+ &
#ifdef EmbReal
                            &ddot(nfs,v2p(1),1,v1p(1),1) &
#else
                            &zdotc(nfs,v2p(1),1,v1p(1),1) &
#endif
                            &*isgn*cm(i,j)
                endif
                if(i==j)cycle
                if(abs(cm(j,i))>1.d-10)then
                    v1p(1:nfs)=v1(ifsl:ifsl+(nfs-1)*nfs_l:nfs_l)
                    v2p(1:nfs)=v2(ifsr:ifsr+(nfs-1)*nfs_l:nfs_l)
                    zes=zes+ &
#ifdef EmbReal
                            &ddot(nfs,v2p(1),1,v1p(1),1) &
#else
                            &zdotc(nfs,v2p(1),1,v1p(1),1) &
#endif                            
                            &*isgn*cm(j,i)
                endif
            enddo; enddo
        enddo
        nbase=nbase+nfs*nfs_l
    enddo
    deallocate(v1p,v2p)
    return

end subroutine expval_gspci_nvij2


! |v2> += \sum_{i,j} (cij f_j f_i^\dagger|v1>
subroutine av1_gspci_nvij(cm,v1,v2)
    use gprec
    use gspci
    use gutil
    implicit none
#ifdef EmbReal
    real(q),intent(in)::cm(dmem%norb,dmem%norb)
    real(q),intent(in)::v1(*)
    real(q),intent(inout)::v2(*)
#else
    complex(q),intent(in)::cm(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(*)
    complex(q),intent(inout)::v2(*)
#endif

    integer ival,istate,i1,i2,jstate,nbase,nnz,inz,nnz2, &
            &itmp,isgn,nfs,nfs_l
    integer,allocatable::is(:),js(:),ibs2(:)
    logical,allocatable::ldiag(:)
#ifdef EmbReal
    real(q),allocatable::cij(:),cji(:)
#else
    complex(q),allocatable::cij(:),cji(:)
#endif

    nnz=0
    do i1=1,dmem%norb; do i2=1,i1
        if(abs(cm(i1,i2))<1.d-10.and.abs(cm(i2,i1))<1.d-10)cycle
        nnz=nnz+1
    enddo; enddo
    allocate(is(nnz),js(nnz),ibs2(nnz),cij(nnz),cji(nnz),ldiag(nnz))
    nnz=0
    do i1=1,dmem%norb; do i2=1,i1
        if(abs(cm(i1,i2))<1.d-10.and.abs(cm(i2,i1))<1.d-10)cycle
        nnz=nnz+1
        is(nnz)=i1; js(nnz)=i2
    enddo; enddo

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
            nnz2=0
            do inz=1,nnz
                itmp=dmem%bs_l(i2)
                isgn=1
                call act_state(itmp,is(inz)-1,.true.,isgn)
                if(isgn==0)cycle
                call act_state(itmp,js(inz)-1,.false.,isgn)
                if(isgn==0)cycle
                if(dmem%ibs_l(itmp)<=0)cycle
                nnz2=nnz2+1
                ibs2(nnz2)=dmem%ibs_l(itmp)
                cij(nnz2)=cm(is(inz),js(inz))*isgn
                cji(nnz2)=cm(js(inz),is(inz))*isgn
                ldiag(nnz2)=is(inz)==js(inz)
            enddo
            do i1=dmem%idx(ival),dmem%idx(ival+1)-1
                istate=nbase+(i1-dmem%idx(ival))*nfs_l+ &
                        &i2-dmem%idx_l(dmem%norb-ival)+1
                if(abs(v1(istate))<1.d-16)cycle
                do inz=1,nnz2
                    jstate=nbase+(i1-dmem%idx(ival))*nfs_l+ &
                            &ibs2(inz)-dmem%idx_l(dmem%norb-ival)+1
                    v2(jstate)=v2(jstate)+cij(inz)*v1(istate)
                    if(ldiag(inz))cycle
                    v2(istate)=v2(istate)+cji(inz)*v1(jstate)
                enddo
            enddo
        enddo
        nbase=nbase+nfs*nfs_l
    enddo
    deallocate(is,js,ibs2,cij,cji,ldiag)
    return

end subroutine av1_gspci_nvij


subroutine calc_ucsr_spci()
    use gprec, only:dp=>q
    use gspci
    use gutil
    implicit none
#ifdef EmbReal
    real(dp) z_row(maxval(dmem%idx(dmem%nval_bot+1:dmem%nval_top+1) &
            &-dmem%idx(dmem%nval_bot:dmem%nval_top)))
#else
    complex(dp) z_row(maxval(dmem%idx(dmem%nval_bot+1:dmem%nval_top+1) &
            &-dmem%idx(dmem%nval_bot:dmem%nval_top)))
#endif
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
                                ! lower trianglar part only
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
    write(0,'(" dim_ucsr =",i0)')dmem%ucsr%i(nstates+1)
    return

end subroutine calc_ucsr_spci


subroutine av1_gspci_dlh(v1,v2)
    use gprec
    use gspci
    use sparse
    use gtime
    implicit none
#ifdef EmbReal
    real(q),intent(in)::v1(*)
    real(q),intent(out)::v2(*)
#else
    complex(q),intent(in)::v1(*)
    complex(q),intent(out)::v2(*)
#endif

    integer i,j

#ifdef DEBUG
    call set_time_point(1,3)
#endif
    ! D_ij
    call av1_gspci_mij(v1,v2)
#ifdef DEBUG
    call set_time_point(2,3)
    call print_time_usage('av1_gspci_mij',3,0)
    call set_time_point(1,3)
#endif
    ! lambda_c (la2)
    call av1_gspci_nvij(dmem%lambdac,v1,v2)
#ifdef DEBUG
    call set_time_point(2,3)
    call print_time_usage('av1_gspci_nvij',3,0)
    call set_time_point(1,3)
#endif
    ! h_loc (ucsr)
    call av1_gspci_ucsr(v1,v2)
#ifdef DEBUG
    call set_time_point(2,3)
    call print_time_usage('av1_gspci_ucsr',3,0)
#endif
    return

end subroutine av1_gspci_dlh


subroutine av1_gspci_ucsr(v1,v2)
    use gprec
    use gspci
    implicit none
#ifdef EmbReal
    real(q),intent(in)::v1(*)
    real(q),intent(inout)::v2(*)
#else
    complex(q),intent(in)::v1(*)
    complex(q),intent(inout)::v2(*)
#endif

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
#ifdef EmbReal
                call daxpy(nfs_l,dmem%ucsr%a(inz), &
#else
                call zaxpy(nfs_l,dmem%ucsr%a(inz), &
#endif
                        &v1(nbase+(j1-dmem%idx(ival))*nfs_l+1),1,&
                        &v2(nbase+(i1-dmem%idx(ival))*nfs_l+1),1)
                if(i1==j1)cycle
#ifdef EmbReal
                call daxpy(nfs_l,dmem%ucsr%a(inz), &
#else
                call zaxpy(nfs_l,conjg(dmem%ucsr%a(inz)), &
#endif
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
#ifdef EmbReal
    real(q) dm(dmem%norb,dmem%norb)
#else
    complex(q) dm(dmem%norb,dmem%norb)
#endif

    ! c_i^\dagger c_j
    dm=0
    call expval_gspci_npij(dmem%v,dm)
    dmem%dm(1:dmem%norb,1:dmem%norb)=dm

    ! f_i^\dagger f_j
    dm=0
    call expval_gspci_nvij(dmem%v,dm)
    dmem%dm(1+dmem%norb:,1+dmem%norb:)=dm

    ! c_j^\dagger f_i
    dm=0
    call expval_gspci_mij(dmem%v,dm)
    dmem%dm(1+dmem%norb:,1:dmem%norb)=dm
#ifdef EmbReal
    dmem%dm(1:dmem%norb,1+dmem%norb:)=transpose(dm)
#else
    dmem%dm(1:dmem%norb,1+dmem%norb:)=transpose(conjg(dm))
#endif
    return

end subroutine calc_dm_spci


subroutine av_gspci(v1,v2,k,primme)
    use gprec
    use gspci
    use gconstant
    implicit none
    integer,intent(in)::k
    integer(8),intent(in)::primme
#ifdef EmbReal
    real(q),intent(in)::v1(*)
    real(q),intent(out)::v2(*)
#else
    complex(q),intent(in)::v1(*)
    complex(q),intent(out)::v2(*)
#endif

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
    use gutil
    implicit none
#ifdef EmbReal
    real(q),intent(in)::cm(dmem%norb,dmem%norb)
    real(q),intent(in)::v1(dmem%nstates)
    real(q),intent(inout)::v2(dmem%nstates)
#else
    complex(q),intent(in)::cm(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(dmem%nstates)
    complex(q),intent(inout)::v2(dmem%nstates)
#endif

    integer i
#ifdef EmbReal
    real(q) cm_(dmem%norb,dmem%norb)
#else
    complex(q) cm_(dmem%norb,dmem%norb)
#endif

    call av1_gspci_npij(cm ,v1,v2)
    cm_=-cm
    call av1_gspci_nvij(cm_,v1,v2)
    v2=v2+v1*trace_a(cm,dmem%norb)
    return

end subroutine av1_gspci_jop


! Expectation value.
! <v2| \sum_{a,b}{cm_{a,b} (nc_{a,b} + nf_{a,b})} |v1>
subroutine v2h_gspci_jop_v1(cm,v1,v2,zes)
    use gprec
    use gspci
    use sparse
    implicit none
#ifdef EmbReal
    real(q),intent(in)::cm(dmem%norb,dmem%norb)
    real(q),intent(in)::v1(dmem%nstates),v2(dmem%nstates)
    real(q),intent(inout)::zes
    real(q),external::ddot
#else
    complex(q),intent(in)::cm(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(dmem%nstates),v2(dmem%nstates)
    complex(q),intent(inout)::zes
    complex(q),external::zdotc
#endif

    call expval_gspci_npij2(cm,v1,v2,zes)
    call expval_gspci_nvij2(cm,v1,v2,zes)
    return

end subroutine v2h_gspci_jop_v1


! The square of one component of angular momentum operator scting on |v1>
! |v2> += (s/l/j)_{x/y/z}^2 |v1>
subroutine av1_gspci_j2op(cm,v1,v2)
    use gprec
    use gspci
    use sparse
    implicit none
#ifdef EmbReal
    real(q),intent(in)::cm(dmem%norb,dmem%norb)
    real(q),intent(in)::v1(dmem%nstates)
    real(q),intent(inout)::v2(dmem%nstates)

    real(q) v1p(dmem%nstates)
#else
    complex(q),intent(in)::cm(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(dmem%nstates)
    complex(q),intent(inout)::v2(dmem%nstates)

    complex(q) v1p(dmem%nstates)
#endif

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
#ifdef EmbReal
    real(q),intent(in)::cm_z(dmem%norb,dmem%norb)
    real(q),intent(in)::v1(dmem%nstates)
    real(q),intent(inout)::v2(dmem%nstates)

    real(q) v1p(dmem%nstates)
    real(q) cm(dmem%norb,dmem%norb)
#else
    complex(q),intent(in)::cm_z(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(dmem%nstates)
    complex(q),intent(inout)::v2(dmem%nstates)

    complex(q) v1p(dmem%nstates)
    complex(q) cm(dmem%norb,dmem%norb)
#endif

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
#ifdef EmbReal
    real(q),intent(in)::cm_n(dmem%norb,dmem%norb)
    real(q),intent(in)::v1(dmem%nstates)
    real(q),intent(inout)::v2(dmem%nstates)

    real(q) v1p(dmem%nstates)
    real(q) cm(dmem%norb,dmem%norb)
#else
    complex(q),intent(in)::cm_n(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(dmem%nstates)
    complex(q),intent(inout)::v2(dmem%nstates)

    complex(q) v1p(dmem%nstates)
    complex(q) cm(dmem%norb,dmem%norb)
#endif

    v1p=0
#ifdef EmbReal
    cm=transpose(cm_n)*sqrt(lambda_j2)
#else
    cm=conjg(transpose(cm_n))*sqrt(lambda_j2)
#endif
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
#ifdef EmbReal
    real(q),intent(in)::cm_vec(dmem%norb,dmem%norb,3)
#else
    complex(q),intent(in)::cm_vec(dmem%norb,dmem%norb,3)
#endif

    integer i
#ifdef EmbReal
    real(q) eval
#else
    complex(q) eval
#endif

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
#ifdef EmbReal
    real(q),intent(in)::cm(dmem%norb,dmem%norb),v1(dmem%nstates)
    real(q),intent(inout)::zes

    real(q) v1p(dmem%nstates)
#else
    complex(q),intent(in)::cm(dmem%norb,dmem%norb),v1(dmem%nstates)
    complex(q),intent(inout)::zes

    complex(q) v1p(dmem%nstates)
#endif

    v1p=0
    call av1_gspci_jop(cm,v1,v1p)
    call v2h_gspci_jop_v1(cm,v1p,v1,zes)
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
#ifdef EmbReal
    real(q),allocatable::rho(:,:)
#else
    complex(q),allocatable::rho(:,:)
#endif
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
    return


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
#ifdef EmbReal
    real(q),intent(out)::rho(n,n)
#else
    complex(q),intent(out)::rho(n,n)
#endif
    
    integer i,i_,j,j_,nbase,nfs_l
#ifdef EmbReal
    real(q),external::ddot
#else
    complex(q),external::zdotc
#endif

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
            rho(i_,j_)= &
#ifdef EmbReal
                    &ddot( &
#else
                    &zdotc( &
#endif
                    &nfs_l,dmem%v(nbase+(j-dmem%idx(ival))*nfs_l+1), &
                    &dmem%v(nbase+(i-dmem%idx(ival))*nfs_l+1),1)
            if(i==j)cycle
#ifdef EmbReal
            rho(j_,i_)=rho(i_,j_)
#else
            rho(j_,i_)=conjg(rho(i_,j_))
#endif
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
#ifdef EmbReal
    real(q),allocatable::phi(:,:)
#else
    complex(q),allocatable::phi(:,:)
#endif

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
    return


end subroutine calc_save_phi_matrix_blks

