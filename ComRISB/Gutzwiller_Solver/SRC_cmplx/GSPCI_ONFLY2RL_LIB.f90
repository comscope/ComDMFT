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

    if(abs(lambda_j2)>1.d-10)then
        call av1_gspci_jz2_pjz(lambda_j2,cm_z,v1,v2)
        call av1_gspci_jnjp(lambda_j2,cm_n,v1,v2)
    endif
    return

end subroutine av1_gspci_j2_npzway


! For physical subspace.
subroutine set_full_fock_states_r_spci()
    use gprec
    use gspci
    implicit none
    integer i

    if(dmem%ngp>0)then
        call set_dmem_phy()
    else
        dmem%bs_r=>dmem%bs; dmem%ibs_r=>dmem%ibs
        dmem%idx_r=>dmem%idx
    endif
    return

end subroutine set_full_fock_states_r_spci


! For variational subspace.
subroutine set_full_fock_states_l_spci()
    use gprec
    use gspci
    implicit none
    integer i

    if(dmem%norb_mott>0.or.dmem%ngp>0)then
        call set_dmem_mott()
    else
        dmem%bs_l=>dmem%bs; dmem%ibs_l=>dmem%ibs
        dmem%idx_l=>dmem%idx
    endif

    dmem%nstates=0
    do i=dmem%nval_bot,dmem%nval_top
        dmem%nstates=dmem%nstates+ &
                &(dmem%idx_r(i+1)-dmem%idx_r(i))* &
                &(dmem%idx_l(dmem%norb-i+1)-dmem%idx_l(dmem%norb-i))
    enddo
    write(0,'(" dim_phi = ",I0)')dmem%nstates
    return

end subroutine set_full_fock_states_l_spci


! setup the "physical" part of the embedding hamiltonian basis.
subroutine set_dmem_phy()
    use gprec
    use gspci
    use gutil
    implicit none
    integer i,n,nbs,isum,sume,ig

    nbs=ishft(1,dmem%norb)
    allocate(dmem%idx_r(0:dmem%norb+1),dmem%bs_r(nbs),dmem%ibs_r(0:nbs-1))
    dmem%ibs_r=0
    dmem%idx_r(0)=1; isum=1
    do n=0,dmem%norb
        fs: do i=dmem%idx(n),dmem%idx(n+1)-1
            do ig=1,dmem%ngp


                call sum_empty_slots_fs(dmem%bs(i),dmem%ngpm(ig), &
                        &dmem%gpm_list(1:dmem%ngpm(ig),ig),sume)
                sume=dmem%ngpm(ig)-sume
                if(sume<dmem%nval_gp(1,1,ig).or. &
                        &sume>dmem%nval_gp(2,1,ig))then
                    cycle fs
                endif
            enddo
            dmem%bs_r(isum)=dmem%bs(i)
            dmem%ibs_r(dmem%bs(i))=isum
            isum=isum+1
        enddo fs
        dmem%idx_r(n+1)=isum
    enddo
    return

end subroutine set_dmem_phy


! setup the "variational" part of the embedding hamiltonian basis.
subroutine set_dmem_mott()
    use gprec
    use gspci
    use gutil
    implicit none
    integer i,n,nbs,isum,sume,ig

    nbs=ishft(1,dmem%norb)
    allocate(dmem%idx_l(0:dmem%norb+1),dmem%bs_l(nbs),dmem%ibs_l(0:nbs-1))
    dmem%ibs_l=0
    dmem%idx_l(0)=1; isum=1
    do n=0,dmem%norb
        if(n<=dmem%norb-dmem%nelect_mott)then
            fs: do i=dmem%idx(n),dmem%idx(n+1)-1
                if(dmem%nelect_mott>0)then
                    call sum_empty_slots_fs(dmem%bs(i),dmem%norb_mott, &
                            &dmem%iorb_mott,sume)
                    if(sume/=dmem%nelect_mott)continue
                endif
                do ig=1,dmem%ngp
                    call sum_empty_slots_fs(dmem%bs(i),dmem%ngpm(ig), &
                            &dmem%gpm_list(1:dmem%ngpm(ig),ig),sume)
                    if(sume<dmem%nval_gp(1,2,ig).or. &
                            &sume>dmem%nval_gp(2,2,ig))then
                        cycle fs
                    endif
                enddo
                dmem%bs_l(isum)=dmem%bs(i)
                dmem%ibs_l(dmem%bs(i))=isum
                isum=isum+1
            enddo fs
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
    complex(q),intent(inout)::npij(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(*)

    integer ival,i1,i,j,nbase,itmp,isgn1,ibs1,nfs_r,nfs_l
    complex(q),external::zdotc

    nbase=0
    do ival=dmem%nval_bot,dmem%nval_top
        nfs_r=dmem%idx_r(ival+1)-dmem%idx_r(ival)
        nfs_l=dmem%idx_l(dmem%norb-ival+1)- &
                &dmem%idx_l(dmem%norb-ival)
        if(nfs_l<=0)cycle
        do i1=dmem%idx_r(ival),dmem%idx_r(ival+1)-1
            do i=1,dmem%norb; do j=1,i
                itmp=dmem%bs_r(i1)
                isgn1=1
                call act_state(itmp,j-1,.false.,isgn1)
                if(isgn1==0)cycle
                call act_state(itmp,i-1,.true.,isgn1)
                if(isgn1==0)cycle
                ibs1=dmem%ibs_r(itmp)
                if(ibs1==0)cycle
                npij(i,j)=npij(i,j)+ &
                        &zdotc(nfs_l, &
                        &v1(nbase+(ibs1-dmem%idx_r(ival))*nfs_l+1), &
                        &1,v1(nbase+(i1-dmem%idx_r(ival))*nfs_l+1),1) &
                        &*isgn1
            enddo; enddo
        enddo
        nbase=nbase+nfs_r*nfs_l
    enddo

    ! Complex conjugate part
    do i=1,dmem%norb; do j=1+i,dmem%norb
        npij(i,j)=conjg(npij(j,i))
    enddo; enddo
    return

end subroutine expval_gspci_npij


! <v2| \sum_ij cij c_i^\dagger c_j + cji c_j^\dagger c_i)|v1>
subroutine expval_gspci_npij2(cm,v1,v2,zes)
    use gprec
    use gspci
    use gutil
    implicit none
    complex(q),intent(in)::cm(dmem%norb,dmem%norb)
    complex(q),intent(inout)::zes
    complex(q),intent(in)::v1(*),v2(*)

    integer ival,i1,i,j,nbase,itmp,isgn1,ibs1,nfs_r,nfs_l
    complex(q),external::zdotc

    nbase=0
    do ival=dmem%nval_bot,dmem%nval_top
        nfs_r=dmem%idx_r(ival+1)-dmem%idx_r(ival)
        nfs_l=dmem%idx_l(dmem%norb-ival+1)- &
                &dmem%idx_l(dmem%norb-ival)
        if(nfs_l<=0)cycle
        do i1=dmem%idx_r(ival),dmem%idx_r(ival+1)-1
            do i=1,dmem%norb; do j=1,i
                if(abs(cm(i,j))<1.d-10.and.abs(cm(j,i))<1.d-10)cycle
                itmp=dmem%bs_r(i1)
                isgn1=1
                call act_state(itmp,j-1,.false.,isgn1)
                if(isgn1==0)cycle
                call act_state(itmp,i-1,.true.,isgn1)
                if(isgn1==0)cycle
                ibs1=dmem%ibs_r(itmp)
                if(ibs1==0)cycle
                if(abs(cm(i,j))>1.d-10)then
                    zes=zes+zdotc(nfs_l, &
                            &v2(nbase+(ibs1-dmem%idx_r(ival))*nfs_l+1), &
                            &1,v1(nbase+(i1-dmem%idx_r(ival))*nfs_l+1),1) &
                            &*isgn1*cm(i,j)
                endif
                if(i==j)cycle
                if(abs(cm(j,i))>1.d-10)then
                    zes=zes+zdotc(nfs_l, &
                            &v2(nbase+(i1-dmem%idx_r(ival))*nfs_l+1), &
                            &1,v1(nbase+(ibs1-dmem%idx_r(ival))*nfs_l+1),1) &
                            &*isgn1*cm(j,i)
                endif
            enddo; enddo
        enddo
        nbase=nbase+nfs_r*nfs_l
    enddo
    return

end subroutine expval_gspci_npij2


! |v2> += (\sum_ij cij c_i^\dagger c_j|v1>
subroutine av1_gspci_npij(cm,v1,v2)
    use gprec
    use gspci
    use gutil
    implicit none
    complex(q),intent(in)::cm(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(*)
    complex(q),intent(inout)::v2(*)

    integer ival,i1,i,j,nbase,itmp,isgn1,ibs1,nfs_r,nfs_l

    nbase=0
    do ival=dmem%nval_bot,dmem%nval_top
        nfs_r=dmem%idx_r(ival+1)-dmem%idx_r(ival)
        nfs_l=dmem%idx_l(dmem%norb-ival+1)- &
                &dmem%idx_l(dmem%norb-ival)
        if(nfs_l<=0)cycle
        do i1=dmem%idx_r(ival),dmem%idx_r(ival+1)-1
            do i=1,dmem%norb; do j=1,i
                if(abs(cm(i,j))<1.d-10.and.abs(cm(j,i))<1.d-10)cycle
                itmp=dmem%bs_r(i1)
                isgn1=1
                call act_state(itmp,j-1,.false.,isgn1)
                if(isgn1==0)cycle
                call act_state(itmp,i-1,.true.,isgn1)
                if(isgn1==0)cycle
                ibs1=dmem%ibs_r(itmp)
                if(ibs1==0)cycle
                if(abs(cm(i,j))>1.d-10)then
                    call zaxpy(nfs_l,cm(i,j)*isgn1, &
                            &v1(nbase+(i1-dmem%idx_r(ival))*nfs_l+1),1,&
                            &v2(nbase+(ibs1-dmem%idx_r(ival))*nfs_l+1),1) 
                endif
                if(i==j)cycle
                if(abs(cm(j,i))>1.d-10)then
                    call zaxpy(nfs_l,cm(j,i)*isgn1, &
                            &v1(nbase+(ibs1-dmem%idx_r(ival))*nfs_l+1),1,&
                            &v2(nbase+(i1-dmem%idx_r(ival))*nfs_l+1),1) 
                endif
            enddo; enddo
        enddo
        nbase=nbase+nfs_r*nfs_l
    enddo
    return

end subroutine av1_gspci_npij


! <v1|f_i^\dagger c_j)|v1>
subroutine expval_gspci_mij(v1,mij)
    use gprec
    use gspci
    use gutil
    implicit none
    complex(q),intent(inout)::mij(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(*)

    integer ival,istate,i1,i2,i,j,jstate,nbase,mbase, &
            &itmp,isgn1,isgn2,ibs1,ibs2, &
            &nfs_r1,nfs_r2,nfs_l1,nfs_l2

    nbase=0
    mbase=(dmem%idx_r(dmem%nval_bot+1)-dmem%idx_r(dmem%nval_bot)) &
            &*(dmem%idx_l(dmem%norb-dmem%nval_bot+1)-&
            & dmem%idx_l(dmem%norb-dmem%nval_bot))
    do ival=dmem%nval_bot+1,dmem%nval_top
        nfs_r1=dmem%idx_r(ival)-dmem%idx_r(ival-1)
        nfs_l1=dmem%idx_l(dmem%norb-ival+2)- &
                &dmem%idx_l(dmem%norb-ival+1)
        nfs_r2=dmem%idx_r(ival+1)-dmem%idx_r(ival)
        nfs_l2=dmem%idx_l(dmem%norb-ival+1)- &
                &dmem%idx_l(dmem%norb-ival)

        if(nfs_l1*nfs_r1<=0)then
            mbase=mbase+nfs_r2*nfs_l2
            cycle
        endif
        if(nfs_l2*nfs_r2<=0)then
            nbase=nbase+nfs_r1*nfs_l1
            cycle
        endif

        do i1=dmem%idx_r(ival),dmem%idx_r(ival+1)-1
            do i=1,dmem%norb; do j=1,dmem%norb
                itmp=dmem%bs_r(i1)
                isgn1=1
                ! c_j |v1>
                call act_state(itmp,j-1,.false.,isgn1)
                if(isgn1==0)cycle
                if(mod(ival-1,2)==1)then
                    isgn1=-isgn1  ! additional sign for f_i^\dagger
                endif
                ibs1=dmem%ibs_r(itmp)
                if(ibs1==0)cycle
                do i2=dmem%idx_l(dmem%norb-ival), &
                        &dmem%idx_l(dmem%norb-ival+1)-1
                    istate=mbase+(i1-dmem%idx_r(ival))*nfs_l2+i2- &
                            &dmem%idx_l(dmem%norb-ival)+1
                    if(abs(v1(istate))<1.d-16)cycle
                    itmp=dmem%bs_l(i2)
                    isgn2=isgn1
                    ! f_i^/dagger c_j |v1>
                    call act_state(itmp,i-1,.true.,isgn2)
                    if(isgn2==0)cycle
                    ibs2=dmem%ibs_l(itmp)
                    if(ibs2==0)cycle
                    jstate=nbase+(ibs1-dmem%idx_r(ival-1))*nfs_l1+ &
                            &ibs2-dmem%idx_l(dmem%norb-ival+1)+1
                    mij(i,j)=mij(i,j)+conjg(v1(jstate))*v1(istate)*isgn2
                enddo
            enddo; enddo
        enddo
        mbase=mbase+nfs_r2*nfs_l2
        nbase=nbase+nfs_r1*nfs_l1
    enddo
    return

end subroutine expval_gspci_mij


! |v2> += \sum_ij (dij^* f_i^\dagger c_j + dij c_j^\dagger f_i)|v1>
subroutine av1_gspci_mij(v1,v2)
    use gprec
    use gspci
    use gutil
    implicit none
    complex(q),intent(in)::v1(*)
    complex(q),intent(inout)::v2(*)

    integer ival,istate,i1,i2,j1,jstate,nbase,nnz,inz,nnz1, &
            &itmp,isgn2,ibs2,nfs_r1,nfs_r2,nfs_l1,nfs_l2
    integer,allocatable::isgn1(:),is(:),js(:),ibs1(:),is1(:)
    complex(q),allocatable::dij(:)

    nnz=count(abs(dmem%daalpha)>=1.d-10)
    allocate(isgn1(nnz),is(nnz),is1(nnz),js(nnz),ibs1(nnz),dij(nnz))

    nnz=0
    do i1=1,dmem%norb; do i2=1,dmem%norb
        if(abs(dmem%daalpha(i1,i2))<1.d-10)cycle
        nnz=nnz+1
        is(nnz)=i1; js(nnz)=i2
    enddo; enddo

    nbase=0
    istate=(dmem%idx_r(dmem%nval_bot+1)-dmem%idx_r(dmem%nval_bot)) &
            &*(dmem%idx_l(dmem%norb-dmem%nval_bot+1)-&
            & dmem%idx_l(dmem%norb-dmem%nval_bot))
    do ival=dmem%nval_bot+1,dmem%nval_top
        nfs_r1=dmem%idx_r(ival)-dmem%idx_r(ival-1)
        nfs_l1=dmem%idx_l(dmem%norb-ival+2)- &
                &dmem%idx_l(dmem%norb-ival+1)
        nfs_r2=dmem%idx_r(ival+1)-dmem%idx_r(ival)
        nfs_l2=dmem%idx_l(dmem%norb-ival+1)- &
                &dmem%idx_l(dmem%norb-ival)
        if(nfs_l1*nfs_r1<=0)then
            istate=istate+nfs_r2*nfs_l2
            cycle
        endif
        if(nfs_l2*nfs_r2<=0)then
            nbase=nbase+nfs_r1*nfs_l1
            cycle
        endif
        do i1=dmem%idx_r(ival),dmem%idx_r(ival+1)-1
            nnz1=1
            do inz=1,nnz
                itmp=dmem%bs_r(i1)
                isgn1(nnz1)=1
                ! c_js |v1>
                call act_state(itmp,js(inz)-1,.false.,isgn1(nnz1))
                if(isgn1(nnz1)==0)cycle
                j1=dmem%ibs_r(itmp)
                if(j1==0)cycle
                ibs1(nnz1)=j1
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
                    ! f_is^\dagger c_js |v1>
                    call act_state(itmp,is1(inz)-1,.true.,isgn2)
                    if(isgn2==0)cycle
                    ibs2=dmem%ibs_l(itmp)
                    if(ibs2==0)cycle
                    jstate=nbase+(ibs1(inz)-dmem%idx_r(ival-1))*nfs_l1+ &
                            &ibs2-dmem%idx_l(dmem%norb-ival+1)+1
                    v2(jstate)=v2(jstate)+conjg(dij(inz))*v1(istate)*isgn2
                    v2(istate)=v2(istate)+dij(inz)*v1(jstate)*isgn2
                enddo
            enddo
        enddo
        nbase=nbase+nfs_r1*nfs_l1
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
    complex(q),intent(in)::v1(*),plist(dmem%norb,dmem%norb,nsym)
    complex(q),intent(inout)::v2(*)

    integer ival,istate,i1,i2,jstate,nbase,itmp,nfs_r,nfs_l,ni_skip,nj_skip, &
            &isym,j1,j2,k,ivalp,jbase
    integer ia(dmem%norb),ib(dmem%norb),ib_skip(dmem%norb), &
            &ja(dmem%norb),jb(dmem%norb),jb_skip(dmem%norb)
    complex(q) coef1,coef2
    real(q) lambda_

    lambda_=lambda_p/nsym
    nbase=0
    do ival=dmem%nval_bot,dmem%nval_top
        ivalp=dmem%norb-ival
        nfs_r=dmem%idx_r(ival+1)-dmem%idx_r(ival)
        nfs_l=dmem%idx_l(ivalp+1)-dmem%idx_l(ivalp)
        if(nfs_l<=0)cycle
        if(ival==0)then
            v2(nbase+1)=v2(nbase+1)+lambda_*v1(nbase+1)*nsym
            nbase=nbase+1
            cycle
        endif
        do i1=dmem%idx_r(ival),dmem%idx_r(ival+1)-1
            itmp=dmem%bs_r(i1)
            call set_occupied_orbitals(itmp,ia,dmem%norb)
            do isym=1,nsym
                call set_skip_orbitals(plist(:,:,isym),dmem%norb, &
                        &ia(1:ival),ival,ib_skip,ni_skip)
                lj1: do j1=dmem%idx_r(ival),dmem%idx_r(ival+1)-1
                    itmp=dmem%bs_r(j1)
                    do k=1,ni_skip
                        if(btest(itmp,ib_skip(k)-1))cycle lj1
                    enddo
                    coef1=1._q
                    call set_occupied_orbitals(itmp,ib,dmem%norb)
                    call calc_fock_state_coef(plist(:,:,isym),dmem%norb, &
                            &ia,ib,ival,coef1)
                    if(abs(coef1)<1.d-12)cycle lj1
                    istate=nbase+(i1-dmem%idx_r(ival))*nfs_l
                    jbase=nbase+(j1-dmem%idx_r(ival))*nfs_l
                    do i2=dmem%idx_l(ivalp),dmem%idx_l(ivalp+1)-1
                        istate=istate+1
                        if(abs(v1(istate))<1.d-16)cycle
                        itmp=dmem%bs_l(i2)
                        call set_occupied_orbitals(itmp,ja,dmem%norb)
                        call set_skip_orbitals(plist(:,:,isym),dmem%norb, &
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
        nbase=nbase+nfs_r*nfs_l
    enddo
    return

end subroutine av1_gspci_pdsym


subroutine expval_gspci_pdsym(plist,nsym,v1,zes)
    use gprec
    use gspci
    use gutil
    implicit none
    integer,intent(in)::nsym
    complex(q),intent(in)::v1(*),plist(dmem%norb,dmem%norb,nsym)
    complex(q),intent(out)::zes

    integer ival,istate,i1,i2,jstate,nbase,itmp, &
            &nfs_r,nfs_l,ni_skip,nj_skip, &
            &isym,j1,j2,k,ivalp,jbase
    integer ia(dmem%norb),ib(dmem%norb),ib_skip(dmem%norb), &
            &ja(dmem%norb),jb(dmem%norb),jb_skip(dmem%norb)
    complex(q) coef1,coef2

    zes=0
    nbase=0
    do ival=dmem%nval_bot,dmem%nval_top
        ivalp=dmem%norb-ival
        nfs_r=dmem%idx_r(ival+1)-dmem%idx_r(ival)
        nfs_l=dmem%idx_l(ivalp+1)-dmem%idx_l(ivalp)
        if(nfs_l<=0)cycle
        if(ival==0)then
            zes=zes+v1(nbase+1)*conjg(v1(nbase+1))*nsym
            nbase=nbase+1
            cycle
        endif
        do i1=dmem%idx_r(ival),dmem%idx_r(ival+1)-1
            itmp=dmem%bs_r(i1)
            call set_occupied_orbitals(itmp,ia,dmem%norb)
            do isym=1,nsym
                call set_skip_orbitals(plist(:,:,isym),dmem%norb, &
                        &ia(1:ival),ival,ib_skip,ni_skip)
                lj1: do j1=dmem%idx_r(ival),dmem%idx_r(ival+1)-1
                    itmp=dmem%bs_r(j1)
                    do k=1,ni_skip
                        if(btest(itmp,ib_skip(k)-1))cycle lj1
                    enddo
                    coef1=1._q
                    call set_occupied_orbitals(itmp,ib,dmem%norb)
                    call calc_fock_state_coef(plist(:,:,isym),dmem%norb, &
                            &ia,ib,ival,coef1)
                    if(abs(coef1)<1.d-12)cycle lj1
                    istate=nbase+(i1-dmem%idx_r(ival))*nfs_l
                    jbase=nbase+(j1-dmem%idx_r(ival))*nfs_l
                    do i2=dmem%idx_l(ivalp),dmem%idx_l(ivalp+1)-1
                        istate=istate+1
                        if(abs(v1(istate))<1.d-16)cycle
                        itmp=dmem%bs_l(i2)
                        call set_occupied_orbitals(itmp,ja,dmem%norb)
                        call set_skip_orbitals(plist(:,:,isym),dmem%norb, &
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
                            zes=zes+coef2*v1(istate)*conjg(v1(jstate))
                        enddo lj2
                    enddo
                enddo lj1
            enddo
        enddo
        nbase=nbase+nfs_r*nfs_l
    enddo
    zes=zes/nsym
    return

end subroutine expval_gspci_pdsym


subroutine chk_expval_gspci_pdsym(plist,nsym,v1)
    use gprec
    use gspci
    implicit none
    integer,intent(in)::nsym
    complex(q),intent(in)::v1(*),plist(dmem%norb,dmem%norb,nsym)

    integer i
    complex(q) ::zes

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
    complex(q),intent(in)::v1(*)
    complex(q),intent(inout)::nvij(dmem%norb,dmem%norb)

    integer ival,i1,i2,i,j,nbase,itmp,isgn,ibs2,nfs_r,nfs_l,ifsr,ifsl
    complex(q),allocatable::v2(:),v2c(:)
    complex(q),external::zdotc

    itmp=maxval(dmem%idx_r(dmem%nval_bot+1:)- &
            &dmem%idx_r(dmem%nval_bot:dmem%nval_top-1))
    allocate(v2(itmp),v2c(itmp))

    nbase=0
    do ival=dmem%nval_bot,dmem%nval_top
        nfs_r=dmem%idx_r(ival+1)-dmem%idx_r(ival)
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
                v2(1:nfs_r)=v1(ifsr:ifsr+(nfs_r-1)*nfs_l:nfs_l)
                if(ibs2==i2)then
                    nvij(i,j)=nvij(i,j)+zdotc(nfs_r,v2(1),1,v2(1),1)*isgn
                else
                    v2c(1:nfs_r)=v1(ifsl:ifsl+(nfs_r-1)*nfs_l:nfs_l)
                    nvij(i,j)=nvij(i,j)+zdotc(nfs_r,v2c(1),1,v2(1),1)*isgn
                endif
            enddo; enddo
        enddo
        nbase=nbase+nfs_r*nfs_l
    enddo
    deallocate(v2,v2c)

    ! Complex conjugate part
    do i=1,dmem%norb; do j=1+i,dmem%norb
        nvij(i,j)=conjg(nvij(j,i))
    enddo; enddo
    return

end subroutine expval_gspci_nvij


! <v2|\sum_ij c_ij f_i^\dagger f_j|v1>
subroutine expval_gspci_nvij2(cm,v1,v2,zes)
    use gprec
    use gspci
    use gutil
    implicit none
    complex(q),intent(in)::cm(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(*),v2(*)
    complex(q),intent(inout)::zes

    integer ival,i1,i2,i,j,nbase,itmp,isgn,ibs2,nfs_r,nfs_l,ifsl,ifsr
    complex(q),allocatable::v1p(:),v2p(:)
    complex(q),external::zdotc

    itmp=maxval(dmem%idx_r(dmem%nval_bot+1:)- &
            &dmem%idx_r(dmem%nval_bot:dmem%nval_top-1))
    allocate(v1p(itmp),v2p(itmp))

    nbase=0
    do ival=dmem%nval_bot,dmem%nval_top
        nfs_r=dmem%idx_r(ival+1)-dmem%idx_r(ival)
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
                    v1p(1:nfs_r)=v1(ifsr:ifsr+(nfs_r-1)*nfs_l:nfs_l)
                    v2p(1:nfs_r)=v2(ifsl:ifsl+(nfs_r-1)*nfs_l:nfs_l)
                    zes=zes+zdotc(nfs_r,v2p(1),1,v1p(1),1)*isgn*cm(i,j)
                endif
                if(i==j)cycle
                if(abs(cm(j,i))>1.d-10)then
                    v1p(1:nfs_r)=v1(ifsl:ifsl+(nfs_r-1)*nfs_l:nfs_l)
                    v2p(1:nfs_r)=v2(ifsr:ifsr+(nfs_r-1)*nfs_l:nfs_l)
                    zes=zes+zdotc(nfs_r,v2p(1),1,v1p(1),1)*isgn*cm(j,i)
                endif
            enddo; enddo
        enddo
        nbase=nbase+nfs_r*nfs_l
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
    complex(q),intent(in)::cm(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(*)
    complex(q),intent(inout)::v2(*)

    integer ival,istate,i1,i2,jstate,nbase,nnz,inz,nnz2, &
            &itmp,isgn,nfs_r,nfs_l
    integer,allocatable::is(:),js(:),ibs2(:)
    logical,allocatable::ldiag(:)
    complex(q),allocatable::cij(:),cji(:)

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
        nfs_r=dmem%idx_r(ival+1)-dmem%idx_r(ival)
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
            do i1=dmem%idx_r(ival),dmem%idx_r(ival+1)-1
                istate=nbase+(i1-dmem%idx_r(ival))*nfs_l+ &
                        &i2-dmem%idx_l(dmem%norb-ival)+1
                if(abs(v1(istate))<1.d-16)cycle
                do inz=1,nnz2
                    jstate=nbase+(i1-dmem%idx_r(ival))*nfs_l+ &
                            &ibs2(inz)-dmem%idx_l(dmem%norb-ival)+1
                    v2(jstate)=v2(jstate)+cij(inz)*v1(istate)
                    if(ldiag(inz))cycle
                    v2(istate)=v2(istate)+cji(inz)*v1(jstate)
                enddo
            enddo
        enddo
        nbase=nbase+nfs_r*nfs_l
    enddo
    deallocate(is,js,ibs2,cij,cji,ldiag)
    return

end subroutine av1_gspci_nvij


subroutine calc_ucsr_spci()
    use gprec, only:dp=>q
    use gspci
    use gutil
    implicit none
    complex(dp) z_row(maxval(dmem%idx_r(dmem%nval_bot+1:dmem%nval_top+1) &
            &-dmem%idx_r(dmem%nval_bot:dmem%nval_top)))
    integer nstates,nnz,irun,ival,istate,i,jstate,j,p,q,r,s, &
            &itmp(4),isgn(4)

    nstates=dmem%idx_r(dmem%nval_top+1)-dmem%idx_r(dmem%nval_bot)
    dmem%ucsr%nrow=nstates
    dmem%ucsr%ncol=nstates
    allocate(dmem%ucsr%i(nstates+1)); dmem%ucsr%i(1)=1

    do irun=1,2
        if(irun==2)then
            allocate(dmem%ucsr%j(nnz),dmem%ucsr%a(nnz))
        endif
        nnz=0; istate=0
        do ival=dmem%nval_bot,dmem%nval_top
            do i=dmem%idx_r(ival),dmem%idx_r(ival+1)-1
                istate=istate+1
                z_row=0
                do p=1,dmem%norb
                    ! <v|p^\dagger
                    isgn(1)=1
                    itmp(1)=dmem%bs_r(i)
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
                        jstate=dmem%ibs_r(itmp(2))
                        if(jstate>i.or.jstate==0)cycle
                        jstate=jstate-dmem%idx_r(ival)+1
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
                                jstate=dmem%ibs_r(itmp(4))
                                if(jstate>i.or.jstate==0)cycle
                                jstate=jstate-dmem%idx_r(ival)+1
                                z_row(jstate)=z_row(jstate)+isgn(4)* &
                                        &dmem%v2e(p,s,q,r)
                            enddo
                        enddo
                    enddo
                enddo
                if(irun==1)then
                    nnz=nnz+count(abs(z_row(1:i-dmem%idx_r(ival)+1))>1.d-10)
                    dmem%ucsr%i(istate+1)=nnz+1
                    cycle
                else
                    do j=dmem%idx_r(ival),i
                        jstate=j-dmem%idx_r(ival)+1
                        if(abs(z_row(jstate))<=1.d-10)cycle
                        nnz=nnz+1
                        dmem%ucsr%j(nnz)=j-dmem%idx_r(dmem%nval_bot)+1
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
    complex(q),intent(inout)::v2(*)

    integer i,j

    ! D_ij
    call av1_gspci_mij(v1,v2)

    ! lambda_c (la2)
    call av1_gspci_nvij(dmem%lambdac,v1,v2)

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

    integer nbase,ival,nfs_r,nfs_l,i1,i1_,j1,inz

    nbase=0
    do ival=dmem%nval_bot,dmem%nval_top
        nfs_r=dmem%idx_r(ival+1)-dmem%idx_r(ival)
        nfs_l=dmem%idx_l(dmem%norb-ival+1)-dmem%idx_l(dmem%norb-ival)
        if(nfs_l<=0)cycle
        do i1=dmem%idx_r(ival),dmem%idx_r(ival+1)-1
            i1_=i1-dmem%idx_r(dmem%nval_bot)+1
            do inz=dmem%ucsr%i(i1_),dmem%ucsr%i(i1_+1)-1
                j1=dmem%ucsr%j(inz)+dmem%idx_r(dmem%nval_bot)-1
                call zaxpy(nfs_l,dmem%ucsr%a(inz), &
                        &v1(nbase+(j1-dmem%idx_r(ival))*nfs_l+1),1,&
                        &v2(nbase+(i1-dmem%idx_r(ival))*nfs_l+1),1)
                if(i1==j1)cycle
                call zaxpy(nfs_l,conjg(dmem%ucsr%a(inz)), &
                        &v1(nbase+(i1-dmem%idx_r(ival))*nfs_l+1),1,&
                        &v2(nbase+(j1-dmem%idx_r(ival))*nfs_l+1),1)
            enddo
        enddo
        nbase=nbase+nfs_r*nfs_l
    enddo
    return

end subroutine av1_gspci_ucsr


subroutine calc_dm_spci()
    use gprec
    use gspci
    use sparse
    implicit none
    complex(q) dm(dmem%norb,dmem%norb)

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
    dmem%dm(1:dmem%norb,1+dmem%norb:)=transpose(conjg(dm))
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
    use gutil
    implicit none
    complex(q),intent(in)::cm(dmem%norb,dmem%norb)
    complex(q),intent(in)::v1(dmem%nstates)
    complex(q),intent(inout)::v2(dmem%nstates)

    integer i
    complex(q) cm_(dmem%norb,dmem%norb)

    call av1_gspci_npij(cm ,v1,v2)
    cm_=-cm
    call av1_gspci_nvij(cm_,v1,v2)
    v2=v2+v1*trace_a(cm,dmem%norb)
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

    call expval_gspci_npij2(cm,v1,v2,zes)
    call expval_gspci_nvij2(cm,v1,v2,zes)
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
        do i1=dmem%idx_r(ival),dmem%idx_r(ival+1)-1
            itmp1=dmem%bs_r(i1)
            isgn=1
            call act_state(itmp1,ia1-1,.false.,isgn)
            if(isgn==0)cycle
            call act_state(itmp1,ia2-1,.true.,isgn)
            if(isgn==0)cycle
            ibs=dmem%ibs_r(itmp1)
            if(ibs==0)cycle
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
        nbase=nbase+(dmem%idx_r(i+1)-dmem%idx_r(i))* &
                &(dmem%idx_l(dmem%norb-i+1)-dmem%idx_l(dmem%norb-i))
    enddo
    do i=dmem%idx_r(ival),dmem%idx_r(ival+1)-1
        i_=dmem%ibs(dmem%bs_r(i))-dmem%idx(ival)+1
        do j=dmem%idx_r(ival),i
            j_=dmem%ibs(dmem%bs_r(j))-dmem%idx(ival)+1
            rho(i_,j_)=zdotc(nfs_l,dmem%v(nbase+(j-dmem%idx(ival))*nfs_l+1),&
                    &1,dmem%v(nbase+(i-dmem%idx(ival))*nfs_l+1),1)
            if(i_==j_)cycle
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

    integer ival,nfs_r,nfs_l,nfs,nbase
    complex(q),allocatable::phi(:,:)

    call gh5_open_w('EMBED_HAMIL_PHIMAT_'//trim(int_to_str(imp))//'.h5',f_id)
    nbase=0
    do ival=dmem%nval_bot,dmem%nval_top
        nfs=dmem%idx(ival+1)-dmem%idx(ival)
        nfs_r=dmem%idx_r(ival+1)-dmem%idx_r(ival)
        nfs_l=dmem%idx_l(dmem%norb-ival+1)-dmem%idx_l(dmem%norb-ival)
        allocate(phi(nfs,nfs))
        call phi_vec_to_mat(dmem%v(nbase+1:nbase+nfs_r*nfs_l),phi,nfs, &
                &dmem%bs_r(dmem%idx_r(ival):dmem%idx_r(ival+1)-1), &
                &nfs_r, &
                &dmem%bs_l(dmem%idx_l(dmem%norb-ival): &
                &dmem%idx_l(dmem%norb-ival+1)-1), &
                &nfs_l,dmem%norb,dmem%ibs,dmem%idx(ival))
        if(maxval(abs(phi))>1.d-16)then
            call gh5_create_group('/valence_block_'// &
                    &trim(int_to_str(ival)),f_id)
            call gh5_write(phi,nfs,nfs,'/valence_block_'// &
                    &trim(int_to_str(ival))//'/PHI',f_id)
        endif
        deallocate(phi)
        nbase=nbase+nfs_r*nfs_l
    enddo
    call gh5_close(f_id)
    return


end subroutine calc_save_phi_matrix_blks

