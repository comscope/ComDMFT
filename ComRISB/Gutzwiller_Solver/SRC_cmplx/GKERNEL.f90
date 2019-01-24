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

module gkernel
    use gprec
    use ghdf5_base
    use gutil
    use bandstru
    use warehouse
    use dcstd
    use gbroyden
    use gparam
    implicit none
    
    type gk_mem
        integer::imix,iembeddiag,nx,nx1,nx2
        integer::iter=0,nbreset=0,max_iter,giter=-1
        real(8)::amix,rtol=1.e-6_q
        real(8)::maxerr=1.e4_q,etot
        real(8),pointer::x(:),fvec(:)
    end type gk_mem

    type(gk_mem)::gkmem

    contains


    subroutine init_gkernel(io)
    integer,intent(in)::io

    call init_g_solver(io)
    if(gkmem%iembeddiag==10)then
        call init_wh_x(io,1._q)
    else
        call init_wh_x(io,.99_q)
    endif
    return

    end subroutine init_gkernel


    subroutine init_newton_rl()
    
    !! initialize solution {R,\lambda} vector
    gkmem%nx1=wh%hm_r%dimhst*wh%r_factor
    gkmem%nx2=wh%hm_l%dimhst
    gkmem%nx=gkmem%nx1+gkmem%nx2
    allocate(gkmem%x(gkmem%nx),gkmem%fvec(gkmem%nx))
    gkmem%x(1:gkmem%nx1)=transfer(wh%r_coef,gkmem%x(1:gkmem%nx1))
    gkmem%x(1+gkmem%nx1:)=wh%la1_coef
    return

    end subroutine init_newton_rl


    subroutine init_newton_n()

    !! initialize solution {n} vector
    gkmem%nx=wh%hm%dimhst
    allocate(gkmem%x(gkmem%nx),gkmem%fvec(gkmem%nx))
    gkmem%x=wh%nks_coef
    return

    end subroutine init_newton_n


    subroutine g_newton_solver(io,av)
    integer,intent(in)::io
    external::av

    integer imore
    real(q),parameter :: rtol=1.e-10_q,epsfcn=1.e-10_q,nitmin=1000

    !! If {n} as the variables like lda+u
    if(gkmem%iembeddiag==10)then
        call init_newton_n()
    else
        call init_newton_rl()
    endif

    if(gkmem%nbreset==0)then
        gkmem%nbreset=gkmem%nx
    endif

    !! For full Mott-localized solution
    if(gkmem%nx<=0)then
        call g_rl_onecycle()
        return
    endif

    if(io>0)then
        write(io,'(" solving nonlinear equations with dim_x = ",I0)')gkmem%nx
    endif

    select case(gkmem%imix)
    case(0)
        do imore=1,3
            call ghybrd( &
                    &av,gkmem%nx,gkmem%x,gkmem%fvec,rtol,epsfcn,io)
            if(gkmem%iter>=nitmin.or.gkmem%maxerr<gkmem%rtol.or. &
                    &gkmem%iter<=3)exit
        enddo
    case(-1,1)
        call broyden_mix( &
                &av,gkmem%nx,gkmem%x,gkmem%fvec,gkmem%amix,gkmem%nbreset, &
                &gkmem%imix)
    case default
        write(0,'(" error in g_newton_solver: undefined imix!")')
        stop
    end select
    deallocate(gkmem%x,gkmem%fvec)
    nullify(gkmem%x,gkmem%fvec)
    return

    end subroutine g_newton_solver


    subroutine init_g_solver(io)
    integer,intent(in)::io

    logical lexist

    call gh5_open_r('GPARAM.h5',f_id)
    call gh5_read(gkmem%imix,'/gimix',f_id)
    call gh5_read(gkmem%iembeddiag,'/giembeddiag',f_id)
    call gh5_read(gkmem%amix,'/gamix',f_id)
    call gh5_read(gkmem%nbreset,'/gnbreset',f_id)
    call gh5_read(gkmem%max_iter,'/gmaxiter',f_id)
    call gh5_close(f_id)

    inquire(file='GMOTT.h5',exist=lexist)
    if(lexist)then
        call gh5_open_r('GMOTT.h5',f_id)
        call gh5_read(gkmem%iembeddiag,'/giembeddiag',f_id)
        call gh5_close(f_id)
    endif

    ! Overwritten by command argument if available
    if(g_maxiter>=0)then
        gkmem%max_iter=g_maxiter
    endif
    if(g_iembeddiag/=-1000)then
        gkmem%iembeddiag=g_iembeddiag
    endif
    if(g_imix/=-100)then
        gkmem%imix=g_imix
    endif
    return

    end subroutine init_g_solver


    subroutine solve_hembed_list(io)
    integer,intent(in)::io

    integer i,imap

    call set_time_point(1,4)

#ifdef mpi_mode
    wh%nc_phy=0; wh%nc_var=0; wh%r0=0; wh%eu2=0
#endif
    !$omp parallel do private(i,imap) schedule(static,1)
    do i=1,wh%num_imp
        if(wh%local_imp(i)==0)cycle
        imap=wh%imap_imp(i)
        if(imap==i)then
            call solve_hembed(i,wh%na2_imp(i)*2,io)
        else
            wh%co(i)%nc_phy=wh%co(imap)%nc_phy
            wh%co(i)%nc_var=wh%co(imap)%nc_var
            wh%co(i)%r0=wh%co(imap)%r0
            wh%eu2(i)=wh%eu2(imap)
        endif
    enddo
    !$omp end parallel do
#ifdef mpi_mode
    call sum_all_mpi(wh%nc_phy,wh%na2112)
    call sum_all_mpi(wh%nc_var,wh%na2112)
    call sum_all_mpi(wh%r0,wh%na2112)
    call sum_all_mpi(wh%eu2,wh%num_imp)
#endif

    call set_time_point(2,4)
    call print_time_usage('solve_hembed',4,io)
    return

    end subroutine solve_hembed_list


    subroutine solve_hembed(i,na4,io)
    integer,intent(in)::i,na4,io

    integer na2,iter,ityp,ierr,ia,ll,stat
    real(q) etot,de,mag_mom
    complex(q) h1e4(na4,na4),dm(na4,na4),h1e2(na4/2,na4/2)
    
    na2=na4/2
    ityp=wh%ityp_imp(i)
    ll = (na2/2-1)/2

    ! remove V2E_imp.inp file in the first iteration.
    if(gkmem%iter<=1)then
        open(unit=123, iostat=stat, file="V2E_"//trim(int_to_str(i))//'.INP', &
                &status='old')
        if(stat==0)then
            close(123, status='delete')
        endif
    endif

    h1e2=wh%co(i)%h1e
    ! adding (orbital dependent) DC term to h1e
    if(associated(wh%co(i)%vext).and.wh%ivext>0)then
        h1e2(:na2,:na2)=h1e2(:na2,:na2)+wh%co(i)%vext
    endif

    call gh5_open_w('EMBED_HAMIL_'//trim(int_to_str(i))//'.h5',f_id) 
    call gh5_write(h1e2,na2,na2,'/H1E',f_id)
    call gh5_write(wh%co(i)%d,na2,na2,'/D',f_id)
    call gh5_write(wh%co(i)%la2,na2,na2,'/LAMBDA',f_id)
    call gh5_write(wh%co(i)%v2e,na2,na2,na2,na2,'/V2E',f_id)
    call gh5_write(na2,'/na2',f_id)
    call gh5_write(wh%nval_bot_ityp(ityp),'/nval_bot',f_id)
    call gh5_write(wh%nval_top_ityp(ityp),'/nval_top',f_id)
    call gh5_write(wh%co(i)%m_struct,na2,na2,'/sigma_struct',f_id)
    call gh5_write(wh%fz(i)%nsorb,'/norb_mott',f_id)
    if(wh%fz(i)%nsorb>0)then
        call gh5_write(wh%fz(i)%nelect,'/nelect_mott',f_id)
        call gh5_write(wh%fz(i)%idx_orb,wh%fz(i)%nsorb,'/iorb_mott',f_id)
    endif
    call gh5_close(f_id)

    de=trace_a(wh%co(i)%la2,na2)
#ifdef flowversion
    if(gkmem%iembeddiag==-1)then
        call system('exe_spci_mott '//int_to_str(i))
    elseif(gkmem%iembeddiag==-2)then
        call system('exe_spci_sjz_mott -sz -i '//int_to_str(i))
    elseif(gkmem%iembeddiag==-3)then
        call system('exe_spci_s2_mott -i '//int_to_str(i))
    elseif(gkmem%iembeddiag==-4)then
        call system('exe_spci_sjz_mott -i '//int_to_str(i))
    elseif(gkmem%iembeddiag==-21)then
        call system('exe_spci_sz_mott '//int_to_str(i))
    elseif(gkmem%iembeddiag==-11)then
        call system('gs_ml.py -i '//int_to_str(i)//' -l '//int_to_str(ll))
    elseif(gkmem%iembeddiag==-12)then
        call system('gs_syten.py -i '//int_to_str(i)//' -l '//int_to_str(ll))
#else
    if(gkmem%iembeddiag==-1)then
        call execute_command_line('exe_spci_mott '//int_to_str(i), exitstat=ierr)
        if(ierr/=0)then
            write(0,'(" Error in running exe_spci_mott!")')
        endif
    elseif(gkmem%iembeddiag==-2)then
        call execute_command_line('exe_spci_sjz_mott -sz -i '//int_to_str(i),&
                &exitstat=ierr)
        if(ierr/=0)then
            write(0,'(" Error in running exe_spci_sjz_mott!")')
        endif
    elseif(gkmem%iembeddiag==-3)then
        call execute_command_line('exe_spci_s2_mott -i '//int_to_str(i),  &
                &exitstat=ierr)
        if(ierr/=0)then
            write(0,'(" Error in running exe_spci_s2_mott!")')
        endif
    elseif(gkmem%iembeddiag==-4)then
        call execute_command_line('exe_spci_sjz_mott -i '// &
                &int_to_str(i),exitstat=ierr)
        if(ierr/=0)then
            write(0,'(" Error in running exe_spci_sjz_mott!")')
        endif
    elseif(gkmem%iembeddiag==-21)then
        call execute_command_line('exe_spci_sz_mott '//int_to_str(i),  &
                &exitstat=ierr)
        if(ierr/=0)then
            write(0,'(" Error in running exe_spci_sz_mott!")')
        endif
    elseif(gkmem%iembeddiag==-11)then
        call execute_command_line(' gs_ml.py -i '//int_to_str(i)//' -l '// &
                &int_to_str(ll), exitstat=ierr)
        if(ierr/=0)then
            write(0,'(" Error in running gs_ml.py!")')
        endif
    elseif(gkmem%iembeddiag==-12)then
        call execute_command_line(' gs_syten.py -i '//int_to_str(i)//' -l '// &
                &int_to_str(ll), exitstat=ierr)
        if(ierr/=0)then
            write(0,'(" Error in running gs_syten.py!")')
        endif
#endif
    else
        write(0,'(" Error in solve_hembed: unsupported embedding solver!")')
        stop 99
    endif

    call gh5_open_r('EMBED_HAMIL_RES_'//trim(int_to_str(i))//'.h5',f_id)
    call gh5_read(etot,'/emol',f_id)
    call gh5_read(dm,na4,na4,'/DM',f_id)
    call gh5_close(f_id)

    if(io>0)then
        write(io,'(" spci e_mol = ",f0.7," e_tot = ",f0.7)')etot,etot+de
    endif

    ! setup h1e4
    h1e4(1:na2,1:na2)=wh%co(i)%h1e
    h1e4(1+na2:,1:na2)=conjg(wh%co(i)%d)
    h1e4(1:na2,1+na2:)=transpose(wh%co(i)%d)
    h1e4(1+na2:,1+na2:)=-wh%co(i)%la2
    wh%eu2(i)=etot-sum(h1e4*dm)

    call get_rn_from_embeddm(dm,wh%co(i)%nc_phy, &
            &wh%co(i)%nc_var,wh%co(i)%r0,na2)
    return

    end subroutine solve_hembed


    subroutine calc_total_energy()
    real(q) etot
    integer io

    io=gp%io
    call calc_et1_list()
    call calc_ebnd()
    call calc_edcla1()
    call calc_edc_list()
    if(gkmem%iembeddiag==10)then
        call calc_eu2_hf_list()
    endif
    etot=sum(bnd%eband)-wh%edcla1-sum(dc%e)+sum(wh%eu2)+sum(wh%et1)
    if(io>0)then
        write(io,'(" e_total_model = ",f0.7)')etot
        if(gkmem%iter>1)then
            write(io,'(" e_total_model diff = ",f0.7)')etot-gkmem%etot
        endif
    endif
    gkmem%etot=etot
    return

    end subroutine calc_total_energy


end module gkernel


subroutine g_fcn(n,x,fvec,iflag)
    use gprec
    use gkernel
    implicit none
    integer,intent(in)::n
    integer,intent(out)::iflag
    real(8),intent(in)::x(n)
    real(8),intent(out)::fvec(n)

    integer io
    real(q) maxerr

    call set_time_point(1,3)

    io=gp%io
    gkmem%iter=gkmem%iter+1

    if(io>0)then
        write(io,'(/,"cygutz iteration = ", i0)')gkmem%iter
        write(io,'("********** wh%x **********")')
        write(io,'(5f14.8)')x
    endif

    if(gkmem%iembeddiag==10)then
        call g_fcn_n(n,x,fvec)
    else
        call g_fcn_rl(n,x,fvec)
    endif

    maxerr=maxval(abs(fvec))
    if(io>0)then
        write(io,'(" iter = ",i0," maxerr = ",f14.8)')gkmem%iter,maxerr
        write(0,'(" iter = ",i0," maxerr = ",f14.8)')gkmem%iter,maxerr
        write(io,'("********** diff_x **********")')
        write(io,'(5f14.8)')fvec
    endif

    if(maxerr<gkmem%maxerr)then
        gkmem%maxerr=maxerr
        if(io>0)then
            call gh5_wrt_wh_rl('WH_RL_BEST.h5')
        endif
    endif

    ! Stop criteria
    if(maxerr<gkmem%rtol)iflag=-1
    ! trade-off
    if(maxerr<gkmem%rtol*10)then
        if(gkmem%giter==-1)then
            gkmem%giter=gkmem%iter
        elseif(gkmem%iter-gkmem%giter>n*3)then
            iflag=-1
        endif
    endif

    if(gkmem%iter>=gkmem%max_iter)iflag=-1

    call set_time_point(2,3)
    call print_time_usage('g_fcn',3,io)
    return

end subroutine g_fcn


subroutine g_fcn_rl(n,x,fvec)
    use gprec
    use gkernel
    implicit none
    integer,intent(in)::n
    real(8),intent(in)::x(n)
    real(8),intent(out)::fvec(n)

    ! x to {r_coef, la1_coef}
    wh%r_coef=transfer(x(1:gkmem%nx1),wh%r_coef)
    wh%la1_coef=x(1+gkmem%nx1:)
    call g_rl_onecycle()
    fvec(1:gkmem%nx1)=transfer(wh%r_coef,fvec(1:gkmem%nx1))-x(1:gkmem%nx1)
    fvec(1+gkmem%nx1:)=wh%nks_coef-wh%ncv_coef
    return

end subroutine g_fcn_rl


subroutine g_rl_onecycle()
    use gprec
    use gmpi
    use gkernel
    use bandstru
    use dcstd
    implicit none
    real(q) etot
    integer io

    io=gp%io
    call hm_expand_all_general(wh%r,wh%r_coef,wh%hm_r,-1,.false.)
    call hm_expand_all_herm(wh%la1,wh%la1_coef,wh%hm_l,-1,.false.)
    call modify_r_la1_frozen(1,30._q)

    if(io>0)then
        call output_matrices('r-in',wh%r,wh%na2112,wh%num_imp, &
                &wh%na2_imp,io,0)
        call output_matrices('la1-in',wh%la1,wh%na2112,wh%num_imp, &
                &wh%na2_imp,io,0)
    endif

    call map_wh_bnd_matrix(wh%r,bnd%r,.false.)
    call map_wh_bnd_matrix(wh%la1,bnd%la1,.false.)
    call calc_band_all(io)
    call gutz_fermi(io)
    call bnd_modify_frozen()
    call modify_r_la1_frozen(1,bnd%ef)
    call calc_mup_dn(io)
    call calc_nks()
    call calc_nks_pp(io)
    call eval_sl_vec_all(1,io)
    call map_wh_bnd_matrix(wh%nks,bnd%nks,.false.)
    call calc_da0()
    call calc_isimix(io)
    call calc_da(io)
    call calc_vdc_list()
    call calc_lambdac(io)
    call solve_hembed_list(io)
    call calc_r01_pp(io)
    call calc_ncvar_pp(io)
    call calc_ncphy_pp(io)
    call eval_sl_vec_all(2,io)
    call calc_total_energy()
    return

end subroutine g_rl_onecycle


subroutine g_fcn_n(n,x,fvec)
    use gprec
    use gkernel
    implicit none
    integer,intent(in)::n
    real(8),intent(in)::x(n)
    real(8),intent(out)::fvec(n)

    integer io

    io=gp%io
    wh%nks_coef=x
    call hm_expand_all_herm(wh%nks,wh%nks_coef,wh%hm,-1,.false.)
    call calc_nks_tot(io)
    if(io>0)then
        call output_matrices('nks-in',wh%nks,wh%na2112,wh%num_imp, &
                &wh%na2_imp,io,0)
    endif
    call calc_la1_hf_list()
    call add_vdc_to_la1_list()
    if(associated(wh%vext).and.wh%ivext>0)then
        wh%la1=wh%la1+wh%vext
    endif
    call calc_la1_pp(io,'la1-inp')
    call map_wh_bnd_matrix(wh%la1,bnd%la1,.false.)
    call calc_band_all(io)
    call gutz_fermi(io)
    call calc_mup_dn(io)
    call calc_nks()
    call calc_nks_pp(io)
    !! True for HF calculation.
    wh%nc_phy=wh%nks
    call eval_sl_vec_all(1,io)
    call map_wh_bnd_matrix(wh%nks,bnd%nks,.false.)
    fvec=wh%nks_coef-x
    call calc_total_energy()
    return

end subroutine g_fcn_n

