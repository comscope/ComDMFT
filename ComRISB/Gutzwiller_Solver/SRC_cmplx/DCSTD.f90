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

module dcstd
    use gprec
    use gconstant
    use warehouse
    use ghdf5_base
    implicit none

    type dc_std
        integer::mode=1
        real(q),allocatable::e(:),vpot(:,:),u_avg(:),j_avg(:),nelf(:,:)
    end type dc_std

    type(dc_std)::dc


    contains


    subroutine init_dc_std(io)
    integer,intent(in)::io

    allocate(dc%vpot(2,wh%num_imp),dc%u_avg(wh%num_imp),dc%j_avg(wh%num_imp),&
            &dc%nelf(2,wh%num_imp),dc%e(wh%num_imp))
    dc%u_avg=0; dc%j_avg=0; dc%nelf=0
    call set_nelf_list(io)
    return

    end subroutine init_dc_std


    subroutine set_nelf_list(io)
    integer,intent(in)::io

    integer i
    logical lexist

    call gh5_open_r('GPARAM.h5',f_id)
    call gh5_read(dc%mode,'/dc_mode',f_id)
    if(dc%mode>0)then
        call gh5_read(dc%u_avg,wh%num_imp,'/dc_u_avg',f_id)
        call gh5_read(dc%j_avg,wh%num_imp,'/dc_j_avg',f_id)
    endif
    if(dc%mode>1)then
        call h5lexists_f(f_id,'/dc_nelf_list',lexist,gh5_err)
        if(lexist)then
            call gh5_read(dc%nelf,2,wh%num_imp,'/dc_nelf_list',f_id)
        else
            do i=1,wh%num_imp
                dc%nelf(:,i)=wh%co(i)%net
            enddo
        endif
    endif
    call gh5_close(f_id)

    if(io>0)then
        select case(dc%mode)
        case(0)
            write(io,'(" no double counting.")')
        case(1)
            write(io,'(" standard fully localized limit (FLL) dc.")')
        case(2)
            write(io,'(" fixed double counting.")')
        case(12)
            write(io,'(" standard FLL dc with dc only updated at the &
                    &electron density cycles.")')
        case default
            stop ' error in GDC_NELF.INP: dc%mode not defined!'
        end select

        if(dc%mode>1)then
            write(io,'(" input nelf:")')
            write(io,'(8x,5(f8.3,2x))')(sum(dc%nelf(:,i)), i=1,wh%num_imp)
        endif

        if(dc%mode>0)then
            write(io,'(" average hubbard u list:")')
            write(io,'(8x,5(f8.3,2x))')dc%u_avg
            write(io,'(" average hund j list:")')
            write(io,'(8x,5(f8.3,2x))')dc%j_avg
        endif
    endif
    return

    end subroutine set_nelf_list


    subroutine update_nelf_list(io)
    integer,intent(in)::io

    integer idx(1),i
    real(q) ndiff(wh%num_imp)
    character(20) fmt

    do i=1,wh%num_imp
        ndiff(i)=sum(wh%co(i)%net-dc%nelf(:,i))
    enddo
    idx=maxloc(abs(ndiff))

    if(io>0)then
        write(io,'(" max nelf diff=",f14.8)')ndiff(idx(1))
    endif

    if(io>0)then
        call gh5_open_w('GDC_NELF_OUT.h5',f_id)
        call gh5_write(dc%nelf,2,wh%num_imp,'/dc_nelf_list_inp',f_id)
        do i=1,wh%num_imp
            dc%nelf(:,i)=wh%co(i)%net
        enddo
        call gh5_write(dc%nelf,2,wh%num_imp,'/dc_nelf_list_out',f_id)
        call gh5_close(f_id)
    endif
    return

    end subroutine update_nelf_list


    subroutine calc_vdc_list()
    integer i,j,na2,na
    real(q) ntot,ntotp

    ntot=1._q; ntotp=1._q
    do i=1,wh%num_imp
        if(dc%mode/=1)then
            ntot=sum(dc%nelf(:,i))
            ntotp=sum(wh%co(i)%net)
        endif
        dc%nelf(:,i)=wh%co(i)%net*ntot/ntotp
    enddo
    do i=1,wh%num_imp
        ntot=sum(dc%nelf(:,i))
        do j=1,2
            dc%vpot(j,i)=dc%u_avg(i)*(ntot-.5_q)- &
                    &dc%j_avg(i)*(dc%nelf(j,i)-0.5_q)
        enddo
    enddo
    ! add to co%la2
    wh%la2=0
    do i=1,wh%num_imp
        na2=wh%na2_imp(i); na=na2/2
        if(.not.associated(wh%db2sab))then
            do j=1,na
                ! spin-index-faster convention
                wh%co(i)%la2(2*j-1,2*j-1)=-dc%vpot(1,i)
                wh%co(i)%la2(2*j,2*j)=-dc%vpot(2,i)
            enddo
        else
             do j=1,na
                ! orbital-index-faster convention
                wh%co(i)%la2(j,j)=-dc%vpot(1,i)
                wh%co(i)%la2(j+na,j+na)=-dc%vpot(2,i)
            enddo
            ! subsequent transformation to the symmetry-adapted basis
            call uhau(wh%co(i)%la2,wh%co(i)%db2sab,na2,na2)
        endif
    enddo
    return

    end subroutine calc_vdc_list


    subroutine add_vdc_to_la1_list()

    call calc_vdc_list()
    wh%la1=wh%la1+wh%la2
    return

    end subroutine add_vdc_to_la1_list


    subroutine calc_edc_list()
    integer i,isp
    real(q) ntot

    if(dc%mode==1)then
        do i=1,wh%num_imp
            dc%nelf(:,i)=wh%co(i)%net
        enddo
    endif
    do i=1,wh%num_imp
        ntot=sum(dc%nelf(:,i))
        dc%e(i)=dc%u_avg(i)*ntot*(ntot-1)/2
        do isp=1,2
            dc%e(i)=dc%e(i)-dc%j_avg(i)*dc%nelf(isp,i)*(dc%nelf(isp,i)-1)/2+ &
                    &+dc%vpot(isp,i)*(wh%co(i)%net(isp)-dc%nelf(isp,i))
        enddo
    enddo
    return
    
    end subroutine calc_edc_list


    subroutine output_energies_dc(io)
    integer,intent(in)::io

    integer i

    write(io,'(" impurity-wise interaction dc energy:")')
    write(io,'(4x,5f14.7)')(dc%e(i), i=1,wh%num_imp)
    write(io,'(" total U-dc energy = ",f0.7)')sum(dc%e)
    return

    end subroutine output_energies_dc



end module dcstd
