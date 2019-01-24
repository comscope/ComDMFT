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

module ghfscf
    use gprec
    use gdiis
    use gutil
    use gparam
    implicit none

    type density_matrix
        complex(q),pointer::a(:,:)=>null()
    end type density_matrix

    type(density_matrix),allocatable::dm_list(:)

    contains

    
    subroutine ghf_init(nhf)
    integer,intent(in)::nhf

    allocate(dm_list(nhf))
    return

    end subroutine ghf_init


    ! with HF specific error vectors.
    subroutine ghf_scf(ihf,norbs,h1e,v2e,dm,e_tot,iter,flevelshift)
    integer,intent(in)::norbs,ihf
    integer,intent(out)::iter
    real(q),intent(out)::e_tot
    real(q),intent(in)::flevelshift
    complex(q),intent(in)::h1e(norbs,norbs),v2e(norbs,norbs,norbs,norbs)
    complex(q),target,intent(out)::dm(norbs,norbs)

    integer,parameter::max_cycle=7000
    real(q),parameter::conv_tol=1.d-5
    integer norb,i
    type(diis_mem) :: dmem
    real(q) mo_energy(norbs),e_tot_last
    logical lconverge
    complex(q) mo_coef(norbs,norbs),dm_last(norbs,norbs),vhf(norbs,norbs),&
            &v_j2e(norbs,norbs,norbs,norbs),errvec(norbs,norbs)

    lconverge=.false.
    norb=norbs/2

    if(.not.associated(dm_list(ihf)%a))then
        ! Get initial guess by 1e
        mo_coef=h1e
        call hermev('v','l',mo_coef,mo_energy,norbs)
    
        if(g_write>=5)then
            write(0,'(" mo_energy:")')
            write(0,'(4x,5f12.6)')mo_energy
        endif
    
        ! Get density matrix, dm_{i,j} = <c_j^\dagger c_i>
        call zgemm('n','c',norbs,norbs,norb,z1,mo_coef,norbs,mo_coef,norbs, &
                &z0,dm,norbs)
    else
        dm = dm_list(ihf)%a
    endif
    call get_coul_exchange(norbs,v2e,v_j2e)
    vhf=0
    call get_veff(norbs,vhf,dm,v_j2e)
    call get_energy_tot(e_tot,norbs,h1e,vhf,dm)
    call diis_vector_init(dmem,norbs*norbs)
    ! begin scf cycles
    iter=1
    do 
        dm_last=dm
        e_tot_last=e_tot        
        mo_coef=h1e+vhf
        if(iter>1.and.(.not.lconverge))then
            ! J. Mol. Struct. 114, 31-34
            ! PCCP, 4, 11
            ! GEDIIS, JCTC, 2, 835
            ! C2DIIS, IJQC, 45, 31
            ! SCF-EDIIS, JCP 116, 8255
            ! error vector = SDF-FDS
            ! = F_ai ~ (S-SDS)*S^{-1}FDS = FDS - SDFDS ~ FDS-SDF in converge
            errvec=dm
            call annxb('n','n',errvec,mo_coef,norbs,1)
            errvec=conjg(transpose(errvec))-errvec
            call gdiis_update(dmem,mo_coef(1,1),errvec(1,1))
        endif

        ! level shift
        if(flevelshift>1.e-4_q.and.(.not.lconverge))then
            mo_coef=mo_coef-dm*flevelshift
            forall(i=1:norbs)mo_coef(i,i)=mo_coef(i,i)+flevelshift
        endif

        call hermev('v','l',mo_coef,mo_energy,norbs)
        ! Get density matrix, dm_{i,j} = <c_j^\dagger c_i>
        call zgemm('n','c',norbs,norbs,norb,z1,mo_coef,norbs,mo_coef, &
                &norbs,z0,dm,norbs)
        call get_veff(norbs,vhf,dm,v_j2e,dm_last)
        call get_energy_tot(e_tot,norbs,h1e,vhf,dm)

        if(g_write>=5)then
            write(0,'(" mo_energy:")')
            write(0,'(4x,5f12.6)')mo_energy
            write(0,'(" Iter = ",I4," HF dm error =",f16.10," etot =",f16.8)') &
                    &iter,maxval(abs(dm-dm_last)),e_tot
        endif

        if(lconverge)exit
        if(abs(e_tot-e_tot_last)<conv_tol.or.iter==max_cycle)then
            lconverge=.true. ! One more iteration without fixing.
        else
            iter=iter+1
        endif
    enddo
    if(iter==max_cycle)then
        stop ' hf does not converge!'
    endif
    if(.not.associated(dm_list(ihf)%a))then
        allocate(dm_list(ihf)%a(norbs,norbs))
    endif
    dm_list(ihf)%a=dm
    ! Transpose to dm_{i,j}=<c_i^\dagger c_j>
    dm=transpose(dm)
    return

    end subroutine ghf_scf


    ! with density matrix difference as the error vector.
    subroutine ghf_scf1(norbs,h1e,h1e_init,v2e,dm,e_tot,iter)
    integer,intent(in)::norbs
    integer,intent(out)::iter
    real(q),intent(out)::e_tot
    complex(q),intent(in)::h1e(norbs,norbs),v2e(norbs,norbs,norbs,norbs), &
            &h1e_init(norbs,norbs)
    complex(q),target,intent(out)::dm(norbs,norbs)

    integer,parameter::max_cycle=7000
    real(q),parameter::conv_tol=1.d-10
    integer norb
    type(diis_mem) :: dmem
    real(q) mo_energy(norbs),e_tot_last
    complex(q) mo_coef(norbs,norbs),dm_last(norbs,norbs),vhf(norbs,norbs),&
            &v_j2e(norbs,norbs,norbs,norbs),errvec(norbs,norbs)

    norb=norbs/2
    ! Get initial guess by 1e
    mo_coef=h1e_init
    call hermev('v','l',mo_coef,mo_energy,norbs)
    ! Get density matrix, dm_{i,j} = <c_j^\dagger c_i>
    call zgemm('n','c',norbs,norbs,norb,z1,mo_coef,norbs,mo_coef,norbs, &
            &z0,dm,norbs)
    call get_coul_exchange(norbs,v2e,v_j2e)
    vhf=0
    call get_veff(norbs,vhf,dm,v_j2e)
    call get_energy_tot(e_tot,norbs,h1e,vhf,dm)
    call diis_vector_init(dmem,norbs*norbs)
    ! begin scf cycles
    do iter=1,max_cycle
        dm_last=dm
        e_tot_last=e_tot        
        mo_coef=h1e+vhf
        call hermev('v','l',mo_coef,mo_energy,norbs)
        ! Get density matrix, dm_{i,j} = <c_j^\dagger c_i>
        call zgemm('n','c',norbs,norbs,norb,z1,mo_coef,norbs,mo_coef, &
                &norbs,z0,dm,norbs)
        errvec=dm-dm_last
        call gdiis_update(dmem,dm(1,1),errvec(1,1))
        call get_veff(norbs,vhf,dm,v_j2e,dm_last)
        call get_energy_tot(e_tot,norbs,h1e,vhf,dm)

        if(g_write>=5)then
            write(0,'(" Iter =",I4," HF dm error =",f16.10," etot =",f16.8)')&
                    &iter,maxval(abs(errvec)),e_tot
        endif

        if(abs(e_tot-e_tot_last)<conv_tol)exit
    enddo
    ! Transpose to dm_{i,j}=<c_i^\dagger c_j>
    dm=transpose(dm)
    return

    end subroutine ghf_scf1


    subroutine get_veff(norbs,vhf,dm,vj2e,dm_last)
    integer,intent(in)::norbs
    complex(q),intent(in)::dm(norbs,norbs),vj2e(norbs,norbs,norbs,norbs)
    complex(q),intent(inout)::vhf(norbs,norbs)
    complex(q),optional,intent(in)::dm_last(norbs,norbs)

    complex(q) dm_t(norbs,norbs)

    if(present(dm_last))then
        dm_t=transpose(dm-dm_last)
    else
        dm_t=transpose(dm)
    endif
    call get_hf_pot(norbs,dm_t,vj2e,vhf)
    return

    end subroutine get_veff


    subroutine get_energy_tot(e_tot,norbs,h1e,vhf,dm)
    real(q),intent(out)::e_tot
    integer,intent(in)::norbs
    complex(q),intent(in)::h1e(norbs,norbs),vhf(norbs,norbs), &
            &dm(norbs,norbs)

    complex(q) dm_t(norbs,norbs),tmp(norbs,norbs)

    dm_t=transpose(dm)
    e_tot=sum(h1e*dm_t)+sum(vhf*dm_t)/2.d0
    return

    end subroutine get_energy_tot


end module ghfscf
