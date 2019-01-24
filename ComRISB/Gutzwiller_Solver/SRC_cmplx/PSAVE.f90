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

module psave
    use gprec
    use gmpi
    use ghdf5_base
    use warehouse
    use bandstru
    use gkernel
    implicit none
    
    contains


    subroutine postcheck(io)
    integer,intent(in)::io


    if(io<0)return
    write(io,'(" maximal local symmetrization error = ", f10.5)')wh%symerr
    if(wh%symerr>0.05_q)then
        write(io,'(&
                &" WARNING: LARGE LOCAL SYMM ERROR!",/, &
                &"   IF YOU DO NOT UNDERSTAND WHAT YOU ARE DOING,",/, &
                &"   REPORT TO YONGXIN YAO AT CYGUTZ@GMAIL.COM!")')
    endif
    return

    end subroutine postcheck


    subroutine postsave()

    call postsave_glog()
    call postsave_gbands()
    return

    end subroutine postsave


    subroutine postsave_glog()

    ! calculate nks,d0 before symmetrization
    call calc_nks()
    call calc_da0()

    ! h1e
    call gh5_open_r('GPARAMBANDS.h5',f_id)
    call gh5_read_wh_h1e(bnd%nspin_in)
    call gh5_close(f_id)
    call rotate_h1e_list()

    if(gp%myrank==gp%master)then
        call gh5_open_w('GLOG.h5',f_id)
        call gh5_create_impurity_groups(f_id)
        call gh5_wrt_wh_matrix_list('/NKS_PRE_SYM',wh%nks)
        call gh5_wrt_wh_matrix_list('/D0_PRE_SYM',wh%d0)
        call gh5_wrt_wh_matrix_list('/H1E_PRE_SYM',wh%h1e)
        call gh5_write(gkmem%maxerr,'/rl_maxerr',f_id)
        call gh5_write(gkmem%etot,'/etot_model',f_id)
        call gh5_write(wh%eu2,wh%num_imp,'/eu2_list',f_id)
        call gh5_write(wh%et1,wh%num_imp,'/et1_list',f_id)
        call gh5_write(bnd%r,wh%nasotot,wh%nasotot,wh%nspin, &
                &"/BND_R",f_id)
        call gh5_write(bnd%la1,wh%nasotot,wh%nasotot,wh%nspin, &
                &"/BND_LAMBDA",f_id)
        call gh5_write(bnd%nrl,wh%nasotot,wh%nasotot,wh%nspin, &
                &"/BND_NRL",f_id)
        call gh5_write(bnd%nc_phy,wh%nasotot,wh%nasotot,wh%nspin, &
                &"/BND_NPHY",f_id)
    endif

    call calc_nks_pp(0)
    call calc_da0_pp(0)
    call calc_herm_matrices_pp(wh%h1e,'h1e',wh%hm,.false.,0,-1)

    if(gp%myrank==gp%master)then
        call gh5_wrt_wh_matrix_list('/NKS_SYM',wh%nks)
        call gh5_wrt_wh_matrix_list('/D0_SYM',wh%d0)
        call gh5_wrt_wh_matrix_list('/H1E_SYM',wh%h1e)
        call gh5_wrt_wh_matrix_list('/LA1_SYM',wh%la1)
        call gh5_wrt_wh_matrix_list('/D_SYM',wh%d)
        call gh5_write(bnd%ef,'/e_fermi',f_id)
        call gh5_write(wh%nspin,'/nspin',f_id)
        call gh5_write(wh%nasotot,'/nasotot',f_id)
        ! Other commonly used data for analysis 
        call gh5_wrt_wh_matrix_list('/R',wh%r)
        call gh5_wrt_wh_matrix_list('/NC_PHY',wh%nc_phy)
        call gh5_close(f_id)
    endif
    return

    end subroutine postsave_glog


    subroutine postsave_gbands()
    integer ivec,isp,ispp,ikp,ikpl,kbase,iks,nbtot,isym,nbands
    complex(q),pointer::p_vk(:,:)
    
    ikpl=0
    do ivec=1,gp%nvec
        kbase=ikpl
        call gh5_open_w(file_name(gp%kvec(ivec,1),'GBANDS'),f_id)
        call gh5_write(gp%kvec(ivec,3)+1,"/IKP_START",f_id)
        call gh5_write(gp%kvec(ivec,3)+gp%kvec(ivec,2),"/IKP_END",f_id)
        do isp=1,wh%nspin
            ispp=min(isp,bnd%nspin_in)
            ikpl=kbase
            call gh5_create_group("/ISPIN_"//trim(int_to_str(isp)),f_id)
            do iks=1,gp%kvec(ivec,2)
                ikp=gp%kvec(ivec,3)+iks
                if(gp%ipar==1)then ! openmp
                    ikpl=ikp
                else
                    ikpl=ikpl+1
                endif
                call gh5_create_group("/ISPIN_"//trim(int_to_str(isp))// &
                        &"/IKP_"//trim(int_to_str(ikp)),f_id)
                nbtot=bnd%ne(1,ikp,ispp)
                call gh5_write(bnd%ek(:,ikp,isp),nbtot,"/ISPIN_"// &
                        &trim(int_to_str(isp))// &
                        &"/IKP_"//trim(int_to_str(ikp))// &
                        &"/ek",f_id)
                nbands=bnd%ne(3,ikp,ispp)-bnd%ne(2,ikp,ispp)+1
                do isym=1,sym%nop
                    call gh5_create_group("/ISPIN_"//trim(int_to_str(isp))//&
                            &"/IKP_"//trim(int_to_str(ikp))// &
                            &"/ISYM_"//trim(int_to_str(isym)),f_id)
                    p_vk(1:nbands,1:nbands)=>bnd%vk(1:nbands**2,isym,ikpl,isp)
                    call gh5_write(p_vk,nbands,nbands,"/ISPIN_"// &
                            &trim(int_to_str(isp))// &
                            &"/IKP_"//trim(int_to_str(ikp))// &
                            &"/ISYM_"//trim(int_to_str(isym))// &
                            &"/EVEC",f_id)
                enddo
            enddo
        enddo
        call gh5_close(f_id)
    enddo
    nullify(p_vk)
    return

    end subroutine postsave_gbands


end module psave
