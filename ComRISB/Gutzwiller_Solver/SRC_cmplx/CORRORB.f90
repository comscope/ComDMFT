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

module corrorb
    use gprec
    use gutil
    use gconstant, only:d1,zi
    implicit none
      
    type corr_orb
        !< number of correlated orbital, spin-orbital, dim*ispin
        integer dim,dim2,dimso
        real(q) :: net(2)
        integer :: dim_hs_l
        complex(q),pointer :: nks(:,:),nc_var(:,:),nc_phy(:,:),isimix(:,:), &
                &r(:,:),r0(:,:),z(:,:),d(:,:),d0(:,:), &
                &la1(:,:),la2(:,:),h1e(:,:),vext(:,:)=>null(), &
                &sx(:,:),sy(:,:),sz(:,:),lx(:,:),ly(:,:),lz(:,:),&
                &hs_l(:,:,:),v2e(:,:,:,:),r_coef(:),v_j2e(:,:,:,:)
        real(q),pointer :: nks_coef(:),la1_coef(:),ncv_coef(:)
        real(q) :: s_val(3,2),l_val(3,2)  !< one-body operator
        complex(q),pointer :: db2sab(:,:)=>null()
        integer,pointer :: m_struct(:,:)=>null()
    end type corr_orb
      
    contains
    

    subroutine eval_co_sl_vec(co,mode)
    type (corr_orb) co
    integer mode
      
    integer i
    complex(q) n_(co%dim2,co%dim2)
      
    if(mode==1)then
        n_=co%nks
    else
        n_=co%nc_phy
    endif
    ! Tr(\rho A) = \sum_{ij}{\rho_ij*A_ji} = \sum_{ij}{\rho^{t}_ji*A_ji}
    co%s_val(1,mode)=sum(co%sx*n_)
    co%s_val(2,mode)=sum(co%sy*n_)
    co%s_val(3,mode)=sum(co%sz*n_)
    co%l_val(1,mode)=sum(co%lx*n_)
    co%l_val(2,mode)=sum(co%ly*n_)
    co%l_val(3,mode)=sum(co%lz*n_)
    return
      
    end subroutine eval_co_sl_vec
      

    subroutine calc_co_lambdac(co)
    type(corr_orb) co
      
    integer i
    real(q) res
    complex(q) rd(co%dim2,co%dim2),pn(co%dim2,co%dim2),h(co%dim2,co%dim2)

    rd=matmul(co%r,transpose(co%d))
    do i=1,co%dim_hs_l
        h=transpose(co%hs_l(:,:,i))
        call pfa_pa(co%nks,pn,h,co%dim2,dsimix,dpsimix) ! p f / p d_n
        res=-real(sum(pn*rd),q)*2
        co%la2=co%la2+co%hs_l(:,:,i)*res
    enddo
    co%la2=co%la2-co%la1
    return
      
    end subroutine calc_co_lambdac
    

    subroutine calc_co_net(co)
    type(corr_orb) co

    integer i
    complex(q) zbuf(co%dim2,co%dim2),sab2db(co%dim2,co%dim2)

    co%net=0
    if(associated(co%db2sab))then
        ! transposition is needed for nks.
        zbuf=transpose(co%nks)
        sab2db=conjg(transpose(co%db2sab))
        call uhau(zbuf,sab2db,co%dim2,co%dim2)
        ! orbital-index faster
        do i=1,co%dim
            co%net(1)=co%net(1)+zbuf(i,i)
            co%net(2)=co%net(2)+zbuf(i+co%dim,i+co%dim)
        enddo
    else
        ! spin-index faster
        do i=1,co%dim
            co%net(1)=co%net(1)+co%nks(2*i-1,2*i-1)
            co%net(2)=co%net(2)+co%nks(2*i,2*i)
        enddo
    endif
    return

    end subroutine calc_co_net


    !*************************************************************************
    ! isimix = (nks(1-nks))^(-1/2)
    !*************************************************************************
    subroutine calc_co_isimix(co)
    type(corr_orb),intent(inout)::co

    complex(q) xn(co%dim2,co%dim2)
      
    xn=co%nks; xn=xn-matmul(xn,xn)
    call atofa(xn,co%isimix,co%dim2,-12,d1,.true.)
    return
      
    end subroutine calc_co_isimix
      
    !*************************************************************************
    ! d = (nks(1-nks))^(-1/2) d0
    !*************************************************************************
    subroutine co_d0_to_d(co)
    type(corr_orb) co
      
    co%d=matmul(co%isimix,co%d0)
    return
      
    end subroutine co_d0_to_d


    subroutine co_nks_patch_order(co,ispo,iso)
    integer,intent(in)::ispo,iso
    type(corr_orb) co


    if(ispo==1)then
        co%nks(1+co%dimso:,1+co%dimso:)=co%nks(1:co%dimso,1:co%dimso)
    endif
    if(iso==1)then
        call orbital_spin_trans(co%nks,co%dim2,.true.)
    endif
    return

    end subroutine co_nks_patch_order


    subroutine calc_co_la1_hf(co)
    type(corr_orb) co

    co%la1=0
    call get_hf_pot(co%dim2,co%nks,co%v_j2e,co%la1)
    co%la1=co%la1+co%h1e
    return

    end subroutine calc_co_la1_hf


    subroutine calc_co_eu2_hf(co, eu2)
    type(corr_orb) co
    real(q),intent(out)::eu2

    complex(q) vhf(co%dim2,co%dim2)

    vhf=0
    call get_hf_pot(co%dim2,co%nks,co%v_j2e,vhf)
    eu2=sum(vhf*co%nks)/2
    return

    end subroutine calc_co_eu2_hf


      
end module corrorb
