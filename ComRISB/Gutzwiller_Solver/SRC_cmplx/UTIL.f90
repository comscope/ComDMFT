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

module gutil
    use gprec
    use gprimme
    use gconstant,only: small,z1,z0,d1,d0,rlbound,rubound
    implicit none
    private
    public :: orbital_spin_trans, a_trans, get_ylm_cr, get_ylm_cr_all,  &
            & inv_aplx, wrt_a, hermev, inv, chk_unitary, &
            & get_hm_expand, uhau, anmxbmm, annxb, trace_a, atofa, pfa_pa, &
            & h_regularize, calc_loewner, dsimix, dpsimix, &
            & file_name, int_to_str, fermi_fun, gauss_fun, set_range, &
            & set_linear_simp, max_interval, lkey_exist, get_eri, &
            & act_state, get_rn_from_embeddm, output_matrices, &
            & chk_eigens_matrix_list, set_binoms, set_fock_state_indices, &
            & primme_diag, dprimme_diag_chk, zprimme_diag_chk,&
            & get_coul_exchange, get_hf_pot, &
            & calc_state_sz, sum_empty_slots_fs,chk_sum_first_empty_slots_fs, &
            & chk_i_in_list, hermev_blk, get_determinant_a, &
            & set_occupied_orbitals, set_skip_orbitals, &
            & calc_fock_state_coef, phi_mat_to_vec, &
            & phi_vec_to_mat
     
    interface primme_diag 
        module procedure dprimme_diag, zprimme_diag
    end interface primme_diag

    interface orbital_spin_trans
        module procedure dspin_blk1_trans, zspin_blk2u_trans, &
                &dspin_blk2_trans, zspin_blk2_trans
    end interface orbital_spin_trans
      
    interface a_trans
        module procedure a4_trans_rc, a4_trans_c, za2_trans, da2_trans
    end interface a_trans
      
    interface inv_aplx
        module procedure zinv_aplx
    end interface inv_aplx
     
    interface get_coul_exchange
        module procedure zget_coul_exchange
    end interface get_coul_exchange

    interface get_hf_pot
        module procedure zget_hf_pot
    end interface get_hf_pot
 
    interface wrt_a
        module procedure dwrt_ann, dwrt_an, zwrt_ann, zwrt_an
    end interface wrt_a
      
    interface hermev
        module procedure zheev_, dsyev_, zheevx_, zheevx_i, dsyevx_, &
                &dsyevx_i, ddiag_jzjp, zdiag_jzjp, zhpev_
    end interface hermev

    interface hermev_blk
        module procedure zheev_blk
    end interface hermev_blk
      
    interface hermevd
        module procedure zheevd_, dsyevd_
    end interface hermevd
      
    interface inv
        module procedure dinv_, zinv_
    end interface inv
      
    interface get_hm_expand
        module procedure get_dhm_expand, get_zhm_expand, &
                &get_dhm_expand_sym, get_zhm_expand_sym
    end interface get_hm_expand
      
    interface uhau
        module procedure duhau, zuhau
    end interface uhau
      
    interface anmxbmm
        module procedure danmxbmm, zanmxbmm
    end interface anmxbmm
      
    interface annxb
        module procedure dannxbnm, zannxbnm, zannxbnn, dannxbnn, &
                &dannxbnmc, zannxbnmc
    end interface annxb
      
    interface get_determinant_a
        module procedure dget_determinant_a, zget_determinant_a
    end interface get_determinant_a

    interface trace_a
        module procedure trace_dann, trace_zann
    end interface trace_a
      
    interface atofa
        module procedure zatofa, datofa, zatofa1, datofa1
    end interface atofa
      
    interface pfa_pa
        module procedure dpfa_pa, zpfa_pa
    end interface pfa_pa
      
    interface h_regularize
        module procedure dh_regularize, zh_regularize
    end interface h_regularize

    interface get_eri
        module procedure dget_eri, zget_eri
    end interface get_eri
      
    interface get_rn_from_embeddm
        module procedure dget_rn_from_embeddm, zget_rn_from_embeddm
    end interface

    interface output_matrices
        module procedure dout_matrix_list,zout_matrix_list,dout_matrix_block,&
                &zout_matrix_block
    end interface output_matrices

    interface chk_eigens_matrix_list
        module procedure chk_eigens_dmatrix_list,chk_eigens_zmatrix_list
    end interface chk_eigens_matrix_list

    interface phi_mat_to_vec
        module procedure dphi_mat_to_vec,dphi_mat_to_vec2, &
                &zphi_mat_to_vec,zphi_mat_to_vec2
    end interface phi_mat_to_vec

    interface phi_vec_to_mat
        module procedure dphi_vec_to_mat,dphi_vec_to_mat2, &
                &zphi_vec_to_mat,zphi_vec_to_mat2
    end interface phi_vec_to_mat

    interface calc_fock_state_coef
        module procedure dcalc_fock_state_coef, zcalc_fock_state_coef
    end interface calc_fock_state_coef

    contains
      
      
    !*************************************************************************
    ! orbital fast <-> spin fast
    !*************************************************************************
    subroutine dspin_blk1_trans(a,n,lsfast)
    integer n
    real(q) a(n)
    logical lsfast
      
    integer i,iadd
    real(q) u(n,n)
      
    u=0; iadd=0
    do i=1,n,2
        iadd=iadd+1
        u(iadd,i)=1._q
    enddo
    do i=2,n,2
        iadd=iadd+1
        u(iadd,i)=1._q
    enddo
    if(lsfast)then; a=matmul(transpose(u),a)
    else          ; a=matmul(u,a)
    endif
    return
      
    end subroutine dspin_blk1_trans
      
      
    !*************************************************************************
    ! orbital fast <-> spin fast
    !*************************************************************************
    subroutine zspin_blk2u_trans(a,n,m,lsfast,luright)
    integer n,m
    complex(q) a(n,m)
    logical lsfast
    logical,optional::luright
      
    integer i,iadd,nm
    logical lur
    complex(q),allocatable::u(:,:)
      
    lur=.false.
    if(present(luright)) lur=luright
    if(lur)then; nm=m; else; nm=n; endif
    allocate(u(nm,nm))
    u=0; iadd=0
    do i=1,nm,2
        iadd=iadd+1
        u(iadd,i)=1._q
    enddo
    do i=2,nm,2
        iadd=iadd+1
        u(iadd,i)=1._q
    enddo
    if(lsfast)then
        if(.not.lur)then
            call zannxbnm('c',u,a,n,m)
        else
            call zanmxbmm('n',a,u,n,m)
        endif
    else
        if(.not.lur)then
            call zannxbnm('n',u,a,n,m)
        else
            call zanmxbmm('c',a,u,n,m)
        endif
    endif
    deallocate(u)
    return
      
    end subroutine zspin_blk2u_trans
      
      
    !*************************************************************************
    ! orbital fast <-> spin fast
    !*************************************************************************
    subroutine dspin_blk2_trans(a,n,lsfast)
    integer n
    real(q) a(n,n)
    logical lsfast
      
    integer i,iadd
    real(q) u(n,n)
      
    u=0; iadd=0
    do i=1,n,2
        iadd=iadd+1
        u(iadd,i)=1._q
    enddo
    do i=2,n,2
        iadd=iadd+1
        u(iadd,i)=1._q
    enddo
    if(lsfast)then
        call duhau(a,u,n,n,trul='c',trur='n')
    else
        call duhau(a,u,n,n,trul='n',trur='c')
    endif
    return
      
    end subroutine dspin_blk2_trans
      
      
    !*************************************************************************
    ! orbital fast <-> spin fast
    !*************************************************************************
    subroutine zspin_blk2_trans(a,n,lsfast)
    integer n
    complex(q) a(n,n)
    logical lsfast
      
    integer i,iadd
    complex(q) u(n,n)
      
    ! u =
    ! 1 0 0 0 0 0 ...
    ! 0 0 1 0 0 0 ...
    ! 0 0 0 0 1 0 ...
    ! .....
    ! 0 1 0 0 0 0 ...
    ! 0 0 0 1 0 0 ...
    ! 0 0 0 0 0 1 ...
    ! .....
    u=0; iadd=0
    do i=1,n,2
        iadd=iadd+1
        u(iadd,i)=1._q
    enddo
    do i=2,n,2
        iadd=iadd+1
        u(iadd,i)=1._q
    enddo
    if(lsfast)then
        call zuhau(a,u,n,n,trul='c',trur='n')
    else
        call zuhau(a,u,n,n,trul='n',trur='c')
    endif
    return
      
    end subroutine zspin_blk2_trans
      
      
    subroutine a4_trans_rc(a,n,u)
    integer n
    real(q) a(n,n,n,n)
    complex(q) u(n,n)
      
    real(q) res
    complex(q),allocatable::b(:,:,:,:)
      
    allocate(b(n,n,n,n)); b=dcmplx(a)
    call a4_trans_c(b,n,u)
    res=maxval(abs(aimag(b)))
    if(res.gt.1.e-6_q)then
        write(0,'(" error in a4_trans_rc: max aimag =",f12.6)')res; stop
    endif
    a=real(b,q)
    return
      
    end subroutine a4_trans_rc
      
      
    subroutine a4_trans_c(a,n,u)
    integer,intent(in)::n
    complex(q),intent(inout)::a(n,n,n,n)
    complex(q),intent(in)::u(n,n)
      
    integer i1,i2
      
    do i1=1,n
        do i2=1,n
            call za2_trans(a(:,:,i1,i2),n,u,1)
        enddo
    enddo
    a=reshape(transpose(reshape(a,(/n*n,n*n/))),(/n,n,n,n/))
    do i1=1,n
        do i2=1,n
            call za2_trans(a(:,:,i1,i2),n,u,1)
        enddo
    enddo
    return
      
    end subroutine a4_trans_c
      
      
    subroutine za2_trans(a,n,u,mode)
    integer n,mode
    complex(q) a(n,n),u(n,n)
      
    if(mode>=0)then
        call zuhau(a,u,n,n)
    else
        call zuhau(a,u,n,n,trul='n',trur='c')
    endif
    return
      
    end subroutine za2_trans
      
      
    subroutine da2_trans(a,n,u,mode)
    integer n,mode
    real(q) a(n,n),u(n,n)
      
    if(mode>=0)then
        call duhau(a,u,n,n)
    else
        call duhau(a,u,n,n,trul='n',trur='c')
    endif
    return
      
    end subroutine da2_trans
      
    
    !*************************************************************************
    ! pre-calculate binomial(n,m) := n!/(m!(n-m)!)
    !*************************************************************************
    subroutine set_binoms(n,binoms)
    integer,intent(in)::n
    integer,intent(out)::binoms(n,0:n)

    real(q) fac(0:n)
    integer i,j,m0

    fac(0)=1._q
    do i=1,n
        fac(i)=i*fac(i-1)
    enddo
    do i=1,n
        do j=0,i
            binoms(i,j)=nint(fac(i)/fac(j)/fac(i-j))
        enddo
    enddo
    return

    end subroutine set_binoms


    subroutine set_fock_state_indices(norb,binom,idx,bs,ibs,bs_sz,sz_orb)
    integer,intent(in)::norb,binom(:,0:)
    integer,pointer,intent(inout)::idx(:),bs(:),ibs(:)
    real(q),pointer,optional,intent(inout)::bs_sz(:)
    real(q),pointer,optional,intent(in)::sz_orb(:)

    integer i,nbs,num1

    nbs=ishft(1,norb)
    allocate(idx(0:norb+1),bs(nbs),ibs(0:nbs-1))
    if(present(bs_sz))then
        allocate(bs_sz(nbs))
    endif
    idx(0:1)=1
    do i=2,norb+1
        idx(i)=idx(i-1)+binom(norb,i-2)
    enddo
    do i=0,nbs-1
        num1=popcnt(i)+1
        bs(idx(num1))=i
        if(present(bs_sz))then
            if(.not.present(sz_orb))then
                call calc_state_sz(i,norb,bs_sz(idx(num1)),0)
            else
                call calc_state_sz2(i,norb,bs_sz(idx(num1)),sz_orb)
            endif
        endif
        idx(num1)=idx(num1)+1
    enddo
    do i=1,nbs
        ibs(bs(i))=i
    enddo
    return

    end subroutine set_fock_state_indices


    !*************************************************************************
    ! complex to real spherical harmonics transformation <c_ylm | r_ylm>
    !*************************************************************************
    subroutine get_ylm_cr(ucr,l)
    integer l
    complex(q) ucr(2*l+1,2*l+1)
      
    integer m,m_
      
    ucr=0
    do m=1,2*l+1
        m_=m-l-1
        if(m_>0)then
            ucr( m_+l+1,m)=(-1)**m_/sqrt(2._q)
            ucr(-m_+l+1,m)=1/sqrt(2._q)
        elseif(m_==0)then
            ucr(l+1,l+1)=1
        else
            ucr( m_+l+1,m)= cmplx(0,1/sqrt(2._q),q)
            ucr(-m_+l+1,m)=-cmplx(0,(-1)**m_/sqrt(2._q),q)
        endif
    enddo
    return
      
    end subroutine get_ylm_cr
      
      
    subroutine get_ylm_cr_all(ylm_cr,lmax)
    integer lmax
    complex(q) ylm_cr(2*lmax+1,2*lmax+1,0:lmax)
      
    integer l,mdim
      
    ylm_cr=0
    do l=0,lmax
        mdim=2*l+1
        call get_ylm_cr(ylm_cr(1:mdim,1:mdim,l),l)
    enddo
    return
      
    end subroutine get_ylm_cr_all
      
      
    !< mode = 0: hermitian matrix; else: non-hermitian, use m*m^\dagger
    subroutine chk_eigens_dmatrix_list(str,am,n112,n2,dim_list,io,mode)
    integer,intent(in) :: n112,n2,io,dim_list(n2),mode
    character(*),intent(in) :: str
    real(q),target,intent(in) :: am(n112)

    integer i,dim2,ip,ibase
    real(q),target::tmp(maxval(dim_list)**2)
    real(q),pointer::p_tmp(:,:)
    real(q) w(maxval(dim_list))

    ibase=0
    do i=1,n2
        dim2=dim_list(i)
        tmp(1:dim2**2)=-am(ibase+1:ibase+dim2**2)
        p_tmp(1:dim2,1:dim2)=>tmp(1:dim2**2)
        if(mode/=0)then
            p_tmp=matmul(p_tmp,transpose(p_tmp))
        endif
        call dsyev_('n','l',p_tmp,w(1:dim2),dim2)
        w=-w
        if(mode/=0)w=sqrt(abs(w))
        if(io>0)then
            write(io,'(" imp=",i3," eigen values of ",a9,":")')i,str
            write(io,'(14f10.4)')w(1:dim2)
        endif
        ibase=ibase+dim2**2
    enddo
    nullify(p_tmp)
    return

    end subroutine chk_eigens_dmatrix_list


    subroutine chk_eigens_zmatrix_list(str,am,n112,n2,dim_list,io,mode)
    integer,intent(in) :: n112,n2,io,dim_list(n2),mode
    character(*),intent(in) :: str
    complex(q),intent(in) :: am(n112)

    integer i,dim2,ip,ibase
    complex(q),target::tmp(maxval(dim_list)**2)
    complex(q),pointer::p_tmp(:,:)
    real(q) w(maxval(dim_list))

    ibase=0
    do i=1,n2
        dim2=dim_list(i)
        tmp(1:dim2**2)=-am(ibase+1:ibase+dim2**2)
        p_tmp(1:dim2,1:dim2)=>tmp(1:dim2**2)
        if(mode/=0)then
            p_tmp=matmul(p_tmp,transpose(p_tmp))
        endif
        call zheev_('n','l',p_tmp,w(1:dim2),dim2)
        w=-w
        if(mode/=0)w=sqrt(abs(w))
        if(io>0)then
            write(io,'(" imp=",i3," eigen values of ",a9,":")')i,str
            write(io,'(14f10.4)')w(1:dim2)
        endif
        ibase=ibase+dim2**2
    enddo
    nullify(p_tmp)
    return

    end subroutine chk_eigens_zmatrix_list


    subroutine dout_matrix_list(str,am,n112,n2,dim_list,io,mode)
    integer,intent(in) :: n112,n2,io,dim_list(n2)
    character(*),intent(in) :: str
    real(q),target,intent(in) :: am(n112)
    integer,intent(in) :: mode

    integer i,dim2,ip,ibase
    real(q) res
    real(q),pointer::p_am(:,:)

    write(io,'("************",2x,a12,2x,"************")')str
    ibase=0
    do i=1,n2
        dim2=dim_list(i)
        write(io,'(" imp=",i3)')i
        p_am(1:dim2,1:dim2)=>am(ibase+1:ibase+dim2*dim2)
        call dwrt_ann(p_am,dim2,io,0,'f9.4')
        if(mode==1)then
            res=0
            do ip=1,dim2
                res=res+p_am(ip,ip)
            enddo
            write(io,'("   sub_tot=",2f10.6)')res
        elseif(mode==-1)then
            write(io,'("   orbital splitting:")')
            write(io,'(14f10.5)')(real(p_am(ip,ip)-p_am(1,1),q), &
                    &ip=1,dim2)
        endif
        ibase=ibase+dim2*dim2
    enddo
    nullify(p_am)
    return

    end subroutine dout_matrix_list


    subroutine zout_matrix_list(str,am,n112,n2,dim_list,io,mode)
    integer,intent(in) :: n112,n2,io,dim_list(n2)
    character(*),intent(in) :: str
    complex(q),target,intent(in) :: am(n112)
    integer,intent(in) :: mode

    integer i,dim2,ip,ibase
    complex(q) res
    complex(q),pointer::p_am(:,:)

    write(io,'("************",2x,a12,2x,"************")')str
    ibase=0
    do i=1,n2
        dim2=dim_list(i)
        write(io,'(" imp=",i3)')i
        p_am(1:dim2,1:dim2)=>am(ibase+1:ibase+dim2*dim2)
        call zwrt_ann(p_am,dim2,io,0,'f9.4')
        if(mode==1)then
            res=0
            do ip=1,dim2
                res=res+p_am(ip,ip)
            enddo
            write(io,'("   sub_tot=",2f10.6)')res
        elseif(mode==-1)then
            write(io,'("   orbital splitting:")')
            write(io,'(14f10.5)')(real(p_am(ip,ip)-p_am(1,1),q), &
                    &ip=1,dim2)
        endif
        ibase=ibase+dim2**2
    enddo
    nullify(p_am)
    return

    end subroutine zout_matrix_list


    subroutine dout_matrix_block(str,am,n12,n2,dim_list,io,mode)
    integer,intent(in)::n12,n2,io,dim_list(n2)
    character(*),intent(in)::str
    real(q),intent(in)::am(n12,n12)
    integer,intent(in)::mode

    integer i,nbase,ibase,dimso,n112
    real(q),allocatable:: am1(:)

    n112=sum(dim_list**2)
    allocate(am1(n112))
    nbase=0
    ibase=0
    do i=1,n2
        dimso=dim_list(i)
        am1(ibase+1:ibase+dimso*dimso)=reshape(am(nbase+1:nbase+dimso, &
                &nbase+1:nbase+dimso),(/dimso*dimso/))
        nbase=nbase+dimso
        ibase=ibase+dimso*dimso
    enddo
    call dout_matrix_list(str,am1,n112,n2,dim_list,io,mode)
    deallocate(am1)
    return

    end subroutine dout_matrix_block


    subroutine zout_matrix_block(str,am,n12,n2,dim_list,io,mode)
    integer,intent(in)::n12,n2,io,dim_list(n2)
    character(*),intent(in)::str
    complex(q),intent(in):: am(n12,n12)
    integer,intent(in)::mode

    integer i,nbase,ibase,dimso,n112
    complex(q),allocatable:: am1(:)

    n112=sum(dim_list**2)
    allocate(am1(n112))
    nbase=0
    ibase=0
    do i=1,n2
        dimso=dim_list(i)
        am1(ibase+1:ibase+dimso*dimso)=reshape(am(nbase+1:nbase+dimso, &
                &nbase+1:nbase+dimso),(/dimso*dimso/))
        nbase=nbase+dimso
        ibase=ibase+dimso*dimso
    enddo
    call zout_matrix_list(str,am1,n112,n2,dim_list,io,mode)
    deallocate(am1)
    return

    end subroutine zout_matrix_block

      
    !*************************************************************************
    ! calculate (a+x*i)^-1 given a^-1; see
    ! http://www.jstor.org/stable/2690437?seq=3
    ! 'on the inverse of the sum of matrices'
    ! much slower than the standard inversion way!!!
    !*************************************************************************
    subroutine zinv_aplx(a,b,n,x)
    integer n
    complex(q) a(n,n),b(n,n)
    real(q) x
      
    integer k
    complex(q) gk
    complex(q) bx(n),v(n)
      
    b=a
    do k=1,n
        gk=-1._q/(1._q+b(k,k)*x)
        bx=b(:,k)*x; v=b(k,:)
        call zgeru(n,n,gk,bx,1,v,1,b,n)
    enddo
    return
      
    end subroutine zinv_aplx


    ! mode>0: add zeros for imaginary part.
    subroutine dwrt_ann(a,n,io,mode,dfmt)
    integer,intent(in)::n,io,mode
    real(q),intent(in)::a(n,n)
    character(*),intent(in)::dfmt

    integer i
    character(20) fmt
      
    write(fmt,'("(",i0,a,")")')n,trim(dfmt)
    write(io,'(" real part")')
    write(io,fmt)a
    if(mode>0)then
        write(io,'(" imag part")')
        write(io,fmt)(0._q, i=1,n*n)
    endif
    return
      
    end subroutine dwrt_ann


    subroutine dwrt_an(a,n,io,mode,dfmt)
    integer,intent(in)::n,io,mode
    real(q),intent(in)::a(n*n)
    character(*),intent(in)::dfmt

    integer i
    character(20) fmt
      
    write(fmt,'("(",i0,a,")")')n,trim(dfmt)
    write(io,'(" real part")')
    write(io,fmt)a
    if(mode>0)then
        write(io,'(" imag part")')
        write(io,fmt)(0._q, i=1,n*n)
    endif
    return
      
    end subroutine dwrt_an


    subroutine zwrt_ann(a,n,io,mode,dfmt)
    integer,intent(in)::n,io,mode
    complex(q),intent(in)::a(n,n)
    character(*),intent(in)::dfmt

    character(20) fmt
      
    write(fmt,'("(",i0,a,")")')n,dfmt
    write(io,'(" real part")')
    write(io,fmt)real(a,q)
    write(io,'(" imag part")')
    write(io,fmt)aimag(a)
    return
      
    end subroutine zwrt_ann
     

    subroutine zwrt_an(a,n,io,mode,dfmt)
    integer,intent(in)::n,io,mode
    complex(q),intent(in)::a(n*n)
    character(*),intent(in)::dfmt

    character(20) fmt
      
    write(fmt,'("(",i0,a,")")')n,dfmt
    write(io,'(" real part")')
    write(io,fmt)real(a,q)
    write(io,'(" imag part")')
    write(io,fmt)aimag(a)
    return
      
    end subroutine zwrt_an


    subroutine zhpev_(jobz,uplo,a,w,z,n)
    integer,intent(in)::n
    character,intent(in)::jobz,uplo
    complex(q),intent(inout)::a(n*(n+1)/2)
    complex(q),intent(out)::z(n,n)
    real(q),intent(out)::w(n)

    integer info
    complex(q) work(2*n)
    real(q) work2(3*n)
    
    call zhpev(jobz,uplo,n,a,w,z,n,work,work2,info)
    if(info.ne.0)then
        write(0,'(" error in zhpev_: info=",i5," n=",i5)')info,n; stop
    endif
    return

    end subroutine zhpev_


    subroutine zheev_blk(uplo,a,w,n)
    integer,intent(in)::n
    character,intent(in)::uplo
    complex(q),intent(inout)::a(n,n)
    real(q),intent(out)::w(n)

    integer i,ihead,nblk


    ihead=1
    ln: do i=1,n
        if(i<n)then
            if(maxval(abs(a(i+1:,ihead:i)))>1.d-7)cycle ln
            if(maxval(abs(a(ihead:i,i+1:)))>1.d-7)cycle ln
        endif
        nblk=i-ihead+1
        if(nblk==1)then
            w(i)=real(a(i,i),q)
            a(i,i)=1._q
        else
            call zheev_('v',uplo,a(ihead:i,ihead:i),w(ihead:i),nblk)
        endif
        ihead=i+1
    enddo ln
    return

    end subroutine zheev_blk

      
    subroutine zheev_(jobz,uplo,a,w,n)
    integer,intent(in)::n
    character,intent(in)::jobz,uplo
    complex(q),intent(inout)::a(n,n)
    real(q),intent(out)::w(n)
      
    integer info,lwork
    real(8),allocatable :: rwork(:)
    complex(8),allocatable :: work(:)
      
    lwork=32*n
    allocate(rwork(max(1,3*n-2)),work(lwork))
    call zheev(jobz,uplo,n,a,n,w,work,lwork,rwork,info)
    if(info.ne.0)then
        write(0,'(" error in zheev_: info=",i5," n=",i5)')info,n; stop
    endif
    deallocate(rwork,work)
    return
      
    end subroutine zheev_
      
      
    subroutine zheevd_(jobz,uplo,a,w,n)
    integer n
    character jobz,uplo
    complex(q) a(n,n)
    real(q)    w(n)
      
    integer info,lwork,lrwork,liwork
    integer,allocatable :: iwork(:)
    real(8),allocatable :: rwork(:)
    complex(8),allocatable :: work(:)
      
    lwork=8*n+n**2; lrwork=1+25*n+2*n**2; liwork=3+25*n
    allocate(rwork(lrwork),work(lwork),iwork(liwork))
    call zheevd(jobz,uplo,n,a,n,w,work,lwork,rwork,lrwork,iwork,liwork,info)
    if(info.ne.0)then
        write(0,'(" error in zheevd_: info=",i5," n=",i5)')info,n; stop
    endif
    deallocate(iwork,rwork,work)
    return
      
    end subroutine zheevd_
      
      
    subroutine zheevx_i(jobz,uplo,a,w,n,il,iu)
    integer,intent(in) :: n, il, iu
    character,intent(in) :: jobz, uplo
    complex(q),intent(inout) :: a(n,n)
    real(q),intent(out) :: w(n)
      
    integer m0
    real(q) vl, vu, abstol
    complex(q),allocatable :: z(:,:)
      
    m0 = iu - il + 1
    abstol = 0._q
    allocate(z(n,m0))
    call zheevx_(jobz,'i',uplo,n,a,vl,vu,il,iu,abstol,m0,w,z)
    a(:, 1:m0) = z
    deallocate(z)
    return
      
    end subroutine zheevx_i
      
      
    subroutine zheevx_(jobz,range,uplo,n,a,vl,vu,il,iu,abstol,m0,w,z)
    integer n,m0,il,iu
    character jobz,range,uplo
    real(q) vl,vu,abstol,w(n)
    complex(q) a(n,n),z(n,m0)
      
    integer lwork,info,m
    integer,allocatable::iwork(:),ifail(:)
    real(q),allocatable::rwork(:)
    complex(q),allocatable::work(:)
      
    lwork=77*n
    allocate(work(lwork),rwork(7*n),iwork(5*n),ifail(n))
    work=0; rwork=0; iwork=0; ifail=0
    call zheevx(jobz,range,uplo,n,a,n,vl,vu,il,iu,abstol,m,w,z,n, &
            &work,lwork,rwork,iwork,ifail,info)
    if(info/=0)then
        write(0,'(" error in zheevx_: info=",i3," ifail=")')info
        write(0,'(10i4)')ifail(1:m); stop
    endif
    deallocate(work,rwork,iwork,ifail)
    return
      
    end subroutine zheevx_
      
      
    subroutine dsyevx_i(jobz,uplo,a,w,n,il,iu)
    integer,intent(in) :: n, il, iu
    character,intent(in) :: jobz, uplo
    real(q),intent(inout) :: a(n,n)
    real(q),intent(out) :: w(n)
      
    integer m0
    real(q) vl, vu, abstol
    real(q),allocatable :: z(:,:)
      
    m0 = iu - il + 1
    abstol = 0._q
    allocate(z(n,m0))
    call dsyevx_(jobz,'i',uplo,n,a,vl,vu,il,iu,abstol,m0,w,z)
    a(:, 1:m0) = z
    deallocate(z)
    return
      
    end subroutine dsyevx_i
      
      
    subroutine dsyevx_(jobz,range,uplo,n,a,vl,vu,il,iu,abstol,m0,w,z)
    integer n,m0,il,iu
    character jobz,range,uplo
    real(q) vl,vu,abstol,w(n)
    real(q) a(n,n),z(n,m0)
      
    integer lwork,info,m
    integer,allocatable::iwork(:),ifail(:)
    real(q),allocatable::work(:)
      
    lwork=77*n
    allocate(work(lwork),iwork(5*n),ifail(n)); work=0; iwork=0; ifail=0
    call dsyevx(jobz,range,uplo,n,a,n,vl,vu,il,iu,abstol,m,w,z,n,work, &
            &lwork,iwork,ifail,info)
    if(info/=0)then
        write(0,'(" error in dstyevx_: info=",i3," ifail=")')info
        write(0,'(10i4)')ifail(1:m); stop
    endif
    deallocate(work,iwork,ifail)
    return
      
    end subroutine dsyevx_
     

    subroutine dprimme_diag(v,n,w,av)
    integer,intent(in) :: n
    real(q),intent(inout) :: v(n)
    real(q),intent(inout) :: w(2)
    external :: av

    integer nest
    integer,parameter::basismax=12,numemax=1,maxmatvecs=3000,printlevel=0
    real(q),parameter::etol=1.e-8_q

    if(abs(w(1)-1._q)<1.e-12_q)then
      nest=1
    else
      nest=0
    endif

    call primme(basismax,numemax,maxmatvecs,etol,printlevel,0,0,n, &
        &nest,w,v,av)
    return

    end subroutine dprimme_diag


    subroutine zprimme_diag(v,n,w,av)
    integer,intent(in) :: n
    complex(q),intent(inout) :: v(n)
    real(q),intent(inout) :: w(2)
    external :: av

    integer nest
    integer,parameter::basismax=12,numemax=1,maxmatvecs=3000,printlevel=0
    real(q),parameter::etol=1.e-8_q

    if(abs(w(1)-1._q)<1.e-12_q)then
      nest=1
    else
      nest=0
    endif

    call primme(basismax,numemax,maxmatvecs,etol,printlevel,0,0,n, &
        &nest,w,v,av)
    return

    end subroutine zprimme_diag


    subroutine dprimme_diag_chk(n,numemax,av,v)
    integer,intent(in) :: n,numemax
    real(q),optional,intent(out)::v(n*numemax)
    external :: av

    integer nest
    integer,parameter::basismax=12,maxmatvecs=3000,printlevel=0
    real(q),parameter::etol=1.e-8_q
    real(q) w(numemax)
    real(q),allocatable::v_(:)

    nest=0
    if(present(v))then
        call primme(basismax,numemax,maxmatvecs,etol,printlevel,0,0,n, &
            &nest,w,v,av)
    else
        allocate(v_(n*numemax))
        call primme(basismax,numemax,maxmatvecs,etol,printlevel,0,0,n, &
            &nest,w,v_,av)
        deallocate(v_)
    endif
    write(0,'(" eigen-values:")')
    write(0,*)w
    return

    end subroutine dprimme_diag_chk


    subroutine zprimme_diag_chk(n,numemax,av,v)
    integer,intent(in) :: n,numemax
    complex(q),optional,intent(out)::v(n*numemax)
    external :: av

    integer nest
    integer,parameter::basismax=12,maxmatvecs=3000,printlevel=0
    real(q),parameter::etol=1.e-8_q
    real(q) w(numemax)
    complex(q),allocatable::v_(:)

    nest=0
    if(present(v))then
        call primme(basismax,numemax,maxmatvecs,etol,printlevel,0,0,n, &
            &nest,w,v,av)
    else
        allocate(v_(n*numemax))
        call primme(basismax,numemax,maxmatvecs,etol,printlevel,0,0,n, &
            &nest,w,v_,av)
        deallocate(v_)
    endif
    write(0,'(" eigen-values:")')
    write(0,*)w
    return

    end subroutine zprimme_diag_chk


    !*************************************************************************
    ! bug in zhevx? for f orbital n_f=8, it gives nonorthogonal eigen vectors
    !*************************************************************************
    subroutine ddiag_jzjp(uplo,rj2,jz,jp,val1,ndim)
    character uplo
    integer ndim
    real(q) jz(ndim,ndim),jp(ndim,ndim)
    real(q) rj2,val1(ndim)
      
    integer nv,mj,iv,iv0,iv1,multj
    real(q) rj,rmj
    real(q) res
      
    rj=sqrt(rj2+.25_q)-0.5_q; multj=nint(2*rj+1)
    nv=ndim/multj
    call dsyev_('v',uplo,jz,val1,ndim)
    rmj=-rj
    do mj=1,multj-1
        res=1/sqrt(rj*(rj+1)-rmj*(rmj+1))
        do iv=1,nv
            iv0=iv+nv*(mj-1); iv1=iv+nv*mj
            call dgemm('n','n',ndim,1,ndim,res,jp,ndim,jz(:,iv0),ndim, &
                    &d0,jz(:,iv1),ndim) ! jp |n,j,mj>
            val1(iv1)=val1(iv0)+1
        enddo
        rmj=rmj+1
    enddo
    return
      
    end subroutine ddiag_jzjp
      
      
    !*************************************************************************
    ! bug in zhevx? for f orbital n_f=8, it gives nonorthogonal eigen vectors
    !*************************************************************************
    subroutine zdiag_jzjp(uplo,rj2,jz,jp,val1,ndim)
    character,intent(in)::uplo
    integer,intent(in)::ndim
    complex(q),intent(inout)::jz(ndim,ndim)
    complex(q),intent(in)::jp(ndim,ndim)
    real(q),intent(in)::rj2
    real(q),intent(inout)::val1(ndim)
      
    integer nv,mj,iv,iv0,iv1,multj
    real(q) rj,rmj
    complex(q) zes
      
    rj=sqrt(rj2+.25_q)-0.5_q; multj=nint(2*rj+1)
    nv=ndim/multj
    call zheev_('v',uplo,jz,val1,ndim)
    rmj=-rj
    do mj=1,multj-1
        zes=1/sqrt(rj*(rj+1)-rmj*(rmj+1))
        do iv=1,nv
            iv0=iv+nv*(mj-1); iv1=iv+nv*mj
            call zgemm('n','n',ndim,1,ndim,zes,jp,ndim,jz(:,iv0),ndim,z0, &
                    &jz(:,iv1),ndim) ! jp |n,j,mj>
            val1(iv1)=val1(iv0)+1
        enddo
        rmj=rmj+1
    enddo
    return
      
    end subroutine zdiag_jzjp
      
      
    !*************************************************************************
    ! not as robust as dsyev, someitmes may fail
    !*************************************************************************
    subroutine dsyevd_(jobz,uplo,a,w,n)
    integer n
    character jobz,uplo
    real(q) a(n,n),w(n)
      
    integer lwork,liwork,info
    integer,allocatable::iwork(:)
    real(q),allocatable::work(:)
      
    lwork=1; liwork=1
    allocate(work(lwork),iwork(liwork))
    call dsyevd(jobz,uplo,n,a,n,w,work,-1,iwork,-1,info )
    lwork=nint(work(1)); liwork=iwork(1)
    deallocate(work,iwork); allocate(work(lwork),iwork(liwork))
    call dsyevd(jobz,uplo,n,a,n,w,work,lwork,iwork,liwork,info )
    deallocate(work,iwork)
    if(info.ne.0)then
        write(0,'(" error in dsyevd_: info=",i5," n=",i5," try dsyev!")')info,n
        call dsyev_(jobz,uplo,a,w,n)
    endif
    return
      
    end subroutine dsyevd_
      
      
    subroutine dsyev_(jobz,uplo,a,w,n)
    integer n
    character jobz,uplo
    real(q) a(n,n),w(n)
      
    integer lwork,info
    real(q),allocatable::work(:)
      
    lwork=1; allocate(work(lwork))
    call dsyev(jobz,uplo,n,a,n,w,work,-1,info )
    lwork=nint(work(1))
    deallocate(work); allocate(work(lwork))
    call dsyev(jobz,uplo,n,a,n,w,work,lwork,info )
    deallocate(work)
    if(info.ne.0)then
        write(0,'(" error in dsyev_: info=",i5," n=",i5)')info,n; stop
    endif
    return
      
    end subroutine dsyev_
      
      
    subroutine zgeev_(jobvr,a,w,n)
    character jobvr
    integer n
    complex(q) a(n,n),w(n)
      
    integer lwork,info
    complex(q) vl(1,1)
    real(q) rwork(2*n)
    complex(q),allocatable::vr(:,:),work(:)
      
    lwork=16*n
    allocate(work(lwork),vr(n,n))
    call zgeev('n',jobvr,n,a,n,w,vl,1,vr,n,work,lwork,rwork,info)
    if(info.ne.0)then
        write(0,'(" error in zgeev_: info=",i5," n=",i5)')info,n; stop
    endif
    a=vr
    return
      
    end subroutine zgeev_
      
      
    !*************************************************************************
    ! inverse of a square matrix
    !*************************************************************************
    subroutine dinv_(a,ndim)
    real(q),intent(inout) :: a(ndim,ndim)
    integer, intent(in)       :: ndim
      
    integer    :: info, lwork, lda
    integer    :: ipiv(ndim)
    real(q):: work(ndim*64)
    lwork = ndim*64
    lda = ndim
      
    call dgetrf( ndim, ndim, a, lda, ipiv, info )
    if (info.ne.0) then
        write(0,*)'dgetrf info=', info; stop
    endif
    call dgetri( ndim, a, lda, ipiv, work, lwork, info )
    if (info.ne.0) then
        write(0,*)'zgetri info=', info; stop
    endif
    return
      
    end subroutine dinv_
      
      
    subroutine zinv_(a,ndim)
    complex(q),intent(inout) :: a(ndim,ndim)
    integer, intent(in)       :: ndim
      
    integer    :: info, lwork, lda
    integer    :: ipiv(ndim)
    complex(q):: work(ndim*64)
    lwork = ndim*64
    lda = ndim
      
    call zgetrf( ndim, ndim, a, lda, ipiv, info )
    if (info.ne.0) then
        write(0,*)'zgetrf info=', info; stop
    endif
    call zgetri( ndim, a, lda, ipiv, work, lwork, info )
    if (info.ne.0) then
        write(0,*)'zgetri info=', info; stop
    endif
    return
      
    end subroutine zinv_
      
      
    subroutine chk_unitary(a,n,maxerr)
    integer n
    complex(q) a(n,n)
    real(q) maxerr
      
    integer i
    complex(q) b(n,n)
      
    b=matmul(a,transpose(conjg(a)))
    do i=1,n
        b(i,i)=b(i,i)-1._q
    enddo
    maxerr=maxval(abs(b))
    return
      
    end subroutine chk_unitary


    !*********************************************************************
    ! h2e &= eri_{pq,rs} p^+ q r^+ s \\
    !     &= (pq|rs) p^+ r^+ s q - (pq|rs) \delta_{qr} p^+ s
    !*********************************************************************
    subroutine dget_eri(v2e,h1e,eri,n)
    integer,intent(in) :: n
    real(q),intent(in) :: v2e(n,n,n,n),h1e(n,n)
    real(q),intent(out) :: eri(n,n,n,n)

    integer i
    complex(q) f1e(n,n)

    f1e=h1e
    do i=1,n
      f1e=f1e-v2e(:,i,i,:)/2.d0
    enddo
    ! Half-filling by definition
    f1e=f1e/dble(n)
    eri=v2e/2.d0
    do i=1,n
        eri(:,:,i,i)=eri(:,:,i,i)+f1e
    enddo
    eri=reshape(transpose(reshape(eri,(/n*n,n*n/))),(/n,n,n,n/))
    do i=1,n
        eri(:,:,i,i)=eri(:,:,i,i)+f1e
    enddo
    return

    end subroutine dget_eri


    subroutine zget_eri(v2e,h1e,eri,n)
    integer,intent(in) :: n
    complex(q),intent(in) :: v2e(n,n,n,n),h1e(n,n)
    complex(q),intent(out) :: eri(n,n,n,n)

    integer i
    complex(q) f1e(n,n)

    f1e=h1e
    do i=1,n
      f1e=f1e-v2e(:,i,i,:)/2.d0
    enddo
    ! Half-filling by definition
    f1e=f1e/dble(n)
    eri=v2e/2.d0
    do i=1,n
        eri(:,:,i,i)=eri(:,:,i,i)+f1e
    enddo
    eri=reshape(transpose(reshape(eri,(/n*n,n*n/))),(/n,n,n,n/))
    do i=1,n
        eri(:,:,i,i)=eri(:,:,i,i)+f1e
    enddo
    return

    end subroutine zget_eri

      
    !*************************************************************************
    ! get hermitian matrix expansion coefficients (real)
    ! mode>0 : forward, get c
    !     <0 : backward, get a
    !     =0 : forward first, then backward (symmetrization)
    !*************************************************************************
    subroutine get_dhm_expand(a,b,n,nb,c,mode,ltrans,lherm)
    integer,intent(in)::n,nb,mode
    real(q),intent(in)::b(n,n,nb)
    real(q),intent(inout)::a(n,n),c(nb)
    logical,intent(in)::ltrans,lherm
      
    integer i
    real(q) b_(n,n)
    real(q),external::ddot
    
    if(mode>=0)then
        do i=1,nb
            if(ltrans)then
                b_=transpose(b(:,:,i))
            else
                b_=b(:,:,i)
            endif
            c(i)=ddot(n*n,b_(1,1),1,a(1,1),1)
        enddo
    endif
    if(mode<=0)then
        a=0
        do i=1,nb
            if(ltrans)then
                b_=transpose(b(:,:,i))
            else 
                b_=b(:,:,i)
            endif
            a=a+b_*c(i)
        enddo
    endif
    return
      
    end subroutine get_dhm_expand
      
      
    subroutine get_zhm_expand(a,b,n,nb,c,mode,ltrans,lherm)
    integer,intent(in)::n,nb,mode
    complex(q),intent(in)::b(n,n,nb)
    complex(q),intent(inout)::a(n,n),c(nb)
    logical,intent(in)::ltrans,lherm
      
    integer i
    complex(q) b_(n,n)
    complex(q),external::zdotc
      
    if(mode>=0)then
        do i=1,nb
            if(ltrans)then
                b_=transpose(b(:,:,i))
            else
                b_=b(:,:,i)
            endif
            c(i)=zdotc(n*n,b_(1,1),1,a(1,1),1)
            if(lherm)c(i)=real(c(i),q)
        enddo
    endif
    if(mode<=0)then
        a=0
        do i=1,nb
            if(ltrans)then
                b_=transpose(b(:,:,i))
            else
                b_=b(:,:,i)
            endif
            a=a+b_*c(i)
        enddo
    endif
    return
      
    end subroutine get_zhm_expand
 

    subroutine get_dhm_expand_sym(a,b,n,nb,ltrans,lherm)
    integer,intent(in)::n,nb
    real(q),intent(in)::b(n,n,nb)
    real(q),intent(inout)::a(n,n)
    logical,intent(in)::ltrans,lherm
      
    integer i
    real(q) b_(n,n),a_(n,n),c
    real(q),external::ddot
      
    a_=a; a=0
    do i=1,nb
        if(ltrans)then
            b_=transpose(b(:,:,i))
        else
            b_=b(:,:,i)
        endif
        c=ddot(n*n,b_(1,1),1,a_(1,1),1)
        a=a+b_*c
    enddo
    return
      
    end subroutine get_dhm_expand_sym
      
    
    ! lherm: hermitian matrix or not.
    subroutine get_zhm_expand_sym(a,b,n,nb,ltrans,lherm)
    integer,intent(in)::n,nb
    complex(q),intent(in)::b(n,n,nb)
    complex(q),intent(inout)::a(n,n)
    logical,intent(in)::ltrans,lherm
      
    integer i
    complex(q) b_(n,n),a_(n,n),c
    complex(q),external::zdotc
      
    a_=a; a=0
    do i=1,nb
        if(ltrans)then
            b_=transpose(b(:,:,i))
        else
            b_=b(:,:,i)
        endif
        c=zdotc(n*n,b_(1,1),1,a_(1,1),1)
        if(lherm)c=real(c,q)
        a=a+b_*c
    enddo
    return
      
    end subroutine get_zhm_expand_sym
      
      
    subroutine zuhau(a,u,n,m,uhau,trul,trur)
    integer,intent(in)::n,m
    complex(q),intent(in)::u(n,m)
    complex(q),intent(inout)::a(n,n)
    complex(q),intent(inout),optional::uhau(m,m)
    character*1,intent(in),optional::trul,trur
      
    character*1 tul,tur
    complex(q),allocatable::a_(:,:),b_(:,:)
      
    if(present(trul))then; tul=trul; else; tul='c'; endif
    if(present(trur))then; tur=trur; else; tur='n'; endif
    allocate(a_(n,m)); a_=0
    call zgemm('n',tur,n,m,n,z1,a,n,u ,n,z0,a_,n)
    allocate(b_(m,m)); b_=0
    call zgemm(tul,'n',m,m,n,z1,u,n,a_,n,z0,b_,m)
    if(present(uhau))then
        uhau=uhau+b_
    else
        if(n/=m)stop ' error in zuhau: n/=m while not preset uhau!'
        a=b_
    endif
    deallocate(a_,b_)
    return
      
    end subroutine zuhau
      
      
    subroutine duhau(a,u,n,m,uhau,trul,trur)
    integer,intent(in)::n,m
    real(q),intent(in)::u(n,m)
    real(q),intent(inout)::a(n,n)
    real(q),intent(inout),optional::uhau(m,m)
    character*1,intent(in),optional::trul,trur
      
    character*1 tul,tur
    real(q),allocatable::a_(:,:),b_(:,:)
      
    if(present(trul))then; tul=trul; else; tul='c'; endif
    if(present(trur))then; tur=trur; else; tur='n'; endif
    allocate(a_(n,m)); a_=0
    call dgemm('n',tur,n,m,n,d1,a,n,u ,n,d0,a_,n)
    allocate(b_(m,m)); b_=0
    call dgemm(tul,'n',m,m,n,d1,u,n,a_,n,d0,b_,m)
    if(present(uhau))then
        uhau=uhau+b_
    else
        if(n/=m)stop ' error in duhau: n/=m while not preset uhau!'
        a=b_
    endif
    deallocate(a_,b_)
    return
      
    end subroutine duhau
      
      
    subroutine zanmxbmm(transb,a,b,n,m)
    character*1 transb
    integer n,m
    complex(q) a(n,m),b(m,m)
      
    complex(q),allocatable::a_(:,:)
      
    allocate(a_(n,m)); a_=0
    call zgemm('n',transb,n,m,m,z1,a,n,b,m,z0,a_,n)
    a=a_
    deallocate(a_)
    return
      
    end subroutine zanmxbmm
      
      
    subroutine danmxbmm(transb,a,b,n,m)
    character*1 transb
    integer n,m
    real(q) a(n,m),b(m,m)
      
    real(q),allocatable::a_(:,:)
      
    allocate(a_(n,m)); a_=0
    call dgemm('n',transb,n,m,m,d1,a,n,b,m,d0,a_,n)
    a=a_
    deallocate(a_)
    return
      
    end subroutine danmxbmm
      
      
    subroutine dannxbnm(transa,a,b,n,m,ab)
    character*1,intent(in)::transa
    integer,intent(in)::n,m
    real(q),intent(in)::a(n,n)
    real(q),intent(inout)::b(n,m)
    real(q),optional,intent(out)::ab(:,:)
      
    real(q),allocatable::b_(:,:)
      
    allocate(b_(n,m)); b_=0
    call dgemm(transa,'n',n,m,n,z1,a,n,b,n,z0,b_,n)
    if(present(ab))then
        ab = b_
    else
        b = b_
    endif
    deallocate(b_)
    return
      
    end subroutine dannxbnm


    subroutine dannxbnmc(transa,a,b,c,n,m)
    character*1,intent(in)::transa
    integer,intent(in)::n,m
    real(q),intent(in)::a(n,n),b(n,m)
    real(q),intent(out)::c(n,m)
      
    call dgemm(transa,'n',n,m,n,z1,a,n,b,n,z0,c,n)
    return
      
    end subroutine dannxbnmc

      
    subroutine zannxbnm(transa,a,b,n,m,ab)
    character*1,intent(in)::transa
    integer,intent(in)::n,m
    complex(q),intent(in)::a(n,n)
    complex(q),intent(inout)::b(n,m)
    complex(q),optional,intent(out)::ab(:,:)
      
    complex(q),allocatable::b_(:,:)
      
    allocate(b_(n,m)); b_=0
    call zgemm(transa,'n',n,m,n,z1,a,n,b,n,z0,b_,n)
    if(present(ab))then
        ab = b_
    else
        b=b_
    endif
    deallocate(b_)
    return
      
    end subroutine zannxbnm


    subroutine zannxbnmc(transa,a,b,c,n,m)
    character*1,intent(in)::transa
    integer,intent(in)::n,m
    complex(q),intent(in)::a(n,n),b(n,m)
    complex(q),intent(out)::c(n,m)
      
    call zgemm(transa,'n',n,m,n,z1,a,n,b,n,z0,c,n)
    return
      
    end subroutine zannxbnmc


    subroutine zannxbnn(transa,transb,a,b,n,mode)
    character*1 transa,transb
    integer n,mode
    complex(q) a(n,n),b(n,n)
      
    complex(q),allocatable::a_(:,:)
      
    allocate(a_(n,n))
    call zgemm(transa,transb,n,n,n,z1,a,n,b,n,z0,a_,n)
    if(mode==1)then
        a=a_
    else
        b=a_
    endif
    deallocate(a_)
    return
      
    end subroutine zannxbnn
      
      
    subroutine dannxbnn(transa,transb,a,b,n,mode)
    character*1 transa,transb
    integer n,mode
    real(q) a(n,n),b(n,n)
      
    real(q),allocatable::a_(:,:)
      
    allocate(a_(n,n))
    call dgemm(transa,transb,n,n,n,z1,a,n,b,n,z0,a_,n)
    if(mode==1)then
        a=a_
    else
        b=a_
    endif
    deallocate(a_)
    return
      
    end subroutine dannxbnn
      

    function trace_dann(a,n) result(tr)
    integer,intent(in)::n
    real(q),intent(in)::a(n,n)
    real(q) tr
      
    integer i
      
    tr=0
    do i=1,n
        tr=tr+a(i,i)
    enddo
    return
      
    end function trace_dann

      
    function trace_zann(a,n) result(tr)
    integer,intent(in)::n
    complex(q),intent(in)::a(n,n)
    complex(q) tr
      
    integer i
      
    tr=0
    do i=1,n
        tr=tr+a(i,i)
    enddo
    return
      
    end function trace_zann
      
      
    subroutine zatofa(a,fa,n,mode,coef,lsym)
    integer,intent(in) :: n,mode
    complex(q),intent(in) :: a(:,:)
    complex(q),intent(out) :: fa(n,n)
    real(q) coef
    logical :: lsym
      
    real(q),allocatable::wr(:)
    complex(q),allocatable::v(:,:),w(:),vinv(:,:)
      
    allocate(w(n),v(n,n)); w=0
    v=a
    if(lsym)then
        allocate(wr(n))
        call zheev_('v','l',v,wr,n)
        w=wr; deallocate(wr)
        call zatofa1(fa,w,v,n,n,mode,coef)
    else
        call zgeev_('v',v,w,n)
        allocate(vinv(n,n)); vinv=v
        call zinv_(vinv,n)
        call zatofa1(fa,w,v,n,n,mode,coef,vinv=vinv)
    endif
    deallocate(w,v)
    return
      
    end subroutine zatofa
      
      
    subroutine datofa(a,fa,n,mode,coef,lsym)
    integer,intent(in) :: n,mode
    real(q),intent(in) :: a(:,:)
    real(q),intent(out) :: fa(n,n)
    real(q) coef
    logical :: lsym
      
    real(q),allocatable::w(:)
    real(q),allocatable::v(:,:),vinv(:,:)
      
    allocate(w(n),v(n,n)); w=0
    v=a
    call dsyev_('v','l',v,w,n)
    if(lsym)then
        call datofa1(fa,w,v,n,n,mode,coef)
    else
        allocate(vinv(n,n)); vinv=v
        call dinv_(vinv,n)
        call datofa1(fa,w,v,n,n,mode,coef,vinv=vinv)
        deallocate(vinv)
    endif
    deallocate(w,v)
    return
      
    end subroutine datofa
      
      
    subroutine datofa1(fa,w,v,n,m,mode,coef,vinv)
    integer,intent(in) :: n,m,mode
    real(q),intent(out) ::  fa(m,m)
    real(q),intent(in) :: w(n),v(m,n),coef
    real(q),intent(in),optional::vinv(n,m)
      
    integer i
    real(q) w_(n),vw(m,n)
      
    w_=w*coef
    select case(mode)
    case(1) ! exponential
        w_=exp(w_)
    case(2) ! log
        do i=1,n
            if(w_(i)<rlbound)then
                w_(i)=rlbound
            endif
            w_(i)=log(w_(i))
        enddo
    case(-1) ! inverse
        do i=1,n
            if(abs(w_(i))<small)then
                ! linear dependent term, ignored
                w_(i)=0._q
            elseif(abs(w_(i))>rubound)then
                w_(i)=1/w_(i)
            endif
        enddo
    case(-12) ! ^(-1/2)
        do i=1,n
            if(w_(i)<rlbound)then
                w_(i)=rlbound
            elseif(w_(i)>rubound)then
                w_(i)=rubound
            endif
            w_(i)=1/sqrt(w_(i))
        enddo
    case default
        stop ' error in datofa: undefined mode!'
    end select
    do i=1,n
        vw(:,i)=v(:,i)*w_(i)
    enddo
    if(present(vinv))then
        call dgemm('n','n',m,m,n,d1,vw,m,vinv,m,d0,fa,m)
    else
        call dgemm('n','c',m,m,n,d1,vw,m,v   ,m,d0,fa,m)
    endif
    return
      
    end subroutine datofa1
      
      
    subroutine zatofa1(fa,w,v,n,m,mode,coef,vinv)
    integer n,m,mode
    complex(q) fa(m,m),v(m,n),w(n)
    real(q) coef
    complex(q),optional::vinv(n,m)
      
    integer i
    complex(q) w_(n)
    complex(q) vw(m,n)
      
    w_=w*coef
    select case(mode)
    case(1) ! exponential
        w_=exp(w_)
    case(2) ! log
        do i=1,n
            if(abs(w_(i))<rlbound)then
                w_(i)=rlbound
            endif
            w_(i)=log(w_(i))
        enddo
    case(-1) ! inverse
        do i=1,n
            if(abs(w_(i))<small)then
                ! linear dependent term, should ignore
                w_(i)=0._q
            else
                w_(i)=1/w_(i)
            endif
        enddo
    case(-12) ! ^(-1/2)
        do i=1,n
            if(abs(w_(i))<rlbound)then
                w_(i)=rlbound
            elseif(abs(w_(i))>rubound)then
                w_(i)=rubound
            endif
            w_(i)=1/sqrt(w_(i))
        enddo
    case default
        stop ' error in datofa: undefined mode!'
    end select
    do i=1,n
        vw(:,i)=v(:,i)*w_(i)
    enddo
    if(present(vinv))then
        call zgemm('n','n',m,m,n,z1,vw,m,vinv,m,z0,fa,m)
    else
        call zgemm('n','c',m,m,n,z1,vw,m,v   ,m,z0,fa,m)
    endif
    return
      
    end subroutine zatofa1
      
      
    !*************************************************************************
    ! a is symmetric, ( \par f(a) / \par d_n --(hermitian component) )
    !*************************************************************************
    subroutine dpfa_pa(a,pa,h,n,f,fp)
    integer n
    real(q) a(n,n),pa(n,n),h(n,n)
    real(q),external::f,fp
      
    real(q) u(n,n),la(n,n)
    real(q) w(n)
      
    u=a
    call dsyev_('v','l',u,w,n)
    call duhau(h,u,n,n) ! in eigenstate reprsentation of a
    call calc_loewner(w,la,n,f,fp)
    pa=la*h
    u=transpose(u)
    call duhau(pa,u,n,n) ! back to original representation
    return
      
    end subroutine dpfa_pa
      
      
    !*************************************************************************
    ! a is hermitian, ( \par f(a) / \par a_ij )
    !*************************************************************************
    subroutine zpfa_pa(a,pa,h,n,f,fp)
    integer n
    complex(q) a(n,n),pa(n,n),h(n,n)
    real(q),external::f,fp
      
    complex(q) u(n,n)
    real(q) w(n),la(n,n)
      
    u=a
    call zheev_('v','l',u,w,n)
    call zuhau(h,u,n,n) ! in eigenstate reprsentation of a
    call calc_loewner(w,la,n,f,fp)
    pa=la*h
    u=transpose(conjg(u))
    call zuhau(pa,u,n,n) ! back to original representation
    return
      
    end subroutine zpfa_pa
      
      
    subroutine dh_regularize(h,n,tol)
    integer n
    real(q) h(n,n),tol
      
    integer i
    real(q) w(n),vw(n,n)
      
    call dsyev_('v','l',h,w,n)
    do i=1,n
        w(i)=max(w(i),tol)
        vw(:,i)=h(:,i)*w(i)
    enddo
    call dannxbnn('n','t',vw,h,n,2)
    return
      
    end subroutine dh_regularize
      
      
    subroutine zh_regularize(h,n,tol)
    integer n
    complex(q) h(n,n)
    real(q) tol
      
    integer i
    real(q) w(n)
    complex(q) vw(n,n)
      
    call zheev_('v','l',h,w,n)
    do i=1,n
        w(i)=max(w(i),tol)
        vw(:,i)=h(:,i)*w(i)
    enddo
    call zannxbnn('n','c',vw,h,n,2)
    return
      
    end subroutine zh_regularize
      
      
    !*************************************************************************
    ! loewner matrix of hermitian a given its eigen values w.
    !*************************************************************************
    subroutine calc_loewner(w,la,n,f,fp)
    integer n
    real(q) w(n),la(n,n)
    real(q),external::f,fp
      
    integer i,j
      
    do i=1,n
        do j=1,n
            if(abs(w(i)-w(j))<1.e-10_q)then
                la(i,j)=fp(w(i))
            else
                la(i,j)=(f(w(i))-f(w(j)))/(w(i)-w(j))
            endif
        enddo
    enddo
    return
      
    end subroutine calc_loewner
      
      
    function dsimix(x)
    real(q),intent(in) :: x
    real(q) :: dsimix
      
    real(q) x_
      
    if(abs(x)<1.e-8_q)then
        x_=1.e-8_q
    elseif(abs(x-1)<1.e-8_q)then
        x_=1-1.e-8_q
    elseif(x<0.or.x>1)then
        write(0,*) " fetal error in dsimix: illegal x=",x; stop
    else
        x_=x
    endif
    dsimix=sqrt(x_*(1-x_))
    return
      
    end function dsimix
      
      
    function dpsimix(x)
    real(q),intent(in) :: x
    real(q) :: dpsimix
      
    real(q) x_
      
    if(abs(x)<1.e-8_q)then
        x_=1.e-8_q
    elseif(abs(x-1)<1.e-8_q)then
        x_=1-1.e-8_q
    elseif(x<0.or.x>1)then
        write(0,*) " fetal error in dsimix: illegal x=",x; stop
    else
        x_=x
    endif
    dpsimix=(.5_q-x)/sqrt(x_*(1-x_))
    return
      
    end function dpsimix
      
      
    subroutine zabs_order(v,e,n)
    integer n
    real(q) e(n)
    complex(q) v(n,n)
      
    integer i,id(n)
    real(q) ebuf(n)
    complex(q),allocatable::zbuf(:,:)
      
    ebuf=abs(e); id=0
    call dsort(n,ebuf,id,.true.)
    allocate(zbuf(n,n))
    ebuf=e; zbuf=v
    do i=1,n
        e(i)=ebuf(id(i))
        v(:,i)=zbuf(:,id(i))
    enddo
    deallocate(zbuf)
    return
      
    end subroutine zabs_order
      
      
    !*************************************************************************
    ! sorts ra in descending/ascending order, and rearanges an index array ib
    !*************************************************************************
    subroutine dsort(n,ra,ib,ldsc)
    integer n
    real(q) ra(n)
    integer ib(n)
    logical ldsc
      
    real(q) rra,fak
    integer iib,l,ir,j,i
      
    do i=1,n
        ib(i)=i
    enddo
    if (n<=1) return
    if(ldsc)then; fak=1._q; else; fak=-1._q; endif
    l=n/2+1; ir=n
    do
        if(l.gt.1)then
            l=l-1
            rra=ra(l); iib=ib(l)
        else
            rra=ra(ir); iib=ib(ir)
            ra(ir)=ra(1); ib(ir)=ib(1)
            ir=ir-1
            if(ir.eq.1)then
                ra(1)=rra; ib(1)=iib
                return
            endif
        endif
        i=l; j=l+l
        do while(j.le.ir)
            if(j.lt.ir)then
                if((ra(j)-ra(j+1))*fak.gt.0)j=j+1
            endif
            if((rra-ra(j))*fak.gt.0)then
                ra(i)=ra(j); ib(i)=ib(j)
                i=j; j=j+j
            else
                j=ir+1
            endif
        enddo
        ra(i)=rra; ib(i)=iib
    enddo
    return
      
    end subroutine dsort
      
      
    subroutine locate_mi(m,n,m1,lm1)
    integer n,m(n),m1,lm1
      
    integer i
      
    lm1=0
    do i=1,n
        if(m1.eq.m(i))then
            lm1=i; exit
        endif
    enddo
    return
      
    end subroutine locate_mi
      
      
    function l22l(l2)
    real(q) l2,l22l
      
    l22l=sqrt(l2+0.25_q)-0.5_q
    return
      
    end function l22l
      
      
    function file_name(i,fname)
    integer i
    character(*) fname
    character*77 file_name
      
    character*7 str
      
    write(str,'(i7)')i
    file_name=fname//'_'//trim(adjustl(str))//'.h5'
    return
      
    end function file_name
      
      
    function int_to_str(i)
    integer i
    character(len = 77) int_to_str
      
    write(int_to_str,*)i
    int_to_str = adjustl(int_to_str)
    return
      
    end function int_to_str
      
      
    function fermi_fun(x)
    real(q) x,fermi_fun
      
    if(x<-200._q)then
        fermi_fun=1._q
    elseif(x<200._q)then
        fermi_fun=1._q/(exp(x)+1._q)
    else
        fermi_fun=0._q
    endif
    return
      
    end function fermi_fun
      
      
    function gauss_fun(x)
    real(q) x,gauss_fun
      
    if(x<-7._q)then
        gauss_fun=2._q
    elseif(x<0._q)then
        gauss_fun=2-erfc(-x)
    elseif(x<7._q)then
        gauss_fun=erfc(x)
    else
        gauss_fun=0._q
    endif
    gauss_fun=gauss_fun/2
    return
      
    end function gauss_fun
      
      
    function max_interval(iarray,n)
    integer n,max_interval,iarray(n)
      
    integer i
      
    max_interval=0
    do i=2,n
        max_interval=max(iarray(i)-iarray(i-1),max_interval)
    enddo
    return
      
    end function max_interval
      
      
    subroutine set_range(x0,x1,x,n)
    integer n
    real(q) x0,x1,x(n)
      
    integer i
    real(q) dx
      
    dx=(x1-x0)/(n-1)
    do i=1,n
        x(i)=x0+(i-1)*dx
    enddo
    return
      
    end subroutine set_range
      
      
    subroutine set_linear_simp(wt,n,dx)
    integer n
    real(8) wt(n),dx
      
    integer i
      
    wt=0
    if(mod(n,2)==0)then
        stop ' error in set_linear_simp: assume odd points!'
    endif
    do i=n,3,-2
        wt(i)  =dx/3._8+wt(i)
        wt(i-1)=dx/3._8*4
        wt(i-2)=dx/3._8
    enddo
    return
      
    end subroutine set_linear_simp
      
      
    function lkey_exist(line, key)
    character(*) line, key
    logical lkey_exist
      
    lkey_exist=.false.
    if(index(line, key)>0.and.(index(line,'#')<=0.or. &
                &index(line,'#')>index(line, key)))then
        lkey_exist=.true.
    endif
    return
      
    end function lkey_exist
      
      
    subroutine get_imin1(i,n,imin)
    integer n,i(n),imin
      
    integer j
      
    imin=100
    do j=1,n
        if(i(j)<=0)cycle
        imin=min(imin,i(j))
    enddo
    return
      
    end subroutine get_imin1
      

    !*************************************************************************
    ! calc. creation a^\dagger or annihilation a^- |  > on exit: integer
    ! isign: sign for this action (antisymmetry)
    !*************************************************************************
    subroutine act_state(state,pos,lact,isign_)
    integer,intent(inout)::state,isign_
    integer,intent(in)::pos
    logical,intent(in)::lact

    logical lbtest
    integer,parameter,dimension(0:1)::iphase=(/1,-1/)

    lbtest=btest(state,pos)
    if(lact.eqv.lbtest) then
        isign_=0
        return
    endif

    isign_=isign_*iphase(poppar(ishft(state,32-pos)))
    if(lact) then
        state=ibset(state,pos)
    else
        state=ibclr(state,pos)
    endif
    return

    end subroutine act_state


    ! Set the indices for occupied orbitals.
    subroutine set_occupied_orbitals(state,ia,n)
    integer,intent(in)::state,n
    integer,intent(out)::ia(n)

    integer i,isum

    ia=0; isum=0
    do i=1,n
        if(btest(state,i-1))then
            isum=isum+1
            ia(isum)=i
        endif
    enddo
    return

    end subroutine set_occupied_orbitals


    subroutine set_skip_orbitals(u,norb,ia,na,iskip,nskip)
    integer,intent(in)::norb,na,ia(na)
    real(q),intent(in)::u(norb,norb)
    integer,intent(out)::nskip,iskip(norb)

    integer i,j

    nskip=0
    l1: do i=1,norb
        l2: do j=1,na
            if(u(i,ia(j))>1.d-10)cycle l1
        enddo l2
        nskip=nskip+1
        iskip(nskip)=i
    enddo l1
    return

    end subroutine set_skip_orbitals


    subroutine dcalc_fock_state_coef(u,n1,ia,ib,n2,coef)
    integer,intent(in)::n1,n2,ia(n2),ib(n2)
    real(q),intent(in)::u(n1,n1)
    real(q),intent(inout)::coef

    integer i,j
    real(q) u_(n2,n2)

    do i=1,n2; do j=1,n2
        u_(i,j)=u(ib(i),ia(j))
    enddo; enddo
    coef=coef*dget_determinant_a(u_,n2)
    return

    end subroutine dcalc_fock_state_coef


    subroutine zcalc_fock_state_coef(u,n1,ia,ib,n2,coef)
    integer,intent(in)::n1,n2,ia(n2),ib(n2)
    complex(q),intent(in)::u(n1,n1)
    complex(q),intent(inout)::coef

    integer i,j
    complex(q) u_(n2,n2)

    do i=1,n2; do j=1,n2
        u_(i,j)=u(ib(i),ia(j))
    enddo; enddo
    coef=coef*zget_determinant_a(u_,n2)
    return

    end subroutine zcalc_fock_state_coef


    subroutine calc_state_sz(bs,norb,sz,m,list)
    integer,intent(in)::bs,norb,m
    integer,optional,intent(in)::list(m)
    real(q),intent(out)::sz

    integer i
    logical linls

    sz=0
    do i=1,norb,2
        if(present(list))then
            call chk_i_in_list(i,m,list,linls)
            if(linls)cycle
        endif
        if(btest(bs,i-1))then ! spin-up
            sz=sz+0.5_q
        endif
    enddo
    do i=2,norb,2
        if(present(list))then
            call chk_i_in_list(i,m,list,linls)
            if(linls)cycle
        endif
        if(btest(bs,i-1))then ! spin-down
            sz=sz-0.5_q
        endif
    enddo
    return

    end subroutine calc_state_sz


    subroutine calc_state_sz2(bs,norb,sz,sz_orb)
    integer,intent(in)::bs,norb
    real(q),intent(in)::sz_orb(:)
    real(q),intent(out)::sz

    integer i

    sz=0
    do i=1,norb
        if(btest(bs,i-1))then 
            sz=sz+sz_orb(i)
        endif
    enddo
    return

    end subroutine calc_state_sz2


    subroutine sum_empty_slots_fs(bs,m,ilist,sume)
    integer,intent(in)::bs,m,ilist(m)
    integer,intent(out)::sume

    integer i

    sume=0
    do i=1,m
        if(.not.btest(bs,ilist(i)-1))then
            sume=sume+1
        endif
    enddo
    return

    end subroutine sum_empty_slots_fs


    subroutine chk_sum_first_empty_slots_fs(bs,m,ilist,mfirste,lfirste)
    integer,intent(in)::bs,m,ilist(m),mfirste
    logical,intent(out)::lfirste

    integer i

    lfirste=.false.
    ! First mfirste slots should be empty
    do i=1,mfirste 
        if(btest(bs,ilist(i)-1))return
    enddo
    ! The remaining should be occupied.
    do i=1+mfirste,m
        if(.not.btest(bs,ilist(i)-1))return
    enddo
    lfirste=.true.
    return

    end subroutine chk_sum_first_empty_slots_fs


    subroutine chk_i_in_list(l,m,list,linls)
    integer,intent(in)::l,m,list(m)
    logical,intent(out)::linls

    integer i

    linls=.false.
    do i=1,m
        if(l==list(i))then
            linls=.true.
            return
        endif
    enddo
    return

    end subroutine chk_i_in_list


    !> Get the impurity 
    !! physical density matrix (ni_{AB} = <c_A^\dagger c_B>), 
    !! bath density matrix (nb_{ab} = <f_b f_a^\dagger>
    !!                              = \delta_{ab}-<f_a^\dagger f_b>)
    !! and r0_{aA} = <c_A^\dagger f_a>.
    !! @param dm complete density matrix <{c_A f_b}^\dagger {c_B f_b}>
    !< Real version.
    subroutine dget_rn_from_embeddm(dm,ni,nb,r0,norb)
    integer,intent(in)::norb !< dimension of {c_A} or {f_a}
    real(q),intent(in)::dm(norb*2,norb*2)
    real(q),intent(out)::ni(norb,norb),nb(norb,norb),r0(norb,norb)

    integer i

    !< Get physical density matrix.
    ni=dm(:norb,:norb)

    !< Bath auxiliary density matrix
    nb=-dm(norb+1:,norb+1:)
    forall(i=1:norb) nb(i,i)=nb(i,i)+1

    !< r0_{aA} = <c_A^\dagger f_a>
    r0=transpose(dm(:norb,norb+1:))
    return

    end subroutine dget_rn_from_embeddm


    !> Get the impurity 
    !! physical density matrix (ni_{AB} = <c_A^\dagger c_B>), 
    !! bath density matrix (nb_{ab} = <f_b f_a^\dagger>
    !!                              = \delta_{ab}-<f_a^\dagger f_b>)
    !! and r0_{aA} = <c_A^\dagger f_a>.
    !! @param dm complete density matrix <{c_A f_b}^\dagger {c_B f_b}>
    !< Complex version.
    subroutine zget_rn_from_embeddm(dm,ni,nb,r0,norb)
    integer,intent(in)::norb !< dimension of {c_A} or {f_a}
    complex(q),intent(in)::dm(norb*2,norb*2)
    complex(q),intent(out)::ni(norb,norb),nb(norb,norb),r0(norb,norb)

    integer i

    !< Get physical density matrix.
    ni=dm(:norb,:norb)

    !< Bath auxiliary density matrix
    nb=-dm(norb+1:,norb+1:)
    forall(i=1:norb) nb(i,i)=nb(i,i)+1

    !< r0_{aA} = <c_A^\dagger f_a>
    r0=transpose(dm(:norb,norb+1:))
    return

    end subroutine zget_rn_from_embeddm


    subroutine zget_coul_exchange(norbs,v2e,v_j2e)
    integer,intent(in)::norbs
    complex(q),intent(in)::v2e(norbs,norbs,norbs,norbs)
    complex(q),intent(out)::v_j2e(norbs,norbs,norbs,norbs)

    integer i,j

    ! swap axes
    do i=1,norbs; do j=1,norbs
        v_j2e(:,i,:,j)=v2e(:,j,:,i)
        enddo; enddo
    v_j2e=v2e-v_j2e
    return

    end subroutine zget_coul_exchange


    subroutine zget_hf_pot(norbs,dm,vj2e,vhf)
    integer,intent(in)::norbs
    complex(q),intent(in)::vj2e(norbs,norbs,norbs,norbs),dm(norbs,norbs)
    complex(q),intent(out)::vhf(norbs,norbs)

    integer i,j

    do i=1,norbs; do j=1,i
        vhf(i,j)=vhf(i,j)+sum(vj2e(:,:,i,j)*dm)
        vhf(j,i)=conjg(vhf(i,j))
    enddo; enddo
    return

    end subroutine zget_hf_pot


    !< Using LU factorization to calculate the determinant of matrix a.
    function dget_determinant_a(a,n) result(zd)
    integer,intent(in)::n
    real(q),intent(in)::a(n,n)
    real(q) zd

    real(q) a_(n,n)
    integer ipiv(n),info,i

    a_=a
    call dgetrf(n,n,a_,n,ipiv,info)
    zd=1._q
    do i=1,n
        zd=zd*a_(i,i)
        if(ipiv(i)>i)then
            zd=-zd
        endif
    enddo
    return

    end function dget_determinant_a


    !< Using LU factorization to calculate the determinant of matrix a.
    function zget_determinant_a(a,n) result(zd)
    integer,intent(in)::n
    complex(q),intent(in)::a(n,n)
    complex(q) zd

    complex(q) a_(n,n)
    integer ipiv(n),info,i

    a_=a
    call zgetrf(n,n,a_,n,ipiv,info)
    zd=1._q
    do i=1,n
        zd=zd*a_(i,i)
        if(ipiv(i)>i)then
            zd=-zd
        endif
    enddo
    return

    end function zget_determinant_a


    ! convert a valence block of phi-vector to phi-matrix.
    subroutine dphi_vec_to_mat(v,a,n,bs_l,nl,norb,ibs,ibs_base)
    integer,intent(in)::n,nl,norb,ibs_base ! nl<n for mott phase
    real(q),intent(in)::v(n*nl)
    integer,intent(in)::bs_l(nl),ibs(0:ishft(1,norb)-1)
    real(q),intent(out)::a(n,n)
    
    integer j,bs,ib,nmax

    a=0 ! \phi_{\Gamma, n}
    nmax=ishft(1,norb)-1
    do j=1,nl
        bs=bs_l(j) 
        bs=nmax-bs ! post particle-hole transformation
        ib=ibs(bs)-ibs_base+1 ! get the correct index for phi-matrix
        a(:,ib)=v(j::nl)
    enddo
    return

    end subroutine dphi_vec_to_mat


    ! convert a valence block of phi-vector to phi-matrix.
    subroutine zphi_vec_to_mat(v,a,n,bs_l,nl,norb,ibs,ibs_base)
    integer,intent(in)::n,nl,norb,ibs_base ! nl<n for mott phase
    complex(q),intent(in)::v(n*nl)
    integer,intent(in)::bs_l(nl),ibs(0:ishft(1,norb)-1)
    complex(q),intent(out)::a(n,n)
    
    integer j,bs,ib,nmax

    a=0 ! \phi_{\Gamma, n}
    nmax=ishft(1,norb)-1
    do j=1,nl
        bs=bs_l(j) 
        bs=nmax-bs ! post particle-hole transformation
        ib=ibs(bs)-ibs_base+1 ! get the correct index for phi-matrix
        a(:,ib)=v(j::nl)
    enddo
    return

    end subroutine zphi_vec_to_mat


    ! convert a valence block of phi-vector to phi-matrix.
    subroutine dphi_mat_to_vec(v,a,n,bs_l,nl,norb,ibs,ibs_base)
    integer,intent(in)::n,nl,norb,ibs_base ! nl<n for mott phase
    real(q),intent(inout)::v(n*nl)
    integer,intent(in)::bs_l(nl),ibs(0:ishft(1,norb)-1)
    real(q),intent(in)::a(n,n)
    
    integer j,bs,ib,nmax

    nmax=ishft(1,norb)-1
    do j=1,nl
        bs=bs_l(j) 
        bs=nmax-bs ! post particle-hole transformation
        ib=ibs(bs)-ibs_base+1 ! get the correct index for phi-matrix
        v(j::nl)=v(j::nl)+a(:,ib)
    enddo
    return

    end subroutine dphi_mat_to_vec


    ! convert a valence block of phi-vector to phi-matrix.
    subroutine zphi_mat_to_vec(v,a,n,bs_l,nl,norb,ibs,ibs_base)
    integer,intent(in)::n,nl,norb,ibs_base ! nl<n for mott phase
    complex(q),intent(inout)::v(n*nl)
    integer,intent(in)::bs_l(nl),ibs(0:ishft(1,norb)-1)
    complex(q),intent(in)::a(n,n)
    
    integer j,bs,ib,nmax

    nmax=ishft(1,norb)-1
    do j=1,nl
        bs=bs_l(j) 
        bs=nmax-bs ! post particle-hole transformation
        ib=ibs(bs)-ibs_base+1 ! get the correct index for phi-matrix
        v(j::nl)=v(j::nl)+a(:,ib)
    enddo
    return

    end subroutine zphi_mat_to_vec


    ! convert a valence block of phi-vector to phi-matrix.
    subroutine dphi_vec_to_mat2(v,a,n,bs_r,nr,bs_l,nl,norb,ibs,ibs_base)
    integer,intent(in)::n,nr,nl,norb,ibs_base ! nl<n for mott phase
    real(q),intent(in)::v(nr*nl)
    integer,intent(in)::bs_r(nr),bs_l(nl),ibs(0:ishft(1,norb)-1)
    real(q),intent(out)::a(n,n)
    
    integer i,j,bs,ibl,ibr,nmax,iv

    a=0 ! \phi_{\Gamma, n}
    iv=0
    nmax=ishft(1,norb)-1
    do i=1,nr
        bs=bs_r(i)
        ibr=ibs(bs)-ibs_base+1 
        do j=1,nl
            bs=bs_l(j)
            bs=nmax-bs ! post particle-hole transformation
            ibl=ibs(bs)-ibs_base+1 ! get the correct index for phi-matrix
            iv=iv+1
            a(ibr,ibl)=v(iv)
        enddo
    enddo
    return

    end subroutine dphi_vec_to_mat2


    ! convert a valence block of phi-vector to phi-matrix.
    subroutine zphi_vec_to_mat2(v,a,n,bs_r,nr,bs_l,nl,norb,ibs,ibs_base)
    integer,intent(in)::n,nr,nl,norb,ibs_base ! nl<n for mott phase
    complex(q),intent(in)::v(nr*nl)
    integer,intent(in)::bs_r(nr),bs_l(nl),ibs(0:ishft(1,norb)-1)
    complex(q),intent(out)::a(n,n)
    
    integer i,j,bs,ibl,ibr,nmax,iv

    a=0 ! \phi_{\Gamma, n}
    iv=0
    nmax=ishft(1,norb)-1
    do i=1,nr
        bs=bs_r(i)
        ibr=ibs(bs)-ibs_base+1 
        do j=1,nl
            bs=bs_l(j)
            bs=nmax-bs ! post particle-hole transformation
            ibl=ibs(bs)-ibs_base+1 ! get the correct index for phi-matrix
            iv=iv+1
            a(ibr,ibl)=v(iv)
        enddo
    enddo
    return

    end subroutine zphi_vec_to_mat2


    ! convert a valence block of phi-vector to phi-matrix.
    subroutine dphi_mat_to_vec2(v,a,n,bs_r,nr,bs_l,nl,norb,ibs,ibs_base)
    integer,intent(in)::n,nr,nl,norb,ibs_base ! nl<n for mott phase
    real(q),intent(inout)::v(nr*nl)
    integer,intent(in)::bs_r(nr),bs_l(nl),ibs(0:ishft(1,norb)-1)
    real(q),intent(in)::a(n,n)
    
    integer i,j,bs,ibl,ibr,nmax,iv

    iv=0
    nmax=ishft(1,norb)-1
    do i=1,nr
        bs=bs_r(i)
        ibr=ibs(bs)-ibs_base+1 
        do j=1,nl
            bs=bs_l(j) 
            bs=nmax-bs ! post particle-hole transformation
            ibl=ibs(bs)-ibs_base+1 ! get the correct index for phi-matrix
            v(iv)=v(iv)+a(ibr,ibl)
        enddo
    enddo
    return

    end subroutine dphi_mat_to_vec2
   

    ! convert a valence block of phi-vector to phi-matrix.
    subroutine zphi_mat_to_vec2(v,a,n,bs_r,nr,bs_l,nl,norb,ibs,ibs_base)
    integer,intent(in)::n,nr,nl,norb,ibs_base ! nl<n for mott phase
    complex(q),intent(inout)::v(nr*nl)
    integer,intent(in)::bs_r(nr),bs_l(nl),ibs(0:ishft(1,norb)-1)
    complex(q),intent(in)::a(n,n)
    
    integer i,j,bs,ibl,ibr,nmax,iv

    iv=0
    nmax=ishft(1,norb)-1
    do i=1,nr
        bs=bs_r(i)
        ibr=ibs(bs)-ibs_base+1 
        do j=1,nl
            bs=bs_l(j) 
            bs=nmax-bs ! post particle-hole transformation
            ibl=ibs(bs)-ibs_base+1 ! get the correct index for phi-matrix
            v(iv)=v(iv)+a(ibr,ibl)
        enddo
    enddo
    return

    end subroutine zphi_mat_to_vec2



end module gutil
  
  
!*****************************************************************************
! down-hill simplex formalism
!*****************************************************************************
subroutine amoeba(p,y,ndim,ftol,funk,itmax,iter)
    implicit real(8)(a-h,o-z)
    integer ndim,itmax,iter
    real(8) p(ndim+1,ndim),y(ndim+1),ftol
    external :: funk
      
    parameter (alpha=1.0_8,beta=0.5_8,gamma=2.0_8)
    real(8),allocatable::pr(:),prr(:),pbar(:)
      
    allocate(pr(ndim),prr(ndim),pbar(ndim)); pr=0; prr=0; pbar=0
    mpts=ndim+1
    iter=0
    do
        ilo=1
        if(y(1).gt.y(2))then
            ihi=1
            inhi=2
        else
            ihi=2
            inhi=1
        endif
        do i=1,mpts
            if(y(i).lt.y(ilo)) ilo=i
            if(y(i).gt.y(ihi))then
                inhi=ihi
                ihi=i
            else if(y(i).gt.y(inhi))then
                if(i.ne.ihi) inhi=i
            endif
        enddo
          
        rtol=abs(y(ihi)-y(ilo))
        if(rtol.lt.ftol)then
            return
        endif
        if(iter.eq.itmax)then
            write(0,'("warning: amoeba exceeding maximum iterations.")')
        endif
        iter=iter+1
        pbar=0.d0
        do i=1,mpts
            if(i.ne.ihi)then
                do j=1,ndim
                    pbar(j)=pbar(j)+p(i,j)
                enddo
            endif
        enddo
        do j=1,ndim
            pbar(j)=pbar(j)/ndim
            pr(j)=(1.+alpha)*pbar(j)-alpha*p(ihi,j)
        enddo
        call funk(pr,ndim,ypr)
        if(ypr.le.y(ilo))then
            do j=1,ndim
                prr(j)=gamma*pr(j)+(1.-gamma)*pbar(j)
            enddo
            call funk(prr,ndim,yprr)
            if(yprr.lt.y(ilo))then
                do j=1,ndim
                    p(ihi,j)=prr(j)
                enddo
                y(ihi)=yprr
            else
                do j=1,ndim
                    p(ihi,j)=pr(j)
                enddo
                y(ihi)=ypr
            endif
        else if(ypr.ge.y(inhi))then
            if(ypr.lt.y(ihi))then
                do j=1,ndim
                    p(ihi,j)=pr(j)
                enddo
                y(ihi)=ypr
            endif
            do j=1,ndim
                prr(j)=beta*p(ihi,j)+(1.-beta)*pbar(j)
            enddo
            call funk(prr,ndim,yprr)
            if(yprr.lt.y(ihi))then
                do j=1,ndim
                    p(ihi,j)=prr(j)
                enddo
                y(ihi)=yprr
            else
                do i=1,mpts
                    if(i.ne.ilo)then
                        do j=1,ndim
                            pr(j)=0.5*(p(i,j)+p(ilo,j))
                            p(i,j)=pr(j)
                        enddo
                        call funk(pr,ndim,y(i))
                    endif
                enddo
            endif
        else
            do j=1,ndim
                p(ihi,j)=pr(j)
            enddo
            y(ihi)=ypr
        endif
    enddo
      
end subroutine amoeba
