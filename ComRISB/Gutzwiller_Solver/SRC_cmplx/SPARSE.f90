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

module sparse
    use gprec
    implicit none

    type dcsr_matrix
        integer :: nrow=0,ncol=0
        real(q),pointer :: a(:)=>null()
        integer,pointer :: i(:)=>null(),j(:)=>null()
    end type dcsr_matrix

    type zcsr_matrix
        integer :: nrow=0,ncol=0
        complex(q),pointer :: a(:)=>null()
        integer,pointer :: i(:)=>null(),j(:)=>null()
    end type zcsr_matrix

    type dcoo_matrix
        integer :: nrow=0,ncol=0,nnz=0
        real(q),pointer :: a(:)=>null()
        integer,pointer :: i(:)=>null(),j(:)=>null()
    end type dcoo_matrix

    type zvector
        complex(q),pointer :: a(:)=>null()
    end type zvector

    type ivector
        integer :: imin=0,imax=0
        integer,pointer :: i(:)=>null()
    end type ivector

    contains
    

    subroutine dcsr_syamuzx_sk(a,x,ax)
    type(dcsr_matrix),intent(in)::a
    complex(q),intent(in)::x(*)
    complex(q),intent(inout)::ax(*)
      
    integer i,k,j
      
    do i=1,a%nrow
        do k=a%i(i),a%i(i+1)-1
            j=a%j(k)
            ax(i)=ax(i)+a%a(k)*x(j)
            ! lower triangular part
            if(j==i)cycle
            ax(j)=ax(j)+a%a(k)*x(i)
        enddo
    enddo
    return
      
    end subroutine dcsr_syamuzx_sk


    subroutine dcsr_sylamuzx_sk(coef,a,x,ax)
    real(q),intent(in)::coef
    type(dcsr_matrix),intent(in)::a
    complex(q),intent(in)::x(*)
    complex(q),intent(inout)::ax(*)
      
    integer i,k,j
      
    do i=1,a%nrow
        do k=a%i(i),a%i(i+1)-1
            j=a%j(k)
            ax(i)=ax(i)+a%a(k)*x(j)*coef
            ! lower triangular part
            if(j==i)cycle
            ax(j)=ax(j)+a%a(k)*x(i)*coef
        enddo
    enddo
    return
      
    end subroutine dcsr_sylamuzx_sk


    subroutine zcsr_syamux_sk(a,x,ax)
    type(zcsr_matrix),intent(in)::a
    complex(q),intent(in)::x(*)
    complex(q),intent(inout)::ax(*)
      
    integer i,k,j
      
    do i=1,a%nrow
        do k=a%i(i),a%i(i+1)-1
            j=a%j(k)
            ax(i)=ax(i)+a%a(k)*x(j)
            ! lower triangular part
            if(j==i)cycle
            ax(j)=ax(j)+conjg(a%a(k))*x(i)
        enddo
    enddo
    return
      
    end subroutine zcsr_syamux_sk
      
      
    subroutine alloc_zcsr(a,nnz,nrow,ncol)
    integer nnz,nrow,ncol
    type(zcsr_matrix)::a
     
    integer ierr

    allocate(a%a(nnz),a%j(nnz),a%i(nrow+1),stat=ierr)
    if(ierr/=0)then
        write(0,'(" Error: allocation in alloc_zcsr with nnz = ",i0)')nnz
        stop 90
    endif
    a%nrow=nrow; a%ncol=ncol
    a%a=0; a%j=0; a%i=0
    return
      
    end subroutine alloc_zcsr


    subroutine copy_zcsr(a,b)
    type(zcsr_matrix)::a,b

    integer nnz

    nnz=a%i(a%nrow+1)-1
    b%nrow=a%nrow; b%ncol=a%ncol
    allocate(b%a(nnz),b%j(nnz),b%i(a%nrow+1))
    b%a=a%a(1:nnz)
    b%j=a%j(1:nnz)
    b%i=a%i
    return

    end subroutine copy_zcsr
    

    subroutine alloc_dcsr(a,nnz,nrow,ncol)
    integer nnz,nrow,ncol
    type(dcsr_matrix)::a
      
    allocate(a%a(nnz),a%j(nnz),a%i(nrow+1))
    a%nrow=nrow; a%ncol=ncol
    a%a=0; a%j=0; a%i=0
    return
      
    end subroutine alloc_dcsr


    subroutine alloc_dcoo(a,nnz,nrow,ncol)
    integer nnz,nrow,ncol
    type(dcoo_matrix)::a
      
    allocate(a%a(nnz),a%j(nnz),a%i(nnz))
    a%nrow=nrow; a%ncol=ncol; a%nnz=nnz
    a%a=0; a%j=0; a%i=0
    return
      
    end subroutine alloc_dcoo


    subroutine dealloc_zcsr(a)
    type(zcsr_matrix)::a
      
    deallocate(a%a,a%j,a%i)
    nullify(a%a,a%j,a%i)
    a%nrow=0; a%ncol=0 
    return
      
    end subroutine dealloc_zcsr


    subroutine dealloc_dcoo(a)
    type(dcoo_matrix)::a
      
    deallocate(a%a,a%j,a%i)
    a%nrow=0; a%ncol=0; a%nnz=0
    return
      
    end subroutine dealloc_dcoo


    subroutine dealloc_dcsr(a)
    type(dcsr_matrix)::a
      
    deallocate(a%a,a%j,a%i)
    nullify(a%a,a%j,a%i)
    a%nrow=0; a%ncol=0 
    return
      
    end subroutine dealloc_dcsr

      
    subroutine zs_dcoomux(zs,dcoo,x,y)
    type(dcoo_matrix),intent(in)::dcoo
    complex(q),intent(in)::zs,x(*)
    complex(q),intent(inout)::y(*)
    
    integer inz,i

    do inz=1,dcoo%nnz
        i=dcoo%i(inz)
        y(i)=y(i)+dcoo%a(inz)*x(dcoo%j(inz))*zs
    enddo
    return

    end subroutine zs_dcoomux


    subroutine zs_dcoohmux(zs,dcoo,x,y)
    type(dcoo_matrix),intent(in)::dcoo
    complex(q),intent(in)::zs,x(*)
    complex(q),intent(inout)::y(*)
    
    integer inz,j

    do inz=1,dcoo%nnz
        j=dcoo%j(inz)
        y(j)=y(j)+dcoo%a(inz)*x(dcoo%i(inz))*zs
    enddo
    return

    end subroutine zs_dcoohmux


    subroutine zvh_dcoo_zv(dcoo,x,y)
    type(dcoo_matrix),intent(in)::dcoo
    complex(q),intent(in)::x(*)
    complex(q),intent(inout)::y

    integer inz,i,j

    do inz=1,dcoo%nnz
        i=dcoo%i(inz)
        j=dcoo%j(inz)
        y=y+conjg(x(i))*dcoo%a(inz)*x(j)
    enddo
    return

    end subroutine zvh_dcoo_zv


    subroutine zv2h_sdcoo_zv1(zs,dcoo,x1,x2,y)
    type(dcoo_matrix),intent(in)::dcoo
    complex(q),intent(in)::zs,x1(*),x2(*)
    complex(q),intent(inout)::y

    integer inz,i,j

    do inz=1,dcoo%nnz
        i=dcoo%i(inz)
        j=dcoo%j(inz)
        y=y+zs*conjg(x2(i))*dcoo%a(inz)*x1(j)
    enddo
    return

    end subroutine zv2h_sdcoo_zv1


    subroutine zv2h_sdcooh_zv1(zs,dcoo,x1,x2,y)
    type(dcoo_matrix),intent(in)::dcoo
    complex(q),intent(in)::zs,x1(*),x2(*)
    complex(q),intent(inout)::y

    integer inz,i,j

    do inz=1,dcoo%nnz
        i=dcoo%i(inz)
        j=dcoo%j(inz)
        y=y+zs*conjg(x2(j))*dcoo%a(inz)*x1(i)
    enddo
    return

    end subroutine zv2h_sdcooh_zv1


    subroutine zvh_dcsr_zv(dcsr,x,y)
    type(dcsr_matrix),intent(in)::dcsr
    complex(q),intent(in)::x(*)
    complex(q),intent(inout)::y

    integer inz,i,j,k

    do i=1,dcsr%nrow; do k=dcsr%i(i),dcsr%i(i+1)-1
        j=dcsr%j(k)
        y=y+conjg(x(i))*dcsr%a(k)*x(j)
    enddo; enddo
    return

    end subroutine zvh_dcsr_zv


    subroutine zvh_sydcsr_zv(dcsr,x,y)
    type(dcsr_matrix),intent(in)::dcsr
    complex(q),intent(in)::x(*)
    complex(q),intent(inout)::y

    integer inz,i,j,k

    do i=1,dcsr%nrow; do k=dcsr%i(i),dcsr%i(i+1)-1
        j=dcsr%j(k)
        y=y+conjg(x(i))*dcsr%a(k)*x(j)
        if(i==j)cycle
        y=y+conjg(x(j))*dcsr%a(k)*x(i)
    enddo; enddo
    return

    end subroutine zvh_sydcsr_zv


    subroutine zvh_zcsr_zv(zcsr,x,y)
    type(zcsr_matrix),intent(in)::zcsr
    complex(q),intent(in)::x(*)
    complex(q),intent(inout)::y

    integer inz,i,j,k

    do i=1,zcsr%nrow; do k=zcsr%i(i),zcsr%i(i+1)-1
        j=zcsr%j(k)
        y=y+conjg(x(i))*zcsr%a(k)*x(j)
    enddo; enddo
    return

    end subroutine zvh_zcsr_zv


    ! zcsr = zs * dcsr
    subroutine zsdcsr(zs,dcsr,zcsr)
    complex(q),intent(in)::zs
    type(dcsr_matrix),intent(in)::dcsr
    type(zcsr_matrix),intent(out)::zcsr

    integer nnz

    nnz=dcsr%i(dcsr%nrow+1)-1
    call alloc_zcsr(zcsr,nnz,dcsr%nrow,dcsr%ncol)
    zcsr%a=zs*dcsr%a(1:nnz)
    zcsr%i=dcsr%i
    zcsr%j=dcsr%j(1:nnz)
    return

    end subroutine zsdcsr


    ! zcsrb = zs * zcsra
    subroutine zszcsr(zs,zcsra,zcsrb)
    complex(q),intent(in)::zs
    type(zcsr_matrix),intent(in)::zcsra
    type(zcsr_matrix),intent(out)::zcsrb

    integer nnz

    nnz=zcsra%i(zcsra%nrow+1)-1
    call alloc_zcsr(zcsrb,nnz,zcsra%nrow,zcsra%ncol)
    zcsrb%a=zs*zcsra%a(1:nnz)
    zcsrb%i=zcsra%i
    zcsrb%j=zcsra%j(1:nnz)
    return

    end subroutine zszcsr


    ! zsum = zcsr + zs * dcsr
    subroutine zcsrplsdcsr1(zs,dcsr,zcsr,zsum)
    complex(q),intent(in)::zs
    type(dcsr_matrix),intent(in)::dcsr
    type(zcsr_matrix),intent(in)::zcsr
    type(zcsr_matrix),intent(out)::zsum
   
    integer i,nnz,k,jpos
    integer iw(zcsr%ncol)

    nnz=dcsr%i(dcsr%nrow+1)+zcsr%i(zcsr%nrow+1)-2
    call alloc_zcsr(zsum,nnz,zcsr%nrow,zcsr%ncol)
    zsum%i(1)=1
    nnz=0
    do i=1,zcsr%nrow
        iw=0
        do k=zcsr%i(i),zcsr%i(i+1)-1
            nnz=nnz+1
            zsum%j(nnz)=zcsr%j(k)
            zsum%a(nnz)=zcsr%a(k)
            iw(zcsr%j(k))=nnz
        enddo
        do k=dcsr%i(i),dcsr%i(i+1)-1
            jpos=iw(dcsr%j(k))
            if(jpos==0)then
                nnz=nnz+1
                zsum%j(nnz)=dcsr%j(k)
                zsum%a(nnz)=zs*dcsr%a(k)
                iw(dcsr%j(k))=nnz
            else
                zsum%a(jpos)=zsum%a(jpos)+zs*dcsr%a(k)
            endif
        enddo
        zsum%i(i+1)=nnz+1
    enddo
    return

    end subroutine zcsrplsdcsr1


    ! zcsr += zs * dcsr
    subroutine zcsrplsdcsr(zs,dcsr,zcsr)
    complex(q),intent(in)::zs
    type(dcsr_matrix),intent(in)::dcsr
    type(zcsr_matrix),intent(inout)::zcsr
   
    type(zcsr_matrix)::zsum

    if(zcsr%nrow==0)then
        call zsdcsr(zs,dcsr,zcsr)
    else
        call zcsrplsdcsr1(zs,dcsr,zcsr,zsum)
        call dealloc_zcsr(zcsr)
        call copy_zcsr(zsum,zcsr)
        call dealloc_zcsr(zsum)
    endif
    return

    end subroutine zcsrplsdcsr


    ! zsum = zs*zcsra + zcsrb
    subroutine zcsraplszcsrb1(zs,zcsra,zcsrb,zsum)
    complex(q),intent(in)::zs
    type(zcsr_matrix),intent(in)::zcsra,zcsrb
    type(zcsr_matrix),intent(out)::zsum
   
    integer i,nnz,k,jpos
    integer iw(zcsra%ncol)

    nnz=zcsra%i(zcsra%nrow+1)+zcsrb%i(zcsrb%nrow+1)-2
    call alloc_zcsr(zsum,nnz,zcsra%nrow,zcsra%ncol)
    zsum%i(1)=1
    nnz=0
    do i=1,zcsrb%nrow
        iw=0
        do k=zcsrb%i(i),zcsrb%i(i+1)-1
            nnz=nnz+1
            zsum%j(nnz)=zcsrb%j(k)
            zsum%a(nnz)=zcsrb%a(k)
            iw(zcsrb%j(k))=nnz
        enddo
        do k=zcsra%i(i),zcsra%i(i+1)-1
            jpos=iw(zcsra%j(k))
            if(jpos==0)then
                nnz=nnz+1
                zsum%j(nnz)=zcsra%j(k)
                zsum%a(nnz)=zs*zcsra%a(k)
                iw(zcsra%j(k))=nnz
            else
                zsum%a(jpos)=zsum%a(jpos)+zs*zcsra%a(k)
            endif
        enddo
        zsum%i(i+1)=nnz+1
    enddo
    return

    end subroutine zcsraplszcsrb1


    ! zsum += zs * zcsr
    subroutine zcsrplszcsra(zs,zcsr,zsum)
    complex(q),intent(in)::zs
    type(zcsr_matrix),intent(in)::zcsr
    type(zcsr_matrix),intent(inout)::zsum
   
    type(zcsr_matrix)::zsum_

    if(zsum%nrow==0)then
        call zszcsr(zs,zcsr,zsum)
    else
        call zcsraplszcsrb1(zs,zcsr,zsum,zsum_)
        call dealloc_zcsr(zsum)
        call copy_zcsr(zsum_,zsum)
        call dealloc_zcsr(zsum_)
    endif
    return

    end subroutine zcsrplszcsra


    subroutine calc_nnz_csrab(nrow,ncol,ia,ja,ib,jb,nnz)
    integer,intent(in)::nrow,ncol,ia(*),ja(*),ib(*),jb(*)
    integer,intent(out)::nnz

    integer i,j,j1,k,k1
    logical lzero(ncol)

    nnz=0
    do i=1,nrow
        lzero=.true.
        do j1=ia(i),ia(i+1)-1
            j=ja(j1)
            do k1=ib(j),ib(j+1)-1
                k=jb(k1)
                if(lzero(k))then
                    nnz=nnz+1
                    lzero(k)=.false.
                endif
            enddo
        enddo
    enddo
    return

    end subroutine calc_nnz_csrab


    ! zcsrab = zcsra * zcsrb
    subroutine zcsramuzcsrb1(zcsra,zcsrb,zcsrab)
    type(zcsr_matrix),intent(in)::zcsra,zcsrb
    type(zcsr_matrix),intent(out)::zcsrab
   
    integer i,nnz,j1,j,k1,k,jpos
    integer iw(zcsrb%ncol)

    call calc_nnz_csrab(zcsra%nrow,zcsrb%ncol,zcsra%i,zcsra%j,zcsrb%i, &
            &zcsrb%j,nnz)
    call alloc_zcsr(zcsrab,nnz,zcsra%nrow,zcsrb%ncol)
    zcsrab%i(1)=1
    nnz=0
    do i=1,zcsra%nrow
        iw=0
        do j1=zcsra%i(i),zcsra%i(i+1)-1
            j=zcsra%j(j1)
            do k1=zcsrb%i(j),zcsrb%i(j+1)-1
                k=zcsrb%j(k1)
                jpos=iw(k)
                if(jpos==0)then
                    nnz=nnz+1
                    zcsrab%j(nnz)=k
                    zcsrab%a(nnz)=zcsra%a(j1)*zcsrb%a(k1)
                    iw(k)=nnz
                else
                    zcsrab%a(jpos)=zcsrab%a(jpos)+zcsra%a(j1)*zcsrb%a(k1)
                endif
            enddo
        enddo
        zcsrab%i(i+1)=nnz+1
    enddo
    return

    end subroutine zcsramuzcsrb1


    ! zsum += s * zcsra * zcsrb
    subroutine zcsradd_szcsramuzcsrb(zs,zcsra,zcsrb,zsum)
    complex(q),intent(in)::zs
    type(zcsr_matrix),intent(in)::zcsra,zcsrb
    type(zcsr_matrix),intent(inout)::zsum
   
    type(zcsr_matrix)::zcsrab

    if(zsum%nrow==0)then
        call zcsramuzcsrb1(zcsra,zcsrb,zsum)
        if(abs(zs-1._q)>1.d-12)then
            zsum%a=zs*zsum%a
        endif
    else
        call zcsramuzcsrb1(zcsra,zcsrb,zcsrab)
        call zcsrplszcsra(zs,zcsrab,zsum) 
        call dealloc_zcsr(zcsrab)
    endif
    return

    end subroutine zcsradd_szcsramuzcsrb


    ! zcsr = s_{a,b}*(dcsra(a,b) + dcsrb(a,b))
    subroutine zsdcsr_22darray_sum(n,s,dcsra,dcsrb,zcsr)
    integer,intent(in)::n
    complex(q),intent(in)::s(n,n)
    type(dcsr_matrix),intent(in)::dcsra(n,n),dcsrb(n,n)
    type(zcsr_matrix),intent(out)::zcsr

    integer i,j

    do i=1,n; do j=1,n
        if(abs(s(i,j))<1.d-12)cycle
        call zcsrplsdcsr(s(i,j),dcsra(i,j),zcsr)
        call zcsrplsdcsr(s(i,j),dcsrb(i,j),zcsr)
    enddo; enddo
    return

    end subroutine zsdcsr_22darray_sum


    subroutine zcsr_array_sqsum(n,zcsr,zsum)
    integer,intent(in)::n
    type(zcsr_matrix),intent(in)::zcsr(n)
    type(zcsr_matrix),intent(out)::zsum

    integer i
    complex(q),parameter::z1=cmplx(1._q,0._q)

    do i=1,n
        call zcsradd_szcsramuzcsrb(z1,zcsr(i),zcsr(i),zsum)
    enddo
    return

    end subroutine zcsr_array_sqsum

      
end module sparse
