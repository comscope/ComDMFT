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

module gmpi
    use gprec
    implicit none
#ifdef mpi_mode
    include "mpif.h"
#endif
    private
    type g_mpi
        integer::master=0
        integer::nprocs=1
        integer::myrank=0
        character*10 cpuid !< myrank in string representation
        integer::ipar=0 !< ipar=0: simple k-parallel; 1: openmp, 2: wien-mpi
        integer::io=-1

        !< k-parallel of wien2k
        integer::nvec=1
        integer,pointer::kvec(:,:)=>null() !< nkp, nkstart0
    end type g_mpi
      
    type(g_mpi) gp
      
    public::gp,init_gmpi

#ifdef mpi_mode
    public::max_master_mpi,min_master_mpi,sum_master_mpi,max_all_mpi, &
            &min_all_mpi,sum_all_mpi
    public::zsum_all_mpi,dsum_all_mpi,gmpi_barrier

    interface max_master_mpi
        module procedure imax1_master_mpi,imax_master_mpi,dmax1_master_mpi, &
                &dmax_master_mpi
    end interface max_master_mpi

    interface min_master_mpi
        module procedure imin1_master_mpi,imin_master_mpi
    end interface min_master_mpi

    interface sum_master_mpi
        module procedure isum1_master_mpi,isum_master_mpi,dsum1_master_mpi, &
                &dsum_master_mpi,zsum1_master_mpi,zsum_master_mpi
    end interface sum_master_mpi

    interface max_all_mpi
        module procedure imax1_all_mpi,imax_all_mpi,dmax1_all_mpi, &
                &dmax_all_mpi
    end interface max_all_mpi

    interface min_all_mpi
        module procedure imin1_all_mpi,imin_all_mpi
    end interface min_all_mpi

    interface sum_all_mpi
        module procedure isum1_all_mpi,isum_all_mpi,dsum1_all_mpi, &
                &dsum_all_mpi,zsum1_all_mpi,zsum_all_mpi
    end interface sum_all_mpi

#endif

    contains
     
    subroutine init_gmpi()
    integer ierr

#ifdef mpi_mode
    call mpi_comm_rank(mpi_comm_world,gp%myrank,ierr)
    call mpi_comm_size(mpi_comm_world,gp%nprocs,ierr)
#endif
    if(gp%myrank==gp%master)then
        gp%io=6
    endif
    return

    end subroutine init_gmpi


#ifdef mpi_mode
      
    !*************************************************************************
    subroutine gmpi_barrier()
    integer ierr
      
    call mpi_barrier(mpi_comm_world,ierr)
      
    end subroutine gmpi_barrier
      
    !*************************************************************************
    subroutine imax1_master_mpi(n)
    integer n
      
    integer i(1)
      
    i(1)=n; call imax_master_mpi(i,1); n=i(1)
    return
      
    end subroutine imax1_master_mpi
      
    !*************************************************************************
    subroutine imax_master_mpi(i,n)
    integer n,i(n)
      
    integer ierr,maxi(n)
      
    call mpi_reduce(i,maxi,n,mpi_integer,mpi_max,gp%master,mpi_comm_world,ierr)
    if(gp%myrank.eq.gp%master)i=maxi
    return
      
    end subroutine imax_master_mpi
      
    !*************************************************************************
    subroutine dmax1_master_mpi(a)
    real(q) a
      
    real(q) b(1)
      
    b(1)=a; call dmax_master_mpi(b,1); a=b(1)
    return
      
    end subroutine dmax1_master_mpi
      
    !*************************************************************************
    subroutine dmax_master_mpi(a,n)
    integer n
    real(q) a(n)
      
    integer ierr
    real(q) maxa(n)
      
    call mpi_reduce(a,maxa,n,mpi_double_precision,mpi_max,gp%master, &
            &mpi_comm_world,ierr)
    if(gp%myrank.eq.gp%master)a=maxa
    return
      
    end subroutine dmax_master_mpi
      
    !*************************************************************************
    subroutine imin1_master_mpi(n)
    integer n
      
    integer m(1)
      
    m(1)=n; call imin_master_mpi(m,1); n=m(1)
    return
      
    end subroutine imin1_master_mpi
      
    !*************************************************************************
    subroutine imin_master_mpi(i,n)
    integer n,i(n)
      
    integer ierr,mini(n)
      
    call mpi_reduce(i,mini,n,mpi_integer,mpi_min,gp%master,mpi_comm_world, &
            &ierr)
    if(gp%myrank.eq.gp%master)i=mini
    return
      
    end subroutine imin_master_mpi
      
    !*************************************************************************
    subroutine isum1_master_mpi(i)
    integer i
      
    integer j(1)
      
    j(1)=i; call isum_master_mpi(j,1); i=j(1)
    return
      
    end subroutine isum1_master_mpi
      
    !*************************************************************************
    subroutine isum_master_mpi(i,n)
    integer n,i(n)
      
    integer ierr
    integer,allocatable::j(:)
      
    allocate(j(n)); j=0
    call mpi_reduce(i,j,n,mpi_integer,mpi_sum,gp%master,mpi_comm_world,ierr)
    if(gp%myrank.eq.gp%master)i=j
    deallocate(j)
    return
      
    end subroutine isum_master_mpi
      
    !*************************************************************************
    subroutine dsum1_master_mpi(a)
    real(q) a
      
    real(q) b(1)
      
    b(1)=a; call dsum_master_mpi(b,1); a=b(1)
    return
      
    end subroutine dsum1_master_mpi
      
    !*************************************************************************
    subroutine dsum_master_mpi(a,n)
    integer n
    real(q) a(n)
      
    integer ierr
    real(q),allocatable::b(:)
      
    allocate(b(n)); b=0
    call mpi_reduce(a,b,n,mpi_double_precision,mpi_sum,gp%master, &
            &mpi_comm_world,ierr)
    if(gp%myrank.eq.gp%master)a=b
    deallocate(b)
    return
      
    end subroutine dsum_master_mpi
      
    !*************************************************************************
    subroutine zsum1_master_mpi(a)
    complex(q) a
      
    complex(q) b(1)
      
    b(1)=a; call zsum_master_mpi(b,1); a=b(1)
    return
      
    end subroutine zsum1_master_mpi
      
    !*************************************************************************
    subroutine zsum_master_mpi(a,n)
    integer n
    complex(q) a(n)
      
    integer ierr
    complex(q),allocatable::b(:)
      
    allocate(b(n)); b=0
    call mpi_reduce(a,b,n,mpi_double_complex,mpi_sum,gp%master, &
            &mpi_comm_world,ierr)
    if(gp%myrank.eq.gp%master)a=b
    deallocate(b)
    return
      
    end subroutine zsum_master_mpi
      
    !*************************************************************************
    subroutine imax1_all_mpi(i)
    integer i
      
    integer j(1)
      
    j(1)=i; call imax_all_mpi(j,1); i=j(1)
    return
      
    end subroutine imax1_all_mpi
      
    !*************************************************************************
    subroutine imax_all_mpi(i,n)
    integer n,i(n)
      
    integer ierr,maxi(n)
      
    call mpi_allreduce(i,maxi,n,mpi_integer,mpi_max,mpi_comm_world,ierr)
    i=maxi
    return
      
    end subroutine imax_all_mpi
      
    !*************************************************************************
    subroutine imin1_all_mpi(i)
    integer i
      
    integer j(1)
      
    j(1)=i; call imin_all_mpi(j,1); i=j(1)
    return
      
    end subroutine imin1_all_mpi
      
    !*************************************************************************
    subroutine imin_all_mpi(i,n)
    integer n,i(n)
      
    integer ierr,mini(n)
      
    call mpi_allreduce(i,mini,n,mpi_integer,mpi_min,mpi_comm_world,ierr)
    i=mini
    return
      
    end subroutine imin_all_mpi
      
    !*************************************************************************
    subroutine dmax1_all_mpi(a)
    real(q) a
      
    real(q) b(1)
      
    b(1)=a; call dmax_all_mpi(b,1); a=b(1)
    return
      
    end subroutine dmax1_all_mpi
      
    !*************************************************************************
    subroutine dmax_all_mpi(a,n)
    integer n
    real(q) a(n)
      
    integer ierr
    real(q) maxa(n)
      
    call mpi_allreduce(a,maxa,n,mpi_double_precision,mpi_max, &
            &mpi_comm_world,ierr)
    a=maxa
    return
      
    end subroutine dmax_all_mpi
      
    !*************************************************************************
    subroutine isum1_all_mpi(i)
    integer i
      
    integer j(1)
      
    j(1)=i; call isum_all_mpi(j,1); i=j(1)
    return
      
    end subroutine isum1_all_mpi
      
    !*************************************************************************
    subroutine isum_all_mpi(i,n)
    integer n,i(n)
      
    integer ierr
    integer,allocatable::j(:)
      
    allocate(j(n)); j=0
    call mpi_allreduce(i,j,n,mpi_integer,mpi_sum,mpi_comm_world,ierr)
    i=j; deallocate(j)
    return
      
    end subroutine isum_all_mpi
      
    !*************************************************************************
    subroutine dsum1_all_mpi(a)
    real(q) a
      
    real(q) b(1)
      
    b(1)=a; call dsum_all_mpi(b,1); a=b(1)
    return
      
    end subroutine dsum1_all_mpi
      
    !*************************************************************************
    subroutine dsum_all_mpi(a,n)
    integer n
    real(q) a(n)
      
    integer ierr
    real(q),allocatable::b(:)
      
    allocate(b(n)); b=0
    call mpi_allreduce(a,b,n,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    a=b; deallocate(b)
    return
      
    end subroutine dsum_all_mpi
      
    !*************************************************************************
    subroutine zsum1_all_mpi(a)
    complex(q) a
      
    complex(q) b(1)
      
    b(1)=a; call zsum_all_mpi(b,1); a=b(1)
    return
      
    end subroutine zsum1_all_mpi
      
    !*************************************************************************
    subroutine zsum_all_mpi(a,n)
    integer n
    complex(q) a(n)
      
    integer ierr
    complex(q),allocatable::b(:)
      
    allocate(b(n)); b=0
    call mpi_allreduce(a,b,n,mpi_double_complex,mpi_sum,mpi_comm_world,ierr)
    a=b; deallocate(b)
    return
      
    end subroutine zsum_all_mpi
      
#endif
      
      
end module gmpi
