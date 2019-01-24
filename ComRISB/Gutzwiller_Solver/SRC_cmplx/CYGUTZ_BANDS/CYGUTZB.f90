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

program cygutz

    use ghdf5_base, only: gh5_init, gh5_end
    use gparam
    use bandstru
    use warehouse
    use gmpi
    use psave
    use gtime
    implicit none
#ifdef mpi_mode
    include "mpif.h"
#endif
    integer ierr
    
    call set_time_point(1,1)
#ifdef mpi_mode
    call mpi_init(ierr)
#endif
    call gh5_init()
    call init_gmpi()
    if(gp%io>0)then
        open(gp%io,file='GUTZ.LOG',status='replace')
#ifdef with_compile_date
        write(gp%io,'(" cygutz (atomic rydberg units) vmpi_5.0 built on", &
                &a12)')__DATE__
#else
        write(gp%io,'(" cygutz (atomic rydberg units) vmpi_5.0")')
#endif
    endif

    call set_gparam(gp%io)
    call init_warehouse(gp%io)
    call set_bnd_info(gp%io)
    call set_gmpi() !< kpt%diml was set here.
    call alloc_bnd()
    call read_bare_hamiltonian()
    call rotate_bare_hamiltonian()

    ! single out the local one-body part.
    call rm_h1e_from_bare_hamiltonian()

    ! check the quasi-particle part with initial guess of {R, lambda}
    call init_wh_x(gp%io, 1._q)
    call map_wh_bnd_matrix(wh%r,bnd%r,.false.)
    call map_wh_bnd_matrix(wh%la1,bnd%la1,.false.)
    call calc_band_all(gp%io)

    ! evaluate fermi energy, possibly more accurate ef for dos calculations.
    call gutz_fermi(gp%io)

    !! Save important data for analysis.
    call postsave_gbands()

    call set_time_point(2,1)
    call print_time_usage('total',1,gp%io)
    if(gp%io>0)then
        close(gp%io)
    endif
    call gh5_end()
#ifdef mpi_mode
    call mpi_finalize(ierr)
#endif
      
end program cygutz
