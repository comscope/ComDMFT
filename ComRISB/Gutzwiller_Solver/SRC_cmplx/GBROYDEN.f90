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
 
module gbroyden
    use gprec
    implicit none
      
    type broyden_data
        integer n,it_reset
        logical :: lfirst=.true.
        real(q) :: alpha_mix
        real(q),allocatable :: x0(:), f0(:), amix(:,:)
    end type broyden_data


    contains
    !*************************************************************************
    subroutine init_broyden_data(br_data,n,alpha_mix,it_reset)
    type(broyden_data) :: br_data
    integer,intent(in) :: n,it_reset
    real(q),intent(in) :: alpha_mix
      
    integer i
      
    if(.not.br_data%lfirst)return
    br_data%n=n
    allocate(br_data%x0(n), br_data%f0(n), br_data%amix(n,n))
    br_data%alpha_mix = alpha_mix
    br_data%it_reset = it_reset
    return
      
    end subroutine init_broyden_data


    subroutine clean_broyden_data(br_data)
    type(broyden_data) :: br_data

    deallocate(br_data%x0, br_data%f0, br_data%amix)

    end subroutine clean_broyden_data

      
    !*************************************************************************
    subroutine inv_broyden_predict(br_data,x,f)
    type(broyden_data) :: br_data
    real(q) x(br_data%n),f(br_data%n)
      
    integer i,j
    real(q) res,vnorm,vfac
    real(q) dx(br_data%n),df(br_data%n),av(br_data%n),va(br_data%n)
    if(br_data%lfirst)then
        br_data%lfirst=.false.
        br_data%amix=0
        do i=1,br_data%n
            br_data%amix(i,i)=br_data%alpha_mix
        enddo
    else
        dx=x-br_data%x0; df=f-br_data%f0
        av=matmul(br_data%amix,df)
        va=df
        vnorm=dot_product(df,df)
        vfac = 1/max(vnorm,1.e-20_q)
        av=av-dx
        do i=1,br_data%n
            do j=1,br_data%n
                br_data%amix(i,j) = br_data%amix(i,j) - vfac*av(i)*va(j)
            enddo
        enddo
    endif
    ! compute correction to vin from inverse jacobian.
    av=matmul(br_data%amix,f)
    ! no too big step
    res=maxval(abs(av))
    if(res.gt.1._q)then
        res=1._q/res
        av=av*res
        br_data%amix=br_data%amix*res
    endif
    ! update x
    br_data%x0=x; br_data%f0=f
    x=x-av
    return
      
    end subroutine inv_broyden_predict
      
    !*************************************************************************
    subroutine broyden_mix(f,n,x,fvec,alpha,ireset,mode)
    integer,intent(in) :: n,mode,ireset
    real(q),intent(inout) :: x(n),fvec(n)
    real(q),intent(in) ::  alpha
    external :: f
     
    type(broyden_data) :: br_data
    integer it,iflag,i_inner,nbadmove
    real(q) alpha_mix
    integer,parameter :: i_inner_max = 27
    real(q) f_max, f_max_best, x_best(n)
      
    f_max_best = 1.e4_q; i_inner=0
    alpha_mix = alpha
    if(mode==-1)then
        call init_broyden_data(br_data,n,alpha,ireset)
    endif
    it=0
    do
        it=it+1
        call f(n,x,fvec,iflag)
        f_max = maxval(abs(fvec))
        if(iflag==-1.or.i_inner>i_inner_max+1.or.alpha_mix<1.e-8_q)then
            if(f_max<f_max_best)then
                exit
            else
                x = x_best
                call f(n,x,fvec,iflag)
                exit
            endif
        endif
        if(f_max<f_max_best)then
            i_inner = 0
            f_max_best = f_max
            x_best = x
        else
            i_inner = i_inner + 1
            if(i_inner>i_inner_max.and.mode==1)then
                alpha_mix = alpha_mix/2
                x = x_best
                i_inner = 0
            endif
        endif
        if(mode==-1)then ! inverse broyden
            ! reset broyden
            if(mod(it,br_data%it_reset)==0)then
                br_data%lfirst=.true.
            endif
            call inv_broyden_predict(br_data,x,fvec)
        elseif(mode==1)then
            x = x+alpha*fvec
        else
            stop " error in broyden mix: unsupported mode!"
        endif
    enddo
    if(mode==-1)then
        call clean_broyden_data(br_data)
    endif
    return
      
    end subroutine broyden_mix
      
end module gbroyden
