module gtime
    use gprec
    implicit none
    type time_tag
        real::ta(2,7) ! Depth of timing tags.
        real(q)::tb(2,7)
    end type time_tag

    type(time_tag)::ttag

    contains


    subroutine set_time_point(it,i)
    integer,intent(in)::i,it

    integer::tib,tirate

    call cpu_time(ttag%ta(it,i))
    call system_clock(tib,tirate)
    ttag%tb(it,i)=dble(tib)/dble(tirate)
    return

    end subroutine set_time_point


    subroutine print_time_usage(task,i,io)
    character(*),intent(in)::task
    integer,intent(in)::i,io

    integer nh,nm
    real(q) time
    real(q),parameter::s_per_m=60._4,s_per_h=3600._4
    character(100) fmt

    if(io<0)return
    time=ttag%ta(2,i)-ttag%ta(1,i)
    nh=int(time/s_per_h); time=mod(time,s_per_h)
    nm=int(time/s_per_m); time=mod(time,s_per_m)
    write(fmt,'("(1x,a",i0,",a6,i5,a3,i2,a3,f5.2,a8)")')len(task)
    write(io,fmt)task," takes",nh," h ",nm," m ",time," s(cpu)."
    time=ttag%tb(2,i)-ttag%tb(1,i)
    nh=int(time/s_per_h); time=mod(time,s_per_h)
    nm=int(time/s_per_m); time=mod(time,s_per_m)
    write(fmt,'("(1x,",i0,"x,6x,i5,a3,i2,a3,f5.2,a9)")')len(task)
    write(io,fmt)nh," h ",nm," m ",time," s(wall)."
    return

    end subroutine print_time_usage


end module gtime
