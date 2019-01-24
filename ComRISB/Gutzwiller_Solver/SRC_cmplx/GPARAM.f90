module gparam
    use gprec
    implicit none

    integer::g_write=-100, g_maxiter=-100, g_iembeddiag=-1000, &
            &g_updaterho=1, g_imix=-100

    contains

    subroutine set_gparam(io)
    integer,intent(in)::io

    integer idx
    character(255) cmd

    call get_command(cmd)
    idx=index(cmd,'-w ')
    if(idx>0)then
        read(cmd(idx+2:),*)g_write
        if(io>0)then
            write(io,'(" command argument g_write = ",i2)')g_write
        endif
    endif

    idx=index(cmd,'-m ')
    if(idx>0)then
        read(cmd(idx+2:),*)g_imix
        if(io>0)then
            write(io,'(" command argument g_imix = ",i2)')g_imix
        endif
    endif

    idx=index(cmd,'-n ')

    if(idx>0)then
        read(cmd(idx+2:),*)g_maxiter
        if(io>0)then
            write(io,'(" command argument g_maxiter = ",i0)')g_maxiter
        endif
    endif

    idx=index(cmd,'-d ')
    if(idx>0)then
        read(cmd(idx+2:),*)g_iembeddiag
        if(io>0)then
            write(io,'(" command argument g_iembeddiag = ",i2)')g_iembeddiag
        endif
    endif

    idx=index(cmd,'-r ')
    if(idx>0)then
        read(cmd(idx+2:),*)g_updaterho
        if(io>0)then
            write(io,'(" command argument g_updaterho = ",i2)')g_updaterho
        endif
    endif
    return

    end subroutine


end module gparam
