subroutine ghybrd(f,n,x,fvec,rtol,epsfcn,io)
    use gprec
    implicit none
    integer,intent(in) :: n,io
    real(q),intent(in) :: rtol,epsfcn
    real(q),intent(inout) :: x(n),fvec(n)
    external :: f
      
    integer info
      
    call hybrd1(f,n,x,fvec,rtol,epsfcn,info)
    select case(info)
    case(0)
        write(io,'(" ghybrd: error! improper input parameters.")')
        write(0,'(" ghybrd: error! improper input parameters.")')
    case(1)
        write(io,'(" ghybrd: success.")')
    case(2)
        write(io,'(" ghybrd: warning! number of calls to fcn has &
                &reached or exceeded 200*(n+1).")')
        write(0,'(" ghybrd: warning! number of calls to fcn has &
                &reached or exceeded 200*(n+1).")')
    case(3)
        write(io,'(" ghybrd: warning! tol is too small. no further &
                &improvement in the approximate solution x is possible.")')
        write(0,'(" ghybrd: warning! tol is too small. no further &
                &improvement in the approximate solution x is possible.")')
    case(4)
        write(io,'(" ghybrd: warning! iteration is not making &
                &good progress.")')
        write(0,'(" ghybrd: warning! iteration is not making &
                &good progress.")')
    end select
    return
      
end subroutine ghybrd
