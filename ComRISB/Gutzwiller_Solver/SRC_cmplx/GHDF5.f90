module ghdf5
    use gprec
    use ghdf5_base
    use sparse
    use hdf5
    implicit none
    private
    public :: gh5_write_compound

    interface gh5_write_compound
      module procedure gh5_write_dcsr_matrix, gh5_write_zcsr_matrix, &
              & gh5_write_dcoo_matrix
    end interface

    contains

    subroutine gh5_write_dcsr_matrix(dcsr, path, f_id)
    type(dcsr_matrix), target, intent(in) :: dcsr
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id

    integer nnz

    call gh5_create_group(path, f_id)
    call gh5_write(dcsr%nrow, path//"/nrow", f_id)
    call gh5_write(dcsr%ncol, path//"/ncol", f_id)
    call gh5_write(1, path//"/base", f_id)
    call gh5_write(dcsr%i, dcsr%nrow + 1, path//"/indptr", f_id)
    nnz = dcsr%i(dcsr%nrow + 1) -1
    call gh5_write(dcsr%j, nnz, path//"/indices", f_id)
    call gh5_write(dcsr%a, nnz, path//"/data", f_id)
    return

    end subroutine gh5_write_dcsr_matrix


    subroutine gh5_write_dcoo_matrix(dcoo, path, f_id)
    type(dcoo_matrix), target, intent(in) :: dcoo
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id

    call gh5_create_group(path, f_id)
    call gh5_write(dcoo%nrow, path//"/nrow", f_id)
    call gh5_write(dcoo%ncol, path//"/ncol", f_id)
    call gh5_write(dcoo%nnz, path//"/nnz", f_id)
    call gh5_write(1, path//"/base", f_id)
    call gh5_write(dcoo%i, dcoo%nnz, path//"/i", f_id)
    call gh5_write(dcoo%j, dcoo%nnz, path//"/j", f_id)
    call gh5_write(dcoo%a, dcoo%nnz, path//"/data", f_id)
    return

    end subroutine gh5_write_dcoo_matrix


    subroutine gh5_write_zcsr_matrix(zcsr, path, f_id)
    type(zcsr_matrix), target, intent(in) :: zcsr
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id

    integer nnz

    call gh5_create_group(path, f_id)
    call gh5_write(zcsr%nrow, path//"/nrow", f_id)
    call gh5_write(zcsr%ncol, path//"/ncol", f_id)
    call gh5_write(1, path//"/base", f_id)
    call gh5_write(zcsr%i, zcsr%nrow + 1, path//"/indptr", f_id)
    nnz = zcsr%i(zcsr%nrow + 1) -1
    call gh5_write(zcsr%j, nnz, path//"/indices", f_id)
    call gh5_write(zcsr%a, nnz, path//"/data", f_id)
    return

    end subroutine gh5_write_zcsr_matrix



end module ghdf5
