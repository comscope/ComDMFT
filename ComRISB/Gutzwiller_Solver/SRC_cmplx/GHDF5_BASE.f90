module ghdf5_base
    use hdf5
    use gprec
    implicit none
    private
    public :: gh5_init, gh5_end, gh5_write, gh5_read, gh5_create_group, &
            & gh5_open_w, gh5_open_r, gh5_open_rw, &
            & gh5_close, f_id, gh5_err
      
    integer(hid_t) :: f_id=-1
    integer(hid_t) :: group_id, dset_id, dspace_id
    integer(hid_t) :: dcomplex_id, string_id
    integer(hsize_t) :: gh5_dims(7), gh5_size, gh5_offset
    type(c_ptr) :: f_ptr
    integer :: gh5_err
    logical :: lgh5_init=.false.
      
    interface gh5_write
        module procedure gh5_write_iscalar, gh5_write_dscalar, &
                & gh5_write_sscalar, &
                & gh5_write_1d_iarray, gh5_write_2d_iarray, &
                & gh5_write_1d_darray, gh5_write_2d_darray, &
                & gh5_write_4d_darray, &
                & gh5_write_1d_zarray, gh5_write_2d_zarray, &
                & gh5_write_4d_zarray, &
                & gh5_write_3d_darray, gh5_write_3d_zarray, &
                & gh5_write_5d_darray, &
                & gh5_write_5d_zarray, &
                & gh5_write_1d_sarray
    end interface
      
    interface gh5_read
        module procedure gh5_read_iscalar, gh5_read_dscalar, &
                & gh5_read_sscalar, &
                & gh5_read_1d_iarray, gh5_read_2d_iarray, &
                & gh5_read_1d_darray, gh5_read_2d_darray, &
                & gh5_read_4d_darray, &
                & gh5_read_1d_zarray, gh5_read_2d_zarray, &
                & gh5_read_4d_zarray, &
                & gh5_read_3d_darray, gh5_read_3d_zarray, &
                & gh5_read_3d_iarray, &
                & gh5_read_5d_darray, gh5_read_5d_zarray, &
                & gh5_read_1d_sarray
    end interface
      
    contains
     

    subroutine gh5_init()
      
    !< initialize fortran interface.
    call h5open_f(gh5_err)
    !< create a new file using default properties.
    call set_dcomplex_id()
    call h5tcopy_f(h5t_fortran_s1, string_id, gh5_err)
    lgh5_init=.true.
    return
      
    end subroutine gh5_init
      

    subroutine gh5_open_w(gh5_file, f_id)
    character(*),intent(in) :: gh5_file
    integer(hid_t),intent(out) :: f_id
      
    if(.not.lgh5_init)call gh5_init()
    !< create a new file using default properties.
    call h5fcreate_f(gh5_file, h5f_acc_trunc_f, f_id, gh5_err)
    return
      
    end subroutine gh5_open_w


    subroutine gh5_open_rw(gh5_file, f_id)
    character(*),intent(in) :: gh5_file
    integer(hid_t),intent(out) :: f_id
      
    if(.not.lgh5_init)call gh5_init()
    !< create a new file using default properties.
    call h5fcreate_f(gh5_file, h5f_acc_rdwr_f, f_id, gh5_err)
    return
      
    end subroutine gh5_open_rw


    subroutine gh5_open_r(gh5_file, f_id)
    character(*),intent(in) :: gh5_file
    integer(hid_t),intent(out) :: f_id

    if(.not.lgh5_init)call gh5_init()
    !< open the file.
    call h5fopen_f(gh5_file, h5f_acc_rdonly_f, f_id, gh5_err)
    return
      
    end subroutine gh5_open_r
      
      
    subroutine gh5_close(f_id)
    integer(hid_t),intent(inout) :: f_id
      
    call h5fclose_f(f_id, gh5_err)
    f_id=-1
    return
      
    end subroutine gh5_close
    

    subroutine gh5_end()
    
    if(.not.lgh5_init)return
    call h5tclose_f(dcomplex_id, gh5_err)
    call h5tclose_f(string_id, gh5_err)
    !< close fortran interface.
    call h5close_f(gh5_err)
    return
      
    end subroutine gh5_end
    

    subroutine set_dcomplex_id()
      
    gh5_size = q*2; gh5_offset = 0
    call h5tcreate_f(h5t_compound_f, gh5_size, dcomplex_id, gh5_err)
    call h5tinsert_f(dcomplex_id, "r", gh5_offset, &
            &h5kind_to_type(q, h5_real_kind), gh5_err)
    gh5_offset = gh5_offset + q
    call h5tinsert_f(dcomplex_id, "i", gh5_offset, &
            &h5kind_to_type(q, h5_real_kind), gh5_err)
    return
      
    end subroutine set_dcomplex_id
    

    subroutine gh5_write_iscalar(i, path, f_id)
    integer, intent(in) :: i
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer ibuf(1), n(1)
      
    ibuf = i; n = 1
    call gh5_write_iarray(ibuf, 1, n, path, f_id)
    return
      
    end subroutine gh5_write_iscalar
    

    subroutine gh5_read_iscalar(i, path, f_id)
    integer, intent(out) :: i
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer ibuf(1), n(1)
      
    n = 1
    call gh5_read_iarray(ibuf, 1, n, path, f_id)
    i = ibuf(1)
    return
      
    end subroutine gh5_read_iscalar
    

    subroutine gh5_write_dscalar(d, path, f_id)
    real(q), intent(in) :: d
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    real(q) dbuf(1)
    integer n(1)
      
    dbuf = d; n = 1
    call gh5_write_darray(dbuf, 1, n, path, f_id)
    return
      
    end subroutine gh5_write_dscalar
    

    subroutine gh5_read_dscalar(d, path, f_id)
    real(q), intent(out) :: d
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    real(q) dbuf(1)
    integer n(1)
      
    n = 1
    call gh5_read_darray(dbuf, 1, n, path, f_id)
    d = dbuf(1)
    return
      
    end subroutine gh5_read_dscalar
    

    subroutine gh5_write_sscalar(sarray, ns, path, f_id)
    integer, intent(in) :: ns
    character(ns), intent(in) :: sarray
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id

    character(ns) sbuf(1)
    integer n(1)

    sbuf = sarray; n = 1
    call gh5_write_sarray(sbuf, ns, 1, n, path, f_id)
    return

    end subroutine gh5_write_sscalar


    subroutine gh5_read_sscalar(sarray, ns, path, f_id)
    integer, intent(in) :: ns
    character(ns), intent(out) :: sarray
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id

    character(ns) sbuf(1)
    integer n(1)

    n = 1
    call gh5_read_sarray(sbuf, ns, 1, n, path, f_id)
    sarray=sbuf(1)
    return

    end subroutine gh5_read_sscalar


    subroutine gh5_write_1d_iarray(iarray, n1, path, f_id)
    integer, intent(in) :: n1, iarray(n1)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer n(1)
      
    n = n1
    call gh5_write_iarray(iarray, 1, n, path, f_id)
    return
      
    end subroutine gh5_write_1d_iarray
    

    subroutine gh5_read_1d_iarray(iarray, n1, path, f_id)
    integer, intent(in) :: n1
    integer, intent(out) :: iarray(n1)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer n(1)
      
    n = n1
    call gh5_read_iarray(iarray, 1, n, path, f_id)
    return
      
    end subroutine gh5_read_1d_iarray
    

    subroutine gh5_write_2d_iarray(iarray, n1, n2, path, f_id)
    integer, intent(in) :: n1, n2, iarray(n1, n2)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer n(2)
      
    n(1) = n1; n(2) = n2
    call gh5_write_iarray(iarray, 2, n, path, f_id)
    return
      
    end subroutine gh5_write_2d_iarray
    

    subroutine gh5_read_2d_iarray(iarray, n1, n2, path, f_id)
    integer, intent(in) :: n1, n2
    integer, intent(out) :: iarray(n1, n2)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer n(2)
      
    n(1) = n1; n(2) = n2
    call gh5_read_iarray(iarray, 2, n, path, f_id)
    return
      
    end subroutine gh5_read_2d_iarray
    

    subroutine gh5_read_3d_iarray(iarray, n1, n2, n3, path, f_id)
    integer, intent(in) :: n1, n2, n3
    integer, intent(out) :: iarray(n1, n2, n3)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer n(3)
      
    n(1) = n1; n(2) = n2; n(3) = n3
    call gh5_read_iarray(iarray, 3, n, path, f_id)
    return
      
    end subroutine gh5_read_3d_iarray


    subroutine gh5_write_1d_darray(darray, n1, path, f_id)
    integer, intent(in) :: n1
    real(q), intent(in) :: darray(n1)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer n(1)
      
    n = n1
    call gh5_write_darray(darray, 1, n, path, f_id)
    return
      
    end subroutine gh5_write_1d_darray
    

    subroutine gh5_read_1d_darray(darray, n1, path, f_id)
    integer, intent(in) :: n1
    real(q), intent(out) :: darray(n1)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer n(1)
      
    n = n1
    call gh5_read_darray(darray, 1, n, path, f_id)
    return
      
    end subroutine gh5_read_1d_darray
    

    subroutine gh5_write_2d_darray(darray, n1, n2, path, f_id)
    integer, intent(in) :: n1, n2
    real(q), intent(in) :: darray(n1, n2)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer n(2)
      
    n(1) = n1; n(2) = n2
    call gh5_write_darray(darray, 2, n, path, f_id)
    return
      
    end subroutine gh5_write_2d_darray
    

    subroutine gh5_read_2d_darray(darray, n1, n2, path, f_id)
    integer, intent(in) :: n1, n2
    real(q), intent(out) :: darray(n1, n2)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer n(2)
      
    n(1) = n1; n(2) = n2
    call gh5_read_darray(darray, 2, n, path, f_id)
    return
      
    end subroutine gh5_read_2d_darray
    

    subroutine gh5_write_3d_darray(darray, n1, n2, n3, path, f_id)
    integer, intent(in) :: n1, n2, n3
    real(q), intent(in) :: darray(n1, n2, n3)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer n(3)
      
    n(1) = n1; n(2) = n2; n(3) = n3
    call gh5_write_darray(darray, 3, n, path, f_id)
    return
      
    end subroutine gh5_write_3d_darray
    

    subroutine gh5_read_3d_darray(darray, n1, n2, n3, path, f_id)
    integer, intent(in) :: n1, n2, n3
    real(q), intent(out) :: darray(n1, n2, n3)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer n(3)
      
    n(1) = n1; n(2) = n2; n(3) = n3
    call gh5_read_darray(darray, 3, n, path, f_id)
    return
      
    end subroutine gh5_read_3d_darray
    

    subroutine gh5_write_4d_darray(darray, n1, n2, n3, n4, path, f_id)
    integer, intent(in) :: n1, n2, n3, n4
    real(q), intent(in) :: darray(n1, n2, n3, n4)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer n(4)
      
    n(1) = n1; n(2) = n2; n(3) = n3; n(4) = n4
    call gh5_write_darray(darray, 4, n, path, f_id)
    return
      
    end subroutine gh5_write_4d_darray
    

    subroutine gh5_read_4d_darray(darray, n1, n2, n3, n4, path, f_id)
    integer, intent(in) :: n1, n2, n3, n4
    real(q), intent(out) :: darray(n1, n2, n3, n4)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer n(4)
      
    n(1) = n1; n(2) = n2; n(3) = n3; n(4) = n4
    call gh5_read_darray(darray, 4, n, path, f_id)
    return
      
    end subroutine gh5_read_4d_darray
    

    subroutine gh5_write_5d_darray(darray, n1, n2, n3, n4, n5, path, f_id)
    integer, intent(in) :: n1, n2, n3, n4, n5
    real(q), intent(in) :: darray(n1, n2, n3, n4, n5)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer n(5)
      
    n(1) = n1; n(2) = n2; n(3) = n3; n(4) = n4; n(5) = n5
    call gh5_write_darray(darray, 5, n, path, f_id)
    return
      
    end subroutine gh5_write_5d_darray
    

    subroutine gh5_read_5d_darray(darray, n1, n2, n3, n4, n5, path, f_id)
    integer, intent(in) :: n1, n2, n3, n4, n5
    real(q), intent(out) :: darray(n1, n2, n3, n4, n5)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer n(5)
      
    n(1) = n1; n(2) = n2; n(3) = n3; n(4) = n4; n(5) = n5
    call gh5_read_darray(darray, 5, n, path, f_id)
    return
      
    end subroutine gh5_read_5d_darray
    

    subroutine gh5_write_1d_zarray(zarray, n1, path, f_id)
    integer, intent(in) :: n1
    complex(q), intent(in) :: zarray(n1)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer n(1)
      
    n = n1
    call gh5_write_zarray(zarray, 1, n, path, f_id)
    return
      
    end subroutine gh5_write_1d_zarray
    

    subroutine gh5_read_1d_zarray(zarray, n1, path, f_id)
    integer, intent(in) :: n1
    complex(q), intent(out) :: zarray(n1)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer n(1)
      
    n = n1
    call gh5_read_zarray(zarray, 1, n, path, f_id)
    return
      
    end subroutine gh5_read_1d_zarray
    

    subroutine gh5_write_1d_sarray(sarray, ns, n1, path, f_id)
    integer, intent(in) :: ns,n1
    character(ns), intent(in) :: sarray(n1)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer n(1)
      
    n = n1
    call gh5_write_sarray(sarray, ns, 1, n, path, f_id)
    return
      
    end subroutine gh5_write_1d_sarray
    

    subroutine gh5_read_1d_sarray(sarray, ns, n1, path, f_id)
    integer, intent(in) :: ns,n1
    character(ns), intent(out) :: sarray(n1)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer n(1)
      
    n = n1
    call gh5_read_sarray(sarray, ns, 1, n, path, f_id)
    return
      
    end subroutine gh5_read_1d_sarray


    subroutine gh5_write_2d_zarray(zarray, n1, n2, path, f_id)
    integer, intent(in) :: n1, n2
    complex(q), intent(in) :: zarray(n1, n2)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer n(2)
      
    n(1) = n1; n(2) = n2
    call gh5_write_zarray(zarray, 2, n, path, f_id)
    return
      
    end subroutine gh5_write_2d_zarray
    

    subroutine gh5_read_2d_zarray(zarray, n1, n2, path, f_id)
    integer, intent(in) :: n1, n2
    complex(q), intent(out) :: zarray(n1, n2)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer n(2)
      
    n(1) = n1; n(2) = n2
    call gh5_read_zarray(zarray, 2, n, path, f_id)
    return
      
    end subroutine gh5_read_2d_zarray
    

    subroutine gh5_write_3d_zarray(zarray, n1, n2, n3, path, f_id)
    integer, intent(in) :: n1, n2, n3
    complex(q), intent(in) :: zarray(n1, n2, n3)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer n(3)
      
    n(1) = n1; n(2) = n2; n(3) = n3
    call gh5_write_zarray(zarray, 3, n, path, f_id)
    return
      
    end subroutine gh5_write_3d_zarray
    

    subroutine gh5_read_3d_zarray(zarray, n1, n2, n3, path, f_id)
    integer, intent(in) :: n1, n2, n3
    complex(q), intent(out) :: zarray(n1, n2, n3)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer n(3)
      
    n(1) = n1; n(2) = n2; n(3) = n3
    call gh5_read_zarray(zarray, 3, n, path, f_id)
    return
      
    end subroutine gh5_read_3d_zarray
    

    subroutine gh5_write_4d_zarray(zarray, n1, n2, n3, n4, path, f_id)
    integer, intent(in) :: n1, n2, n3, n4
    complex(q), intent(in) :: zarray(n1, n2, n3, n4)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer n(4)
      
    n(1) = n1; n(2) = n2; n(3) = n3; n(4) = n4
    call gh5_write_zarray(zarray, 4, n, path, f_id)
    return
      
    end subroutine gh5_write_4d_zarray
    

    subroutine gh5_read_4d_zarray(zarray, n1, n2, n3, n4, path, f_id)
    integer, intent(in) :: n1, n2, n3, n4
    complex(q), intent(out) :: zarray(n1, n2, n3, n4)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer n(4)
      
    n(1) = n1; n(2) = n2; n(3) = n3; n(4) = n4
    call gh5_read_zarray(zarray, 4, n, path, f_id)
    return
      
    end subroutine gh5_read_4d_zarray
    

    subroutine gh5_write_5d_zarray(zarray, n1, n2, n3, n4, n5, path, f_id)
    integer, intent(in) :: n1, n2, n3, n4, n5
    complex(q), intent(in) :: zarray(n1, n2, n3, n4, n5)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer n(5)
      
    n(1) = n1; n(2) = n2; n(3) = n3; n(4) = n4; n(5) = n5
    call gh5_write_zarray(zarray, 5, n, path, f_id)
    return
      
    end subroutine gh5_write_5d_zarray
    

    subroutine gh5_read_5d_zarray(zarray, n1, n2, n3, n4, n5, path, f_id)
    integer, intent(in) :: n1, n2, n3, n4, n5
    complex(q), intent(out) :: zarray(n1, n2, n3, n4, n5)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    integer n(5)
      
    n(1) = n1; n(2) = n2; n(3) = n3; n(4) = n4; n(5) = n5
    call gh5_read_zarray(zarray, 5, n, path, f_id)
    return
      
    end subroutine gh5_read_5d_zarray
    

    subroutine gh5_write_iarray(iarray, nd, n, path, f_id)
    integer, intent(in) :: nd, n(nd), iarray(*)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    gh5_dims(1:nd) = n
    ! create the dataspace.
    call h5screate_simple_f(nd, gh5_dims(1:nd), dspace_id, gh5_err)
    ! create the dataset with default properties.
    call h5dcreate_f(f_id, path, h5t_native_integer, dspace_id, &
            &dset_id, gh5_err)
    call h5dwrite_f(dset_id, h5t_native_integer, iarray, gh5_dims(1:nd), &
            &gh5_err)
    ! end access to the dataset and release resources used by it.
    call h5dclose_f(dset_id, gh5_err)
    ! terminate access to the data space.
    call h5sclose_f(dspace_id, gh5_err)
    return
      
    end subroutine gh5_write_iarray
    

    subroutine gh5_read_iarray(iarray, nd, n, path, f_id)
    integer, intent(in) :: nd, n(nd)
    integer, intent(out) :: iarray(*)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    gh5_dims(1:nd) = n
    ! open path using the default propertie
    call h5dopen_f(f_id, path, dset_id, gh5_err)
    ! read the data using the default properties.
    call h5dread_f(dset_id, h5t_native_integer, iarray, gh5_dims(1:nd), &
            &gh5_err)
    ! end access to the dataset and release resources used by it.
    call h5dclose_f(dset_id, gh5_err)
    return
      
    end subroutine gh5_read_iarray
    

    subroutine gh5_write_darray(darray, nd, n, path, f_id)
    integer, intent(in) :: nd, n(nd)
    real(q), intent(in) :: darray(*)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    gh5_dims(1:nd) = n
    ! create the dataspace.
    call h5screate_simple_f(nd, gh5_dims(1:nd), dspace_id, gh5_err)
    ! create the dataset with default properties.
    call h5dcreate_f(f_id, path, h5t_native_double, dspace_id, dset_id, &
            & gh5_err)
    call h5dwrite_f(dset_id, h5t_native_double, darray, gh5_dims(1:nd), &
            & gh5_err)
    ! end access to the dataset and release resources used by it.
    call h5dclose_f(dset_id, gh5_err)
    ! terminate access to the data space.
    call h5sclose_f(dspace_id, gh5_err)
    return
      
    end subroutine gh5_write_darray
    

    subroutine gh5_read_darray(darray, nd, n, path, f_id)
    integer, intent(in) :: nd, n(nd)
    real(q), intent(out) :: darray(*)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    gh5_dims(1:nd) = n
    ! open path using the default propertie
    call h5dopen_f(f_id, path, dset_id, gh5_err)
    ! read the data using the default properties.
    call h5dread_f(dset_id, h5t_native_double, darray, gh5_dims(1:nd), &
            &gh5_err)
    ! end access to the dataset and release resources used by it.
    call h5dclose_f(dset_id, gh5_err)
    return
      
    end subroutine gh5_read_darray
    

    subroutine gh5_write_zarray(zarray, nd, n, path, f_id)
    integer, intent(in) :: nd, n(nd)
    complex(q), target, intent(in) :: zarray(*)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    gh5_dims(1:nd) = n
    ! create the dataspace.
    call h5screate_simple_f(nd, gh5_dims(1:nd), dspace_id, gh5_err)
    ! create the dataset with default properties.
    call h5dcreate_f(f_id, path, dcomplex_id, dspace_id, dset_id, gh5_err)
    f_ptr = c_loc(zarray)
    call h5dwrite_f(dset_id, dcomplex_id, f_ptr, gh5_err)
    ! end access to the dataset and release resources used by it.
    call h5dclose_f(dset_id, gh5_err)
    ! terminate access to the data space.
    call h5sclose_f(dspace_id, gh5_err)
    return
      
    end subroutine gh5_write_zarray
    

    subroutine gh5_read_zarray(zarray, nd, n, path, f_id)
    integer, intent(in) :: nd, n(nd)
    complex(q), target, intent(out) :: zarray(*)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    gh5_dims(1:nd) = n
    ! open path using the default propertie
    call h5dopen_f(f_id, path, dset_id, gh5_err)
    f_ptr = c_loc(zarray)
    call h5dread_f(dset_id, dcomplex_id, f_ptr, gh5_err)
    ! end access to the dataset and release resources used by it.
    call h5dclose_f(dset_id, gh5_err)
    return
      
    end subroutine gh5_read_zarray
    

    subroutine gh5_write_sarray(sarray, ns, nd, n, path, f_id)
    integer, intent(in) :: ns, nd, n(nd)
    character(len=ns), target, intent(in) :: sarray(*)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    gh5_size = ns
    call h5tset_size_f(string_id, gh5_size, gh5_err)
    gh5_dims(1:nd) = n
    ! create the dataspace.
    call h5screate_simple_f(nd, gh5_dims(1:nd), dspace_id, gh5_err)
    ! create the dataset with default properties.
    call h5dcreate_f(f_id, path, string_id, dspace_id, dset_id, gh5_err)
    f_ptr = c_loc(sarray(1)(1:1))
    call h5dwrite_f(dset_id, string_id, f_ptr, gh5_err)
    ! end access to the dataset and release resources used by it.
    call h5dclose_f(dset_id, gh5_err)
    ! terminate access to the data space.
    call h5sclose_f(dspace_id, gh5_err)
    return
      
    end subroutine gh5_write_sarray


    subroutine gh5_read_sarray(sarray, ns, nd, n, path, f_id)
    integer, intent(in) :: ns, nd, n(nd)
    character(len=ns), target, intent(out) :: sarray(*)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    gh5_size = ns
    call h5tset_size_f(string_id, gh5_size, gh5_err)
    gh5_dims(1:nd) = n
    ! open path using the default propertie
    call h5dopen_f(f_id, path, dset_id, gh5_err)
    f_ptr = c_loc(sarray(1)(1:1))
    call h5dread_f(dset_id, string_id, f_ptr, gh5_err)
    ! end access to the dataset and release resources used by it.
    call h5dclose_f(dset_id, gh5_err)
    return
      
    end subroutine gh5_read_sarray


    subroutine gh5_create_group(path, f_id)
    character(len=*), intent(in) :: path
    integer(hid_t), intent(in) :: f_id
      
    call h5gcreate_f(f_id, path, group_id, gh5_err)
    call h5gclose_f(group_id, gh5_err)
    return
      
    end subroutine gh5_create_group
      
      
end module ghdf5_base
