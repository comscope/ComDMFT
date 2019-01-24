program exespci
    use gspci
    use gtime
    implicit none
    
    call set_time_point(1,1)
    call gspci_gh5exe()
    call set_time_point(2,1)
    call print_time_usage('exe_spci',1,0)

end program exespci
