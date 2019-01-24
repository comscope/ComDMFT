      subroutine output_boson_q_nu(bos,ch)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      double precision, intent(in) :: bos(0:n_nu,n_symb)
      character*4, intent(in) :: ch
      integer :: k,i_len,i_nu,m
      m=min(n_symb,5)
      if(maswrk) then
        i_len=len_trim(allfile)
        if(ubi=='dft') open(3,file=allfile(1:i_len)//'_NU_Qsym_DFT'//ch)
        if(ubi=='dif') open(3,file=allfile(1:i_len)//'_NU_Qsym_DIF'//ch)
        if(ubi==' hf') open(3,file=allfile(1:i_len)//'_NU_Qsym_HF'//ch)
        if(ubi==' gw') open(3,file=allfile(1:i_len)//'_NU_Qsym_GW'//ch)
        if(ubi==' qp') open(3,file=allfile(1:i_len)//'_NU_Qsym_QP'//ch)
        if(ubi=='gwg') open(3,file=allfile(1:i_len)//'_NU_Qsym_GWG'//ch)
        if(ubi=='dmf') open(3,file=allfile(1:i_len)//'_NU_Qsym_DMF'//ch)
        do i_nu=0,n_nu
          write(3,'(f11.5,5(1x,e12.5))')w_nu(i_nu)*evolt/2,
     &      (bos(i_nu,k),k=1,m)
        enddo
        close(3)
      endif
      end
