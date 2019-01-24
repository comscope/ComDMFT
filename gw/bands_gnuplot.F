      subroutine bands_gnuplot(nbn,nbn0,epl)
	use manager_mod
	use parallel_mod
	use solid_mod
	use units_mod
      implicit none
      integer, intent(in) :: nbn,nbn0
      real*8, intent(in) :: epl(nbn0,num_k_plot,nspin)
	integer :: i_len,ispin,k,i
	character*72 :: f_bnd
	character*3, allocatable :: glabel(:)
	real*8, allocatable :: e(:)
	allocate(e(2:n_k_div))
	if(.not.maswrk) return
	i_len=len_trim(allfile)
	if(ubi=='dft') then
	  f_bnd=allfile(1:i_len)//'_dft_band_LAPW.dat'
	else if(ubi==' hf') then
	  f_bnd=allfile(1:i_len)//'_hf_band.dat'
	else if(ubi==' qp') then
	  f_bnd=allfile(1:i_len)//'_qp_band.dat'
	else if(ubi==' gw') then
	  f_bnd=allfile(1:i_len)//'_gw_band.dat'
	else if(ubi=='gwg') then
	  f_bnd=allfile(1:i_len)//'_gwg_band.dat'
	endif
	open(3,file=f_bnd)
c ----------------------------------------------------------
      do ispin=1,nspin
        do i=1,nbn
          do k=1,num_k_plot
            write(3,'(2E16.8)')k_dist(k),
     &			               (epl(i,k,ispin)-chem_pot)*evolt/2
	    enddo
          write(3,*) ' '
        enddo
      enddo
      close(3)
c --------------------------------------------------------------
	if(ubi=='dft') then
	  open(3,file=allfile(1:i_len)//'_dft_band_LAPW.gnu')
	else if(ubi==' hf') then
	  open(3,file=allfile(1:i_len)//'_hf_band.gnu')
	else if(ubi==' qp') then
	  open(3,file=allfile(1:i_len)//'_qp_band.gnu')
	else if(ubi==' gw') then
	  open(3,file=allfile(1:i_len)//'_gw_band.gnu')
	else if(ubi=='gwg') then
	  open(3,file=allfile(1:i_len)//'gwg_band.gnu')
	endif
c ----------------------------------------------------------------
	allocate(glabel(num_k_path+1))
      glabel(1)=' '//k_label(1,1)//' '
      do i=2,num_k_path
        if(k_label(2,i-1)/=k_label(1,i)) then
          glabel(i)=k_label(2,i-1)//'/'//k_label(1,i)
        else
          glabel(i)=' '//k_label(2,i-1)//' '
        end if
      end do
      glabel(num_k_path+1)=' '//k_label(2,num_k_path)//' '
c ---------------------------------------------------------------------
      write(3,*)
     & 'set loadpath "C:/P4Win/Gnuplot/gnuplot/bin/share/PostScript/"'
      write(3,*)'set terminal postscript "fontsize" 26'
	if(ubi=='dft') then
        write(3,*)'set output "'//allfile(1:i_len)//
     &                       '_lda_bnd_LAPW.eps"'
	else if(ubi==' hf') then
        write(3,*)'set output "'//allfile(1:i_len)//'_hf_bnd.eps"'
	else if(ubi==' qp') then
        write(3,*)'set output "'//allfile(1:i_len)//'_qp_bnd.eps"'
	else if(ubi==' gw') then
        write(3,*)'set output "'//allfile(1:i_len)//'_gw_bnd.eps"'
	else if(ubi=='gwg') then
        write(3,*)'set output "'//allfile(1:i_len)//'_gwg_bnd.eps"'
      endif
      write(3,701)k_dist(num_k_plot),emindos*evolt/2,emaxdos*evolt/2
      do i = 1, num_k_path-1
        write(3,705)sum(kpath_len(1:i)),emindos*evolt/2,
     &	            sum(kpath_len(1:i)),emaxdos*evolt/2
      enddo
      write(3,702, advance="no")glabel(1),0.d0,
     &	            (glabel(i+1),sum(kpath_len(1:i)),i=1,num_k_path-1)
      write(3,703)glabel(num_k_path+1),sum(kpath_len(:))
      write(3,*)'plot ','"'//trim(f_bnd)//'" with lines lw 3 lt 1 lc 1'
      write(3,*)'unset output'
      write(3,*)'unset terminal'
      write(3,*)'unset loadpath'   
      close(3)
      deallocate(glabel,e)
701   format('set data style lines',/,'set nokey',/,
     &       'set xrange [0:',F8.5,']',/,
     &       'set yrange [',F10.5,' :',F10.5,']')
702   format('set xtics (',:20('"',A3,'" ',F8.5,','))
703   format(A3,'" ',F8.5,')')
705   format('set arrow from ',F8.5,',',F10.5,' to ',F8.5,',',F10.5,
     &       ' nohead')
      end 
