      subroutine properties(iexit1)
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      implicit none
      integer, intent(inout) :: iexit1
c     --- Numbers of Occupation ----------------------------------------
      call occ_number
c     --- Output GW Energy ---------------------------------------------
      call etot_gw_2
c     --- Density of States (DOS) --------------------------------------
      if (dosplot) then
        if(iexit1==1.or.ubi/='dft') call output_dos
      endif
c     --- Bands plotting -----------------------------------------------
      if (bandplot) then
        if(iexit1==1) then
          if(ubi=='dft') call plot_bands_lapw   !! LAPW bands
        endif
      endif
      end
