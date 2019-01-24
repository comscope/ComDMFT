      subroutine init_vh
************************************************************
*   The initiation of different arrays for smultro         *
************************************************************
	use atom_mod
	use manager_mod
	use solid_mod
      implicit none
      allocate(indfgtl(limlpb,nsort))
      call fgtlind
      allocate(fgtlr(maxfgtl,nsort,nplwrop))
      allocate(fgtli(maxfgtl,nsort,nplwrop))
      allocate(agtlr(maxfgtl,nsort,nplwrop))
      allocate(agtli(maxfgtl,nsort,nplwrop))
      allocate(bgtlr(nplwrop,maxfgtl,nsort))
      allocate(bgtli(nplwrop,maxfgtl,nsort))
      allocate(cgtlr(nplwrop,maxfgtl,nsort))
      allocate(cgtli(nplwrop,maxfgtl,nsort))
      call setfgtl
      call timel('***** setfgtl finished *************')
      end