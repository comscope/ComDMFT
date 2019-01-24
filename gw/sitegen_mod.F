      module sitegen_mod
      public
	integer, allocatable :: lastg_site(:),
     &	                        iref_vgen(:)
      real*8, allocatable :: vect_site(:,:),
     &	                     rco_site(:),
     &                       tauc_site(:,:),
     &                       vec_vgen(:,:),
     &                       rw_vgen(:)
      public lastg_site
      end
