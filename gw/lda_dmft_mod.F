      module lda_dmft_mod
      logical :: is_lda_dmft
      integer :: n_dmft_bnd_max,n_dmft_wan_max
      real*8 :: dmft_temperature
      integer, allocatable :: dmft_bnd_index(:,:),
     $  n_dmft_bnd(:),n_dmft_wan(:)
c$$$  $  ,is_dmft_bnd_included(:)
      real*8, allocatable ::
     $  gfun_dmft_0(:,:),gfund_dmft_0(:,:)
      complex*16, allocatable ::
     $  disentangle_mat(:,:,:),
     $  dmft_den_mat(:,:,:),
     $  z_bnd_dmft_0(:,:,:,:),
     $  ev_bnd_dmft_0(:,:,:,:)
      end
      
