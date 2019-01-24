      module wannier90_mod

      character*300 :: seed_name
      character*20 :: checkpoint
      character*20, allocatable :: atom_symbols(:)

      integer ::
     $  nntot,
     $  num_bands,
     $  num_wann,
     $  num_nnmax,
     $  nrpts,
     $  nexband

      logical :: gamma_only,
     $  have_disentangled

      integer, allocatable ::
     $  nnlist(:,:),
     $  nncell(:,:,:),
     $  proj_l(:),
     $  proj_m(:),
     $  proj_radial(:),
     $  exclude_bands(:),
     $  include_bands(:),
     $  exclude_bands_table(:),
     $  proj_s(:),
     $  atom_wan_initial(:),
     $  atom_wan(:),
     $  ldim_wan(:),
     $  rdim_wan(:,:),
     $  lval_wan(:,:),
     $  ndegen(:),
     $  irvec(:,:),
     $  ndimwin(:),
     $  alimj(:,:),
     $  wan_correlated(:)

      
      double precision :: 
     $  spread(3),
     $  dis_win_min,
     $  dis_win_max,
     $  omega_invariant


      double precision,allocatable ::
     $  kpt_latt(:,:),
     $  proj_site(:,:),
     $  proj_x(:,:),
     $  proj_z(:,:),
     $  proj_zona(:),
     $  eigenvalues(:,:),
     $  wann_centers(:,:),
     $  wann_spreads(:),
     $  proj_s_qaxisx(:,:),
     $  atom_wan_distance(:,:)



      complex*16,allocatable :: 
     $  m_matrix(:,:,:,:),
     $  a_matrix(:,:,:),
     $  u_matrix(:,:,:),
     $  u_matrix_opt(:,:,:),
c$$$     $  wan_g_overlap(:,:,:,:),
     $  wan_g_overlap(:,:),
     $  v_matrix(:,:,:),
     $  z_wan_bnd(:,:,:),
     $  ev_wan_bnd(:,:,:),
     $  radial_coeff_wan(:,:,:),
     $  hfmat(:,:,:)

      
      logical, allocatable :: 
     $  lwindow(:,:)
      end
