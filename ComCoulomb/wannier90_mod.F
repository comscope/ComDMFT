      module wannier90_mod

      character*300 :: seed_name
      character*20 :: checkpoint

      character*1000 :: 
     &  rstpath,
     &  wanpath

      integer ::
     $  num_alimj,
     $  max_nalimj,
     $  nntot,
     $  num_bands,
     $  num_wann,
     $  nrpts,
     $  num_orb_cor,
     $  num_orb_iter,
     $  num_orb_rest,
     $  num_orb_proj,
     $  pol_iter_for_u,
     $  natom_ciw,
     $  n_ciw_max,
     $  nwpb_max,
     $  nwpb_all,
     $  num_kpts


      
      logical ::
     $  have_disentangled

      
      integer, allocatable ::
     $  mt2wan(:,:,:),
     &  corr_iter_rel(:,:),
     $  alimj(:,:,:),
     $  nalimj(:),
     $  exclude_bands(:),
     $  include_bands(:),
     $  ndegen(:),
     $  irvec(:,:),
     $  ndimwin(:),
     $  wan_correlated(:),
     $  wan_itinerant(:),
     $  mt2ciw(:,:,:),
     $  mt2ciw_rel_l(:,:,:),      
     &  corr_iter(:,:),
     $  at_bonding(:,:),
     $  at_bonding_rel(:,:),
     $  ind_bands_foru(:,:,:),
     $  include_bands_foru(:,:,:,:,:),
     $  n_ciw_atom(:),
     $  iatom_ciw(:),
     $  ciw_ind(:,:),
     $  at_orb_ciw(:,:),
     $  nwpb(:)

      
      double precision :: 
     $  omega_invariant



      double precision,allocatable ::
     $  alimj_xaxis(:,:,:),
     $  alimj_zaxis(:,:,:),
     $  kpt_latt(:,:),
     $  eigenvalues(:,:),
     $  wann_centers(:,:),
     $  wann_spreads(:),
     $  eigenvalues_cor(:,:,:),
     $  eigenvalues_iter(:,:,:),
     $  eigenvalues_rest(:,:,:)


      complex*16,allocatable :: 
     $  alimj_coeff(:,:),
     $  m_matrix(:,:,:,:),
     $  u_matrix(:,:,:),
     $  u_matrix_opt(:,:,:),
     $  v_matrix(:,:,:),
     &  m_pi(:,:,:,:,:),
     &  pi_ww(:,:,:,:,:),
     $  v_matrix_cor(:,:,:,:),
     $  v_matrix_iter(:,:,:,:),      
     $  v_matrix_rest(:,:,:,:),
     $  prod_wan2(:,:,:,:,:,:),


     $  w_k_den(:,:,:,:,:,:,:),
     $  w_loc(:,:,:,:,:,:,:),      
     $  epsinv_k_den(:,:,:,:,:,:,:),
     $  epsinv_loc(:,:,:,:,:,:,:),      
     $  eps_k_den(:,:,:,:,:,:,:),
     $  eps_loc(:,:,:,:,:,:,:),      
     $  u_k_den(:,:,:,:,:,:,:),
     $  v_k_den(:,:,:,:,:,:),      
     $  w_wpb_k(:,:,:,:,:,:,:),
     $  w_wpb_local(:,:,:,:,:),
     $  chi0_wpb_k(:,:,:,:,:,:,:),
     $  chi0_wpb_local(:,:,:,:,:),                  
     $  u_loc(:,:,:,:,:,:,:),
     $  v_loc(:,:,:,:,:,:),      
     $  epsinv_high_k_den(:,:,:,:,:,:,:),
     $  epsinv_high_loc(:,:,:,:,:,:,:),
     $  eps_high_k_den(:,:,:,:,:,:,:),
     $  eps_high_loc(:,:,:,:,:,:,:),   
     $  chi0_low_k_den(:,:,:,:,:,:,:),
     $  chi0_low_loc(:,:,:,:,:,:,:),   
     $  chi0_high_k_den(:,:,:,:,:,:,:),
     $  chi0_high_loc(:,:,:,:,:,:,:),      
     $  chi0_k_den(:,:,:,:,:,:,:),
     $  chi0_loc(:,:,:,:,:,:,:),      

     $  chi_high_k_den(:,:,:,:,:,:,:),
     $  chi_high_loc(:,:,:,:,:,:,:),      
     $  chi_k_den(:,:,:,:,:,:,:),
     $  chi_loc(:,:,:,:,:,:,:),      

     $  chi_high_k_den_tau(:,:,:,:,:,:,:),
     $  chi_k_den_tau(:,:,:,:,:,:,:),
     $  wc_k_den_tau(:,:,:,:,:,:,:),
     $  proj_renorm(:,:,:,:)
      
      logical, allocatable ::
     $  natom_cor(:),
     $  natom_iter(:),
     $  lwindow(:,:)
      end
      
      
      




      


      
