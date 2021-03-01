      module comlowh_mod


      character*300 :: seed_name

      character*3 :: den_mode      

      integer ::
     $  cal_mode,
     $  num_bands,
     $  num_wann,
     $  num_wancell,
     $  num_proj,
     $  num_atom,
     $  num_orb_max,
     $  num_rest,
     $  nrpts,
     $  ndiv(3),
     $  ndiv_fi(3),
     $  nqdiv,
     $  nqdiv_fi,
     $  n_omega,
     $  ndim_kk_fi_max,
     $  ndim_kk_co_max,      
     $  is_matsubara,
     $  is_recal_ef,
     $  n_i,
     $  n_ioac,      
     $  n_iac,
     $  n_iac_nm,
     $  n_iac_mat,
     $  n_iac_mat_i,
     $  n_iac_mat_f,
     $  ndeg_ioac_max,
     $  iun,
     $  num_cell,
     $  is_magnetic,
     $  is_kpath,
     $  is_spinorbit,
     $  is_metal,
     $  num_alimj,
     $  max_nalimj,
     $  is_qp


      integer, allocatable ::
     $  iiiio(:,:),
     $  num_orb(:),
     $  norb_iac(:),
     $  include_bands(:),
     $  ndegen(:),
     $  irvec(:,:),
     $  iac(:,:),
     $  ndeg_iac(:),
     $  ndeg_ioac(:,:),
     $  ndeg_itot(:,:),            
     $  wan_correlated(:,:),
     $  wan_cor_index(:,:,:),      
     $  ndim_kk_fi(:),
     $  ndim_kk_co(:),      
     $  n_mpi_kk_fi(:),
     $  n_mpi_kk_co(:),      
     $  ioac(:,:,:),
     $  iio_diagonal(:,:),
     $  nalimj(:),
     $  alimj(:,:,:),
     $  bnd_in_win(:,:),
     $  nbnd_in_win(:),
     $  bnd_in_win_co(:,:),
     $  nbnd_in_win_co(:),
     $  to_orig_bnd(:,:)


      double precision :: 
     $  ewin_min,
     $  ewin_max,
     $  gamma_selfenergy,
     $  gamma_real,
     $  ef,
     $  dc_shift,
     $  real_lattice(3,3),
     $  recip_lattice(3,3),
     $  recip_lattice_conv(3,3),      
     $  beta,
     $  delta_charge,
     $  nelectron,
     $  temperature


      


      double precision,allocatable ::
     $  kpt_latt(:,:),
     $  kpt_latt_fi(:,:),
     $  eigenvalues(:,:),
     $  eigenvalues_co(:,:),
     $  eigenvalues_co_rest(:,:),      
     $  eigenvalues_co_all(:,:),            
     $  eigenvalues_fi(:,:),
     $  frequency_dble(:),
     $  lpoint(:,:),
     $  alimj_xaxis(:,:,:),
     $  alimj_zaxis(:,:,:)



      complex*16,allocatable :: 
     $  vw_matrix_co(:,:,:),      
     $  v_matrix(:,:,:),
     $  w_matrix_co(:,:,:),            
     $  w_matrix_fi(:,:,:),
     $  frequency(:),
     $  e_imp(:,:),
     $  sig(:,:,:),
     $  dc(:,:),
     $  zinv_m1(:,:),      
     $  delta_mat(:,:,:,:),      
     $  projected_g(:,:,:,:,:),
     $  h_k_wan_fi(:,:,:),
     $  proj(:,:,:,:,:),
     $  proj_renorm(:,:,:,:,:),
     $  proj_co(:,:,:,:,:),      
     $  proj_renorm_co(:,:,:,:,:),            
     $  gloc(:,:,:,:),
     $  g_lat_inv(:,:,:,:),
     $  trans_basis(:,:,:),            
     $  e_projected_mat(:,:,:),
     $  alimj_coeff(:,:),
     $  h_nonint(:,:,:),
     $  h_nonint_co(:,:,:),
     $  n_matrix(:,:,:)

      integer :: mpierr,me, me1,master,nproc,irank,isize
      logical :: maswrk,zinv_m1_exist,is_trans_basis


      end
      
      
      




      


      
