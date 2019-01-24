      module comwann_mod

      integer ::
     &  irel,
     &  iun,
     &  limlb,
     &  maxia,
     &  maxib,
     &  maxic,
     &  magn,
     &  maxb,
     &  maxel,
     &  maxlfun,
     &  maxnrad,
     &  maxntle,
     &  maxpb,
     &  maxplw,
     &  maxwig,
     &  maxwig_j,
     &  maxwf,
     &  maxylm,
     &  master,
     &  me,
     &  mpierr,
     &  natom,
     &  nbasmpw,
     &  nbndf,
     &  nbndtm,
     &  ndiv(3),
     &  nfun,
     &  ngroup,
     &  nplwbas,
     &  nproc,
     &  nplwdif,
     &  npnt,
     &  nqdiv,
     &  nrdiv(3),
     &  nrel,
     &  nsort,
     &  nspin,
     &  nspin_0,
     &  numplw,
     &  n_tau,
     &  ubi_num,
     $  ind_low,
     $  ind_high

      integer, allocatable ::
     &  gbs_sym(:,:),
     &  i_kref(:),
     &  igbs(:,:),
     &  indgb(:,:),
     &  indplw(:,:,:),
     &  ip(:,:),
     &  iplf_bk(:,:),
     &  k_group(:),
     &  ind_wf(:,:),
     &  indxcg(:),
     &  indxcgr(:,:,:,:),
     &  indbas(:,:,:),
     &  indbasa(:,:,:,:),
     &  io_lem(:),
     &  lf_isz(:,:),
     &  li_val(:,:),
     &  is(:),
     &  jcg(:),
     &  jcgr(:),
     &  lfun(:),
     &  lfunm(:),
     &  lmb(:),
     &  nbask(:),
     &  n_mpi_kk(:),
     &  n_mpi_k(:),
     &  n_bnd(:,:),
     &  nrad(:),
     &  ntle(:,:),
     &  ndim_k(:),
     &  ndim_kk(:)


      double precision ::
     &  amega,
     &  b_a,
     &  c_a,
     &  clight,
     &  chem_pot,
     &  chem_pot_qp,
     &  evolt,
     &  gbas(3,3),
     &  par,
     &  pi,
     &  rbas(3,3),
     &  rydberg,
     &  bohr,
     $  cut_froz,
     $  cut_low,
     $  cut_total

      double precision, allocatable ::
     &  bett(:),
     &  cg(:),
     &  dr(:,:),
     &  e_bnd(:,:,:),
     &  e_qp(:,:,:),
     &  ffsmt(:,:,:,:,:,:,:),
     &  gbs(:,:),
     &  gfun(:,:),
     &  gfund(:,:),
     &  h(:),
     &  pnt(:,:),
     &  r(:,:),
     &  shift(:,:),
     &  sovr(:),
     &  sovi(:),
     &  tau(:,:),
     &  tshift(:,:,:),
     &  u(:,:),
     &  z(:)

      logical ::
     &  complex_ro,
     &  goparr,
     &  maswrk

      character*3 :: 
     &  ubi_0
     
      character*72 :: 
     &  allfile

      character*1000 :: 
     & rstpath      

      character*3, allocatable :: 
     $  augm(:,:,:)

      character*4, allocatable :: 
     $  txtel(:)

      complex*16, allocatable ::
     $  cgr(:),
     &  ev_bnd(:,:,:,:),
     &  ev_qp(:,:,:,:),
     &  uj(:,:),
     &  z_bnd(:,:,:,:),
     &  z_qp(:,:,:,:)

      parameter (pi=3.141592653589793238d0)
      parameter (bohr=0.52917720859d0)
      parameter (rydberg=13.605691752d0)
      parameter (evolt=27.2107d0)

      character*300 ::
     $ seed_name
      character*20, allocatable ::
     $ atom_symbols(:)

      integer ::
     $  nntot,
     $  num_bands,
     $  num_wann,
     $  num_nnmax,
     $  nrpts,
     $  nexband,
     $  max_nalimj,
     $  num_alimj,
!     $  num_kpts,
!     $  mp_grid(3),
     $  writewan,
!     $  write_wan_pdos,
!     $  cal_photoemission_input,
!     $  cal_optics_input,
     $  nbnd_in_froz,
     $  num_kpts_fi,
     $  mp_grid_fi(3),
     $  num_iter,
     $  dis_num_iter      
!     $  ndim_mp_grid_fi_max



      logical ::
     $  gamma_only


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
     $  atom_wan(:),
!     $  ldim_wan(:),
!     $  rdim_wan(:,:),
!     $  lval_wan(:,:),
     $  ndegen(:),
     $  irvec(:,:),
     $  alimj(:,:,:),
     $  nalimj(:),
     $  ndim_mp_grid_fi(:)
!     $  n_mpi_mp_grid_fi(:)


      double precision ::
     $  spread(3),
     $  dis_win_min,
     $  dis_win_max,
     $  dis_froz_min,
     $  dis_froz_max,
     $  photoemission_emin,
     $  photoemission_emax


      double precision,allocatable ::
     $  kpt_latt(:,:),
!     $  kpt_latt_fi(:,:),
     $  proj_site(:,:),
     $  proj_x(:,:),
     $  proj_z(:,:),
     $  proj_zona(:),
     $  eigenvalues(:,:),
     $  eigenvalues_fi(:,:),
     $  wann_centers(:,:),
     $  wann_spreads(:),
     $  proj_s_qaxisx(:,:),
     $  atom_wan_distance(:,:),
     $  alimj_xaxis(:,:,:),
     $  alimj_zaxis(:,:,:),
c$$$     $  weight_in_froz_win(:,:,:),
     $  atom_xaxis(:,:),
     $  atom_zaxis(:,:)
!     $  wb(:,:),
!     $  bk(:,:,:)


      complex*16,allocatable ::
     $  m_matrix(:,:,:,:),
     $  a_matrix(:,:,:),
     $  u_matrix(:,:,:),
     $  w_matrix_fi(:,:,:),
     $  u_matrix_opt(:,:,:),
     $  wan_g_overlap(:,:),
     $  v_matrix(:,:,:),
     $  z_wan_bnd(:,:,:),
     $  ev_wan_bnd(:,:,:),
!     $  radial_coeff_wan(:,:,:),
     $  alimj_coeff(:,:),




     $  lte_wan(:,:,:),
     $  gks_wan(:,:,:),
     $  lte_bnd(:,:,:),
     $  gks_bnd(:,:,:),
     $  momentum(:,:,:,:),
     $  momentum_wan(:,:,:,:)

      logical, allocatable ::
     $  lwindow(:,:)
      end



