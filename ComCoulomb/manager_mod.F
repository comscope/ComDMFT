      module manager_mod
      logical :: complex_ro,
     ^  el_gas,
     &  we_need_u,
     $  u_linear,
     $  p_local_k,
     $  p_local_r      
      character*3 :: ubi,
     &  v_repr
      character*4 :: u_mode
      character*40 :: mem_where

      character*8 :: 
     &  obietto
      character*72 :: allfile
      integer :: 
     $  interp_nu_d,
     S  dosplotflag,
     &  iun,
     &  irel,
     &  itest,
     &  magn,
     &  maxylm,
     &  n_nu,
     &  n_rfreq,
     &  n_nu_asy,
     &  n_nu_d,
     &  n_nu_exa,
     &  n_omega,
     &  n_omega_asy,
     &  n_omega_d,
     &  n_omega_exa,
     &  n_tau,
     &  ndiv_r1,
     &  nrel,
     &  nrel_core,
     &  nspin,
     &  nspin_0,
     &  nspin_1,
     &  ubi_num,
     &  p_low_mode
      double precision ::
     &  dist_rc,
     &  betta_t,
     &  chem_pot,
     &  chem_pot_qp,
     &  eps_pb,
     &  eps_v_coul,
     &  projemin,
     &  projemax,
     &  enpj_emin,
     &  enpj_emax,
     &  r_ref(3),
     &  mem_max,
     &  mem_using,
     &  r1_begin(3),
     &  r1_end(3),
     &  maxrealfreq
      integer, allocatable :: 
     &  ind_nu_int(:),
     &  indxcg(:),
     &  indxcgr(:,:,:,:),
     &  jcg(:),
     &  jcgr(:)
      double precision, allocatable :: ab_nu_from_tau(:,:,:),
     &  bett(:),
     &  et_sum_nu(:),
     &  cg(:),
     &  tau_from_nu(:,:,:),
     &  tau_mesh(:),
     &  w_nu(:),
     &  w_nu_d(:)
      complex*16, allocatable :: cgr(:)
      end
