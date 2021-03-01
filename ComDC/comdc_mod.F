      module comdc_mod


      integer ::
     $  n_omega,
     $  n_tau,
     $  n_nu,
     $  iun,
     $  ndim_omega_max,
     $  ndim_tau_max,
     $  ndim_nu_max,            
     $  num_orb,
     $  is_magnetic,
     $  is_spinorbit,
     $  n_magn


      logical ::
     $  is_dynamical_f0,
     $  is_trans_dc

      integer, allocatable ::
     $  ndim_omega(:),
     $  ndim_tau(:),
     $  ndim_nu(:),
     $  n_mpi_omega(:),
     $  n_mpi_tau(:),
     $  n_mpi_nu(:)



      double precision :: 
     $  beta

      double precision,allocatable ::
     $  tau(:),
     $  dynamical_f0(:),
     $  ff(:)


      complex*16,allocatable :: 
     $  gimp(:,:,:,:),
     $  omega(:),
     $  nu(:),
     $  trans_dc(:,:,:)

      
      integer ::
     $  mpierr,me, me1,master,nproc,irank,isize
      logical ::
     $  maswrk

      end
      
      
      




      


      
