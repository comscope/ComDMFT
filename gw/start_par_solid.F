      subroutine start_par_solid
	use atom_mod
	use manager_mod
	use parallel_mod
	use solid_mod
      implicit none
	allocate(ndim_k(nproc))
	allocate(n_mpi_k(nproc))
	allocate(ndim_kk(nproc))
	allocate(n_mpi_kk(nproc))
	allocate(ndim3_kk(nproc_k))
	allocate(n3_mpi_kk(nproc_k))
	allocate(ndim3_kkt(nproc_tau))
	allocate(n3_mpi_kkt(nproc_tau))
	allocate(ndim3_k(nproc_k))
	allocate(n3_mpi_k(nproc_k))
c ---- Memory distribution for K-points ---------------------------------
	call size_shift_par(npnt,nproc,ndim_k,n_mpi_k)
	call size_shift_par(npnt,nproc_k,ndim3_k,n3_mpi_k)
c ---- Memory distribution for K_All-points -----------------------------
	call size_shift_par(nqdiv,nproc,ndim_kk,n_mpi_kk)
	call size_shift_par(nqdiv,nproc_k,ndim3_kk,n3_mpi_kk)
	call size_shift_par(nqdiv,nproc_tau,ndim3_kkt,n3_mpi_kkt)
	allocate(ndim3_nr_red(nproc_k))
	allocate(n3_mpi_nr_red(nproc_k))
	call size_shift_par(nr_red,nproc_k,ndim3_nr_red,n3_mpi_nr_red)
	allocate(ndimt_nr_red(nproc_tau))
	allocate(nt_mpi_nr_red(nproc_tau))
	call size_shift_par(nr_red,nproc_tau,ndimt_nr_red,nt_mpi_nr_red)
	allocate(ndimt_nr(nproc_tau))
	allocate(nt_mpi_nr(nproc_tau))
	call size_shift_par(nr_full,nproc_tau,ndimt_nr,nt_mpi_nr)
	ndim4_kk=nqdiv
c ----------------------------------------------------------------------
	allocate(ndim_nfun(nproc_k))
	allocate(n_mpi_nfun(nproc_k))
	call size_shift_par(nfun,nproc_k,ndim_nfun,n_mpi_nfun)
	allocate(ndim_nrr(nproc_k))
	allocate(n_mpi_nrr(nproc_k))
	call size_shift_par(nrr_max,nproc_k,ndim_nrr,n_mpi_nrr)
	allocate(ndim_nrr_cor(nproc_k))
	allocate(n_mpi_nrr_cor(nproc_k))
	call size_shift_par(nrr_max_cor,nproc_k,ndim_nrr_cor,
     &	                  n_mpi_nrr_cor)
	allocate(ndim4_pbr(1,npnt))
	allocate(n4_mpi_pbr(1,npnt))
      end
