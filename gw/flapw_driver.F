      module flapw_driver
! Copyright (c) 2022-2023, Vincent Sacksteder


      implicit none
      private


      public :: print_qp_Z_Sigma
      public :: print_Sigma

      contains




      subroutine print_qp_Z_Sigma()
      use atom_mod
      use etot_mod
      use manager_mod
!      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none
      integer :: jband,iband,j0,i0,global_k
      complex*16 :: cc
      real*8, allocatable :: eigs(:)
      complex*16, allocatable :: sig_0(:,:),zzz(:,:)
      complex*16, allocatable :: zhalf(:,:)
      integer :: ispin,ind_k,this_n_bnd,ibnd,jbnd
      complex*16, allocatable :: tmpmat(:,:),qqp(:,:)
      complex*16, allocatable :: tmpmat2(:,:),qpHam(:,:)
      complex*16, allocatable :: zzz_global(:,:),
     c                             sig_0_global(:,:),
     c                             qpHam_global(:,:)
      integer :: ind_omega,i_omega
      complex*16, allocatable :: sig_c_omega0(:,:,:),sig_c_diag(:,:)

        if((ubi==' hf').or.(ubi=='dft')) then
           return
        end if

      if(maswrk) then
        if(ubi==' gw') then
          open(3,file='qpZ.gw')
          open(4,file='qpSigma.gw')
          open(5,file='qpHam.gw')
          open(6,file='qpEig.gw')
        else if(ubi=='gwg') then
          open(3,file='qpZ.gwg')
          open(4,file='qpSigma.gwg')
          open(5,file='qpHam.gwg')
          open(6,file='qpEig.gwg')
        else if(ubi==' qp') then
          open(3,file='qpZ.qp')
          open(4,file='qpSigma.qp')
          open(5,file='qpHam.qp')
          open(6,file='qpEig.qp')
        endif
      endif

      allocate(sig_0_global(nbndf,nbndf))
      allocate(zzz_global(nbndf,nbndf))
      allocate(qpHam_global(nbndf,nbndf))

! Almost all of this routine is inside these loops, except for deallocations and file closes.
      do ispin=1,nspin
        do global_k = 1,npnt

! prep the variables for the mpi gather
          zzz_global=(0.d0,0.d0)
          sig_0_global=(0.d0,0.d0)
          qpHam_global=(0.d0,0.d0)
          this_n_bnd=n_bnd(global_k,ispin)

! do the on-node loop over k's to figure out whether the local node owns the k specified by global_k.
          do ind_k=1,ndim3_k(me3_k+1)
          if((n_mpi_kk(me+1)+ind_k).eq.global_k)then


          allocate(tmpmat(this_n_bnd,this_n_bnd))
          allocate(tmpmat2(this_n_bnd,this_n_bnd))
          allocate(qpHam(this_n_bnd,this_n_bnd))
          allocate(sig_0(this_n_bnd,this_n_bnd))
          allocate(zzz(this_n_bnd,this_n_bnd))
          allocate(zhalf(this_n_bnd,this_n_bnd))
          allocate(qqp(this_n_bnd,this_n_bnd))
          allocate(eigs(this_n_bnd))



          ! get qqp, the matrix that relates qsgw/linearized bands to Hartree-Fock bands. q_qp is the version of this that has been parallelized.
          ! Copy from q_qp to qqp.
            qqp=(0.d0,0.d0)
            do ibnd=1,this_n_bnd
               do jbnd=1,this_n_bnd
                     qqp(ibnd,jbnd)=q_qp(ibnd,jbnd,ind_k,ispin)
               enddo ! jbnd=1,this_n_bnd
            enddo ! ibnd=1,this_n_bnd
            ! Done getting qqp



      allocate(sig_c_omega0(this_n_bnd,this_n_bnd,2))
      allocate(sig_c_diag(this_n_bnd,0:n_omega))
      sig_c_omega0=(0.d0,0.d0)
      sig_c_diag=(0.d0,0.d0)

      do ind_omega=1,ndim3_omega
        i_omega=me3_tau*ndim3_omega+ind_omega-1


! fill sig_c_diag with the diagonal elements of sig_c_omega
        do jbnd=1,this_n_bnd
          sig_c_diag(jbnd,i_omega)=dcmplx(sig_c_omega(jbnd,jbnd,1,
     c      ind_omega,ind_k,ispin),
     c      sig_c_omega(jbnd,jbnd,2,ind_omega,ind_k,
     c      ispin))
        enddo

      enddo  !! over ind_omega

      do ind_omega=1,ndim3_omega
        i_omega=me3_tau*ndim3_omega+ind_omega-1

! fill sig_c_omega0(:,:,1,2) with the i_omega=0,1 elements of sig_c_omega
! we use nbndf because sig_c_omega is of size nbndf.
        if(i_omega<=1) then
          call ferm_unpack_omega(sig_c_omega0(1,1,i_omega+1),
     c      sig_c_omega(1,1,1,ind_omega,ind_k,
     c      ispin),this_n_bnd,nbndf,this_n_bnd)
        endif

      enddo  !! over ind_omega

! gather over frequencies
      if(nproc_tau/=1) then
        call dgop(sig_c_diag,2*this_n_bnd*(n_omega+1),'  +',comm_pnk)
        call dgop(sig_c_omega0,2*this_n_bnd*this_n_bnd*2,'  +',comm_pnk)
      endif

! About sigma_0_z_ren:
! as input, sig_c_omega0 is the i_omega=0,1 elements of sig_c_omega
! as output, contains the omega=0 value of the self-energy, and the omega-derivative of the self-energy
! The sig_c_omega0(:,:,1) elements contain sigma(omega=0)
! The sig_c_omega0(:,:,2) elements contain the derivative dsigma/domega
! sig_c_diag is the diagonal elements of sig_c_omega, and is unchanged
!      call sigma_0_z_ren(this_n_bnd,this_n_bnd,sig_c_omega0,sig_c_diag)

!  move sig_c_omega0(:,:,1) to sig_0 and sig_c_omega0(:,:,2) to zzz.
      do ibnd=1,this_n_bnd
        do jbnd=1,this_n_bnd
              sig_0(ibnd,jbnd)=sig_c_omega0(ibnd,jbnd,1)
              zzz(ibnd,jbnd)=sig_c_omega0(ibnd,jbnd,2)
        enddo ! jbnd=1,this_n_bnd
      enddo ! ibnd=1,this_n_bnd


      deallocate(sig_c_omega0)
      deallocate(sig_c_diag)


!  zzz = Z^-1, so invert Z^-1 to get Z. (This is not the best way to invert anything! Better to use zgetrf/zgetri.)
          tmpmat=zzz


! Do the diagonalization only if this node  is the root of times/omega's.
      if(me_tau_kk_pbr==0) then
      call eig_val_solver(this_n_bnd,
     c         this_n_bnd,tmpmat,eigs)
      end if

! distribute the results of the diagonalization to all processors used for parallelizing tau
      if(nproc_tau/=1) then
        call brdcst(msgdbl,eigs,8*this_n_bnd,0,comm_tau_kk_pbr)
        call brdcst(msgdbl,tmpmat,16*this_n_bnd*this_n_bnd,0,
     c          comm_tau_kk_pbr)
      endif

! complain if the first eigenvalue of Z is negative
      if(eigs(1)<=0.d0) then
        print *,'me : ',me,' problem in get_Z_Sigma_qp'
        print *,' D(1):',eigs(1)
      endif


! Build Z.
! tmpmat contains the unitary matrix containing the eigenvectors of Z^{-1}.
      zzz=(0.d0,0.d0)
      do jbnd=1,this_n_bnd
        do ibnd=1,this_n_bnd
          cc=conjg(tmpmat(jbnd,ibnd))/eigs(ibnd)
          call zaxpy(this_n_bnd,cc,
     c         tmpmat(1,ibnd),1,zzz(1,jbnd),1)
        enddo
      enddo
! Done building Z.

! Build sqrt(Z).
! tmpmat contains the unitary matrix containing the eigenvectors of Z^{-1}.
      zhalf=(0.d0,0.d0)
      do jbnd=1,this_n_bnd
        do ibnd=1,this_n_bnd
          cc=conjg(tmpmat(jbnd,ibnd))/sqrt(eigs(ibnd))
          call zaxpy(this_n_bnd,cc,
     c         tmpmat(1,ibnd),1,zhalf(1,jbnd),1)
        enddo
      enddo
! Done building zhalf = sqrt(Z).


       ! add the Hartree-Fock energy to the self-energy
       ! G^{-1} = w+chem_pot-e_bnd - Sigma
       do ibnd=1,this_n_bnd
          sig_0(ibnd,ibnd) = sig_0(ibnd,ibnd)
     c       + e_bnd(ibnd,n_mpi_kk(me+1)+ind_k,ispin) - chem_pot
       end do

! Build the quasiparticle Hamiltonian.
      tmpmat2=(0.d0,0.d0)
      qpHam=(0.d0,0.d0)
      call zgemm('n','n',this_n_bnd,this_n_bnd,this_n_bnd,
     c   (1.d0,0.d0),sig_0,this_n_bnd,zhalf,
     c   this_n_bnd,(0.d0,0.d0),tmpmat2,this_n_bnd)
      call zgemm('n','n',this_n_bnd,this_n_bnd,this_n_bnd,
     c       (1.d0,0.d0),zhalf,this_n_bnd,tmpmat2,this_n_bnd,
     c       (0.d0,0.d0),qpHam,this_n_bnd)
! qpHam = sqrt(Z) * (Sigma + e_bnd - chem_pot) * sqrt(Z). q_qp and e_qp are the eigenvectors and eigenvalues of qpHam, up to a shift by chem_pot.
! Done building the quasiparticle Hamiltonian.


! Transform the self-energy and Z factor to the quasiparticle basis.
! sigma = qqp * sigma * qqp^\dagger.  
! zzz = qqp * zzz * qqp^\dagger.
      call zgemm('n','c',this_n_bnd,
     c             this_n_bnd,this_n_bnd,
     c             (1.d0,0.d0),zzz,
     c             this_n_bnd,qqp,
     c             this_n_bnd,(0.d0,0.d0),
     c             tmpmat,this_n_bnd)

      call zgemm('n','n',this_n_bnd,
     c             this_n_bnd,this_n_bnd,
     c             (1.d0,0.d0),qqp,
     c             this_n_bnd,tmpmat,
     c             this_n_bnd,(0.d0,0.d0),
     c             zzz,this_n_bnd)

      call zgemm('n','c',this_n_bnd,
     c             this_n_bnd,this_n_bnd,
     c             (1.d0,0.d0),sig_0,
     c             this_n_bnd,qqp,
     c             this_n_bnd,(0.d0,0.d0),
     c             tmpmat,this_n_bnd)

      call zgemm('n','n',this_n_bnd,
     c             this_n_bnd,this_n_bnd,
     c             (1.d0,0.d0),qqp,
     c             this_n_bnd,tmpmat,
     c             this_n_bnd,(0.d0,0.d0),
     c             sig_0,this_n_bnd)
! Done transforming to the quasiparticle basis.



! Copy sig_0 to sig_0_global, zzz to zzz_global, qpHam to qpHam_global.
! (No reason to have this duplication, except if one is thinking of changing the code to collect all k-points in zzz, sig_0, qpHam.)
      do ibnd=1,this_n_bnd
      do jbnd=1,this_n_bnd
         zzz_global(ibnd,jbnd) = zzz(ibnd,jbnd)
         sig_0_global(ibnd,jbnd) = sig_0(ibnd,jbnd)
         qpHam_global(ibnd,jbnd) = qpHam(ibnd,jbnd)
      enddo ! jbnd=1,this_n_bnd
      enddo ! ibnd=1,this_n_bnd

      deallocate(eigs)
      deallocate(tmpmat)
      deallocate(tmpmat2)
      deallocate(qpHam)
      deallocate(zzz)
      deallocate(zhalf)
      deallocate(sig_0)
      deallocate(qqp)

! end of  the on-node loop over k's to figure out whether the local node owns the k specified by global_k.
         endif ! ((n_mpi_kk(me+1)+ind_k).eq.global_k)
         enddo ! ind_k=1,ndim3_k(me3_k+1)

! Gather data distributed across k-points.
      if(nproc_k/=1) call dgop(zzz_global,2*nbndf*nbndf,'  +',comm_k)
      if(nproc_k/=1) call dgop(sig_0_global,2*nbndf*nbndf,'  +',comm_k)
      if(nproc_k/=1) call dgop(qpHam_global,2*nbndf*nbndf,'  +',comm_k)


! Print out Z and Sigma and qp Hamiltonian and eigenvalues to disk.
      if(maswrk) then
      do ibnd=1,this_n_bnd
      do jbnd=1,this_n_bnd
         write(3,'(i5,i5,i5,i5,1(1x,f12.7))')ispin,
     c      global_k,ibnd,jbnd,(zzz_global(ibnd,jbnd))

         write(4,'(i5,i5,i5,i5,1(1x,f12.7))')ispin,
     c      global_k,ibnd,jbnd,(sig_0_global(ibnd,jbnd))

         write(5,'(i5,i5,i5,i5,1(1x,f12.7))')ispin,
     c      global_k,ibnd,jbnd,(qpHam_global(ibnd,jbnd))

         if(ibnd .eq. jbnd) then
         write(6,'(i5,i5,i5,2(1x,f12.7))')ispin,
     c      global_k,ibnd,(e_bnd(ibnd,global_k,ispin)),
     c      (e_qp(ibnd,global_k,ispin))
         end if ! (ibnd .eq. jbnd) then

      enddo ! jbnd=1,this_n_bnd
      enddo ! ibnd=1,this_n_bnd
      endif ! maswrk

      enddo ! global_k = 1,npnt
      enddo   !! over ispin

      if(maswrk) then
        close(3)
        close(4)
        close(5)
        close(6)
      endif

      deallocate(sig_0_global)
      deallocate(zzz_global)
      deallocate(qpHam_global)

      end subroutine ! print_qp_Z_Sigma


      subroutine print_Sigma()
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      implicit none

      integer :: ispin, ind_k, iband, jband, this_n_bnd, ind_dos,
     c           iedos,  ind_omega, i_omega, global_k,
     c           i0,j0
      real*8 :: deltae,om0

      complex*16, allocatable :: sig_om(:,:,:)
      real*8, allocatable :: weigd(:)
      complex*16, allocatable :: achain(:,:,:)
      real*8, allocatable :: dossp(:)
      complex*16, allocatable :: sig_omega_global(:,:,:)
      complex*16, allocatable :: sig_freq_global(:,:,:)
      complex*16, allocatable :: qqp(:,:)
      real*8, allocatable :: s0(:,:,:,:),s0_cheb(:,:,:,:)
      integer, allocatable :: nd_dos(:),nmpi_dos(:)
      complex*16, allocatable :: tmpmat(:,:), tmpmat2(:,:)
      complex*16, allocatable :: green_om(:,:)
      complex*16 :: omi,gg
      integer :: info
      integer, allocatable :: ipiv(:)


      if((ubi==' hf').or.(ubi=='dft')) then
         return
      end if

      allocate(dossp(0:ndos))
      allocate(sig_omega_global(nbndf,nbndf,0:n_omega))

      if(ubi==' gw') then
         allocate(sig_freq_global(nbndf,nbndf,0:ndos))
      endif ! (ubi==' gw')

      dossp=0.d0

      allocate(nd_dos(nproc_tau))
      allocate(nmpi_dos(nproc_tau))
      call size_shift_par(ndos+1,nproc_tau,nd_dos,nmpi_dos)

      if(maswrk) then
        if(ubi==' gw') then
          open(3,file='qpSigmaImagFreq.gw')
          open(4,file='qpSigmaRealFreq.gw')
          open(5,file='qpDOS.gw')
        else if(ubi=='gwg') then
          open(3,file='qpSigmaImagFreq.gwg')
          open(5,file='qpDOS.gwg')
        else if(ubi==' qp') then
          open(3,file='qpSigmaImagFreq.qp')
          open(5,file='qpDOS.qp')
        endif

      endif

      deltae=(emaxdos-emindos)/ndos

! This loop over spin and k contains the rest of the routine, except for file closes and deallocations.
      do ispin=1,nspin
      do global_k = 1,npnt

          sig_omega_global=(0.d0,0.d0)

          if(ubi==' gw') then
             sig_freq_global=(0.d0,0.d0)
          endif ! (ubi==' gw')

          this_n_bnd=n_bnd(global_k,ispin)

          ! parallelize over k
          do ind_k=1,ndim3_k(me3_k+1)
          if((n3_mpi_k(me3_k+1)+ind_k).eq.global_k)then

            allocate(ipiv(this_n_bnd))
            allocate(weigd(this_n_bnd))
            allocate(sig_om(this_n_bnd,this_n_bnd,0:n_omega))
            allocate(qqp(this_n_bnd,this_n_bnd))
            if(ubi==' gw') then
             allocate(achain(0:n_cont_frac,this_n_bnd,this_n_bnd))
            endif
            allocate(tmpmat(this_n_bnd,this_n_bnd))
            allocate(tmpmat2(this_n_bnd,this_n_bnd))
            allocate(green_om(this_n_bnd,this_n_bnd))


            ! get qqp, the matrix that relates qsgw/linearized bands to Hartree-Fock bands. q_qp is the version of this that has been parallelized.
            ! Copy from q_qp to qqp.
              qqp=(0.d0,0.d0)
             do iband=1,this_n_bnd
                do jband=1,this_n_bnd
                   qqp(iband,jband)=q_qp(iband,jband,ind_k,ispin)
                enddo
             enddo


! Begin preparing to analytically continue the frequency-dependent part of the self-energy.
          if(ubi==' qp'.or.ubi==' gw') then

! Get sig_om from sig_c_omega.
! We use nbndf because sig_c_omega is of size nbndf.
          sig_om=(0.d0,0.d0)
          do ind_omega=1,ndim3_omega
             i_omega=me3_tau*ndim3_omega+ind_omega-1
             call ferm_unpack_omega(sig_om(1,1,i_omega),
     c         sig_c_omega(1,1,1,ind_omega,ind_k,
     c         ispin),this_n_bnd,nbndf,this_n_bnd)
          enddo

! Gather sig_om over frequencies.
          if(nproc_tau/=1) call dgop(sig_om,
     c          2*this_n_bnd*this_n_bnd*(n_omega+1),'  +',comm_pnk)
          end if ! (ubi==' qp'.or.ubi==' gw')) then


! Copy sig_om into sig_omega_global, and add e_bnd - chem_pot.
! sig_omega_global = sig_omega + e_bnd - chem_pot
! sig_omega_global will be printed out in qpSigmaImagFreq.*
          do ind_omega=0,n_omega

          do iband=1,this_n_bnd
          do jband=1,this_n_bnd

          tmpmat=sig_om(iband,jband,ind_omega)

          if(iband .eq. jband) then
             tmpmat(iband,jband) =
     c          tmpmat(iband,jband) +
     c          e_bnd(iband,global_k,ispin) - chem_pot
          end if ! (iband .eq. jband)

          end do ! jband=1,this_n_bnd
          end do ! iband=1,this_n_bnd

! Transform sig_omega_global to the quasiparticle basis.
! sigma = qqp * sig_omega_global * qqp^\dagger
! sigma = qqp * (sig_omega + e_bnd - chem_pot) * qqp^\dagger
          call zgemm('n','c',this_n_bnd,
     c             this_n_bnd,this_n_bnd,
     c             (1.d0,0.d0),tmpmat,
     c             this_n_bnd,qqp,
     c             this_n_bnd,(0.d0,0.d0),
     c             tmpmat2,this_n_bnd)

          call zgemm('n','n',this_n_bnd,
     c             this_n_bnd,this_n_bnd,
     c             (1.d0,0.d0),qqp,
     c             this_n_bnd,tmpmat2,
     c             this_n_bnd,(0.d0,0.d0),
     c             tmpmat,this_n_bnd)

          sig_omega_global = (0.d0,0.d0)
          do iband=1,this_n_bnd
          do jband=1,this_n_bnd
             sig_omega_global(iband,jband,ind_omega)=tmpmat(iband,jband)
          end do ! jband=1,this_n_bnd
          end do ! iband=1,this_n_bnd

          end do  ! ind_omega=0,n_omega


          if(ubi==' gw') then

          if(n_cont_frac<=n_omega_exa) then

c --- Another variant of AC: Continued Fractions approximation ------
            call pade_fermion_chain(sig_om,this_n_bnd,this_n_bnd,
     c                              achain(0,1,1),
     c                              n_cont_frac)

          else ! not(n_cont_frac<=n_omega_exa)

          allocate(s0(this_n_bnd,this_n_bnd,2,0:n_omega))
          allocate(s0_cheb(this_n_bnd,this_n_bnd,2,0:n_omega_asy-2))
          s0=0.d0

            do ind_omega=1,ndim3_omega
            i_omega=me3_tau*ndim3_omega+ind_omega-1

! populate s0 from sig_c_omega
              s0(:,:,:,i_omega)=
     c           sig_c_omega(1:this_n_bnd,1:this_n_bnd,:,
     c           ind_omega,ind_k,ispin)
            enddo

! gather s0 over freqencies
          if(nproc_tau/=1) call dgop(s0,
     c       2*this_n_bnd*this_n_bnd*(n_omega+1),'  +',comm_pnk)

! s0_cheb = s0 * cheb_omega.
          call dgemm('n','t',2*this_n_bnd*this_n_bnd,
     c          n_omega_asy-1,n_omega_asy-1,1.d0,
     c          s0(1,1,1,n_omega_exa+1),2*this_n_bnd*this_n_bnd,
     c          cheb_omega(0,1),n_omega_asy-1,0.d0,s0_cheb(1,1,1,0),
     c                 2*this_n_bnd*this_n_bnd)

! Prepare for analytic continuation and put the results in achain, which is this_n_bnd.
            call ac_fermion_chain(s0,s0_cheb,this_n_bnd,this_n_bnd,
     c              achain(0,1,1),n_cont_frac)

            deallocate(s0)
            deallocate(s0_cheb)

          endif ! not(n_cont_frac<=n_omega_exa)

          endif ! (ubi==' gw') ! Done with preparing to analytically continue the frequency-dependent part of the self-energy.



      do ind_dos=1,nd_dos(me3_tau+1)

        iedos=nmpi_dos(me3_tau+1)+ind_dos-1

        if(ubi=='dft'.or.ubi==' hf'.or.ubi==' qp') then
        om0=dcmplx(emindos+iedos*deltae,0.d0)
        omi=dcmplx(0.d0,e_small)

        ! Evaluate weigd using the Cauchy formula.
            weigd=0.d0

            if(ubi=='dft'.or.ubi==' hf') then

! this calculates all the poles which will contribute to the dos
              do iband=1,this_n_bnd
                 weigd(iband)=-wgt(global_k)*imag((1.d0,0.d0)/(om0+omi
     c              +chem_pot-e_bnd(iband,global_k,ispin)))/pi
              enddo ! iband

            else if(ubi==' qp') then

              do iband=1,this_n_bnd
                gg=imag((1.d0,0.d0)/(om0+omi
     c             +chem_pot_qp-e_qp(iband,global_k,ispin)))

! this prepares for calculating the trace, performed in the quasiparticle basis.
                do jband=1,this_n_bnd
                weigd(jband)=weigd(jband)
     c             -abs(q_qp(iband,jband,ind_k,ispin))**2
     c                          *wgt(global_k)*gg/pi
                enddo ! jband

              enddo ! iband

          endif ! (ubi==' qp')


        else if(ubi==' gw') then
! For gw, calculate the Green's function and take its trace.


! finish analytically continuing
          call pade_fermion_chain_out(this_n_bnd,this_n_bnd,
     c                achain(0,1,1),om0,tmpmat2,n_cont_frac)

          green_om=(0.d0,0.d0)

         ! sig_freq_global = tmpmat2 + e_bnd - chem_pot
         ! tmpmat2 = -tmpmat2 - e_bnd + chem_pot + om0
         ! green_om = 1
          do jband=1,this_n_bnd
            do iband=1,this_n_bnd

              tmpmat(iband,jband) = tmpmat2(iband,jband)
              tmpmat2(iband,jband)=-tmpmat2(iband,jband)


              if(jband==iband) then
                tmpmat2(iband,iband)=tmpmat2(iband,iband)
     c               +om0+chem_pot-e_bnd(iband,global_k,ispin)

                 tmpmat(iband,jband) = tmpmat(iband,jband)
     c              -chem_pot+e_bnd(iband,global_k,ispin)

                 green_om(iband,iband)=(1.d0,0.d0)

              endif


            enddo ! iband=1,this_n_bnd
          enddo ! jband=1,this_n_bnd

! zgetrf/zgetrs inverts tmpmat2 and puts the result in green_om
          CALL ZGETRF(this_n_bnd,this_n_bnd,
     c         tmpmat2,this_n_bnd,ipiv,INFO)
          CALL ZGETRS('n',this_n_bnd,this_n_bnd,
     c         tmpmat2,this_n_bnd,ipiv,green_om,this_n_bnd,info)

! prep for calculating the trace of the Green's function green_om.
          do iband=1,this_n_bnd
             weigd(iband)=-imag(tmpmat2(iband,iband))/pi*wgt(global_k)
          enddo


! Transform sig_freq_global to the quasiparticle basis.
! sig_freq_global = qqp * sig_freq_global * qqp^\dagger.
! sig_freq_global = qqp * (tmpmat2 + e_bnd - chem_pot) * qqp^\dagger.
      call zgemm('n','c',this_n_bnd,
     c             this_n_bnd,this_n_bnd,
     c             (1.d0,0.d0),tmpmat,
     c             this_n_bnd,qqp,
     c             this_n_bnd,(0.d0,0.d0),
     c             tmpmat2,this_n_bnd)

      call zgemm('n','n',this_n_bnd,
     c             this_n_bnd,this_n_bnd,
     c             (1.d0,0.d0),qqp,
     c             this_n_bnd,tmpmat2,
     c             this_n_bnd,(0.d0,0.d0),
     c             tmpmat,this_n_bnd)

      sig_freq_global = (0.d0,0.d0)
      do iband=1,this_n_bnd
      do jband=1,this_n_bnd
         sig_freq_global(iband,jband,iedos) = tmpmat(iband,jband)
      end do ! jband=1,this_n_bnd
      end do ! iband=1,this_n_bnd


           endif ! (ubi==' gw')

! finish taking the trace
          dossp(iedos)=0.d0
          do iband=1,this_n_bnd
            dossp(iedos)=dossp(iedos)
     c        +(weigd(iband)*2.d0/dfloat(nspin)/nrel)
          enddo               !!! over bands


          enddo ! ind_dos=1,nd_dos(me3_tau+1)

! gather information distributed over tau/omega
          if(nproc_tau/=1) then
             call dgop(dossp,(1+ndos),'  +',comm_pnk)

             if(ubi==' gw') then
                call dgop(sig_freq_global,
     c              2*nbndf*nbndf*(1+ndos),'  +',comm_pnk)
             endif ! (ubi==' gw')

          endif

            deallocate(ipiv)
            deallocate(weigd)
            deallocate(sig_om)
            deallocate(qqp)
            deallocate(tmpmat)
            deallocate(tmpmat2)
            deallocate(green_om)

            if(ubi==' gw') then
              deallocate(achain)
            endif


            endif ! ((n_mpi_kk(me+1)+ind_k).eq.global_k)
            enddo ! ind_k=1,nd_k

       ! Gather data distributed across k-points, to the root.
       if(nproc_k/=1) call dgop(sig_omega_global,
     c    2*nbndf*nbndf*(1+n_omega),'  +',comm_k)
       if(nproc_k/=1) call dgop(sig_freq_global,
     c    2*nbndf*nbndf*(1+n_omega),'  +',comm_k)
       if(nproc_k/=1) call dgop(dossp,(1+ndos),'  +',comm_k)

       ! Print to disk.
       if(maswrk) then

       do iedos=0,ndos

       write(5,'(i5,i5,i5,1(1x,f12.7))')ispin,
     c      global_k,iedos,dossp(iedos)


       if(ubi==' gw') then

       do iband=1,this_n_bnd
       do jband=1,this_n_bnd
          write(4,'(i5,i5,i5,i5,i5,1(1x,f12.7))')ispin,
     c      global_k,iedos,iband,jband,
     c      (sig_freq_global(iband,jband,iedos))
       enddo ! jbnd=1,this_n_bnd
       enddo ! ibnd=1,this_n_bnd

       endif ! (ubi==' gw')

       end do ! iedos=0,ndos

       do iedos=0,n_omega
       do iband=1,this_n_bnd
       do jband=1,this_n_bnd
          write(3,'(i5,i5,i5,i5,i5,1(1x,f12.7))')ispin,
     c      global_k,iedos,iband,jband,
     c      (sig_omega_global(iband,jband,iedos))
       enddo ! jbnd=1,this_n_bnd
       enddo ! ibnd=1,this_n_bnd
       end do ! iedos=0,n_omega

       endif ! maswrk
       ! Done printing

      enddo ! global_k = 1,npnt
      enddo   !! over ispin

      if(maswrk) then
        close(3)

        if(ubi==' gw') then
           close(4)
        endif ! (ubi==' gw')

        close(5)
      endif

      deallocate(dossp)
      deallocate(sig_omega_global)

      if(ubi==' gw') then
         deallocate(sig_freq_global)
      endif ! (ubi==' gw')

      deallocate(nd_dos)
      deallocate(nmpi_dos)

      end subroutine ! print_Sigma


      end module flapw_driver
