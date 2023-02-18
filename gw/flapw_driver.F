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
      integer :: global_n_bnd
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

        if(ubi==' hf') then
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
          global_n_bnd = 0

! do the on-node loop over k's to figure out whether the local node owns the k specified by global_k.
          do ind_k=1,ndim3_k(me3_k+1)
          if((n_mpi_kk(me+1)+ind_k).eq.global_k)then

          this_n_bnd=n_bnd(global_k,ispin)

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



! todotodowhy nbndf?
      allocate(sig_c_omega0(nbndf,nbndf,2))
      allocate(sig_c_diag(nbndf,0:n_omega))
      sig_c_omega0=(0.d0,0.d0)
      sig_c_diag=(0.d0,0.d0)

      do ind_omega=1,ndim3_omega
        i_omega=me3_tau*ndim3_omega+ind_omega-1


! fill sig_c_diag with the diagonal elements of sig_c_omega
        do jbnd=1,this_n_bnd
          sig_c_diag(jbnd,i_omega)=dcmplx(sig_c_omega(jbnd,jbnd,1,
     &      ind_omega,ind_k,ispin),
     &      sig_c_omega(jbnd,jbnd,2,ind_omega,ind_k,
     &      ispin))
        enddo

      enddo  !! over ind_omega

      do ind_omega=1,ndim3_omega
        i_omega=me3_tau*ndim3_omega+ind_omega-1

! fill sig_c_omega0(:,:,1,2) with the i_omega=0,1 elements of sig_c_omega
        if(i_omega<=1) then
          call ferm_unpack_omega(sig_c_omega0(1,1,i_omega+1),
     &      sig_c_omega(1,1,1,ind_omega,ind_k,
     &      ispin),this_n_bnd,nbndf,nbndf)
        endif

      enddo  !! over ind_omega

! gather over frequencies
      if(nproc_tau/=1) then
        call dgop(sig_c_diag,2*nbndf*(n_omega+1),'  +',comm_pnk)
        call dgop(sig_c_omega0,2*nbndf*nbndf*2,'  +',comm_pnk)
      endif

! as input, sig_c_omega0 is the i_omega=0,1 elements of sig_c_omega
! as output, contains the omega=0 value of the self-energy, and the omega-derivative of the self-energy
! The sig_c_omega0(:,:,1) elements contain sigma(omega=0)
! The sig_c_omega0(:,:,2) elements contain the derivative dsigma/domega
! sig_c_diag is the diagonal elements of sig_c_omega, and is unchanged
! The only matrix elements that are touched are between (1:this_n_bnd), not (1:nbndf)
      call sigma_0_z_ren(this_n_bnd,nbndf,sig_c_omega0,sig_c_diag)

!  move sig_c_omega0(:,:,1) to sig_0 and sig_c_omega0(:,:,2) to zzz.
      do ibnd=1,this_n_bnd
        do jbnd=1,this_n_bnd
              sig_0(ibnd,jbnd)=sig_c_omega0(ibnd,jbnd,1)
              zzz(ibnd,jbnd)=sig_c_omega0(ibnd,jbnd,2)
        enddo ! jbnd=1,this_n_bnd
      enddo ! ibnd=1,this_n_bnd


      deallocate(sig_c_omega0)
      deallocate(sig_c_diag)


!  zzz = Z^-1, so invert Z^-1 to get Z (this is not the best way to invert anything! Better to use zgetrf/zgetri)
          tmpmat=zzz


! Do the diagonalization only if this node  is the root of times/omega's.
      if(me_tau_kk_pbr==0) then
      call eig_val_solver(this_n_bnd,
     c         this_n_bnd,tmpmat,eigs)
      end if

! distribute the results of the diagonalization to all processors used for parallelizing tau
      if(nproc_tau/=1) then
        call brdcst(msgdbl,eigs,8*this_n_bnd,0,comm_tau_kk_pbr)
        call brdcst(msgdbl,tmpmat,16*this_n_bnd*this_n_bnd,0,comm_tau_kk_pbr)
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



! Copy sig_0 to sig_0_global, zzz to zzz_global, qpHam to qpHam_global, this_n_bnd to global_n_bnd.
! (No reason to have this duplication, except if one is thinking of changing the code to collect all k-points in zzz, sig_0, qpHam.)
      global_n_bnd = this_n_bnd
      do ibnd=1,global_n_bnd
      do jbnd=1,global_n_bnd
         zzz_global(ibnd,jbnd) = zzz(ibnd,jbnd)
         sig_0_global(ibnd,jbnd) = sig_0(ibnd,jbnd)
         qpHam_global(ibnd,jbnd) = qpHam(ibnd,jbnd)
      enddo ! jbnd=1,global_n_bnd
      enddo ! ibnd=1,global_n_bnd

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
      if(nproc_k/=1) call dgop(zzz_global,nbndf*nbndf,'  +',comm_k)
      if(nproc_k/=1) call dgop(sig_0_global,nbndf*nbndf,'  +',comm_k)
      if(nproc_k/=1) call dgop(qpHam_global,nbndf*nbndf,'  +',comm_k)
      if(nproc_k/=1)  call igop(global_n_bnd,1,'  +',comm_k) ! todotodo I think that global_n_bnd is not necessary

! Print out Z and Sigma to disk.
      if(maswrk) then
      do ibnd=1,global_n_bnd
      do jbnd=1,global_n_bnd
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

      enddo ! jbnd=1,global_n_bnd
      enddo ! ibnd=1,global_n_bnd
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

      allocate(dossp(0:ndos))
      allocate(sig_omega_global(nbndf,nbndf,0:n_omega))
      allocate(sig_freq_global(nbndf,nbndf,0:ndos))

      if(ubi==' hf') then
         return
      end if

      allocate(nd_dos(nproc_tau))
      allocate(nmpi_dos(nproc_tau))
      call size_shift_par(ndos+1,nproc_tau,nd_dos,nmpi_dos)
      dossp=0.d0
      dos_k=0.d0

      if(maswrk) then
        if(ubi==' gw') then
          open(3,file='qpSigmaImagFreq.gw')
          open(4,file='qpSigmaRealFreq.gw')
          open(5,file='qpDOS.gw')
        else if(ubi=='gwg') then
          open(3,file='qpSigmaImagFreq.gwg')
          open(4,file='qpSigmaRealFreq.gwg')
          open(5,file='qpDOS.gwg')
        else if(ubi==' qp') then
          open(3,file='qpSigmaImagFreq.qp')
          open(4,file='qpSigmaRealFreq.qp')
          open(5,file='qpDOS.qp')
        endif
      endif

      deltae=(emaxdos-emindos)/ndos

! This loop over spin and k contains the rest of the routine, except for file closes and deallocations.
      do ispin=1,nspin
      do global_k = 1,npnt

          dossp = 0.d0
          sig_omega_global=(0.d0,0.d0)
          sig_freq_global=(0.d0,0.d0)
          this_n_bnd=n_bnd(global_k,ispin)

          ! parallelize over k
          do ind_k=1,ndim3_k(me3_k+1)
          if((n3_mpi_k(me3_k+1)+ind_k).eq.global_k)then
 
      
            allocate(weigd(this_n_bnd))
            allocate(sig_om(this_n_bnd,this_n_bnd,0:n_omega))
            allocate(qqp(this_n_bnd,this_n_bnd))
            if(ubi==' gw') then
              allocate(achain(0:n_cont_frac,this_n_bnd,this_n_bnd))
            endif


            ! get qqp, the matrix that relates qsgw/linearized bands to Hartree-Fock bands. q_qp is the version of this that has been parallelized.
            ! Copy from q_qp to qqp.
              qqp=(0.d0,0.d0)
             do iband=1,this_n_bnd
                do jband=1,this_n_bnd
                   qqp(iband,jband)=q_qp(iband,jband,ind_k,ispin)
                 enddo
              enddo

! todotodo use qqp

! Begin preparing to analytically continue the frequency-dependent part of the self-energy.
          if(ubi==' qp'.or.ubi==' gw')) then

! Get sig_om from sig_c_omega.
          sig_om=(0.d0,0.d0)
          do ind_omega=1,ndim3_omega
          i_omega=me3_tau*ndim3_omega+ind_omega-1
            call ferm_unpack_omega(sig_om(1,1,i_omega),
     &         sig_c_omega(1,1,1,ind_omega,ind_k,
     &         ispin),this_n_bnd,this_n_bnd,nbndf)
          enddo
! Gather sig_om over frequencies.
          if(nproc_tau/=1) call dgop(sig_om,2*this_n_bnd*this_n_bnd*(n_omega+1),'  +',
     &                                 comm_pnk)

          end if ! (ubi==' qp'.or.ubi==' gw')) then


! Copy sig_om into sig_omega_global, and add e_bnd - chem_pot.
! sig_omega_global = sig_omega + e_bnd - chem_pot
! sig_omega_global will be printed out in qpSigmaImagFreq.*
          do ind_omega=0,n_omega
          do iband=1,this_n_bnd
          do jband=1,this_n_bnd
          sig_omega_global(iband,jband,ind_omega)=
     c             sig_om(iband,jband,ind_omega)

          if(iband .eq. jband) then
          sig_omega_global(iband,jband,ind_omega) =
     c          sig_omega_global(iband,jband,ind_omega) +
     c          e_bnd(iband,n0_k+ind_k,ispin) - chem_pot
          end if ! (iband .eq. jband)

          end do ! jband=1,this_n_bnd
          end do ! iband=1,this_n_bnd
          end do  ! ind_omega=0,n_omega


          if(ubi==' gw')) then

          if(n_cont_frac<=n_omega_exa) then

c --- Another variant of AC: Continued Fractions approximation ------
            call pade_fermion_chain(sig_om,this_n_bnd,this_n_bnd,
     &                              achain(0,1,1),
     &                              n_cont_frac)

          else ! not(n_cont_frac<=n_omega_exa)

          allocate(s0(this_n_bnd,this_n_bnd,2,0:n_omega))
          allocate(s0_cheb(this_n_bnd,this_n_bnd,2,0:n_omega_asy-2))
          s0=0.d0

            do ind_omega=1,ndim3_omega
            i_omega=me3_tau*ndim3_omega+ind_omega-1

! populate s0 from sig_c_omega
              s0(:,:,:,i_omega)=
     &      sig_c_omega(1:this_n_bnd,1:this_n_bnd,:,ind_omega,ind_k,ispin)
            enddo

! gather s0 over freqencies
          if(nproc_tau/=1) call dgop(s0,2*this_n_bnd*this_n_bnd*(n_omega+1),'  +',
     &                                   comm_pnk)

! s0_cheb = s0 * cheb_omega.
          call dgemm('n','t',2*this_n_bnd*this_n_bnd,n_omega_asy-1,n_omega_asy-1,1.d0,
     &                 s0(1,1,1,n_omega_exa+1),2*this_n_bnd*this_n_bnd,cheb_omega(0,1),
     &                 n_omega_asy-1,0.d0,s0_cheb(1,1,1,0),2*this_n_bnd*this_n_bnd)

! Prepare for analytic continuation and put the results in achain, which is this_n_bnd.
            call ac_fermion_chain(s0,s0_cheb,this_n_bnd,this_n_bnd,
     &                            achain(0,1,1),
     &                            n_cont_frac)

            deallocate(s0)
            deallocate(s0_cheb)

          endif ! not(n_cont_frac<=n_omega_exa)

          endif ! (ubi==' gw') ! Done with preparing to analytically continue the frequency-dependent part of the self-energy.



      do ind_dos=1,nd_dos(me3_tau+1)
        iedos=nmpi_dos(me3_tau+1)+ind_dos-1
        om0=emindos+iedos*deltae

        if(ubi=='dft'.or.ubi==' hf'.or.ubi==' qp') then
! Evaluate weigd using the Cauchy formula.
            weigd=0.d0
            om0=dcmplx(e,0.d0)
            omi=dcmplx(0.d0,e_small)

            if(ubi=='dft'.or.ubi==' hf') then

              do iband=1,this_n_bnd
                gg=imag((1.d0,0.d0)/(om0+omi+chem_pot-e_bnd(iband,global_k,ispin)))
                  weigd(iband)=-wgt(k)*gg/pi
              enddo ! iband

            else if(ubi==' qp') then

              do iband=1,this_n_bnd
                gg=imag((1.d0,0.d0)/(om0+omi+chem_pot_qp-e_qp(iband,global_k,ispin)))

                do jband=1,this_n_bnd
                weigd(jband)=weigd(jband)
     c     -wgt(k)*gg*abs(q_qp(iband,jband,ind_k,ispin))**2
     c                          /pi
                enddo ! jband

              enddo ! iband

          endif ! (ubi==' qp')


! For gw, take the trace of the Green's function.
        else if(ubi==' gw') then

! todotodo setup inv_green,green_om  should be of size (this_n_bnd,this_n_bnd)

! finish analytically continuing, do something about nbndf to change to this_n_bnd
          call pade_fermion_chain_out(this_n_bnd,nbndf,achain(0,1,1),
     &                                  om0,inv_green,n_cont_frac)
          green_om=(0.d0,0.d0)
          do jband=1,this_n_bnd
            do iband=1,this_n_bnd
              inv_green(iband,jband)=-inv_green(iband,jband)
              sig_freq_global(iband,jband,iedos) = inv_green(iband,jband)
              if(jband==iband) then
                inv_green(iband,iband)=inv_green(iband,iband)+om0+chem_pot-e_bnd(iband,k,ispin)
                green_om(iband,iband)=(1.d0,0.d0)
              endif
            enddo
          enddo


! zgetrf/zgetrs inverts a and puts the result in green_om
          CALL ZGETRF(this_n_bnd,this_n_bnd,inv_green,this_n_bnd,ipi,INFO)
          CALL ZGETRS('n',this_n_bnd,this_n_bnd,inv_green,this_n_bnd,ipi,green_om,this_n_bnd,info)

! calculate the trace of the Green's function green_om.
          do iband=1,this_n_bnd
             weigd(iband)=-imag(green_om(iband,iband))/pi*wgt(k)
          enddo

           endif ! (ubi==' gw')

          dossp(iedos)=0.d0
          do ibnd=1,ndim
            dossp(iedos)=dossp(iedos)+(weigd(ibnd)*2.d0/dfloat(nspin)/nrel)
          enddo               !!! over bands


          enddo ! ind_dos=1,nd_dos(me3_tau+1)

          if(nproc_tau/=1) then
             call dgop(dossp,(1+ndos),'  +',comm_pnk)
             call dgop(sig_freq_global,2*nbndf*nbndf*(1+ndos),'  +',comm_pnk)
          endif

            deallocate(weigd)
            deallocate(sig_om)
            deallocate(qqp)
            if(ubi==' gw') then
              allocate(achain(0:n_cont_frac,this_n_bnd,this_n_bnd))
            endif


            endif ! ((n_mpi_kk(me+1)+ind_k).eq.global_k)
            enddo ! ind_k=1,nd_k

! todotodo gather over k, to the root

      ! Print to disk.
       if(maswrk) then

       do iedos=0,ndos

       write(5,'(i5,i5,i5,1(1x,f12.7))')ispin,
     c      global_k,iedos,dossp(iedos)

       do iband=1,this_n_bnd
       do jband=1,this_n_bnd
          write(4,'(i5,i5,i5,i5,i5,1(1x,f12.7))')ispin,
     c      global_k,iedos,iband,jband,
     c      (sig_freq_global(iband,jband,iedos))
       enddo ! jbnd=1,this_n_bnd
       enddo ! ibnd=1,this_n_bnd

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
        close(4)
        close(5)
      endif


      deallocate(dossp)
      deallocate(sig_omega_global)
      deallocate(sig_freq_global)
      deallocate(nd_dos)
      deallocate(nmpi_dos)

      end subroutine ! print_Sigma


      end module flapw_driver
