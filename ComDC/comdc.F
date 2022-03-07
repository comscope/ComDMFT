
      program ComDC


      use comdc_mod

      implicit none
      include "mpif.h"
      integer ::
     $  im,ii,jj,kk,iomega,itau,inu,io
      double precision ::pi
      character*30 :: ch
      double precision, allocatable :: 
     $  tempdble1(:,:),tempdble2(:,:,:,:),
     $  matsubara(:)
      complex*16 :: ai


      CALL MPI_INIT(mpierr)
      
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRANK,MPIERR)
      ME = irank
      me1=me+1
      MASTER = 0                                                        
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,MPIERR)      


      NPROC = isize
      MASWRK = ME .EQ. MASTER

      iun=993

      pi=datan2(1.0d0, 1.0d0)*4.0d0
      ai=dcmplx(0.0d0, 1.0d0)

      if (me .lt. 10) then
        write(ch, '(a, i1, a)') "check_comdc_", me, ".dat"
      elseif (me .lt. 100) then
        write(ch, '(a, i2, a)') "check_comdc_", me, ".dat"
      elseif (me .lt. 1000) then
        write(ch, '(a, i3, a)') "check_comdc_", me, ".dat"
      elseif (me .lt. 10000) then
        write(ch, '(a, i4, a)') "check_comdc_", me, ".dat"
      elseif (me .lt. 100000) then
        write(ch, '(a, i5, a)') "check_comdc_", me, ".dat"
      elseif (me .lt. 1000000) then
        write(ch, '(a, i6, a)') "check_comdc_", me, ".dat"        
      endif

c$$$      open(179, file=trim(ch))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     1. read parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (maswrk) then

        is_magnetic=0
        open(unit=10, file='comdc.ini', status='old')
        read(10,*) num_orb
        read(10,*) is_spinorbit
        if (is_spinorbit .eq. 0) then
          read(10,*) is_magnetic
        endif
        close(10)
        if (is_magnetic .eq. 1) then
          n_magn=2
        else
          n_magn=1
        endif

          

        inquire(file="trans_dc.dat", exist=is_trans_dc)
        
        if (is_trans_dc) then
          allocate(trans_dc(num_orb, num_orb, n_magn))
          trans_dc=0.0d0
          
          allocate(tempdble1(2,num_orb))          
          open(unit=10, file='trans_dc.dat', status='old')
          do im=1, n_magn
            do ii=1, num_orb
              read(10,*) ((tempdble1(kk,jj), kk=1, 2), jj=1, num_orb)
              do jj=1, num_orb
                trans_dc(ii,jj,im)
     $            =dcmplx(tempdble1(1,jj),tempdble1(2,jj))
              enddo
            enddo
          enddo
          close(10)
          deallocate(tempdble1)
        endif


        n_omega = 0
        open(unit=10, file='g_loc.dat', status='old')        
        DO
          READ(10,*,iostat=io)
          IF (io.ne.0) EXIT
          n_omega = n_omega + 1
        END DO
        CLOSE (10)
        write(*,*) 'n_omega', n_omega
        
        allocate(matsubara(n_omega))

        allocate(gimp(num_orb,num_orb,n_omega,n_magn))
        
        allocate(tempdble2(2,num_orb,num_orb,n_magn))
        
        open(unit=10, file='g_loc.dat', status='old')
        do iomega=1, n_omega
          read(10,*) matsubara(iomega),
     $      ((((tempdble2(kk,ii,jj,im), kk=1, 2), ii=1, num_orb),
     $      jj=1, num_orb), im=1, n_magn)
          write(*,*) iomega, matsubara(iomega)
          call flush(6)
          do im=1, n_magn
            do ii=1, num_orb
              do jj=1, num_orb
                gimp(ii,jj,iomega,im)
     $            =dcmplx(tempdble2(1,ii,jj,im),tempdble2(2,ii,jj,im))
              enddo
            enddo
          enddo
        enddo
        close(10)

        deallocate(tempdble2)        

        beta=pi/matsubara(1)
        deallocate(matsubara)


        allocate(ff(0:3))        
        open(unit=10, file='slater.dat', status='old')
        if (
     $    ((num_orb .eq. 1).and.(is_spinorbit .eq. 0))
     $    .or.
     $    ((num_orb .eq. 2).and.(is_spinorbit .eq. 1))        
     $    ) then        
          read(10,*) ff(0)
        elseif (
     $    ((num_orb .eq. 3).and.(is_spinorbit .eq. 0))
     $    .or.
     $    ((num_orb .eq. 6).and.(is_spinorbit .eq. 1))        
     $    ) then        
          read(10,*) (ff(ii), ii=0, 1)          
        elseif (
     $    ((num_orb .eq. 5).and.(is_spinorbit .eq. 0))
     $    .or.
     $    ((num_orb .eq. 10).and.(is_spinorbit .eq. 1))        
     $    ) then        
          read(10,*) (ff(ii), ii=0, 2)
        elseif (
     $      (num_orb .eq. 7) .and. (is_spinorbit .eq. 0)
     $      .or.
     $      (num_orb .eq. 14) .and. (is_spinorbit .eq. 1)
     $      ) then          
          read(10,*) (ff(ii), ii=0, 3)
        endif
        close(10)
      endif

      call mpi_bcast(num_orb,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(is_spinorbit,1,mpi_integer,0,MPI_COMM_WORLD,
     $  mpierr)
      call mpi_bcast(is_trans_dc,1,mpi_logical,0,MPI_COMM_WORLD,
     $  mpierr)      
      call mpi_bcast(is_magnetic,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(n_magn,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)      
      call mpi_bcast(n_omega,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)
      call mpi_bcast(beta,1,mpi_double_precision,0,
     $  MPI_COMM_WORLD,mpierr)      
      
      if (.not. maswrk) then
        allocate(ff(0:3))
        
        allocate(gimp(num_orb,num_orb,n_omega,n_magn))
        if (is_trans_dc) then
          allocate(trans_dc(num_orb, num_orb, n_magn))
          trans_dc=0.0d0
        endif
      endif
      
      if (is_trans_dc) then      
        call mpi_bcast(trans_dc,num_orb**2*n_magn,mpi_double_complex,0,
     $    MPI_COMM_WORLD,mpierr)
      endif
      call mpi_bcast(ff(0),4,mpi_double_precision,0,
     $  MPI_COMM_WORLD,mpierr)
      call mpi_bcast(gimp,num_orb**2*n_omega*n_magn,mpi_double_complex,
     $  0,MPI_COMM_WORLD,mpierr)
      

      if (maswrk) then

        
        inquire(file="dynamical_f0.dat", exist=is_dynamical_f0)
        if (is_dynamical_f0) then        
          n_nu = -1
          open(unit=10, file='dynamical_f0.dat', status='old')        
          DO
            READ(10,*,iostat=io)
            IF (io.ne.0) EXIT
            n_nu = n_nu + 1
          END DO
          CLOSE (10)
          allocate(dynamical_f0(0:n_nu))
          open(unit=10, file='dynamical_f0.dat', status='old')
          do inu=0, n_nu
            read(10,*) dynamical_f0(inu)
          enddo
          close(10)
        else
          n_nu=n_omega-1
        endif
        write(*,*) 'n_nu', n_nu                
      endif
      
      call mpi_bcast(is_dynamical_f0,1,mpi_logical,0,MPI_COMM_WORLD,
     $  mpierr)
      call mpi_bcast(n_nu,1,mpi_integer,0,MPI_COMM_WORLD,
     $  mpierr)            
      
      if (is_dynamical_f0) then
        if (.not. maswrk) then
          allocate(dynamical_f0(0:n_nu))
        endif
        call mpi_bcast(dynamical_f0(0),n_nu+1,mpi_double_precision,0,
     $    MPI_COMM_WORLD,mpierr)        
      endif

      allocate(ndim_omega(nproc)) !
      ndim_omega=0
      allocate(n_mpi_omega(nproc)) !
      n_mpi_omega=0
      
      call size_shift_par(n_omega,
     $  nproc,ndim_omega,n_mpi_omega)
      ndim_omega_max=maxval(ndim_omega)
      
      allocate(ndim_nu(nproc))  !
      ndim_nu=0
      allocate(n_mpi_nu(nproc)) !
      n_mpi_nu=0
      
      call size_shift_par
     $  (n_nu+1,nproc,ndim_nu,n_mpi_nu)
      ndim_nu_max=maxval(ndim_nu)      
      
      allocate(ndim_tau(nproc)) !
      ndim_tau=0
      allocate(n_mpi_tau(nproc)) !
      n_mpi_tau=0
      
      n_tau=n_nu*4          ! should be an even number
        
      call size_shift_par
     $  (n_tau+1,nproc,ndim_tau,n_mpi_tau) 
      ndim_tau_max=maxval(ndim_tau)
      
      allocate(nu(0:n_nu))
      nu=0.0d0
      allocate(omega(n_omega))
      omega=0.0d0
      allocate(tau(0:n_tau))
      tau=0.0d0

      do itau=0, n_tau          
        if (itau .le. n_tau/2) then
          tau(itau)=((dble(itau))/dble(n_tau/2))**3/2.0d0*beta
        else
          tau(itau)=(beta-tau(n_tau-itau))
        endif          
      enddo

      do inu=0, n_nu
        nu(inu)=inu*2.0d0*pi/beta*ai
      enddo
      
      do iomega=1, n_omega
        omega(iomega)=(2*iomega-1)*pi/beta*ai
      enddo
      
      
      call gw_double_counting

      if (is_trans_dc) then
        deallocate(trans_dc)
      endif
      deallocate(ff)
      deallocate(gimp)

      deallocate(ndim_omega)
      deallocate(n_mpi_omega)
      deallocate(ndim_nu)
      deallocate(n_mpi_nu)
      deallocate(ndim_tau)
      deallocate(n_mpi_tau)
      deallocate(nu)
      deallocate(omega)
      deallocate(tau)
      if (is_dynamical_f0) then      
        deallocate(dynamical_f0)
      endif      

c$$$      close(179)
      call mpi_finalize(mpierr)
      end

      subroutine gw_double_counting
! coulomb matrix 1*(r1) 2(r1) 3*(r2) 4(r2)
      
      use comdc_mod
      implicit none
      include "mpif.h"
      
      integer :: iorb,jorb,korb,lorb,
     $  im,ind_tau,itau,iomega,ind_omega,
     $  tau1,tau2,tau3,tau4,mtau,inu,ind_nu,
     $  wcmat0_irank,
     $  ii,jj,kk,ll

      double precision :: pi,slaterf(0:3)
      complex*16 ::
     $  ai,
     $  tempdcmplx,
     $  one_matrix_f(num_orb,num_orb),
     $  g_tau(num_orb,num_orb,0:n_tau,n_magn),
     $  g2mat(num_orb, num_orb,n_magn),
     $  g3mat(num_orb, num_orb,n_magn),
     $  gtemp(num_orb, num_orb),      
     $  v_h(num_orb,num_orb,n_magn),
     $  sig_x(num_orb,num_orb,n_magn),
     $  pola_tau(0:n_tau),
     $  pola_nu_temp1(0:n_tau),
     $  pola_nu_temp2(0:n_tau),            
     $  pola_nu(num_orb,num_orb,num_orb,num_orb,
     $  ndim_nu_max),
     $  epsilonmat(num_orb,num_orb,num_orb,num_orb,
     $  ndim_nu_max),
     $  epsiloninv(num_orb,num_orb,num_orb,num_orb,
     $  ndim_nu_max),
     $  umatrix(num_orb,num_orb,num_orb,num_orb),
     $  wmat(num_orb,num_orb,num_orb,num_orb,
     $  ndim_nu_max),
     $  wcmat(num_orb,num_orb,num_orb,num_orb,
     $  ndim_nu_max),
     $  wcmat0(num_orb,num_orb,num_orb,num_orb),
     $  wc_tau(0:n_tau),
     $  wc2mat(num_orb,num_orb,num_orb,num_orb),
     $  gwc_tau(num_orb, num_orb,0:n_tau,n_magn),
     $  gwc_omega_temp1(0:n_tau),
     $  gwc_omega_temp2(0:n_tau),      
     $  gwc_omega(num_orb,num_orb,n_omega,n_magn),
     $  sig_omega(num_orb,num_orb,n_omega,n_magn),
     $  gwc0(num_orb,num_orb,n_magn),
     $  gwc0_deriv(num_orb,num_orb,n_magn)
c$$$     $  exterp_g_corr1(num_orb,num_orb,0:n_tau,n_magn),
c$$$     $  exterp_g_corr2(num_orb,num_orb,0:n_tau,n_magn),
c$$$     $  exterp_g_corr3(num_orb,num_orb,0:n_tau,n_magn),
c$$$     $  exterp_wc_corr2(num_orb,num_orb,num_orb,num_orb,
c$$$     $  0:n_tau),
c$$$     $  exterp_wc_corr4(num_orb,num_orb,num_orb,num_orb,
c$$$     $  0:n_tau)
      

!!!!! ucal

c$$$      if (maswrk) then
c$$$        open(unit=111, file='umat.dat')
c$$$        do ii=n_iac_mat_i, n_iac_mat_f
c$$$          if ((is_magnetic .eq. 1) .and. (ii .eq. 0)) cycle
c$$$          do iorb=1, norb_iac(ii)
c$$$            do jorb=1, norb_iac(ii)
c$$$              do korb=1, norb_iac(ii)
c$$$                do lorb=1, norb_iac(ii)
c$$$                  write(111, '(5i5, 2f12.6)') 
c$$$     $              ii,iorb, jorb, korb, lorb,
c$$$     $              umatrix(iorb,jorb,korb,lorb,ii)
c$$$                enddo
c$$$              enddo
c$$$            enddo
c$$$          enddo
c$$$        enddo
c$$$        close(111)
c$$$      endif
      
      pi=4.0d0*datan2(1.0d0, 1.0d0)
      
      ai=dcmplx(0.0d0, 1.0d0)

      
      g2mat=0.0d0
      g3mat=0.0d0
      g_tau=0.0d0
      v_h=0.0d0
      sig_x=0.0d0
      pola_tau=0.0d0
      pola_nu=0.0d0
      epsilonmat=0.0d0
      epsiloninv=0.0d0
      wmat=0.0d0
      wcmat=0.0d0
      wc2mat=0.0d0
      wc_tau=0.0d0
      gwc_tau=0.0d0
      gwc_omega=0.0d0
      sig_omega=0.0d0      

c$$$      exterp_g_corr1=0.0d0
c$$$      exterp_g_corr2=0.0d0
c$$$      exterp_g_corr3=0.0d0
c$$$      
c$$$      exterp_wc_corr2=0.0d0
c$$$      exterp_wc_corr4=0.0d0

      one_matrix_f=0.0d0
      do iorb=1, num_orb
        one_matrix_f(iorb,iorb)=1.0d0
      enddo
      
      do im=1,n_magn
        g2mat(:,:,im)=(gimp(:,:,n_omega,im)
     $    +transpose(dconjg(gimp(:,:,n_omega,im))))
     $    /2.0d0*omega(n_omega)**2
        g3mat(:,:,im)=(gimp(:,:,n_omega,im)
     $    -transpose(dconjg(gimp(:,:,n_omega,im)))
     $    -one_matrix_f*2.0/omega(n_omega))
     $    /2.0d0*omega(n_omega)**3        
c$$$  
        gtemp=0.0d0
        do iorb=1, num_orb
          do jorb=iorb, num_orb        
            gtemp(iorb,jorb)=
     $        one_matrix_f(iorb,jorb)/omega(n_omega)
     $        +g2mat(iorb,jorb,im)/omega(n_omega)**2
     $        +g3mat(iorb,jorb,im)/omega(n_omega)**3
          enddo
        enddo
      
        do iorb=1, num_orb
          do jorb=iorb, num_orb
            if (cdabs(gtemp(iorb,jorb)-gimp(iorb,jorb,n_omega,im))
     $        .gt. 1.0d-6) then
c$$$              write(179,*) 'ERROR in g_imp extrapolation'
              if (maswrk) then
                write(iun,*) 'ERROR in g_imp extrapolation'
              endif
              call ending
            endif
          enddo
        enddo
        
        do iorb=1, num_orb
          do jorb=1, num_orb
            do ind_tau=1, ndim_tau(me1)
              itau=n_mpi_tau(me1)+ind_tau-1
!     exact interval
              do iomega=1, n_omega
                g_tau(iorb,jorb,itau,im)
     $            =g_tau(iorb,jorb,itau,im)
     $            +1.0d0/beta*cdexp(-tau(itau)*omega(iomega))
     $            *(
     $            gimp(iorb,jorb,iomega,im)
     $            -one_matrix_f(iorb,jorb)/omega(iomega)
     $            -g2mat(iorb,jorb,im)/(omega(iomega))**2
     $            -g3mat(iorb,jorb,im)/(omega(iomega))**3                  
     $            )
     $            +1.0d0/beta*cdexp(tau(itau)*omega(iomega))
     $            *(
     $            dconjg(gimp(jorb,iorb,iomega,im))
     $            +one_matrix_f(iorb,jorb)/omega(iomega)
     $            -g2mat(iorb,jorb,im)/(omega(iomega))**2
     $            +g3mat(iorb,jorb,im)/(omega(iomega))**3                 
     $            )
                
c$$$                exterp_g_corr1(iorb,jorb,itau,im)
c$$$     $            =exterp_g_corr1(iorb,jorb,itau,im)
c$$$     $            +1.0d0/beta*cdexp(-tau(itau)*omega(iomega))
c$$$     $            *(
c$$$     $            -one_matrix_f(iorb,jorb)/omega(iomega)
c$$$     $            )
c$$$     $            +1.0d0/beta*cdexp(tau(itau)*omega(iomega))
c$$$     $            *(
c$$$     $            +one_matrix_f(iorb,jorb)/omega(iomega)
c$$$     $            )
c$$$                
c$$$                exterp_g_corr2(iorb,jorb,itau,im)
c$$$     $            =exterp_g_corr2(iorb,jorb,itau,im)
c$$$     $            +1.0d0/beta*cdexp(-tau(itau)*omega(iomega))
c$$$     $            *(
c$$$     $            -g2mat(iorb,jorb,im)/(omega(iomega))**2
c$$$     $            )
c$$$     $            +1.0d0/beta*cdexp(tau(itau)*omega(iomega))
c$$$     $            *(
c$$$     $            -g2mat(iorb,jorb,im)/(omega(iomega))**2
c$$$     $            )
c$$$                
c$$$                exterp_g_corr3(iorb,jorb,itau,im)
c$$$     $            =exterp_g_corr3(iorb,jorb,itau,im)
c$$$     $            +1.0d0/beta*cdexp(-tau(itau)*omega(iomega))
c$$$     $            *(
c$$$     $            -g3mat(iorb,jorb,im)/(omega(iomega))**3                  
c$$$     $            )
c$$$     $            +1.0d0/beta*cdexp(tau(itau)*omega(iomega))
c$$$     $            *(
c$$$     $            +g3mat(iorb,jorb,im)/(omega(iomega))**3                 
c$$$     $            )
              enddo
!     asymptotic interval
              g_tau(iorb,jorb,itau,im)
     $          =g_tau(iorb,jorb,itau,im)
     $          -one_matrix_f(iorb,jorb)/2.0d0
     $          +g2mat(iorb,jorb,im)*beta/2.0d0
     $          *(tau(itau)/beta-1.0d0/2.0d0)
     $          -g3mat(iorb,jorb,im)*beta**2/4.0d0
     $          *((tau(itau)/beta)**2-(tau(itau)/beta))
              
c$$$              exterp_g_corr1(iorb,jorb,itau,im)
c$$$     $          =exterp_g_corr1(iorb,jorb,itau,im)
c$$$     $          -one_matrix_f(iorb,jorb)/2.0d0                
c$$$              
c$$$              exterp_g_corr2(iorb,jorb,itau,im)
c$$$     $          =exterp_g_corr2(iorb,jorb,itau,im)
c$$$     $          +g2mat(iorb,jorb,im)*beta/2.0d0
c$$$     $          *(tau(itau)/beta-1.0d0/2.0d0)                
c$$$              
c$$$              exterp_g_corr3(iorb,jorb,itau,im)
c$$$     $          =exterp_g_corr3(iorb,jorb,itau,im)
c$$$     $          -g3mat(iorb,jorb,im)*beta**2/4.0d0
c$$$     $          *((tau(itau)/beta)**2-(tau(itau)/beta))                
            enddo
          enddo
        enddo
      enddo
      call mpi_allreduce_dcmplx(g_tau,
     $  num_orb**2*(n_tau+1)*n_magn,mpi_sum,mpi_comm_world)

c$$$      call mpi_allreduce_dcmplx(exterp_g_corr1,
c$$$     $  num_orb**2*(n_tau+1)*n_magn,mpi_sum,mpi_comm_world)
c$$$      call mpi_allreduce_dcmplx(exterp_g_corr2,
c$$$     $  num_orb**2*(n_tau+1)*n_magn,mpi_sum,mpi_comm_world)
c$$$      call mpi_allreduce_dcmplx(exterp_g_corr3,
c$$$     $  num_orb**2*(n_tau+1)*n_magn,mpi_sum,mpi_comm_world)

c$$$      if (maxval(cdabs(exterp_g_corr3)) .gt. 1.0d-4) then
c$$$        write(179,*) maxval(cdabs(exterp_g_corr2)),
c$$$     $    maxval(cdabs(exterp_g_corr3))
c$$$        write(179,*) 'higher order correction to the g is necessary'
c$$$        if (maswrk) then
c$$$          write(iun,*) maxval(cdabs(exterp_g_corr2)),
c$$$     $      maxval(cdabs(exterp_g_corr3))        
c$$$          write(iun,*) 'higher order correction to the g is necessary'
c$$$        endif
c$$$        call ending
c$$$      endif

      if (is_dynamical_f0) then
        slaterf=ff+(/dynamical_f0(0),0.0d0, 0.0d0,0.0d0/)
      else
        slaterf=ff
      endif
      
      call cal_u_matrix(num_orb,slaterf,umatrix)

      do im=1, n_magn
!     hartree
        do iorb=1,num_orb
          do jorb=1,num_orb
            do korb=1,num_orb
              do lorb=1,num_orb
                if (is_spinorbit .eq. 0) then
                  if (is_magnetic .eq. 0) then
                    tempdcmplx=g_tau(lorb,korb,n_tau,im)*2                    
                  else
                    tempdcmplx
     $                =g_tau(lorb,korb,n_tau,1)
     $                +g_tau(lorb,korb,n_tau,2)
                  endif
                else
                  tempdcmplx
     $              =g_tau(lorb,korb,n_tau,im)
                endif
                v_h(iorb,jorb,im)=v_h(iorb,jorb,im)
     &            -umatrix(iorb,jorb,korb,lorb)
     &            *tempdcmplx
              enddo
            enddo
          enddo
        enddo
      enddo

      slaterf=ff
      call cal_u_matrix(num_orb,slaterf,umatrix)
      
      do im=1, n_magn
c$$$  !     fock
        do lorb=1,num_orb        
          do korb=1,num_orb
            do jorb=1,num_orb
              do iorb=1,num_orb
                sig_x(iorb,jorb,im)=sig_x(iorb,jorb,im)
     $            +umatrix(iorb,lorb,korb,jorb)
     &            *g_tau(lorb,korb,n_tau,im)
              enddo
            enddo
          enddo
        enddo
      enddo
      
!     polarizability 
      
c$$$      if (maswrk) write(iun,'(a, i5, a)') 'impurity', im,
c$$$     $  '   polarizability (tau) and (nu)'
      
!     polarizability at tau
!     tau1, tau2, tau3, tau4
!     n,    p,    q,    m        
      do tau4=1,num_orb
        do tau1=1,num_orb
          do tau2=1,num_orb
            do tau3=1,num_orb
              pola_tau=0.0d0
              do itau=0, n_tau
                mtau=n_tau-itau                  
                if (is_spinorbit .eq. 0) then
                  if (is_magnetic .eq. 0) then
                    pola_tau(itau)
     $                =pola_tau(itau) ! tau4,norb, tau2,tau3 order
     $                -g_tau(tau1,tau2,itau,1) !n,p
     $                *g_tau(tau3,tau4,mtau,1) !q,m
     $                *2.0d0
                  else
                    pola_tau(itau)
     $                = pola_tau(itau)
     $                -g_tau(tau1,tau2,itau,1)
     $                *g_tau(tau3,tau4,mtau,1)
     $                -g_tau(tau1,tau2,itau,2)
     $                *g_tau(tau3,tau4,mtau,2)
                  endif
                else
                  pola_tau(itau)
     $              =pola_tau(itau)
     $              -g_tau(tau1,tau2,itau,1)
     $              *g_tau(tau3,tau4,mtau,1)                    
                endif                  
              enddo
              do ind_nu=1, ndim_nu(me1)
                inu=n_mpi_nu(me1)+ind_nu-1
                pola_nu_temp1=0.0d0
                pola_nu_temp2=0.0d0                                    
                do itau=0, n_tau
                  pola_nu_temp1(itau) 
     $              =pola_nu_temp1(itau)
     $              +pola_tau(itau)*cdexp(nu(inu)*tau(itau))
                enddo
                call fderiv_dcmplx(-1,n_tau+1,tau(0),pola_nu_temp1(0),
     $            pola_nu_temp2(0))
                pola_nu(tau4,tau1,tau2,tau3,ind_nu)
     $            =pola_nu_temp2(n_tau)
              enddo
            enddo
          enddo
        enddo
      enddo
      
      wcmat0=0.0d0
      wcmat0_irank=0

      do ind_nu=1, ndim_nu(me1)
        inu=n_mpi_nu(me1)+ind_nu-1

        if (is_dynamical_f0) then        
          slaterf=ff+(/dynamical_f0(inu),0.0d0, 0.0d0,0.0d0/)
        else
          slaterf=ff
        endif
        call cal_u_matrix(num_orb,slaterf,umatrix)        

        call zgemm('n','n',num_orb**2,num_orb**2,
     $    num_orb**2,(-1.0d0,0.0d0),
     $    umatrix(1,1,1,1),num_orb**2,
     $    pola_nu(1,1,1,1,ind_nu),num_orb**2,
     $    (0.0d0,0.0d0),epsilonmat(1,1,1,1,ind_nu),num_orb**2)
        do tau1=1, num_orb
          do tau2=1,num_orb             
            epsilonmat(tau1,tau2,tau1,tau2,ind_nu)
     $        =epsilonmat(tau1,tau2,tau1,tau2,ind_nu)+1.0d0
          enddo
        enddo          
        call dcmplx_matinv(epsilonmat(1,1,1,1,ind_nu),
     $    epsiloninv(1,1,1,1,ind_nu),
     $    num_orb**2, num_orb**2)
        
        call zgemm('n','n',num_orb**2,num_orb**2,
     $    num_orb**2,(1.0d0,0.0d0),
     $    epsiloninv(1,1,1,1,ind_nu),num_orb**2,
     $    umatrix(1,1,1,1),num_orb**2,
     $    (0.0d0,0.0d0),wmat(1,1,1,1,ind_nu),num_orb**2)

        slaterf=ff
        call cal_u_matrix(num_orb,slaterf,umatrix)
        
        wcmat(:,:,:,:,ind_nu)
     $    =wmat(:,:,:,:,ind_nu)-umatrix(:,:,:,:)
        
        if (inu .eq. n_nu) then
          wcmat0=wcmat(:,:,:,:,ind_nu)
          wcmat0_irank=irank
        endif
      enddo
      
      call mpi_allreduce_int(wcmat0_irank,
     $  1,mpi_sum,mpi_comm_world)
      call mpi_bcast(wcmat0,num_orb**4,
     $  mpi_double_complex,wcmat0_irank,MPI_COMM_WORLD,mpierr)
      
      wc2mat=nu(n_nu)**2*wcmat0
      
      do tau1=1,num_orb
        do tau2=1,num_orb
          do tau3=1,num_orb
            do tau4=1,num_orb        
              wc_tau=0.0d0
              
              do itau=0, n_tau
                
!     exact 
c$$$  do inu=2, n_nu
                do ind_nu=1, ndim_nu(me1)
                  inu=n_mpi_nu(me1)+ind_nu-1
                  
                  if (inu .eq. 0) then
!     nu=0 component
                    wc_tau(itau)
     $                =wc_tau(itau)
     $                +1.0d0/beta
     $                *wcmat(tau1,tau4,tau3,tau2,ind_nu)
                  else
                    wc_tau(itau)
     $                =wc_tau(itau)
     $                +2.0d0/beta*
     $                (
     $                wcmat(tau1,tau4,tau3,tau2,ind_nu)
     $                -wc2mat(tau1,tau4,tau3,tau2)/(nu(inu))**2
     $                )
     $                *dcos(dimag(nu(inu))*tau(itau))
c$$$                    exterp_wc_corr2(tau1,tau4,tau3,tau2,itau)
c$$$     $                =exterp_wc_corr2(tau1,tau4,tau3,tau2,itau)
c$$$     $                +2.0d0/beta*
c$$$     $                (
c$$$     $                -wc2mat(tau1,tau4,tau3,tau2)/(nu(inu))**2
c$$$     $                )
c$$$     $                *dcos(dimag(nu(inu))*tau(itau))
                  endif
                enddo
!     asymptotic
                if (irank .eq. wcmat0_irank) then
                  wc_tau(itau)
     $              =wc_tau(itau)
     $              -wc2mat(tau1,tau4,tau3,tau2)
     $              *beta/2.0d0
     $              *((tau(itau)/beta)**2
     $              -(tau(itau)/beta)+1.0d0/6.0d0)
c$$$                  exterp_wc_corr2(tau1,tau4,tau3,tau2,itau)
c$$$     $              =exterp_wc_corr2(tau1,tau4,tau3,tau2,itau)
c$$$     $              -wc2mat(tau1,tau4,tau3,tau2)                    
c$$$     $              *beta/2.0d0
c$$$     $              *((tau(itau)/beta)**2
c$$$     $              -(tau(itau)/beta)+1.0d0/6.0d0)
                endif                  
              enddo
              
              call mpi_allreduce_dcmplx(wc_tau,
     $          n_tau+1,mpi_sum,mpi_comm_world)
              do im=1, n_magn
                do itau=0, n_tau
                  gwc_tau(tau1,tau2,itau,im)
     $              =gwc_tau(tau1,tau2,itau,im)
     $              -g_tau(tau4,tau3,itau,im)
     $              *wc_tau(itau)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

c$$$      call mpi_allreduce_dcmplx(exterp_wc_corr2(1,1,1,1,1),
c$$$     $  num_orb**4*(n_tau+1),mpi_sum,mpi_comm_world)

c$$$      gwc0=0.0d0
c$$$      gwc0_deriv=0.0d0
      
      do im=1, n_magn
        do tau1=1,num_orb
          do tau2=1,num_orb
            gwc_omega_temp1=0.0d0
            gwc_omega_temp2=0.0d0                                                
c$$$            do itau=0, n_tau
c$$$              gwc_omega_temp1(itau)
c$$$     $          =gwc_tau(tau1,tau2,itau,im)
c$$$            enddo
c$$$            call fderiv_dcmplx(-1,n_tau+1,tau(0),gwc_omega_temp1(0),
c$$$     $        gwc_omega_temp2(0))
c$$$            gwc0(tau1,tau2,im)=gwc_omega_temp2(n_tau)
c$$$
c$$$            gwc_omega_temp1=0.0d0
c$$$            gwc_omega_temp2=0.0d0                                                
c$$$            do itau=0, n_tau
c$$$              gwc_omega_temp1(itau)
c$$$     $          =gwc_tau(tau1,tau2,itau,im)*tau(itau)
c$$$            enddo
c$$$            call fderiv_dcmplx(-1,n_tau+1,tau(0),gwc_omega_temp1(0),
c$$$     $        gwc_omega_temp2(0))
c$$$            gwc0_deriv(tau1,tau2,im)=gwc_omega_temp2(n_tau)            
            
            do ind_omega=1, ndim_omega(me1)
              iomega=n_mpi_omega(me1)+ind_omega          
              
              gwc_omega_temp1=0.0d0
              gwc_omega_temp2=0.0d0                                    
              do itau=0, n_tau
                gwc_omega_temp1(itau) 
     $            =gwc_omega_temp1(itau)
     $            +gwc_tau(tau1,tau2,itau,im)
     $            *cdexp(omega(iomega)*tau(itau))
              enddo
              call fderiv_dcmplx(-1,n_tau+1,tau(0),gwc_omega_temp1(0),
     $          gwc_omega_temp2(0))              
              gwc_omega(tau1,tau2,iomega,im)=gwc_omega_temp2(n_tau)
            enddo
          enddo
        enddo
      enddo
      
      call mpi_allreduce_dcmplx(gwc_omega,
     $  num_orb**2*n_omega*n_magn,mpi_sum,mpi_comm_world)

      
      do im=1, n_magn
        do tau1=1,num_orb
          do tau2=1,num_orb
            do iomega=1, n_omega

              sig_omega(tau1,tau2,iomega,im)
     $          =v_h(tau1,tau2,im)
     $          +sig_x(tau1,tau2,im)
     $          +gwc_omega(tau1,tau2,iomega,im)
            enddo
          enddo
        enddo
      enddo

      if (maswrk) then

        slaterf=ff
        call cal_u_matrix(num_orb,slaterf,umatrix)
        
        open(unit=111, file='vmat.dat')
        do iorb=1, num_orb
          do jorb=1, num_orb
            do korb=1, num_orb
              do lorb=1, num_orb
                write(111, '(4i5, 2f12.6)') 
     $            iorb, jorb, korb, lorb,
     $            umatrix(iorb,jorb,korb,lorb)
              enddo
            enddo
          enddo
        enddo
        close(111)

        slaterf=ff+(/dynamical_f0(0),0.0d0, 0.0d0,0.0d0/)
        call cal_u_matrix(num_orb,slaterf,umatrix)
        
        open(unit=111, file='u0mat.dat')
        do iorb=1, num_orb
          do jorb=1, num_orb
            do korb=1, num_orb
              do lorb=1, num_orb
                write(111, '(4i5, 2f12.6)') 
     $            iorb, jorb, korb, lorb,
     $            umatrix(iorb,jorb,korb,lorb)
              enddo
            enddo
          enddo
        enddo
        close(111)        

        open(unit=111, file='wcmat0.dat')
        do iorb=1, num_orb
          do jorb=1, num_orb
            do korb=1, num_orb
              do lorb=1, num_orb
                write(111, '(4i5, 2f12.6)') 
     $            iorb, jorb, korb, lorb,
     $            wcmat0(iorb,jorb,korb,lorb)
              enddo
            enddo
          enddo
        enddo
        close(111)

        open(unit=111, file='nimp.dat')
        do im=1, n_magn
          do tau1=1, num_orb
            do tau2=1, num_orb
              write(111,'(3i5, 2f12.6)')
     $          im,tau1, tau2, -g_tau(tau1,tau2,n_tau,im)
            enddo
          enddo
        enddo
        close(111)        

        open(unit=111, file='hartree.dat')
        do im=1, n_magn
          do tau2=1, num_orb          
            do tau1=1, num_orb
              write(111,'(3i5, 2f12.6)')
     $          im,tau1, tau2, v_h(tau1,tau2,im)
            enddo
          enddo
        enddo
        close(111)
        open(unit=111, file='exchange.dat')
        do im=1, n_magn
          do tau2=1, num_orb          
            do tau1=1, num_orb
              write(111,'(3i5, 2f12.6)')
     $          im,tau1, tau2, sig_x(tau1,tau2,im)
            enddo
          enddo
        enddo
        close(111)

c$$$        open(unit=111, file='sig0.dat')
c$$$        do im=1, n_magn
c$$$          do tau1=1, num_orb
c$$$            do tau2=1, num_orb
c$$$              write(111,'(3i5, 2f12.6)')
c$$$     $          im,tau1, tau2, gwc0(tau1,tau2,im)
c$$$c$$$     $          +v_h(tau1,tau2,im)+sig_x(tau1,tau2,im)
c$$$            enddo
c$$$          enddo
c$$$        enddo
c$$$        close(111)
c$$$
c$$$        open(unit=111, file='sig_deriv.dat')
c$$$        do im=1, n_magn
c$$$          do tau1=1, num_orb
c$$$            do tau2=1, num_orb
c$$$              write(111,'(3i5, 2f12.6)')
c$$$     $          im,tau1, tau2, gwc0_deriv(tau1,tau2,im)
c$$$            enddo
c$$$          enddo
c$$$        enddo
c$$$        close(111)                
c$$$
c$$$        open(unit=111, file='gwc_omega.dat')
c$$$        do im=1, n_magn
c$$$          do iomega=1, n_omega
c$$$            do iorb=1, num_orb
c$$$              do jorb=1, num_orb
c$$$                write(111, '(4(i5,2x), 2f12.6)')
c$$$     $            im,iomega,iorb,jorb,gwc_omega(iorb,jorb,iomega,im)
c$$$              enddo
c$$$            enddo
c$$$          enddo
c$$$        enddo
c$$$        close(111)        
        
        
        open(unit=10, file='sig_mat.dat')
        do iomega=1, n_omega
          write(10,'(f20.12, 5x, 10000(2f20.12, 2x))')
     $      dimag(omega(iomega)),
     $      (((sig_omega(tau1,tau2,iomega,im), tau1=1, num_orb),
     $      tau2=1, num_orb),
     $      im=1, n_magn)
        enddo
        close(10)


        open(unit=10, file='sig_gw_mat.dat')
        do iomega=1, n_omega
          write(10,'(f20.12, 5x, 10000(2f20.12, 2x))')
     $      dimag(omega(iomega)),
     $      (((sig_omega(tau1,tau2,iomega,im)
     $      -v_h(tau1,tau2,im), tau1=1, num_orb),
     $      tau2=1, num_orb),
     $      im=1, n_magn)
        enddo
        close(10)
        open(unit=10, file='sig_gwc_mat.dat')
        do iomega=1, n_omega
          write(10,'(f20.12, 5x, 10000(2f20.12, 2x))')
     $      dimag(omega(iomega)),
     $      (((gwc_omega(tau1,tau2,iomega,im),
     $      tau1=1, num_orb),
     $      tau2=1, num_orb),
     $      im=1, n_magn)
        enddo
        close(10)                

        
      endif
        
      
      
      end      



      subroutine size_shift_par(n,npr,ndim_pr,n_mpi_pr)

c$$$  parallization of processor
c$$$  among the "npr" processor, we allocate "n" jobs
c$$$  to illustrate if we have 5 processor and 8 jobs
c$$$  the first 3 processors work on 2 jobs and the last work on 1 job      
c$$$  ndim_pr: the number of jobs in each processor: ex) (2, 2, 2, 1, 1)
c$$$  n_mpi_pr: (the index of first jobs-1) in each processor ex) (0, 2, 4, 6, 7)

      implicit none
      integer, intent(in) :: n,npr
      integer, intent(out) :: ndim_pr(npr),n_mpi_pr(npr)
      integer :: ish,irest,ind,i
      ish=n/npr
      irest=n-ish*npr
      ndim_pr(1:irest)=ish+1
      ndim_pr(irest+1:npr)=ish
      ind=0
      do i=1,npr
        n_mpi_pr(i)=ind
        ind=ind+ndim_pr(i)
      enddo
      end

      

      SUBROUTINE mpi_allreduce_int(X,LENX,OP,comm)
      implicit none
      INCLUDE 'mpif.h'
      integer,intent(in) :: lenx,op,comm
      integer, intent(inout) :: x(lenx)
      integer :: mxbuff
      integer,allocatable :: buff(:)
      integer :: i,j,npass,length,locx,ierror

      mxbuff=8192
      allocate(buff(mxbuff))
      NPASS = (LENX-1)/MXBUFF + 1
      LENGTH = MXBUFF
      LOCX = 1
      DO I=1,NPASS
        IF(I.EQ.NPASS) LENGTH = LENX - MXBUFF*(NPASS-1)
        CALL MPI_ALLREDUCE(X(LOCX),BUFF,LENGTH,MPI_INTEGER,
     $    OP,comm,IERROR)
c$$$        if (ierror .ne. MPI_SUCCESS) then
c$$$          write(179,*) 'mpi_allreduce_int fail',ierror
c$$$          call flush(179)
c$$$        endif
c$$$  CALL ICOPY(LENGTH,BUFF,1,X(LOCX),1)
        do j=1,length
          x(locx+j-1)=buff(j)
        enddo
        LOCX = LOCX + MXBUFF
      enddo
      deallocate(buff)
      END

      SUBROUTINE mpi_allreduce_dble(X,LENX,OP,comm)
      implicit none
      INCLUDE 'mpif.h'
      integer,intent(in) :: lenx,op,comm
      double precision, intent(inout) :: x(lenx)
      integer :: mxbuff
      double precision, allocatable :: buff(:)
      integer :: i,j,npass,length,locx,ierror

      mxbuff=4096
      allocate(buff(mxbuff))
      NPASS = (LENX-1)/MXBUFF + 1
      LENGTH = MXBUFF
      LOCX = 1
      DO I=1,NPASS
        IF(I.EQ.NPASS) LENGTH = LENX - MXBUFF*(NPASS-1)
        CALL MPI_ALLREDUCE(X(LOCX),BUFF,LENGTH,MPI_double_precision,
     $    OP,comm,IERROR)
c$$$        if (ierror .ne. MPI_SUCCESS) then
c$$$          write(179,*) 'mpi_allreduce_dble fail',ierror
c$$$          call flush(179)
c$$$        endif
        CALL DCOPY(LENGTH,BUFF,1,X(LOCX),1)
c$$$  do j=1,length
c$$$  x(locx+j-1)=buff(j)
c$$$  enddo
        LOCX = LOCX + MXBUFF
      enddo
      deallocate(buff)
      END

      SUBROUTINE mpi_allreduce_dcmplx(X,LENX,OP,comm)
      implicit none
      INCLUDE 'mpif.h'
      integer,intent(in) :: lenx,op,comm
      complex*16, intent(inout) :: x(lenx)
      integer :: mxbuff
      complex*16,allocatable :: buff(:)
      integer :: i,j,npass,length,locx,ierror


      mxbuff=2048
      allocate(buff(mxbuff))
      NPASS = (LENX-1)/MXBUFF + 1
      LENGTH = MXBUFF
      LOCX = 1
      DO I=1,NPASS
        IF(I.EQ.NPASS) LENGTH = LENX - MXBUFF*(NPASS-1)
        CALL MPI_ALLREDUCE(X(LOCX),BUFF,LENGTH,MPI_double_complex,
     $    OP,comm,IERROR)
c$$$        if (ierror .ne. MPI_SUCCESS) then
c$$$          write(179,*) 'mpi_allreduce_dble fail',ierror
c$$$          call flush(179)
c$$$        endif
        CALL ZCOPY(LENGTH,BUFF,1,X(LOCX),1)
c$$$  do j=1,length
c$$$  x(locx+j-1)=buff(j)
c$$$  enddo
        LOCX = LOCX + MXBUFF
      enddo
      deallocate(buff)
      END


      subroutine dcmplx_matinv(mat, invmat, dim,dimmax)
      implicit none
      integer, intent(in) :: dim, dimmax
      complex*16, intent(in) :: mat(dimmax, dimmax)
      complex*16, intent(out) :: invmat(dimmax,dimmax)
      
      integer :: ipiv(dim),info
      complex*16 :: wrk(dim**2)

      invmat=mat
      call zgetrf(dim,dim,invmat,
     $  dimmax,ipiv,info)
c$$$      if (info .ne. 0) then
c$$$        write(179, *) 'error in zgetrf', info
c$$$      endif
      call zgetri(dim,invmat,dimmax,
     $  ipiv,wrk,dim**2,info)          
c$$$      if (info .ne. 0) then
c$$$        write(179, *) 'error in zgetri', info
c$$$      endif
      end

      DOUBLE PRECISION function iFactorial(j)
      IMPLICIT NONE
      INTEGER, intent(in) :: j
      INTEGER :: i
      DOUBLE PRECISION :: x
      if (j<0) print *,
     $  "iFactorial defined only for non-negative numbers!"
      x=1
      iFactorial = x
      if (j.eq.1) return
      DO i=2,j
        x = x*i
      END DO
      iFactorial = x
      return
      end function iFactorial
c$$$
      DOUBLE PRECISION function dFactorial(x)
      IMPLICIT NONE
      DOUBLE PRECISION, intent(in) :: x
      DOUBLE PRECISION, PARAMETER :: spi2 = 0.8862269254527579
      DOUBLE PRECISION :: y, r
      r=1
      y=x
      DO WHILE(y.gt.1.0)
        r= r * y
        y= y -1.
      ENDDO
      IF (abs(y-0.5).LT.1e-10) r = r*spi2
      dFactorial = r
      return
      END function dFactorial
c$$$  
      DOUBLE PRECISION function mone(i)
      INTEGER, intent(in) :: i
      mone = 1 - 2*MOD(abs(i),2)
      return
      end function mone
c$$$
      DOUBLE PRECISION function Ddelta(j1, j2, j)
      IMPLICIT NONE
      DOUBLE PRECISION, intent(in) :: j1, j2, j
!     function calls
      DOUBLE PRECISION :: dFactorial
      Ddelta = sqrt(dFactorial(j1+j2-j)*dFactorial(j1-j2+j)
     $  *dFactorial(-j1+j2+j)/dFactorial(j1+j2+j+1))
      return
      END function Ddelta

      DOUBLE PRECISION function f3j_int(j1, m1, j2, m2, j3, m3)
c$$$  wigner 3j symbol 
      IMPLICIT NONE
      integer, intent(in) :: j1, j2, j3, m1, m2, m3
      INTEGER            :: tmin, tmax, t
      DOUBLE PRECISION             :: sum, v1, v2, dn,
     $  j1d,j2d,j3d,m1d,m2d,m3d
!     function calls
      DOUBLE PRECISION,external :: dFactorial
      DOUBLE PRECISION,external :: iFactorial
      DOUBLE PRECISION,external :: Ddelta
      DOUBLE PRECISION,external :: mone
      f3j_int=0
      j1d=j1*1.0d0
      j2d=j2*1.0d0
      j3d=j3*1.0d0
      m1d=m1*1.0d0
      m2d=m2*1.0d0
      m3d=m3*1.0d0
      IF (abs(m1d+m2d+m3d) .GT. 1e-10) return
      IF (abs(j1d-j2d)-1e-14.GT.j3d.OR.j3d.GT.j1d+j2d+1e-14) return
      if (abs(m1d).GT.j1d.OR.abs(m2d).GT.j2d.OR.abs(m3d).GT.j3d)
     $  return
      tmin = INT(max(max(0.0,j2d-j3d-m1d),j1d-j3d+m2d)+1e-14)
      tmax = INT(min(min(j1d+j2d-j3d,j1d-m1d),j2d+m2d)+1e-14)
      sum=0
      DO t=tmin, tmax
        v1 = dFactorial(j3d-j2d+m1d+t)*dFactorial(j3d-j1d-m2d+t)
        v2 = dFactorial(j1d+j2d-j3d-t)*dFactorial(j1d-m1d-t)
     $    *dFactorial(j2d+m2d-t)
        sum = sum + mone(t)/(iFactorial(t)*v1*v2)
      END DO
      dn = dFactorial(j1d+m1d)*dFactorial(j1d-m1d)*dFactorial(j2d+m2d)
     $  *dFactorial(j2d-m2d)*dFactorial(j3d+m3d)*dFactorial(j3d-m3d)
      f3j_int = mone(INT(j1d-j2d-m3d))*Ddelta(j1d,j2d,j3d)*sqrt(dn)*sum
      return
      END function f3j_int      

      subroutine cal_coulomb_matrix_from_slater(lval, ff, mat)
      implicit none
      integer, intent(in) :: lval
      double precision, intent(in) :: ff(0:3)
      complex*16, intent(out) :: mat
     $  (2*lval+1,2*lval+1,2*lval+1,2*lval+1)      

      integer :: i1,i2,i3,i4,kk0,kk,mi1,mi3,qq,mq
      double precision, external :: f3j_int

      mat=0.0d0

      do kk=0,2*lval,2
        do qq=-kk, kk
          do i1=-lval, lval
            do i2=-lval, lval
              do i3=-lval, lval
                do i4=-lval, lval                  
                  mq=-qq
                  mi1=-i1
                  mi3=-i3
                  kk0=kk/2
                  mat(i1+lval+1,i2+lval+1,i3+lval+1,i4+lval+1)
     &              =mat(i1+lval+1,i2+lval+1,i3+lval+1,i4+lval+1)
     $              +(2.0d0*lval+1.0d0)**2*ff(kk0)
     $              *f3j_int(lval,0,kk,0,lval,0)**2
     $              *f3j_int(lval,mi1,kk,qq,lval,i2)
     $              *f3j_int(lval,mi3,kk,mq,lval,i4)
     $              *(-1)**(i1+i4)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      end


!     Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
!     This file is modified to be suitable for complex variables by sangkook choi
!     distributed under the terms of the GNU Lesser General Public
!     License. See the file COPYING for license details.
      
!     BOP
!     !ROUTINE: fderiv_dcmplx
!     !INTERFACE:
      subroutine fderiv_dcmplx(m,n,x,f,g)
!     !INPUT/OUTPUT PARAMETERS:
!     m : order of derivative (in,integer)
!     n : number of points (in,integer)
!     x : abscissa array (in,real(n))
!     f : function array (in,complex(n))
!     g : (anti-)derivative of f (out,complex(n))
!     !DESCRIPTION:
!     Given function $f$ defined on a set of points $x_i$ then if $m\ge 0$ this
!     routine computes the $m$th derivative of $f$ at each point. If $m<0$ the
!     anti-derivative of $f$ given by
!     $$ g(x_i)=\int_{x_1}^{x_i} f(x)\,dx $$
!     is calculated. If $m=-1$ then an accurate integral is computed by fitting
!     the function to a clamped cubic spline_dcmplx. When $m=-3$ the fast but low
!     accuracy trapezoidal integration method is used. Simpson's integration,
!     which is slower but more accurate than the trapezoidal method, is used if
!     $m=-2$.
!     
!     !REVISION HISTORY:
!     Created May 2002 (JKD)
!     EOP
!     BOC
      implicit none
!     arguments
      integer, intent(in) :: m,n
      double precision, intent(in) :: x(n)
      complex*16, intent(in) :: f(n)
      complex*16, intent(out) :: g(n)
!     local variables
      integer :: i
      double precision ::  x0,x1,x2,dx
!     automatic arrays
      complex*16 :: cf(3,n)
c$$$  if (n.le.0) then
c$$$  write(*,*)
c$$$  write(*,'("Error(fderiv_dcmplx): invalid number of points : ",I8)') n
c$$$  write(*,*)
c$$$  stop
c$$$  end if
      select case(m)
      case(-3)
!     low accuracy trapezoidal integration
        g(1)=0.d0
        do i=1,n-1
          g(i+1)=g(i)+0.5d0*(x(i+1)-x(i))*(f(i+1)+f(i))
        end do
        return
      case(-2)
!     medium accuracy Simpson integration
        g(1)=0.d0
        do i=1,n-2
          x0=x(i)
          x1=x(i+1)
          x2=x(i+2)
          g(i+1)=g(i)+(x0-x1)*(f(i+2)*(x0-x1)**2
     $      +f(i+1)*(x2-x0)*(x0+2.d0*x1-3.d0*x2)
     $      +f(i)*(x2-x1)*(2.d0*x0+x1-3.d0*x2))/(6.d0*(x0-x2)*(x1-x2))
        end do
        x0=x(n)
        x1=x(n-1)
        x2=x(n-2)
        g(n)=g(n-1)+(x1-x0)*(f(n-2)*(x1-x0)**2
     $    +f(n)*(x1-x2)*(3.d0*x2-x1-2.d0*x0)
     $    +f(n-1)*(x0-x2)*(3.d0*x2-2.d0*x1-x0))/(6.d0*(x2-x1)*(x2-x0))
        return
      case(0)
        g(:)=f(:)
        return
      case(4:)
        g(:)=0.d0
        return
      end select
!     high accuracy integration/differentiation from spline_dcmplx interpolation
      call spline_dcmplx(n,x,f,cf)
      select case(m)
      case(:-1)
        g(1)=0.d0
        do i=1,n-1
          dx=x(i+1)-x(i)
          g(i+1)=g(i)+(((0.25d0*cf(3,i)*dx
     $      +0.3333333333333333333d0*cf(2,i))*dx
     $      +0.5d0*cf(1,i))*dx+f(i))*dx
        end do
      case(1)
        g(:)=cf(1,:)
      case(2)
        g(:)=2.d0*cf(2,:)
      case(3)
        g(:)=6.d0*cf(3,:)
      end select
      return
      end subroutine
      
      
!     Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
!     This file is distributed under the terms of the GNU Lesser General Public
!     This file is modified to be suitable for complex variables by sangkook choi  !     License. See the file COPYING for license details.
      
!     BOP
!     !ROUTINE: spline_dcmplx
!     !INTERFACE:
      subroutine spline_dcmplx(n,x,f,cf)
!     !INPUT/OUTPUT PARAMETERS:
!     n  : number of points (in,integer)
!     x  : abscissa array (in,real(n))
!     f  : input data array (in,complex*16(n))
!     cf : cubic spline_dcmplx coefficients (out,complex*16(3,n))
!     !DESCRIPTION:
!     Calculates the coefficients of a cubic spline_dcmplx fitted to input data. In other
!     words, given a set of data points $f_i$ defined at $x_i$, where
!     $i=1\ldots n$, the coefficients $c_j^i$ are determined such that
!     $$ y_i(x)=f_i+c_1^i(x-x_i)+c_2^i(x-x_i)^2+c_3^i(x-x_i)^3, $$
!     is the interpolating function for $x\in[x_i,x_{i+1})$. The coefficients are
!     determined piecewise by fitting a cubic polynomial to adjacent points.
!     
!     !REVISION HISTORY:
!     Created November 2011 (JKD)
!     EOP
!     BOC
      implicit none
!     arguments
      integer, intent(in) :: n
      double precision, intent(in) :: x(n)
      complex*16, intent(in) :: f(n)
      complex*16, intent(out) :: cf(3,n)
!     local variables
      integer ::  i
      double precision :: x0,x1,x2,x3
      complex*16 :: y0,y1,y2,y3,c1,c2,c3,t0,t1,t2,t3,t4,t5,t6
c$$$  if (n.le.0) then
c$$$  write(*,*)
c$$$  write(*,'("Error(spline_dcmplx): n <= 0 : ",I8)') n
c$$$  write(*,*)
c$$$  stop
c$$$  end if
      if (n.eq.1) then
        cf(:,1)=0.d0
        return
      end if
      if (n.eq.2) then
        cf(1,1)=(f(2)-f(1))/(x(2)-x(1))
        cf(2:3,1)=0.d0
        cf(1,2)=cf(1,1)
        cf(2:3,2)=0.d0
        return
      end if
      if (n.eq.3) then
        x0=x(1)
        x1=x(2)-x0
        x2=x(3)-x0
        y0=f(1)
        y1=f(2)-y0
        y2=f(3)-y0
        t0=1.d0/(x1*x2*(x2-x1))
        t1=x1*y2
        t2=x2*y1
        c1=t0*(x2*t2-x1*t1)
        c2=t0*(t1-t2)
        cf(1,1)=c1
        cf(2,1)=c2
        cf(3,1)=0.d0
        t3=2.d0*c2
        cf(1,2)=c1+t3*x1
        cf(2,2)=c2
        cf(3,2)=0.d0
        cf(1,3)=c1+t3*x2
        cf(2,3)=c2
        cf(3,3)=0.d0
        return
      end if
      y0=f(1)
      y1=f(2)-y0
      y2=f(3)-y0
      y3=f(4)-y0
      x0=x(1)
      x1=x(2)-x0
      x2=x(3)-x0
      x3=x(4)-x0
      t0=1.d0/(x1*x2*x3*(x1-x2)*(x1-x3)*(x2-x3))
      t1=x1*x2*y3
      t2=x2*x3*y1
      t3=x3*x1*y2
      t4=x1**2
      t5=x2**2
      t6=x3**2
      y1=t3*t6-t1*t5
      y3=t2*t5-t3*t4
      y2=t1*t4-t2*t6
      c1=t0*(x1*y1+x2*y2+x3*y3)
      c2=-t0*(y1+y2+y3)
      c3=t0*(t1*(x1-x2)+t2*(x2-x3)+t3*(x3-x1))
      cf(1,1)=c1
      cf(2,1)=c2
      cf(3,1)=c3
      cf(1,2)=c1+2.d0*c2*x1+3.d0*c3*t4
      cf(2,2)=c2+3.d0*c3*x1
      cf(3,2)=c3
      if (n.eq.4) then
        cf(1,3)=c1+2.d0*c2*x2+3.d0*c3*t5
        cf(2,3)=c2+3.d0*c3*x2
        cf(3,3)=c3
        cf(1,4)=c1+2.d0*c2*x3+3.d0*c3*t6
        cf(2,4)=c2+3.d0*c3*x3
        cf(3,4)=c3
        return
      end if
      do i=3,n-2
        y0=f(i)
        y1=f(i-1)-y0
        y2=f(i+1)-y0
        y3=f(i+2)-y0
        x0=x(i)
        x1=x(i-1)-x0
        x2=x(i+1)-x0
        x3=x(i+2)-x0
        t1=x1*x2*y3
        t2=x2*x3*y1
        t3=x3*x1*y2
        t0=1.d0/(x1*x2*x3*(x1-x2)*(x1-x3)*(x2-x3))
        c3=t0*(t1*(x1-x2)+t2*(x2-x3)+t3*(x3-x1))
        t4=x1**2
        t5=x2**2
        t6=x3**2
        y1=t3*t6-t1*t5
        y2=t1*t4-t2*t6
        y3=t2*t5-t3*t4
        cf(1,i)=t0*(x1*y1+x2*y2+x3*y3)
        cf(2,i)=-t0*(y1+y2+y3)
        cf(3,i)=c3
      end do
      c1=cf(1,n-2)
      c2=cf(2,n-2)
      c3=cf(3,n-2)
      cf(1,n-1)=c1+2.d0*c2*x2+3.d0*c3*t5
      cf(2,n-1)=c2+3.d0*c3*x2
      cf(3,n-1)=c3
      cf(1,n)=c1+2.d0*c2*x3+3.d0*c3*t6
      cf(2,n)=c2+3.d0*c3*x3
      cf(3,n)=c3
      return
      end subroutine

      subroutine cal_rotmat_mlms2jmj(ll, rotmat)
      
      implicit none
      integer, intent(in) :: ll
      complex*16, intent(out) :: rotmat(4*ll+2,4*ll+2)
      
      integer:: jj2,indms,indml,indmlms,indi,indmj,indjmj,
     $  iorb,jorb
      double precision :: jj,ulmu,fac1,fac2

!     default order:
!     mj=> fastest index, starting from smallest,
!     i=>next fastest, starting from smallest,

!     ml=> fastest index, starting from smallest,
!     ms=>next fastest, starting from smallest,            
      
      do indms=-1,1,2      
        do indml=-ll, ll
          indmlms=indml+ll+1+(indms+1)/2*(2*ll+1)
          do indi=-1,1,2
            jj2=2*ll+indi
            do indmj=-jj2, jj2, 2
              indjmj=(indmj+jj2)/2+1+(indi+1)/2*(jj2-1)
              ulmu=dble(indmj)/2.0d0/(dble(ll)+0.5d0)
              if (indi .eq. indms) then
                fac1=1.0d0
              else
                fac1=-1.0d0
              endif
              fac2=1.0d0
              if (indi .eq. -1 .and. indms .eq. 1) then
                fac2=-1.0d0
              endif
              if (indmj .eq. indms+indml*2) then
                rotmat(indjmj,indmlms)=1.0d0/dsqrt(2.0d0)
     $            *dsqrt(1+ulmu*fac1)*fac2
              endif
            enddo
          enddo
        enddo
      enddo


      end      


      subroutine cal_rotmat_cmplx2real(norb, rotmat)
      implicit none
      integer, intent(in) :: norb
      complex*16, intent(out) :: rotmat(norb,norb)

      integer :: iorb,jorb

      rotmat=0.0d0

      rotmat(norb/2+1,norb/2+1)=1.0d0       
      do iorb=1, norb/2
        rotmat(iorb+norb/2+1,iorb+norb/2+1)
     $    =1.0d0/dsqrt(2.0d0)*(-1)**iorb
        rotmat(iorb+norb/2+1,-iorb+norb/2+1)
     $    =1.0d0/dsqrt(2.0d0)
        rotmat(-iorb+norb/2+1,iorb+norb/2+1)
     $    =1.0d0/dsqrt(2.0d0)*dcmplx(0.0d0, 1.0d0)*(-1)**(iorb+1)
        rotmat(-iorb+norb/2+1,-iorb+norb/2+1)
     $    =1.0d0/dsqrt(2.0d0)*dcmplx(0.0d0, 1.0d0)
      enddo
      end            


      subroutine rotate_umatrix(norb,umatrix, rotmat,daggerflag)
      implicit none

      integer, intent(in) :: norb,daggerflag
      complex*16, intent(inout) :: umatrix(norb,norb,norb,norb)
      complex*16, intent(in) :: rotmat(norb,norb)


      integer :: iorb1,iorb2,jorb1,jorb2,korb1,korb2,lorb1,lorb2,
     $  iorb,jorb,korb,lorb 
      complex*16 :: umatrix_temp(norb,norb,norb,norb),
     $  rotmat_temp(norb,norb)
      
      if (daggerflag .eq. 1) then
        rotmat_temp=dconjg(transpose(rotmat))
      else
        rotmat_temp=rotmat
      endif
      
      
      umatrix_temp=umatrix
      umatrix=0.0d0

      do lorb2=1, norb
        do korb2=1, norb
          do jorb2=1, norb
            do iorb2=1, norb            
              do iorb1=1, norb
                umatrix(iorb1,jorb2,korb2,lorb2)
     $            =umatrix(iorb1,jorb2,korb2,lorb2)
     $            +rotmat_temp(iorb1,iorb2)
     $            *umatrix_temp(iorb2,jorb2,korb2,lorb2)
              enddo
              
            enddo
          enddo
        enddo
      enddo

      umatrix_temp=umatrix
      umatrix=0.0d0

      do lorb2=1, norb
        do korb2=1, norb
          do jorb2=1, norb
            do iorb2=1, norb            
              do jorb1=1, norb            
                umatrix(iorb2,jorb1,korb2,lorb2)
     $            =umatrix(iorb2,jorb1,korb2,lorb2)
     $            +dconjg(rotmat_temp(jorb1,jorb2))
     $            *umatrix_temp(iorb2,jorb2,korb2,lorb2)
              enddo

            enddo
          enddo
        enddo
      enddo


      umatrix_temp=umatrix      
      umatrix=0.0d0
      do lorb2=1, norb
        do korb2=1, norb
          do jorb2=1, norb
            do iorb2=1, norb            
              do korb1=1, norb                        
                umatrix(iorb2,jorb2,korb1,lorb2)
     $            =umatrix(iorb2,jorb2,korb1,lorb2)
     $            +rotmat_temp(korb1,korb2)
     $            *umatrix_temp(iorb2,jorb2,korb2,lorb2)
              enddo

            enddo
          enddo
        enddo
      enddo


      umatrix_temp=umatrix      
      umatrix=0.0d0
      do lorb2=1, norb
        do korb2=1, norb
          do jorb2=1, norb
            do iorb2=1, norb            
              do lorb1=1, norb                        
                umatrix(iorb2,jorb2,korb2,lorb1)
     $            =umatrix(iorb2,jorb2,korb2,lorb1)
     $            +dconjg(rotmat_temp(lorb1,lorb2))
     $            *umatrix_temp(iorb2,jorb2,korb2,lorb2)
              enddo

            enddo
          enddo
        enddo
      enddo                  


      end

      subroutine cal_u_matrix(norb,slaterf,umatrix)
      use comdc_mod      
      implicit none
      integer, intent(in) :: norb
      double precision,intent(in) :: slaterf(0:3)
      complex*16,intent(out) :: umatrix(norb,norb,norb,norb)
      

      integer :: kk, ll, lval,norb_2
      complex*16, allocatable :: 
     $  rotmat_cmplx2real(:,:),
     $  rotmat_new(:,:),
     $  umatrix_temp(:,:,:,:),
     $  rotmat_mlms2jmj(:,:)      

      umatrix=0.0d0            

      if  (is_spinorbit .eq. 0)then
        lval=(norb-1)/2
        call cal_coulomb_matrix_from_slater(lval,slaterf,umatrix)
          
        allocate(rotmat_cmplx2real(norb,norb))
        rotmat_cmplx2real=0.0d0
        call cal_rotmat_cmplx2real(norb,rotmat_cmplx2real)
          

        if (is_trans_dc) then
          allocate(rotmat_new(norb,norb))
          rotmat_new=0.0d0
          call zgemm('n','n',norb,norb,norb,(1.d0,0.d0),
     $      trans_dc,norb,rotmat_cmplx2real,norb,
     $      (0.d0,0.d0),rotmat_new,norb)
          call rotate_umatrix(norb,umatrix,rotmat_new,0)
          deallocate(rotmat_new)          
        else
          call rotate_umatrix(norb,umatrix,rotmat_cmplx2real,0)
        endif

          
        deallocate(rotmat_cmplx2real)

      else
!     for wannier90 convension
c$$$  !        j*2     5 5 5 5 5 5  7 7 7 7  7  7  7  7
c$$$  !        mj*2   -5-3-1 1 3 5 -7-5-3-1  1  3  5  7
c$$$  ! default index  1 2 3 4 5 6  7 8 9 10 11 12 13 14        
        lval=(norb/2-1)/2
        norb_2=norb/2
        allocate(umatrix_temp(norb_2,norb_2,norb_2,norb_2))
        call cal_coulomb_matrix_from_slater(lval,slaterf,umatrix_temp)
        do kk=1,2
          do ll=1,2
            umatrix(
     $        (1+(kk-1)*norb_2):(norb_2+(kk-1)*norb_2),
     $        (1+(kk-1)*norb_2):(norb_2+(kk-1)*norb_2),
     $        (1+(ll-1)*norb_2):(norb_2+(ll-1)*norb_2),
     $        (1+(ll-1)*norb_2):(norb_2+(ll-1)*norb_2)
     $        )
     $        =umatrix_temp
          enddo
        enddo          
        
        allocate(rotmat_mlms2jmj(norb,norb))
        rotmat_mlms2jmj=0.0d0
        
        call cal_rotmat_mlms2jmj(lval, rotmat_mlms2jmj)
        if (is_trans_dc) then        
          allocate(rotmat_new(norb,norb))
          rotmat_new=0.0d0
          call zgemm('n','n',norb,norb,
     $      norb,(1.d0,0.d0),
     $      trans_dc,norb,
     $      rotmat_mlms2jmj,norb,
     $      (0.d0,0.d0),
     $      rotmat_new,norb)        
          call rotate_umatrix(norb, umatrix, rotmat_new,0)
          deallocate(rotmat_new)
        else
          call rotate_umatrix(norb, umatrix, rotmat_mlms2jmj,0)
        endif
        deallocate(umatrix_temp)
        deallocate(rotmat_mlms2jmj)

        
      endif
      end


      SUBROUTINE ENDING
      implicit none
      include "mpif.h"
      integer :: ierror

c$$$      write(179,*) 'call ending'
c$$$      call flush(179)
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
      CALL MPI_FINALIZE(IERROR)
      stop      
      END                                                                     
