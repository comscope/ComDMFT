      subroutine wannier_interface
      
      use comwann_mod      

      implicit none             
      include "mpif.h"

      integer :: ispin,ii,jj,kk,ib,indexb,
     $  tempint,ind_k,ibnd,k,jbnd,ind,iq0,irk,iwan,
     $  k0,i,aa,bb,cc,
     $  indloc1(1),indloc2(1),indloc3(1),indloc4(1),indloc5(1),
     $  indtmp,j,m,tempintarr(3),
     $  jwan,ndiv_fi(3),
     $  me1
!     $  wanntot
!     $  ndim,nb,
!     $  jbnd1,jbnd2,kvec(3),
!     $  ndim_arpes_kpath(nproc),
!     $  ndim_arpes_kpath_max,
!     $  ir,
!     $  ndim_kk_fi(nproc),ndim_kk_fi_max,
!     $  tempk(3)
      double precision :: bvec(3),tempvec(3),tempvec2(3),
!     $  kqvec(3),
     $  distance,temp_distance,
     $  nelectron
!     $  rho_projected
!      complex*16 ::htemp
!      character*300 ::  kstar_ch  
!      logical :: file_exist
!      character*20 :: ch_atom 
!      integer, allocatable :: ipiv(:)      
      double precision, allocatable :: tempdble(:),
     $  tempdble2(:)
!     $  kpath(:,:),
!     $  kgrid_fi(:,:)
!     $  dummy_dble1(:,:),dummy_dble2(:,:)
!     $  arpes_kpath(:,:)

c$$$      complex*16, allocatable :: 
!     $  g_out4(:,:), g_out5(:,:,:),
!     $  rhomat(:,:),
!     $  zmat(:,:,:),zmat0(:,:,:),
!     $  work(:,:),
!     $  dummy_dcmplx(:,:,:),
c$$$     $  wan_rgrid(:,:,:,:,:)
!     $  wan_rgrid_fine(:,:,:,:)
!     $  momentum_tmp(:,:,:)
      character*300 :: ch

!---------------------------------------------------------------------
!      ai=dcmplx(0.0d0,1.0d0)
      me1=me+1
      if (magn .gt. 1) then
c        write(179,*) 'ComWann has not been tested for '//
c     $    'the case of magn==2'        
        if (maswrk) then
          write(993,*) 'ComWann has not been tested for '//
     $      'the case of magn==2'
        endif
        call ending
      endif

      if (maswrk) then
        if (nspin .eq. 1) then
          write(iun,*) "!-----------------------------------!"        
          write(iun,*) "       Total number of spin: 1"
          write(iun,*) "!-----------------------------------!"
        elseif (nspin .eq. 2) then
          write(iun,*) "!-----------------------------------!"        
          write(iun,*) "       Total number of spin: 2"
          write(iun,*) "!-----------------------------------!"
        end if
      endif



      do ispin=1, nspin

        if (nspin .eq. 1) then
          write(seed_name,'(a)') 'wannier'
        elseif (nspin .eq. 2) then
          if (ispin .eq. 1) then
            if (maswrk) then
              write(iun,*) "!-----------------------------------!"
              write(iun,*) "Information for the up spin electron"
              write(iun,*) "!-----------------------------------!"
            endif
            write(seed_name,'(a)') 'wannier_up'
            
          else
            if (maswrk) then
              write(iun,*) "!-------------------------------------!"
              write(iun,*) "Information for the down spin electron"
              write(iun,*) "!-------------------------------------!"
            endif
            write(seed_name,'(a)') 'wannier_dn'
          end if
        endif
        
        num_nnmax=12

        
        allocate(proj_site(3,nbndf))
        proj_site=0.0d0
        allocate(proj_z(3,nbndf))
        proj_z=0.0d0
        allocate(proj_x(3,nbndf))
        proj_x=0.0d0
        allocate(proj_zona(nbndf))
        proj_zona=0.0d0
        
        allocate(nnlist(nqdiv,num_nnmax))
        nnlist=0
        allocate(nncell(3,nqdiv,num_nnmax))
        nncell=0
        allocate(proj_l(nbndf))
        proj_l=0
        allocate(proj_m(nbndf))
        proj_m=0
        allocate(proj_radial(nbndf))
        proj_radial=0
        allocate(exclude_bands(nbndf))
        exclude_bands=0
        allocate(proj_s(nbndf))
        proj_s=0
        allocate(proj_s_qaxisx(3,nbndf))
        proj_s_qaxisx=0.0d0
        
        allocate(exclude_bands_table(nbndf))
        exclude_bands_table=0

        allocate(atom_xaxis(3,natom))
        allocate(atom_zaxis(3,natom))        
        
        call read_local_axis
        
        call cal_num_bands(ispin)

        allocate(eigenvalues(num_bands,nqdiv))
        eigenvalues=0.0d0

        do irk = 1, nqdiv
          k0=i_kref(irk)
          do ibnd = 1, num_bands
            jbnd=ibnd+ind_low
            if (ubi_num .ge. 4) then
              eigenvalues(ibnd,irk)=
     $          (e_qp(jbnd,k0,ispin)-chem_pot)*evolt/2.0d0
            else
              eigenvalues(ibnd,irk)=
     $          (e_bnd(jbnd,k0,ispin)-chem_pot)*evolt/2.0d0
            endif

          enddo
        enddo

        write(ch, '(a, a)') trim(seed_name), '.eig'

        if (maswrk) then
          open (unit = 990, file = trim(ch), form = 'formatted')
          do irk = 1, nqdiv
            do ibnd = 1, num_bands
              write(990,'(2i8,2x,F17.12)') ibnd, irk,
     $          eigenvalues(ibnd,irk)
            enddo
          enddo
          close(unit = 990)
        endif        

        allocate(z_wan_bnd(nfun,num_bands,npnt))
        allocate(ev_wan_bnd(nbasmpw,num_bands,npnt))
        z_wan_bnd=(0.d0,0.d0)
        ev_wan_bnd=(0.d0,0.d0)
        

        do ind_k=1,ndim_k(me1)
          k=n_mpi_k(me1)+ind_k
          do ibnd=1, num_bands
            jbnd=ibnd+ind_low
            if (ubi_num .eq. 4) then
              z_wan_bnd(:,ibnd,k)=z_qp(:,jbnd,ind_k,ispin)
              ev_wan_bnd(:,ibnd,k)=ev_qp(:,jbnd,ind_k,ispin)
!            elseif (ubi_num .gt. 4) then
!              z_wan_bnd(:,ibnd,k)=z_qp_2(:,jbnd,ind_k,ispin)
!              ev_wan_bnd(:,ibnd,k)=ev_qp_2(:,jbnd,ind_k,ispin)
            else
              z_wan_bnd(:,ibnd,k)=z_bnd(:,jbnd,ind_k,ispin)
              ev_wan_bnd(:,ibnd,k)=ev_bnd(:,jbnd,ind_k,ispin)
            endif
          enddo
        enddo
        if(nproc/=1) then
          call mpi_allreduce_dcmplx
     $      (z_wan_bnd,nfun*num_bands*npnt,mpi_sum,mpi_comm_world)
          call mpi_allreduce_dcmplx
     $      (ev_wan_bnd,nbasmpw*num_bands*npnt,mpi_sum,mpi_comm_world)
        endif        

        call count_nbnd_in_froz_win        

        call wannier_orbital(ispin)        


        call read_wann_inip
        call write_wannier_win(ispin)

        if (maswrk) then
          write(iun,*) "Calling wannier_setup"
        endif
        
        gamma_only = .false.
        
        if (maswrk) then
          call wannier_setup(trim(seed_name),ndiv,nqdiv,
     $      transpose(rbas)*par*bohr,transpose(gbas)*2.0d0*pi/par/bohr,
     $      kpt_latt,nbndf,
     $      natom,atom_symbols,tau*par*bohr,gamma_only,.false.,
     $      nntot,nnlist,nncell,num_bands,num_wann,
     $      proj_site,proj_l,proj_m,proj_radial,proj_z,
     $      proj_x,proj_zona,exclude_bands,proj_s,proj_s_qaxisx)
        endif

        call mpi_bcast(nntot,1,mpi_integer,0,mpi_comm_world,mpierr)
        call mpi_bcast(nnlist,nqdiv*num_nnmax,mpi_integer,0,
     $    mpi_comm_world,mpierr)
        call mpi_bcast(nncell,3*nqdiv*num_nnmax,mpi_integer,0,
     $    mpi_comm_world,mpierr)
        call mpi_bcast(num_bands,1,mpi_integer,0,mpi_comm_world,mpierr)
        call mpi_bcast(num_wann,1,mpi_integer,0,mpi_comm_world,mpierr)
        call mpi_bcast(proj_site,3*nbndf,mpi_double_precision,0,
     $    mpi_comm_world,mpierr)
        call mpi_bcast(proj_l,nbndf,mpi_integer,0,mpi_comm_world,mpierr)
        call mpi_bcast(proj_m,nbndf,mpi_integer,0,mpi_comm_world,mpierr)
        call mpi_bcast(proj_radial,nbndf,mpi_integer,0,
     $    mpi_comm_world,mpierr)
        call mpi_bcast(proj_z,3*nbndf,mpi_double_precision,0,
     $    mpi_comm_world,mpierr)
        call mpi_bcast(proj_x,3*nbndf,mpi_double_precision,0,
     $    mpi_comm_world,mpierr)
        call mpi_bcast(proj_zona,nbndf,mpi_double_precision,0,
     $    mpi_comm_world,mpierr)
        call mpi_bcast(exclude_bands,nbndf,mpi_integer,0,
     $    mpi_comm_world,mpierr)
        call mpi_bcast(proj_s,nbndf,mpi_integer,0,mpi_comm_world,mpierr)
        call mpi_bcast(proj_s_qaxisx,3*nbndf,mpi_double_precision,0,
     $    mpi_comm_world,mpierr)

        if ((mod(num_bands, 2) .eq. 1) .and. (irel .eq. 2)) then
c          write(179,*) 'num_bands should be even number if irel .eq. 2'
          write(993,*) 'num_bands should be even number if irel .eq. 2'
          call ending
        endif

!        allocate(ldim_wan(num_wann))
!        allocate(rdim_wan(3,num_wann))
!        allocate(lval_wan(3,num_wann))
!        allocate(radial_coeff_wan(100,3,num_wann))

        allocate(tempdble(num_wann))
        allocate(tempdble2(num_wann))
        
        allocate(include_bands(num_bands)) 
        include_bands=0

        allocate(a_matrix(num_bands,num_wann,nqdiv))
        a_matrix=0.0d0
        allocate(M_Matrix(num_bands,num_bands,nntot,nqdiv))
        m_matrix=0.0d0
        allocate(atom_wan(num_wann))
        atom_wan=0
        allocate(atom_wan_distance(3,num_wann))
        atom_wan_distance=0.0d0



        allocate(U_matrix(num_wann,num_wann,nqdiv))
        u_matrix=0.0d0
        allocate(v_matrix(num_bands,num_wann,nqdiv))
        v_matrix=0.0d0
        allocate(U_matrix_opt(num_bands,num_wann,nqdiv))
        u_matrix_opt=0.0d0
        allocate(wan_g_overlap(num_wann,num_wann))
        wan_g_overlap=0.0d0
        allocate(lwindow(num_bands,nqdiv))
        lwindow=.false.
        allocate(wann_centers(3,num_wann))
        wann_centers=0.0d0
        allocate(wann_spreads(num_wann))
        wann_spreads=0.0d0

        nexband=0
        do ib=1,nbndf
          indexb=exclude_bands(ib)
          if (indexb .gt. nbndf .or. indexb .lt. 0) then
            write(*,*) "something wrong in exclude_bands", indexb, ib
          elseif (indexb .ne. 0) then
            nexband=nexband+1
            exclude_bands_table(indexb)=1
          endif
        enddo
        

        
        if (nbndf  .eq. nexband+num_bands) then
          tempint=0
          do ib=1,nbndf
            if (exclude_bands_table(ib) .eq. 0) then
              tempint=tempint+1
              include_bands(tempint)=ib
            end if
          enddo
        end if

        if (maswrk) then
          write(iun,*) 'include_bands', num_bands
          do ibnd=1, num_bands
            write(iun,*) include_bands(ibnd)
          enddo
        endif

        


        if (maswrk) then
          write(iun,*) "nntot",nntot
          write(iun,*) "kvec and bvec in gbas"
          do ind=1, nqdiv
            do iq0=1, nntot
              bvec=
     $          kpt_latt(:,nnlist(ind,iq0))+nncell(:,ind,iq0)
     $          -kpt_latt(:,ind)
              write(iun,'(6f12.6)') kpt_latt(:,ind), bvec
            enddo
          enddo
          write(iun,*) "num_bands",num_bands
          write(iun,*) "num_wann",num_wann
        endif


        call poorman_wan_proj(ispin)

        
!     Allocate M_Matrix
        
        if (maswrk) write(iun,*) "Calculating overlaps"
        call wannier_mmn(ispin)


        call timel("Calling wannier_run")

        if (maswrk) then
           write(iun,*) 'num_wann', num_wann
           write(iun,*) 'num_bands', num_bands
           write(iun,*) 'eigenvalues', eigenvalues
          call wannier_run(trim(seed_name),ndiv,nqdiv,
     $      transpose(rbas)*par*bohr,transpose(gbas)*2.0d0*pi/par/bohr,
     $      kpt_latt,num_bands,num_wann,nntot,natom,atom_symbols,
     $      tau*par*bohr,gamma_only,M_matrix,A_Matrix,eigenvalues,
     $      U_matrix,U_matrix_opt,lwindow,
     $      wann_centers,wann_spreads,spread)
        endif

        call mpi_bcast(lwindow,num_bands*nqdiv,
     $    mpi_logical,0, mpi_comm_world,mpierr)
        call mpi_bcast(wann_centers,3*num_wann,
     $    mpi_double_precision,0, mpi_comm_world,mpierr)
        call mpi_bcast(wann_spreads,num_wann,
     $    mpi_double_precision,0, mpi_comm_world,mpierr)
        call mpi_bcast(spread,3,
     $    mpi_double_precision,0, mpi_comm_world,mpierr)


        do ii=1, num_wann        
          distance=100.0d0
          do jj=1, natom
            do aa=-ndiv(1)+1, ndiv(1)-1
              do bb=-ndiv(2)+1, ndiv(2)-1
                do cc=-ndiv(3)+1, ndiv(3)-1
                  tempvec=(tau(:,jj)+
     $              rbas(:, 1)*aa+rbas(:,2)*bb+rbas(:,3)*cc)*par*bohr
     $              -wann_centers(:,ii)
                  temp_distance=dsqrt(sum(tempvec*tempvec))
                  if (distance .gt. temp_distance) then
                    distance=temp_distance
                    atom_wan(ii)=jj
                    atom_wan_distance(:,ii)=tempvec
                  endif
                enddo
              enddo
            enddo
          enddo
        enddo

        if (maswrk) then
          write(iun,*) 'wannier_center and atom'
          do ii=1, num_wann        
            tempvec=matmul(transpose(gbas),
     $        wann_centers(:,ii)/(par*bohr))
            tempvec2=matmul(transpose(gbas),
     $        tau(:,atom_wan(ii))
     $        -wann_centers(:,ii)/(par*bohr))
            distance=dsqrt(sum(atom_wan_distance(:,ii)**2))
            write(iun,'(2i5, 4f12.6,4x,2(3f12.6,2x))') ii,
     $        atom_wan(ii),atom_wan_distance(:,ii),
     $        distance,tempvec,tempvec2
          enddo
          write(iun,*) 'wannier_function centered at each atom'
          do ii=1, natom
            do jj=1, num_wann
              if (atom_wan(jj) .eq. ii) then
                write(iun,*) ii,jj
              endif
            enddo
          enddo

          write(iun,*) 'nearest neighbors and shift vectors'
          do iwan=1, num_wann        
            do jwan=1, num_wann
              distance=10000.0d0
              do ii=-5, 5
                do jj=-5, 5
                  do kk=-5, 5
                    tempvec=
     $                wann_centers(:,iwan)
     $                -wann_centers(:,jwan)
     $                -matmul(rbas, (/dble(ii), dble(jj), dble(kk)/))
     $                *(par*bohr)
                    temp_distance=dsqrt(sum(tempvec**2))
                    if (distance .gt. temp_distance) then
                      distance=temp_distance
                      tempintarr=(/ii,jj,kk/)
                    endif
                  enddo
                enddo
              enddo
              write(iun,'(7i5,f12.6)') iwan, atom_wan(iwan),
     $          jwan, atom_wan(jwan), tempintarr, distance
            enddo
          enddo


          write(iun,'(a)') 'calculate v_matrix'
          
          v_matrix=0.0d0
          do irk=1,nqdiv
            do j=1,num_wann
              tempvec2=matmul(transpose(gbas),
     $          wann_centers(:,j)/(par*bohr)-tau(:,atom_wan(j)))
!              kdottau=2.0d0*pi*
!     $          sum(tempvec2*kpt_latt(:,irk))
              do i=1,num_wann
                indtmp=0
                do m=1,num_bands            
                  if (lwindow(m,irk)) then
                    indtmp=indtmp+1
                    v_matrix(m,j,irk)=
     $                v_matrix(m,j,irk)
     $                +u_matrix_opt(indtmp,i,irk)
     $                *u_matrix(i,j,irk)
                  endif
                enddo
              enddo
            enddo
          enddo

          do ii=1,num_wann      !wan
            do jj=1,num_wann    !g
              do irk=1, nqdiv
                do ibnd=1, num_bands
                  if (lwindow(ibnd,irk)) then
                    wan_g_overlap(ii,jj)=
     $                wan_g_overlap(ii,jj)
     $                +1.0d0/dble(nqdiv)
     $                *dconjg(v_matrix(ibnd,ii,irk))
     $                *a_matrix(ibnd,jj,irk)
                  endif
                enddo
              enddo
            enddo
          enddo
          

        endif
        
        call mpi_bcast(u_matrix,num_wann*num_wann*nqdiv,
     $    mpi_double_complex,0, mpi_comm_world,mpierr)
        call mpi_bcast(wan_g_overlap,num_wann*num_wann,
     $    mpi_double_complex,0, mpi_comm_world,mpierr)
        call mpi_bcast(u_matrix_opt,num_bands*num_wann*nqdiv,
     $    mpi_double_complex,0, mpi_comm_world,mpierr)
        call mpi_bcast(v_matrix,num_bands*num_wann*nqdiv,
     $    mpi_double_complex,0, mpi_comm_world,mpierr)



        if (maswrk) then          
          write(iun,*) '|<w|g>|'
          write(iun,*) 'where g is initial guess and'
          write(iun,*) 'w is wannier function'
          do ii=1, num_wann     ! wan 
            
            tempdble=0.0d0
            tempdble=cdabs(wan_g_overlap(ii,:))
            tempdble2=tempdble
            indloc1=maxloc(tempdble)
            tempdble(indloc1(1))=0.0d0
            indloc2=maxloc(tempdble)
            tempdble(indloc2(1))=0.0d0            
            indloc3=maxloc(tempdble)
            tempdble(indloc3(1))=0.0d0            
            indloc4=maxloc(tempdble)
            tempdble(indloc4(1))=0.0d0            
            indloc5=maxloc(tempdble)

            write(iun,'(2i5,5(2i5,f12.6))') ii, 
     $        atom_wan(ii),
     $        atom_wan(indloc1(1)),indloc1,
     $        tempdble2(indloc1(1)),
     $        atom_wan(indloc2(1)),indloc2,
     $        tempdble2(indloc2(1)),
     $        atom_wan(indloc3(1)),indloc3,
     $        tempdble2(indloc3(1)),
     $        atom_wan(indloc4(1)),indloc4,
     $        tempdble2(indloc4(1)),
     $        atom_wan(indloc5(1)),indloc5,
     $        tempdble2(indloc5(1))            

          enddo
        endif

        allocate(lte_wan(nfun,num_wann,nqdiv))
        lte_wan=0.0d0
        allocate(gks_wan(nbasmpw,num_wann,nqdiv))
        gks_wan=0.0d0

        allocate(lte_bnd(nfun,num_bands,nqdiv))
        lte_bnd=0.0d0
        allocate(gks_bnd(nbasmpw,num_bands,nqdiv))
        gks_bnd=0.0d0        
        
        call za_wan_for_u(ispin)


        allocate(ndegen(8*nqdiv)) !
        ndegen=0
        allocate(irvec(3,8*nqdiv)) !
        irvec=0
        
        call wigner_seitz(nrpts, ndegen, irvec)

#ifdef USE_HDF5
        call write_wannier_dat_hdf5
#else
        if (maswrk) then
          call write_wannier_dat
        endif
#endif        
          
        if (writewan .eq. 1) then
c$$$          if (irel.eq.2) then
c$$$            if (maswrk) then
c$$$              write(iun,*)
c$$$     $          'printing out wannier function '//
c$$$     $          'for the fully-relativistic case '//
c$$$     $          'is not possible yet'
c$$$            endif
c$$$            call ending
c$$$
c$$$          endif

c$$$          allocate(wan_rgrid(nrdiv(1)*ndiv(1),
c$$$     $      nrdiv(2)*ndiv(2),nrdiv(3)*ndiv(3),nrel,num_wann)) !
c$$$          wan_rgrid=0.0d0

          call wannier_realgrid(nrdiv,ispin)
c$$$          deallocate(

        endif

c$$$        if (maswrk) then
c$$$          hmat=0.0d0
c$$$          do irk=1, nqdiv
c$$$            do iwan=1, num_wann
c$$$              do jwan=1, num_wann
c$$$                do ibnd=1, num_bands
c$$$                  hmat(iwan,jwan)=hwam(iwan,jwan)
c$$$     $              +eigenvalues(ibnd,irk)*v_matrix(ibnd,jwan,irk)
c$$$     $              *dconjg(v_matrix(ibnd,iwan,irk))/dble(nqdiv)
c$$$                enddo
c$$$              enddo
c$$$            enddo
c$$$          enddo
c$$$          
c$$$          write(iun,*) 'hmat'
c$$$          do iwan=1, num_wann
c$$$            write(iun,'(1000(2f12.6, 2x))')
c$$$     $        (hmat(iwan,jwan), jwan=1, num_wann)
c$$$          enddo
c$$$        endif
                  


        deallocate(ndegen)
        deallocate(irvec)
        deallocate(atom_xaxis)
        deallocate(atom_zaxis)        
        
!        deallocate(ldim_wan)
!        deallocate(rdim_wan)
!        deallocate(lval_wan)
!        deallocate(radial_coeff_wan)
        
        deallocate(tempdble)
        deallocate(tempdble2)
        deallocate(proj_site)
        deallocate(proj_z)
        deallocate(proj_x)
        deallocate(proj_zona)
        deallocate(nnlist)
        deallocate(nncell)
        deallocate(proj_l)
        deallocate(proj_m)
        deallocate(proj_radial)
        deallocate(exclude_bands)
        deallocate(proj_s)
        deallocate(proj_s_qaxisx)
        deallocate(exclude_bands_table)
        deallocate(include_bands)
        deallocate(A_Matrix)
        deallocate(M_Matrix)
        deallocate(atom_wan)
        deallocate(atom_wan_distance)
        deallocate(eigenvalues)
        deallocate(z_wan_bnd)
        deallocate(ev_wan_bnd)
        deallocate(U_matrix)
        deallocate(wan_g_overlap)
        deallocate(v_matrix)
        deallocate(U_matrix_opt)
        deallocate(lwindow)
        deallocate(wann_centers)
        deallocate(wann_spreads)
        deallocate(alimj)
        deallocate(nalimj)
        deallocate(alimj_coeff)
        deallocate(alimj_xaxis)
        deallocate(alimj_zaxis)        


      enddo


      end


!================================!
      subroutine wigner_seitz(nrptss,ndegenn, irvecc)
!================================!
!     Calculates a grid of lattice vectors r that fall inside (and eventually  !
!     on the surface of) the Wigner-Seitz supercell centered on the            ! 
!     origin of the Bravais superlattice with primitive translations           !
!     ndiv(1)*a_1, ndiv(2)*a_2, and ndiv(3)*a_3                       !
!==========================================================================!
      use comwann_mod
      implicit none



      integer, intent(out) :: nrptss, ndegenn(8*nqdiv),
     $  irvecc(3, 8*nqdiv)

      integer       :: ndiff(3)
      double precision :: dist(125),tot,dist_min,real_metric(3,3)
      integer       :: n1,n2,n3,i1,i2,i3,icnt,i,j,ir

      ndegenn=0
      irvecc=0

c$$$  if (timing_level>1.and.on_root)
c$$$  $  call io_stopwatch('postw90_common: wigner_seitz',1)

!     The Wannier functions live in a periodic supercell of the real space unit 
!     cell. This supercell is ndiv(i) unit cells long along each primitive
!     translation vector a_i of the unit cell
!     
!     We loop over grid points r on a cell that is approx. 8 times
!     larger than this "primitive supercell."
!     
!     One of these points is in the W-S supercell if it is closer to R=0 than any
!     of the other points R (where R are the translation vectors of the 
!     supercell). In practice it is sufficient to inspect only 125 R-points.

!     In the end, nrpts contains the total number of grid points that have been 
!     found in the Wigner-Seitz cell
      real_metric=matmul(transpose(rbas), rbas)*(par*bohr)**2

      nrptss = 0  
      do n1 = -ndiv(1) , ndiv(1)  
        do n2 = -ndiv(2), ndiv(2)  
          do n3 = -ndiv(3),  ndiv(3)  
!     Loop over the 125 points R. R=0 corresponds to i1=i2=i3=0, 
!     or icnt=63
            icnt = 0  
            do i1 = -2, 2  
              do i2 = -2, 2  
                do i3 = -2, 2  
                  icnt = icnt + 1  
!     Calculate distance squared |r-R|^2
                  ndiff(1) = n1 - i1 * ndiv(1)  
                  ndiff(2) = n2 - i2 * ndiv(2)  
                  ndiff(3) = n3 - i3 * ndiv(3)  
                  dist(icnt) = 0.0d0
                  do i = 1, 3  
                    do j = 1, 3  
                      dist(icnt)=dist(icnt)+
     $                  dble(ndiff(i))*real_metric(i,j)*dble(ndiff(j))
                    enddo
                  enddo
                enddo
              enddo
            enddo
            dist_min=minval(dist)
            if (abs(dist(63) - dist_min ) .lt.1.0d-7) then ! R=0
              nrptss = nrptss + 1  
              ndegenn(nrptss)=0
              do i=1,125
                if(abs(dist(i)-dist_min).lt.1.0d-7)
     $            ndegenn(nrptss)=ndegenn(nrptss)+1
              end do
              irvecc(1, nrptss) = n1  
              irvecc(2, nrptss) = n2   
              irvecc(3, nrptss) = n3   
!     
!     Remember which grid point r is at the origin
!     
!              if (n1==0 .and. n2==0 .and. n3==0) rpt_origin=nrptss
c$$$  endif
            end if
!     n3
          enddo
!     n2
        enddo
!     n1
      enddo
!     

      tot = 0.0d0
      do ir = 1, nrptss  
!     
!     Corrects weights in Fourier sums for R-vectors on the boundary of the 
!     W-S supercell 
!     
        tot=tot+1.0d0/dble(ndegenn(ir))
      enddo
      if (abs(tot-dble(ndiv(1)*ndiv(2)*ndiv(3))) > 1.0d-8)
     $  call ending
     $  ('ERROR in wigner_seitz: error in finding Wigner-Seitz points')

      end subroutine wigner_seitz







      subroutine wannier_interpolation
!     be careful,
!     input w.r.t. band basis
!     output w.r.t. wannier function if do_diagonalize==1
!     diag_coarse(num_bands,nqdiv): double precision
!     mat_coarse(num_bands,num_bands,nqdiv): complex*16
!     mat_fine(num_wann,num_wann,ndim_kk_fi_max): complex*16
!     eig_fine(num_wann,ndim_kk_fin_max): double precision
     $  (is_diagonal,do_diagonalize,
     $  mat_coarse,diag_coarse,
     $  mat_fine,eig_fine,
     $  nqdiv_fi,ndim_kk_fi,n_mpi_kk_fi,ndim_kk_fi_max,kvec)

      use comwann_mod

      implicit none
      include 'mpif.h'
      integer, intent(in) :: nqdiv_fi,
     $  ndim_kk_fi(nproc),n_mpi_kk_fi(nproc),ndim_kk_fi_max,
     $  do_diagonalize,is_diagonal
      double precision, intent(in) ::
     $  diag_coarse(num_bands,nqdiv),kvec(3,nqdiv_fi)
      complex*16,intent(in) :: 
     $  mat_coarse(num_bands,num_bands,nqdiv)                  
      double precision, intent(out) ::
     $  eig_fine(num_wann,ndim_kk_fi_max)
      complex*16,intent(out) :: 
     $  mat_fine(num_wann,num_wann,ndim_kk_fi_max)
      
      integer :: me1,iwan,jwan,irk,iband,jband,
     $ i1,i2,i3,j1,ik,jj
      complex*16 :: h_k_wan_co(nqdiv),ai,h_r_wan(nrpts)

      ai=dcmplx(0.0d0, 1.0d0)
      me1=me+1

      mat_fine=0.0d0
      eig_fine=0.0d0
c      write(179,*) 'in wannier_interpolation'
      call flush(179)
      
      do iwan=1, num_wann
        do jwan=1, num_wann

          h_r_wan=0.0d0
          h_k_wan_co=0.0d0
          do irk=1, nqdiv
            do iband=1, num_bands
              if (is_diagonal .eq. 1) then
                h_k_wan_co(irk)=h_k_wan_co(irk)+
     $            diag_coarse(iband,irk)*v_matrix(iband,jwan,irk)
     $            *dconjg(v_matrix(iband,iwan,irk))
              else
                do jband=1, num_bands                  
                  h_k_wan_co(irk)=h_k_wan_co(irk)+
     $              mat_coarse(iband,jband,irk)*v_matrix(jband,jwan,irk)
     $              *dconjg(v_matrix(iband,iwan,irk))
                enddo
              endif
            enddo
          enddo
          
!     fourier transform to r space
          
          do i1=1, nrpts
            do j1=1, nqdiv
              h_r_wan(i1)=h_r_wan(i1)
     $          +h_k_wan_co(j1)*
     $          cdexp(-2.0*ai*pi*
     $          sum(irvec(:,i1)*kpt_latt(:,j1)))
     $          /dble(nqdiv)
            enddo
          enddo
          
          do ik=1, ndim_kk_fi(me1)
            irk=n_mpi_kk_fi(me1)+ik
            do j1=1, nrpts
              mat_fine(iwan,jwan,ik)
     $          =mat_fine(iwan,jwan,ik)
     $          +h_r_wan(j1)*
     $          cdexp(2.0*ai*pi*
     $          sum(kvec(:,irk)
     $          *irvec(:,j1))
     $          )/dble(ndegen(j1))
            enddo
          enddo
        enddo
      enddo

      if (do_diagonalize .eq. 1) then
        do ik=1, ndim_kk_fi(me1)
          call hermitianeigen_cmplxdouble(num_wann,eig_fine(1,ik),
     $      mat_fine(1,1,ik))
        enddo
      endif
c      write(179,*) 'end of wannier_interpolation'      
      end



      subroutine za_wan_for_u(ispin)

      use comwann_mod
      implicit none
      include "mpif.h"
      
      integer, intent(in) :: ispin
      integer :: ind_k, k,
     $  k0,iq
      
      lte_bnd=0.0d0
      gks_bnd=0.0d0      

      do ind_k=1,ndim_kk(me+1)  ! k vector
        k=n_mpi_kk(me+1)+ind_k
        k0=i_kref(k)
        call sym_z_0(lte_bnd(1,1,k),k,z_wan_bnd(1,1,k0),
     &    num_bands,k_group(k),pnt(1,k))
        call sym_a_2(gks_bnd(1,1,k),k,ev_wan_bnd(1,1,k0),
     &    num_bands,k_group(k),k0)
      enddo


      if(nproc/=1) then
        call mpi_allreduce_dcmplx
     $    (lte_bnd,nfun*num_bands*nqdiv,
     $    mpi_sum,mpi_comm_world)
        call mpi_allreduce_dcmplx
     $    (gks_bnd,nbasmpw*num_bands*nqdiv,
     $    mpi_sum,mpi_comm_world)
      endif

      do iq=1, nqdiv
        call zgemm('n','n',nfun,num_wann,num_bands,(1.0d0,0.0d0),
     $    lte_bnd(1,1,iq),nfun,v_matrix(1,1,iq),num_bands,
     $    (0.0d0,0.0d0),lte_wan(1,1,iq),nfun)
        call zgemm('n','n',nbasmpw,num_wann,num_bands,(1.0d0,0.0d0),
     $    gks_bnd(1,1,iq),nbasmpw,v_matrix(1,1,iq),num_bands,
     $    (0.0d0,0.0d0),gks_wan(1,1,iq),nbasmpw)
      enddo

c$$$      if (maswrk) then
c$$$        open(unit=334, file=trim(seed_name)//'_z_wan_band.dat',
c$$$     $    form='unformatted')
c$$$        write(334) lte_bnd
c$$$        close(334)
c$$$        open(unit=334, file=trim(seed_name)//'_a_wan_band.dat',
c$$$     $    form='unformatted')
c$$$        write(334) gks_bnd
c$$$        close(334)
c$$$      endif
      
      end
      


      subroutine read_wann_inip

      use comwann_mod            

      implicit none
      include "mpif.h"      
      integer :: iwan,iorb,ii,jorb,iatom,lval,mval1
      double precision :: tempvec1(3),tempvec2(3),mval2,ival

      if (maswrk) then

        
        if (irel .le. 1) then
          num_alimj=3
        else
          num_alimj=4
        endif

        open(unit=111, file=trim(seed_name)//'.inip') ! initial projector
        allocate(nalimj(num_wann))
        do iwan=1, num_wann
          nalimj(iwan)=1
        enddo
        max_nalimj=maxval(nalimj)              
        allocate(alimj(num_alimj,max_nalimj,num_wann))
        alimj=0.0d0
        allocate(alimj_coeff(max_nalimj,num_wann))
        alimj_coeff=0.0d0
        allocate(alimj_xaxis(3,max_nalimj,num_wann))
        alimj_coeff=0.0d0
        allocate(alimj_zaxis(3,max_nalimj,num_wann))
        alimj_coeff=0.0d0                        

        do iwan=1, num_wann
          if (irel .le. 1) then
            read(111, *) 
     $        iatom, lval, mval1,
     $        (alimj_xaxis(ii, 1,iwan), ii=1, 3),
     $        (alimj_zaxis(ii, 1,iwan), ii=1, 3)
            alimj(1,1,iwan)=iatom
            alimj(2,1,iwan)=lval
            alimj(3,1,iwan)=mval1+lval+1
          else
            read(111, *) 
     $        iatom, lval, ival,mval2,
     $        (alimj_xaxis(ii, 1,iwan), ii=1, 3),
     $        (alimj_zaxis(ii, 1,iwan), ii=1, 3)
            alimj(1,1,iwan)=iatom
            alimj(2,1,iwan)=lval
            alimj(3,1,iwan)=nint(ival*2.0d0)
            alimj(4,1,iwan)=nint(mval2*2.0d0)
          endif          
          do iorb=1, nalimj(iwan)
            alimj_coeff(iorb,iwan)=dcmplx(1.0d0, 0.0d0)
          enddo
        enddo
        close(111)

      endif

      call mpi_bcast(num_alimj,1,mpi_integer,
     $  0,mpi_comm_world,mpierr)

      if (.not. maswrk) allocate(nalimj(num_wann))      
      call mpi_bcast(nalimj,num_wann,mpi_integer,
     $  0,mpi_comm_world,mpierr)
      
      if (.not. maswrk) then

        max_nalimj=maxval(nalimj)            
        allocate(alimj(num_alimj,max_nalimj,num_wann))
        alimj=0.0d0
        allocate(alimj_coeff(max_nalimj,num_wann))
        alimj_coeff=0.0d0
        allocate(alimj_xaxis(3,max_nalimj,num_wann))
        alimj_xaxis=0.0d0
        allocate(alimj_zaxis(3,max_nalimj,num_wann))
        alimj_zaxis=0.0d0                        
      endif
      call mpi_bcast(num_wann,1,mpi_integer,
     $  0,mpi_comm_world,mpierr)
      call mpi_bcast(alimj,num_alimj*max_nalimj*num_wann,mpi_integer,
     $  0,mpi_comm_world,mpierr)
      call mpi_bcast(alimj_coeff,max_nalimj*num_wann,
     $  mpi_double_complex,0,mpi_comm_world,mpierr)
      call mpi_bcast(alimj_xaxis,3*max_nalimj*num_wann,
     $  mpi_double_precision,0,mpi_comm_world,mpierr)
      call mpi_bcast(alimj_zaxis,3*max_nalimj*num_wann,
     $  mpi_double_precision,0,mpi_comm_world,mpierr)

      end







c      subroutine find_fine_grid
c      use comwann_mod      
c      implicit none
c
c      integer :: ii,jj,ik,tempk(3)
c      double precision :: fbz_vol, g_abs(3),fac
c      double precision, external :: dcrossdot      
c
c!      pi=4.0d0*datan2(1.0d0,1.0d0)
c
c      fbz_vol=dcrossdot(gbas(:,1),gbas(:,2),gbas(:,3))
c     $  *(2.0d0*pi/par)**3/0.001d0
c
c      do ii=1, 3
c        g_abs(ii)=dsqrt(sum(gbas(:,ii)**2))
c      enddo
c
c      fac=(fbz_vol/g_abs(1)/g_abs(2)/g_abs(3))**(1.0d0/3.0d0)
c      
c      do ii=1, 3
c        mp_grid_fi(ii)=ceiling(fac*g_abs(ii))
c      enddo
c
c      num_kpts_fi=mp_grid_fi(1)*mp_grid_fi(2)*mp_grid_fi(3)
c
c      if (maswrk) then
c        write(iun,'(a,i5)') 'num_kpts_fi=', num_kpts_fi
c        write(iun,'(a,3i5)') 'mp_grid_fi=', mp_grid_fi
c      endif
c
c      
c      end
      

      subroutine wannier_interpolation_hamiltonian
     $  (nqdiv_fi,ndim_kk_fi,n_mpi_kk_fi,ndim_kk_fi_max,kvec,
     $  mat_fine,eig_fine)      

      use comwann_mod      

      implicit none
      include 'mpif.h'

      integer, intent(in) :: nqdiv_fi,
     $  ndim_kk_fi(nproc),n_mpi_kk_fi(nproc),ndim_kk_fi_max
      double precision, intent(in) ::
     $  kvec(3,nqdiv_fi)      

      double precision :: eig_fine(num_wann,ndim_kk_fi_max)
      complex*16, intent(out) ::
     $  mat_fine(num_wann,num_wann,ndim_kk_fi_max)

      complex*16 :: dummy_mat(num_bands,num_bands,nqdiv)
      
      call wannier_interpolation(1,1,dummy_mat,eigenvalues,
     $  mat_fine,eig_fine,
     $  nqdiv_fi,ndim_kk_fi,n_mpi_kk_fi,ndim_kk_fi_max,kvec)            
      end





c      subroutine cal_wan_pdos
c
c      use comwann_mod      
c
c      implicit none
c      include 'mpif.h'
c
c      integer :: ii,ibnd,iwan,ik
c      
c      double precision :: 
c     $  ee,pdos(800,num_wann)
c
c      pdos=0.0d0
c      
c      do ii=1, 800
c        ee=(ii-400)*0.1d0
c        do ibnd=1, num_wann
c          do iwan=1, num_wann
c            do ik=1, ndim_mp_grid_fi(me+1)            
c              pdos(ii, iwan)=pdos(ii, iwan)
c     $          +1.0d0/0.1d0/dsqrt(2.0d0*pi)
c     $          *dexp(-(ee-eigenvalues_fi(ibnd,ik))**2/2.0d0/0.1d0**2)
c     $          *cdabs(w_matrix_fi(iwan,ibnd,ik))**2/dble(num_kpts_fi)
c            enddo
c          enddo
c        enddo
c      enddo
c
c      call mpi_allreduce_dble(pdos,800*num_wann,mpi_sum,mpi_comm_world)
c      
c      if (maswrk) then      
c        open (unit = 990, file = adjustl(trim(seed_name))//'.pdos', 
c     $    form = 'formatted')
c        write(990,'(2(a,i5))') '#n_omega=', 800, ', num_wann=', num_wann
c        do ii=1, 800
c          ee=(ii-400)*0.1d0          
c          write(990,'(f12.6,5x,1000f12.6)')
c     $      ee,(pdos(ii,iwan), iwan=1, num_wann)
c        enddo
c      endif
c      end                  

c      subroutine prepare_optics_input(ispin)
c
c      use comwann_mod      
c
c      implicit none
c      include 'mpif.h'
c
c      integer, intent(in) :: ispin
c
c      integer :: i1,i2,i3,ig1,ig2,ig3,j1,
c     $  ik,irk,ibnd,jbnd,iwan,ii,
c     $  rspace_dim,rspace_grid(3),
c     $  nqvec_fi(3),nqvec_fi_shifted(3),rvec(3),rvec_shifted(3),
c     $  gvec(3),gvec_shifted(3),idir,me1,
c     $  proc_num,indout,jwan,iband,jband
c
c
c      complex*16 :: ai
c      integer :: irequest,istatus(mpi_status_size)      
c      logical :: mpiflag
c
c      complex*16 :: momentum_wan_r(num_wann,num_wann,nrpts,3)
c      complex*16, allocatable :: tempmat(:,:,:,:),tempmat1(:,:),
c     $  tempmat2(:,:,:)
c
c
c      if (irel .eq. 2) return
c      
c      ai=dcmplx(0.0d0, 1.0d0)
c      me1=me+1
c      if (irel .le. 1) then
c        allocate(momentum(num_bands,num_bands,nqdiv,3))
c        momentum=0.0d0
c        call moment_operator_nospinorbit(ispin,momentum)
c        
c        if (maswrk) then
c
c          open(unit=55, file='momentum_optics_bnd.dat')
c          do idir=1,3
c            do irk=1, nqdiv            
c              do iwan=1, num_bands
c                do jwan=1, num_bands
c                  write(55,'(6i5,2e20.12)') iwan,jwan,
c     $              nint(kpt_latt(:,irk)*ndiv),
c     $              idir,momentum(iwan,jwan,irk,idir)
c                enddo
c              enddo
c            enddo
c          enddo
c          close(55)
c          
c          allocate(momentum_wan(num_wann,num_wann,nqdiv,3))
c          momentum_wan=0.0d0
c          do idir=1,3
c            do irk=1, nqdiv            
c              do iwan=1, num_wann
c                do jwan=1, num_wann
c                  do iband=1, num_bands
c                    do jband=1, num_bands                  
c                      momentum_wan(iwan,jwan,irk,idir)
c     $                  =momentum_wan(iwan,jwan,irk,idir)+
c     $                  momentum(iband,jband,irk,idir)
c     $                  *v_matrix(jband,jwan,irk)
c     $                  *dconjg(v_matrix(iband,iwan,irk))
c                    enddo
c                  enddo
c                enddo
c              enddo
c            enddo
c          enddo
c          open(unit=55, file='momentum_optics_wan.dat')
c          do idir=1,3
c            do irk=1, nqdiv            
c              do iwan=1, num_wann
c                do jwan=1, num_wann          
c                  write(55,'(6i5,2e20.12)') iwan,jwan,
c     $              nint(kpt_latt(:,irk)*ndiv),
c     $              idir,momentum_wan(iwan,jwan,irk,idir)
c                enddo
c              enddo
c            enddo
c          enddo
c          close(55)
c
c          momentum_wan_r=0.0d0
c
c          do idir=1,3
c            do iwan=1, num_wann
c              do jwan=1, num_wann          
c                do i1=1, nrpts
c                  do j1=1, nqdiv
c                    momentum_wan_r(iwan,jwan,i1,idir)
c     $                =momentum_wan_r(iwan,jwan,i1,idir)
c     $                +momentum_wan(iwan,jwan,j1,idir)*
c     $                cdexp(-2.0*ai*pi*
c     $                sum(irvec(:,i1)*kpt_latt(:,j1)))
c     $                /dble(nqdiv)
c                  enddo
c                enddo
c
c              enddo
c            enddo
c          enddo
c
c          open(unit=55, file='momentum_optics_wan_r.dat')
c          do idir=1,3
c            do irk=1, nrpts
c              do iwan=1, num_wann
c                do jwan=1, num_wann          
c                  write(55,'(6i5,2e20.12)') iwan,jwan,irvec(:,irk),idir,
c     $              momentum_wan_r(iwan,jwan,irk,idir)
c                enddo
c              enddo
c            enddo
c          enddo
c          close(55)
c          
c        endif
c        
c        
c      endif
c      
c      deallocate(momentum)
c      if (maswrk) deallocate(momentum_wan)      
c      end


c      subroutine prepare_photoemission_input
c     $  (nqdiv_fi,kvec,ndim_kk_fi,n_mpi_kk_fi,ndim_kk_fi_max,
c     $  wan_rgrid,ispin)
c
c      use comwann_mod      
c
c      implicit none
c      include 'mpif.h'
c
c      integer, intent(in) :: nqdiv_fi,
c     $  ndim_kk_fi(nproc),n_mpi_kk_fi(nproc),ndim_kk_fi_max,ispin
c      double precision,intent(in) :: 
c     $  kvec(3,nqdiv_fi)
c
c      complex*16,intent(in) :: 
c     $  wan_rgrid(nrdiv(1)*ndiv(1),nrdiv(2)*ndiv(2),
c     $  nrdiv(3)*ndiv(3),num_wann)
c
c      integer :: i1,i2,i3,ig1,ig2,ig3,
c     $  ik,irk,ibnd,jbnd,iwan,ii,jj,kk,rspace_dim,rspace_grid(3),
c     $  nqvec_fi(3),nqvec_fi_shifted(3),rvec(3),rvec_shifted(3),
c     $  gvec(3),gvec_shifted(3),idir,me1,
c     $  ind_k,cnt,proc_num,indout
c
c      double precision :: 
c     $  kpath_pw(3,nqdiv_fi),
c     $  eig_fine(num_wann,ndim_kk_fi_max),
c     $  eig_fine_allk(num_wann,nqdiv_fi),
c     $  energy,e_pw
c
c      complex*16 ::
c!     $  ai,
c     $  pw_wan(3,num_wann,ndim_kk_fi_max),
c     $  pw_wan_tmp(3,num_wann),
c     $  mat_fine(num_wann,num_wann,ndim_kk_fi_max),
c     $  mat_fine_tmp(num_wann,num_wann)
c      
c
c      integer :: irequest,istatus(mpi_status_size)      
c      logical :: mpiflag
c      
c!      ai=dcmplx(0.0d0, 1.0d0)
c      me1=me+1
c
c      call wannier_interpolation_hamiltonian(
c     $  nqdiv_fi,ndim_kk_fi,n_mpi_kk_fi,ndim_kk_fi_max,
c     $  kvec,mat_fine,eig_fine)
c
c      eig_fine_allk=0.0d0
c      do ik=1, ndim_kk_fi(me1)
c        irk=n_mpi_kk_fi(me1)+ik        
c        eig_fine_allk(:,irk)=eig_fine(:,ik)
c      enddo
c      
c      call mpi_allreduce_dble
c     $  (eig_fine_allk,nqdiv_fi*num_wann,mpi_sum,mpi_comm_world)
c
c
c      if (maswrk) then
c        if (ispin .eq. 1) then
c          write(119, '(4(a13, i5))')
c     $      '# is_soc   : ', 0,
c     $      ', num_spin : ', nspin,
c     $      ', num_k    :', nqdiv_fi,
c     $      ', num_bands: ', num_wann
c        endif
c        do ik=1, nqdiv_fi
c          do ibnd=1, num_wann
c            write(119, '(f20.12)')
c     $        eig_fine_allk(ibnd,ik)
c          enddo
c        enddo
c      endif
c      
c      if (maswrk) then
c
c
c        cnt=0
c        do ii=-nrdiv(1)/2, nrdiv(1)/2
c          do jj=-nrdiv(2)/2, nrdiv(2)/2
c            do kk=-nrdiv(3)/2, nrdiv(3)/2        
c              energy=sum((matmul(gbas,(/ii,jj,kk/))*2.0d0*pi/par)**2) ! energy in rydberg
c              if(
c     $          (energy*rydberg .gt. photoemission_emin)
c     $          .and.
c     $          (energy*rydberg .lt. photoemission_emax)
c     $          ) then
c                cnt=cnt+1
c              endif
c            enddo
c          enddo
c        enddo
c
c        
c        if (ispin .eq. 1) then
c          if (irel .eq. 2) then
c            write(129) 1,nspin,cnt,nqdiv_fi,num_wann
c          else
c            write(129) 0,nspin,cnt,nqdiv_fi,num_wann
c          endif
c          
c          if (irel .eq. 2) then
c            write(139) 1,nspin,nqdiv_fi,num_wann
c          else
c            write(139) 0,nspin,nqdiv_fi,num_wann
c          endif          
c          
c          do ii=-nrdiv(1)/2, nrdiv(1)/2
c            do jj=-nrdiv(2)/2, nrdiv(2)/2
c              do kk=-nrdiv(3)/2, nrdiv(3)/2        
c                energy=sum((matmul(gbas,(/ii,jj,kk/))*2.0d0*pi/par)**2) ! energy in rydberg
c                if(
c     $            (energy*rydberg .gt. photoemission_emin)
c     $            .and.
c     $            (energy*rydberg .lt. photoemission_emax)
c     $            ) then
c                  write(129) ii,jj,kk
c                endif
c              enddo
c            enddo
c          enddo
c        endif
c      endif
c
c!     moment
c      if (.not. maswrk) then
c        do ik=1, ndim_kk_fi(me1)
c          irk=n_mpi_kk_fi(me1)+ik
c          call mpi_isend(mat_fine(1,1,ik),
c     $      num_wann**2,mpi_double_complex,0,irk,
c     $      mpi_comm_world,irequest,mpierr)
c          call mpi_wait(irequest,istatus,mpierr)
c          call mpi_test(irequest,mpiflag,istatus,mpierr)
c          if (.not.mpiflag)
c     $      call mpi_abort(mpi_comm_world,mpierr)
c        enddo
c      else
c        do irk=1, nqdiv_fi
c          call find_proc_ind(ndim_kk_fi, n_mpi_kk_fi,nproc,irk,
c     $      proc_num,indout)
c          if (proc_num .eq. 0)then
c            write(139) mat_fine(:,:,indout)
c          else
c            mat_fine_tmp=0.0d0
c            call mpi_irecv(mat_fine_tmp,
c     $        num_wann**2,mpi_double_complex,proc_num,
c     $        irk,mpi_comm_world,irequest,mpierr)
c            call mpi_wait(irequest,istatus,mpierr)
c            call mpi_test(irequest,mpiflag,istatus,mpierr)
c            if (.not.mpiflag)
c     $        call mpi_abort(mpi_comm_world,mpierr)
c            write(139) mat_fine_tmp
c          endif
c        enddo
c      endif                     
c      
c
c      if (maswrk) then
c        open(unit=118, file='matrix_element')
c      endif
c      do ii=-nrdiv(1)/2, nrdiv(1)/2
c        do jj=-nrdiv(2)/2, nrdiv(2)/2
c          do kk=-nrdiv(3)/2, nrdiv(3)/2
c            if (maswrk) then
c              write(iun,'(a, 3i5)') 'g vector scan', ii,jj,kk
c            endif
c            energy=sum((matmul(gbas,(/ii,jj,kk/))*2.0d0*pi/par)**2) ! energy in rydberg
c            if(
c     $        (energy*rydberg .gt. photoemission_emin)
c     $        .and.
c     $        (energy*rydberg .lt. photoemission_emax)
c     $        ) then
c              
c              
c              kpath_pw=0.0d0
c              
c              do ik=1,nqdiv_fi
c                kpath_pw(:,ik)=kvec(:,ik)+(/ii,jj,kk/)
c              enddo
c
c              e_pw=sum((matmul(gbas,kpath_pw)*2.0d0*pi/par)**2)*rydberg
c              
c              call cal_pw_wan(nqdiv_fi,ndim_kk_fi,n_mpi_kk_fi,
c     $          ndim_kk_fi_max,kpath_pw,wan_rgrid,pw_wan)
c
c
c              if (.not. maswrk) then
c                do ik=1, ndim_kk_fi(me1)
c                  irk=n_mpi_kk_fi(me1)+ik
c                  call mpi_isend(pw_wan(1,1,ik),
c     $              3*num_wann,mpi_double_complex,0,irk,
c     $              mpi_comm_world,irequest,mpierr)
c                  call mpi_wait(irequest,istatus,mpierr)
c                  call mpi_test(irequest,mpiflag,istatus,mpierr)
c                  if (.not.mpiflag)
c     $              call mpi_abort(mpi_comm_world,mpierr)
c                enddo
c              else
c                do irk=1, nqdiv_fi
c                  call find_proc_ind(ndim_kk_fi, n_mpi_kk_fi,nproc,irk,
c     $              proc_num,indout)
c                  if (proc_num .eq. 0)then
c                    write(129) pw_wan(:,:,indout)
c                    if ((ii.eq.2).and.(jj.eq.2).and.(kk.eq.0)) then
c                      do ibnd=1,num_wann
c                        do idir=1, 3
c                          write(118,'(2f20.12)')
c     $                      pw_wan(idir,ibnd,indout)
c                        enddo
c                      enddo
c                    endif
c                  else
c                    pw_wan_tmp=0.0d0
c                    call mpi_irecv(pw_wan_tmp,
c     $                3*num_wann,mpi_double_complex,proc_num,
c     $                irk,mpi_comm_world,irequest,mpierr)
c                    call mpi_wait(irequest,istatus,mpierr)
c                    call mpi_test(irequest,mpiflag,istatus,mpierr)
c                    if (.not.mpiflag)
c     $                call mpi_abort(mpi_comm_world,mpierr)
c                    write(129) pw_wan_tmp
c
c                    if ((ii.eq.2).and.(jj.eq.2).and.(kk.eq.0)) then
c                      do ibnd=1,num_wann
c                        do idir=1, 3
c                          write(118,'(2f20.12)')
c     $                      pw_wan_tmp(idir,ibnd)
c                        enddo
c                      enddo
c                    endif
c                    
c                  endif
c                enddo
c              endif
c              call mpi_barrier(mpi_comm_world,mpierr)
c            endif
c          enddo
c        enddo
c      enddo
c      if (maswrk)then
c        close(118)
c      endif
c      
c      end      











      subroutine find_proc_ind(ndim_proc, n_mpi_proc, nproc,indin,
     $  proc_num,indout)
!     indin starting from 1
!     proc_num starting from 0
!     indout starting from 1      

      implicit none
      integer, intent(in) :: nproc, ndim_proc(nproc),n_mpi_proc(nproc),
     $  indin
      integer, intent(out) :: proc_num,indout
      integer :: ip
      
      proc_num=0
      indout=0
      do ip=0, nproc-1
        if (indin .ge. n_mpi_proc(ip+1)+1) then
          proc_num=ip
        endif
      enddo

      indout=indin-n_mpi_proc(proc_num+1)
      end      


      subroutine write_wannier_dat
      use comwann_mod

      implicit none

      integer num_kpts
      logical have_disentangled

      open(unit=10, file=adjustl(trim(seed_name))//'.dat',
     $  form='unformatted')

      write(10) rbas*par*bohr
      write(10) gbas*2.0d0*pi/par/bohr
      write(10) num_bands
      write(10) num_wann
      write(10) ndiv 
      write(10) kpt_latt
      write(10) include_bands
      write(10) v_matrix
      write(10) eigenvalues
      write(10) nrpts
      write(10) irvec
      write(10) ndegen
      write(10) u_matrix      
      have_disentangled=.true.
      write(10) have_disentangled
      if(have_disentangled) then
        write(10) lwindow
        write(10) u_matrix_opt
      endif

      close(10)                                                                                       
      end                                           
      
