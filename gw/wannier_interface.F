      subroutine wannier_interface
      
      use wannier90_mod
      use parallel_mod      
      use manager_mod      
      use solid_mod      
      use atom_mod      
      use units_mod      


      implicit none             
      include "mpif.h"
      integer :: ik,ispin,ii,jj,kk,isort,ib,indexb,
     $  tempint,ind_k,ibnd,k,jbnd,ind,iq0,iatom,irk,iwan,iband,n,
     $  iband1,iband2,k0,i,kq,kq0,k1,k2,k3,aa,bb,cc,
     $  indloc1(1),indloc2(1),indloc3(1),indloc4(1),indloc5(1),
     $  indtmp,j,m,tempintarr(3),
     $  tauindex,qindex,jwan,kwan,ndiv_fi(3),nqdiv_fi,
     $  me1,ind_wann,i1,i2,i3,j1,j2,j3,
     $  wanntot,iatom_c,atom_count,info,ind_omega,
     $  i_omega, ndim,ikstar, kstar

c$$$      integer, allocatable :: wan_correlated(:,:)
      double precision :: bvec(3),tempvec(3),tempvec2(3),kqvec(3),
     $  distance,temp_distance,kdottau,dblesum,
     $  dos_correlated(1600,5,2), nelectron,
     $  e_projected,rho_projected
      double precision, allocatable :: tempdble(:),tempdble2(:)
c$$$  $,  kpath(:,:),kpath_point(:,:),tempval(:)
      complex*16 ::ai,htemp
!     tempmatrix1(:,:),tempmatrix2(:,:), tempdcmplx(:), 
      character*300 :: dummy,ch, kstar_ch, i_omega_ch
      complex*16, allocatable :: g_out1(:,:), g_out2(:,:), g_out3(:,:),
     $  g_out4(:,:), g_out5(:,:),hmat(:,:), rhomat(:,:)

      
!---------------------------------------------------------------------
      ai=dcmplx(0.0d0,1.0d0)
      me1=me+1

      if (magn .gt. 1) then
        write(179,*) 'wannier interface only works if magn<=1'
        write(993,*) 'wannier interface only works if magn<=1'
        call ending
      endif
      
      if (maswrk) then
        write(iun,*) ""
        write(iun,*) ""
        write(iun,*) ""
        write(iun,*) "!----------------------------------------------!"
        write(iun,*) "!----------------------------------------------!"
        write(iun,*) "!                  wannier                     !"
        write(iun,*) "!----------------------------------------------!"
        write(iun,*) "!                  Warning                     !"
        write(iun,*) "!   (1) you should set 'gamma_only' '.false.'  !"
        write(iun,*) "!    in your win file                          !"
        write(iun,*) "!    even for gamma only calculation           !"
        write(iun,*) "!   (2) it's working only if magn<=1       .'  !"
        write(iun,*) "!----------------------------------------------!"

        write(iun,*) ""
        write(iun,*) "Information for the wannier code calls"
        write(iun,*) ""
        write(iun, *) 'atom_symbols kpt_latt are defined in solid.f'
        write(iun,*) 

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

c$$$        allocate(wan_correlated(maxel_loc, natom_loc))
c$$$        wan_correlated=0
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

        if (maswrk) then
          write(iun,*) "Calling wannier_setup"
        endif
        
        gamma_only = .false.
        if (maswrk) then
          write(iun,*) "atom_symbols"
          do iatom=1, natom
            write(iun,*) atom_symbols(iatom)
          enddo
        endif
        num_wann=0
        
        if (maswrk) then
          call wannier_setup(trim(seed_name),ndiv,nqdiv,
     $      transpose(rbas)*par*bohr,transpose(gbas)*2.0d0*pi/par/bohr,
     $      kpt_latt,nbndf,
     $      natom,atom_symbols,tau*par*bohr,gamma_only,.false.,
     $      nntot,nnlist,nncell,num_bands,num_wann,
     $      proj_site,proj_l,proj_m,proj_radial,proj_z,
     $      proj_x,proj_zona,exclude_bands,proj_s,proj_s_qaxisx,
     $      dis_win_min,dis_win_max)

          write(iun,*)
          write(iun,*) 'dis_win_min', dis_win_min
          write(iun,*) 'dis_win_max', dis_win_max
          write(iun,*) 'kindex test'

          do k=1, nqdiv
            k0=i_kref(k)
            write(iun,'(i5, 3f12.6)') k, kpt_latt(:,k)
            do iq0=1,nntot      ! b vector
              bvec=(kpt_latt(:,nnlist(k,iq0))+nncell(:,k,iq0)
     $          -kpt_latt(:,k))
              call zone1_number(matmul(gbas, kpt_latt(:,nnlist(k,iq0))),
     $          rb0,ndiv,jj)
              kk=index_k1(jj)
              write(iun,'(i5, 3f12.6, i5, 3f12.6)') iq0,
     $          bvec,nnlist(k,iq0)-kk,
     $          kpt_latt(:,nnlist(k,iq0))-
     $          matmul(transpose(rbas),pnt(:,kk))
            enddo
            write(iun,*)
            write(iun,*)
            write(iun,*)
          enddo

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
        call mpi_bcast(dis_win_min,1,mpi_double_precision,0,
     $    mpi_comm_world,mpierr)
        call mpi_bcast(dis_win_max,1,mpi_double_precision,0,
     $    mpi_comm_world,mpierr)

        if ((mod(num_bands, 2) .eq. 1) .and. (irel .eq. 2)) then
          write(179,*) 'num_bands should be even number if irel .eq. 2'
          write(993,*) 'num_bands should be even number if irel .eq. 2'
          call ending
        endif

        allocate(ldim_wan(num_wann))
        allocate(rdim_wan(3,num_wann))
        allocate(lval_wan(3,num_wann))
        allocate(radial_coeff_wan(100,3,num_wann))

        allocate(tempdble(num_wann))
        allocate(tempdble2(num_wann))
c$$$        allocate(tempdcmplx(num_wann))
c$$$        allocate(tempmatrix1(num_wann,num_wann))
c$$$        allocate(tempmatrix2(num_bands,num_bands))
        
        allocate(include_bands(num_bands)) 
        include_bands=0

        allocate(a_matrix(num_bands,num_wann,nqdiv))
        a_matrix=0.0d0
        allocate(M_Matrix(num_bands,num_bands,nntot,nqdiv))
        m_matrix=0.0d0
        allocate(atom_wan_initial(num_wann))
        atom_wan_initial=0
        allocate(atom_wan(num_wann))
        atom_wan=0
        allocate(atom_wan_distance(3,num_wann))
        atom_wan_distance=0.0d0
        allocate(eigenvalues(num_bands,nqdiv))
        eigenvalues=0.0d0

        allocate(z_wan_bnd(nfun,num_bands,npnt))
        allocate(ev_wan_bnd(nbasmpw,num_bands,npnt))
        z_wan_bnd=(0.d0,0.d0)
        ev_wan_bnd=(0.d0,0.d0)

c$$$        allocate(z_wan_all(nfun,num_wann,nqdiv))
c$$$        allocate(ev_wan_all(nbasmpw,num_wann,nqdiv))
c$$$        z_wan_all=(0.d0,0.d0)
c$$$        ev_wan_all=(0.d0,0.d0)        

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
c$$$        allocate(unitarymat(num_wann,num_wann))
c$$$        unitarymat=0.0d0
c$$$        do iwan=1, num_wann
c$$$          unitarymat(iwan,iwan)=1.0d0
c$$$        enddo

        if (maswrk) then
          write(iun,*) 'atom_wan_initial'
        endif

        do ii=1, num_wann        
          do jj=1, natom
            tempvec=matmul(transpose(gbas), tau(:,jj))-proj_site(:,ii)
            if (maswrk) then
              write(iun,'(i5, 3f12.6,i5,12f12.6)') ii,proj_site(:,ii),
     $          jj, tau(:,jj),matmul(transpose(gbas), tau(:,jj)),
     $          tempvec,
     $          dmod(tempvec(1), 1.0d0),
     $          dmod(tempvec(2), 1.0d0),
     $          dmod(tempvec(3), 1.0d0)
            endif

            if ((abs(dmod(tempvec(1), 1.0d0)) .lt. 1.0d-4) .and.
     $        (abs(dmod(tempvec(2), 1.0d0)) .lt. 1.0d-4) .and.
     $        (abs(dmod(tempvec(3), 1.0d0)) .lt. 1.0d-4)) then
              atom_wan_initial(ii)=jj
            endif
          enddo
        enddo

        
        do ii=1, num_wann
          write(179, '(i5, 3f12.6)') 
     $      , atom_wan_initial(ii), proj_site(:, ii)
        enddo

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

        
        do ind_k=1,ndim3_k(me3_k+1)
          k=n3_mpi_k(me3_k+1)+ind_k
          do ibnd=1, num_bands
            jbnd=include_bands(ibnd)
            if (ubi_num .eq. 4) then
              z_wan_bnd(:,ibnd,k)=z_qp(:,jbnd,ind_k,ispin)
              ev_wan_bnd(:,ibnd,k)=ev_qp(:,jbnd,ind_k,ispin)
            elseif (ubi_num .gt. 4) then
              z_wan_bnd(:,ibnd,k)=z_qp_2(:,jbnd,ind_k,ispin)
              ev_wan_bnd(:,ibnd,k)=ev_qp_2(:,jbnd,ind_k,ispin)
            else
              z_wan_bnd(:,ibnd,k)=z_bnd(:,jbnd,ind_k,ispin)
              ev_wan_bnd(:,ibnd,k)=ev_bnd(:,jbnd,ind_k,ispin)
            endif
          enddo
        enddo
        if(nproc_k/=1) then
          call mpi_allreduce_dcmplx
     $      (z_wan_bnd,nfun*num_bands*npnt,mpi_sum,comm_tau)
          call mpi_allreduce_dcmplx
     $      (ev_wan_bnd,nbasmpw*num_bands*npnt,mpi_sum,comm_tau)
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
          write(iun,*) "atom_wan_initial,proj_site,proj_z,proj_x,"
     $      //"proj_l,proj_m"
          do ii=1, num_wann
            write(iun,'(i5, 4x,3f12.6,4x,3f12.6,4x,3f12.6,4x,i5, i5 )')
     $        atom_wan_initial(ii), proj_site(:,ii),proj_z(:,ii),
     $        proj_x(:,ii),proj_l(ii),proj_m(ii)
          enddo
        endif

        do irk = 1, nqdiv
          k0=i_kref(irk)
          do ibnd = 1, num_bands
            jbnd=include_bands(ibnd)
            if (ubi_num .ge. 4) then
              eigenvalues(ibnd,irk)=
     $          (e_qp(jbnd,k0,ispin)-chem_pot)*evolt/2.0d0
            else
              eigenvalues(ibnd,irk)=
     $          (e_bnd(jbnd,k0,ispin)-chem_pot)*evolt/2.0d0
            endif

          enddo
        enddo


        if (maswrk)
     $    write(iun,*) "Calculating projections on to trial functions"
        if (irel .le. 1) then
          call wannier_amn_from_win(ispin)
          
!     Write A_Matrix to file!! do it later
          if (maswrk) then
            open(unit = 992, file = trim(seed_name)//'.amn')
            write(992,'(3i8)')
            write(992,'(3i8)') num_bands, nqdiv, num_wann
            do irk = 1, nqdiv
              do iwan = 1, num_wann
                do iband = 1, num_bands
                  write(992,'(3i8,2x,2F17.12)') iband, iwan, irk,
     $              A_Matrix(iband,iwan,irk)
                enddo
              enddo
            enddo
            close(unit = 992)
          endif
        else
          call cal_mtorbital_bandprojection(ispin)

          allocate(alimj(5,num_wann))
          alimj=0
          if (maswrk) then
            open(unit=111, file='wannier.irel2')
            do iwan=1, num_wann
              read(111, *) (alimj(ii,iwan),ii=1, 5)
            enddo
          endif
          call mpi_bcast(alimj,num_wann*5,mpi_integer,0,
     $      mpi_comm_world,mpierr)          
          do iwan=1, num_wann
            proj_site(:,iwan)
     $        =matmul(transpose(gbas), tau(:,alimj(1,iwan)))
          enddo
!     renormalize and diagonalization
          ! renormalize
          call poorman_wan_proj(ispin)
        endif

!     Allocate M_Matrix


        if (maswrk) write(iun,*) "Calculating overlaps"
        call wannier_mmn(ispin)
        write(ch, '(a, a)') trim(seed_name), '.mmn'
        call timel("Writing mmn to "//trim(ch)//' in eV')
        
        if (maswrk) then
          open(unit = 990, file = trim(ch), form = 'formatted')
          write(990,*) 
          write(990,'(3i8)') num_bands, nqdiv, nntot
          do irk = 1, nqdiv
            do n = 1, nntot
              write(990,'(5i8)') irk, nnlist(irk,n), nncell(:,irk,n)
              do iband2 = 1, num_bands
                do iband1 = 1, num_bands
                  write(990,'(2F17.12)') M_matrix(iband1,iband2,n,irk)
                enddo
              enddo
            enddo
          enddo
          close(990)

        endif

        

        write(ch, '(a, a)') trim(seed_name), '.eig'
        call timel("Writing eigenvalues to "//trim(ch)//' in eV')

        if (maswrk) then
          open (unit = 990, file = trim(ch), form = 'formatted')
          do irk = 1, nqdiv
            do iband = 1, num_bands
              write(990,'(2i8,2x,F17.12)') iband, irk,
     $          eigenvalues(iband,irk)
            enddo
          enddo
          close(unit = 990)
        endif

        call timel("Calling wannier_run")

        if (maswrk) then
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
            write(iun,'(3i5, 4f12.6,4x,2(3f12.6,2x))') ii,
     $        atom_wan_initial(ii),atom_wan(ii),atom_wan_distance(:,ii),
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

          write(iun,*) 'location in the wigner_seitz cell'
          do iwan=1, num_wann        
            distance=10000.0d0
            do ii=-15, 15
              do jj=-15, 15
                do kk=-15, 15
                  tempvec=
     $              wann_centers(:,iwan)
     $              -matmul(rbas, (/dble(ii), dble(jj), dble(kk)/))
     $              *(par*bohr)
                  temp_distance=dsqrt(sum(tempvec**2))
                  if (distance .gt. temp_distance) then
                    distance=temp_distance
                    tempintarr=(/ii,jj,kk/)
                  endif
                enddo
              enddo
            enddo
            write(iun,'(4i5,f12.6)') iwan, tempintarr, distance
          enddo

          write(iun,'(a)') 'calculate v_matrix'
          
          v_matrix=0.0d0
          do irk=1,nqdiv
            do j=1,num_wann
              tempvec2=matmul(transpose(gbas),
     $          wann_centers(:,j)/(par*bohr)-tau(:,atom_wan(j)))
              kdottau=2.0d0*pi*
     $          sum(tempvec2*kpt_latt(:,irk))
              do i=1,num_wann
                indtmp=0
                do m=1,num_bands            
                  if (lwindow(m,irk) .eq. .true.) then
                    indtmp=indtmp+1
                    v_matrix(m,j,irk)=
     $                v_matrix(m,j,irk)
     $                +u_matrix_opt(indtmp,i,irk)
c$$$     $                *u_matrix(i,j,irk)*cdexp(ai*kdottau)
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
                  if (lwindow(ibnd,irk) .eq. .true.) then
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
          
c$$$  write(iun,*) 'wan_g_overlap test'

c$$$  do ii=1, num_wann
c$$$  dblesum=0.0d0
c$$$  do irk=1,nqdiv
c$$$  do ibnd=1, num_bands
c$$$  if (lwindow(ibnd,irk) .eq. .true.) then                
c$$$  dblesum=dblesum+cdabs(v_matrix(ibnd,ii,irk))**2
c$$$  endif
c$$$  enddo
c$$$  enddo
c$$$  write(iun,*) ii, sum(cdabs(wan_g_overlap(ii,:))**2),
c$$$  $        dblesum/dble(nqdiv)
c$$$  enddo

c$$$  ! v_matrix_test
c$$$  write(iun,*) 'v_matrix_test'
c$$$  tempmatrix1=
c$$$  $      matmul(dconjg(transpose(v_matrix(:,:,1))), v_matrix(:,:,1))
c$$$  tempmatrix2=
c$$$  $      matmul(v_matrix(:,:,1),dconjg(transpose(v_matrix(:,:,1))))
c$$$  write(iun,*) '<v+|v>'
c$$$  do ii=1, num_wann
c$$$  write(iun,'(100(2f10.5, 2x))')
c$$$  $        (tempmatrix1(ii,jj), jj=1, num_wann)
c$$$  enddo
c$$$  
c$$$  write(iun,*) '<v+|v>'
c$$$  do ii=1, num_bands
c$$$  write(iun,'(100(2f10.5, 2x))')
c$$$  $        (tempmatrix2(ii,jj), jj=1, num_bands)
c$$$  enddo
c$$$  
c$$$  write(iun,*) '|v>'
c$$$  do ii=1, num_wann
c$$$  write(iun,*) ii, sum(cdabs(v_matrix(:,ii,1))**2)
c$$$  enddo
c$$$  do ii=1, num_bands
c$$$  write(iun,*) ii, sum(cdabs(v_matrix(ii,:,1))**2)
c$$$  enddo

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
     $        atom_wan_initial(indloc1(1)),indloc1,
     $        tempdble2(indloc1(1)),
     $        atom_wan_initial(indloc2(1)),indloc2,
     $        tempdble2(indloc2(1)),
     $        atom_wan_initial(indloc3(1)),indloc3,
     $        tempdble2(indloc3(1)),
     $        atom_wan_initial(indloc4(1)),indloc4,
     $        tempdble2(indloc4(1)),
     $        atom_wan_initial(indloc5(1)),indloc5,
     $        tempdble2(indloc5(1))
          enddo
        endif

c$$$        
c$$$        do iatom_c=1, natom_loc
c$$$          iatom=iat_loc(iatom_c)        
c$$$          isort=is(iatom)
c$$$          atom_count=0
c$$$          do iwan=1, num_wann
c$$$            if (atom_wan(iwan) .eq. iatom) then
c$$$              atom_count=atom_count+1
c$$$              wan_correlated(atom_count,iatom_c)=iwan
c$$$            endif
c$$$          enddo
c$$$        enddo


c$$$        call dmft_projector(wan_correlated,ispin)
        
        call za_wan_for_u(ispin)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


c$$$!     cal nelectron
c$$$
c$$$        if (maswrk) then
c$$$          nelectron=0.0
c$$$          do ik=1, nqdiv
c$$$            do ibnd=1, num_bands
c$$$              do iwan=1, num_wann
c$$$                if (eigenvalues(ibnd,ik) .lt. 0.0d0) then
c$$$                  nelectron=nelectron+1.0d0/dble(nqdiv)*
c$$$     $              cdabs(v_matrix(ibnd,iwan,ik))**2
c$$$                endif
c$$$              enddo
c$$$            enddo
c$$$          enddo
c$$$          write(iun,*) 'nelectron=', nelectron
c$$$        endif
c$$$        
c$$$        call mpi_bcast(nelectron,1,mpi_double_precision,
c$$$     $    0,mpi_comm_world,mpierr)
        
        if (maswrk) then
          write(ch, '(a, a)') trim(seed_name), '.umatrix'
          open(unit = 990, file = trim(ch), form = 'formatted')
          write(990,*) 
          write(990,'(3i8)') num_wann, num_wann, nqdiv
          do irk = 1, nqdiv
            do iband1 = 1, num_wann
              do iband2 = 1, num_wann
                write(990,'(2F17.12)') u_matrix(iband1,iband2,irk)
              enddo
            enddo
          enddo
          close(990)
          write(ch, '(a, a)') trim(seed_name), '.umatrix_opt'
          open(unit = 990, file = trim(ch), form = 'formatted')
          write(990,*) 
          write(990,'(3i8)') num_bands, num_wann, nqdiv
          do irk = 1, nqdiv
            do iband1 = 1, num_bands
              do iband2 = 1, num_wann
                write(990,'(2F17.12)') u_matrix_opt(iband1,iband2,irk)
              enddo
            enddo
          enddo
          close(990)
          write(ch, '(a, a)') trim(seed_name), '.vmatrix'
          open(unit = 990, file = trim(ch), form = 'formatted')
          write(990,*) 
          write(990,'(3i8)') num_bands, num_wann, nqdiv
          do irk = 1, nqdiv
            do iband1 = 1, num_bands
              do iband2 = 1, num_wann
                write(990,'(2F17.12)') v_matrix(iband1,iband2,irk)
              enddo
            enddo
          enddo
          close(990)

          open(unit=334, file=trim(seed_name)//'_include_bands.dat')
          do ibnd=1, num_bands
            write(334,*) include_bands(ibnd)
          enddo
          close(334)          
          
c$$$          write(ch, '(a, a)') trim(seed_name), '.dos_correlated'
c$$$          open(unit = 990, file = trim(ch), form = 'formatted')
c$$$          do ii = 1, 2
c$$$            do iwan = 1, 5
c$$$              do iband = 1, 1600
c$$$                write(990,'(2F17.12)') iband*0.01d0-10.0d0,
c$$$     $            dos_correlated(iband,iwan,ii)
c$$$              enddo
c$$$            enddo
c$$$          enddo
c$$$          close(990)

          open(unit=990, file='gw+dmft.ini')
          write(990,'(3i5)') ndiv_fi
          write(990,'(i5)') 
          write(990,'(i5)') maxel_loc
          write(990,'(i5)') natom_loc
c$$$          write(990,'(100i5)') wan_correlated
c$$$          write(990,'(i5)') num_wann
c$$$          write(990,'(f12.6)') nelectron
          close(990)

          
!     cal occupation
c$$$
          
          write(993,*) "e_projected and rho_projected from coarse grid"
          do iwan=1, num_wann
            e_projected=0.0d0
            rho_projected=0.0d0
            do ik=1, nqdiv
              do ibnd=1, num_bands
                e_projected=e_projected
     $            +eigenvalues(ibnd,ik)
     $            *cdabs(v_matrix(ibnd,iwan,ik))**2
     $            /dble(nqdiv)
                if (eigenvalues(ibnd,ik) .lt. 0.0d0) then
                  rho_projected=rho_projected
     $              +cdabs(v_matrix(ibnd,iwan,ik))**2
     $              /dble(nqdiv)
                endif
              enddo
            enddo
            write(993,*) e_projected, rho_projected
          enddo


!     cal hmatrix and rhomat
          allocate(hmat(num_wann,num_wann))
          hmat=0.0d0
          allocate(rhomat(num_wann,num_wann))
          rhomat=0.0d0

          do iwan=1, num_wann
            do jwan=1, num_wann            
              do ik=1, nqdiv
                do ibnd=1, num_bands
                  hmat(iwan,jwan)=hmat(iwan,jwan)
     $              +eigenvalues(ibnd,ik)
     $              *dconjg(v_matrix(ibnd,iwan,ik))
     $              *v_matrix(ibnd,jwan,ik)                  
     $              /dble(nqdiv)
                  if (eigenvalues(ibnd,ik) .lt. 0.0d0) then
                    rho_projected=rho_projected
     $                +dconjg(v_matrix(ibnd,iwan,ik))
     $                *v_matrix(ibnd,jwan,ik)                       
     $                /dble(nqdiv)
                  endif
                enddo
              enddo
            enddo
          enddo
          write(993,*) "hmat from coarse grid"          
          do iwan=1, num_wann
            write(993, '(1000(2f12.6, 4x))')
     $        (hmat(iwan,jwan), jwan=1, num_wann)
          enddo
          write(993,*) "rhomat from coarse grid"          
          do iwan=1, num_wann
            write(993, '(1000(2f12.6, 4x))')
     $        (rhomat(iwan,jwan), jwan=1, num_wann)
          enddo

          write(993,*) "hmat and rhomat test"
          hmat=0.0d0
          rhomat=0.0d0
          do iwan=1, num_wann
            do jwan=1, num_wann            
              do ik=1, nqdiv
                do ibnd=1, num_bands
                  hmat(iwan,jwan)=hmat(iwan,jwan)
     $              +eigenvalues(ibnd,ik)
     $              *(v_matrix(ibnd,iwan,ik))
     $              *dconjg(v_matrix(ibnd,jwan,ik))                  
     $              /dble(nqdiv)
                  if (eigenvalues(ibnd,ik) .lt. 0.0d0) then
                    rho_projected=rho_projected
     $                +(v_matrix(ibnd,iwan,ik))
     $                *dconjg(v_matrix(ibnd,jwan,ik))
     $                /dble(nqdiv)
                  endif
                enddo
              enddo
            enddo
          enddo
          write(993,*) "hmat from coarse grid"          
          do iwan=1, num_wann
            write(993, '(1000(2f12.6, 4x))')
     $        (hmat(iwan,jwan), jwan=1, num_wann)
          enddo
          write(993,*) "rhomat from coarse grid"          
          do iwan=1, num_wann
            write(993, '(1000(2f12.6, 4x))')
     $        (rhomat(iwan,jwan), jwan=1, num_wann)
          enddo
          deallocate(hmat,rhomat)
        endif
        

        if ((ubi_num .ge. 5) .and. (postprocessflag .eq. 0)) then
          allocate(g_out1(nbndf, nbndf))
          allocate(g_out2(nbndf, nbndf))
          allocate(g_out3(num_bands, num_bands))
          allocate(g_out4(num_bands, num_wann))
          allocate(g_out5(num_wann, num_wann))

          do ind_k=1,ndim3_k(me3_k+1)
            k=n3_mpi_k(me3_k+1)+ind_k
            ndim=n_bnd(k,ispin)              

            do ind_omega=1,ndim3_omega
              i_omega=me3_tau*ndim3_omega+ind_omega-1


              if (i_omega .lt. 10) then
                write(i_omega_ch, '(i1)') i_omega
              elseif (i_omega .lt. 100) then
                write(i_omega_ch, '(i2)') i_omega
              elseif (i_omega .lt. 1000) then
                write(i_omega_ch, '(i3)') i_omega
              elseif (i_omega .lt. 10000) then
                write(i_omega_ch, '(i4)') i_omega
              elseif (i_omega .lt. 100000) then
                write(i_omega_ch, '(i5)') i_omega
              endif
              g_out1=0.0d0
              g_out2=0.0d0
              g_out1=
     $          dcmplx(g_full_omega(:,:,1,ind_omega,ind_k,ispin),
     $          g_full_omega(:,:,2,ind_omega,ind_k,ispin))

!     green's function in qp basis
c$$$              if (ubi_num .eq. 4) then
c$$$                call zgemm('n','n',ndim,ndim,
c$$$     $            ndim,(1.d0,0.d0),
c$$$     &            g_out1,nbndf,
c$$$     $            q_qp(1,1,ind_k,ispin),nbndf,(0.d0,0.d0),
c$$$     $            g_out2,nbndf)              
c$$$                g_out1=0.0d0
c$$$                call zgemm('c','n',ndim,ndim,
c$$$     $            ndim,(1.d0,0.d0),
c$$$     &            q_qp(1,1,ind_k,ispin), nbndf, g_out2,nbndf,
c$$$     $            (0.d0,0.d0),
c$$$     $            g_out1,nbndf) 
c$$$              elseif (ubi_num .gt. 4) then             
              call zgemm('n','n',ndim,ndim,
     $          ndim,(1.d0,0.d0),
     &          g_out1,nbndf,
     $          q_qp_2(1,1,ind_k,ispin),nbndf,(0.d0,0.d0),
     $          g_out2,nbndf)              
              g_out1=0.0d0
              call zgemm('c','n',ndim,ndim,
     $          ndim,(1.d0,0.d0),
     &          q_qp_2(1,1,ind_k,ispin), nbndf, g_out2,nbndf,
     $          (0.d0,0.d0),
     $          g_out1,nbndf) 
c$$$  endif
              
!     green's function in selected qp basis
              g_out3=0.0d0

              do ibnd=1, num_bands
                ii=include_bands(ibnd)
                do jbnd=1, num_bands
                  jj=include_bands(jbnd)
                  g_out3(ibnd,jbnd)=g_out1(ii,jj)
                enddo
              enddo

!     green's function in tau basis for all k in the first BZ

              do ikstar=1, k_star(k)
                kstar=k_list(ikstar,k)
                
                if (kstar .lt. 10) then
                  write(kstar_ch, '(i1)') kstar
                elseif (kstar .lt. 100) then
                  write(kstar_ch, '(i2)') kstar
                elseif (kstar .lt. 1000) then
                  write(kstar_ch, '(i3)') kstar
                elseif (kstar .lt. 10000) then
                  write(kstar_ch, '(i4)') kstar
                elseif (kstar .lt. 100000) then
                  write(kstar_ch, '(i5)') kstar
                endif

                g_out4=0.0d0                
                g_out5=0.0d0
                call zgemm('n','n',num_bands,num_wann,
     $            num_bands,(1.d0,0.d0),
     $            g_out3,num_bands,
     $            v_matrix(1,1,kstar),num_bands,(0.d0,0.d0),
     $            g_out4,num_bands)
                call zgemm('c','n',num_wann,num_wann,
     $            num_bands,(1.d0,0.d0),
     $            v_matrix(1,1,kstar),num_bands,
     $            g_out4,num_bands,(0.d0,0.d0),
     $            g_out5,num_wann)
                
                open(4,file=trim(seed_name)//'_gfull_'
     $            //trim(kstar_ch)
     $            //'_'//trim(i_omega_ch)//'.dat')
                do ii=1, num_wann
                  do jj=1, num_wann
                    write(4,'(2i5, 2f20.12)') ii, jj, 
     $                g_out5(ii, jj)/rydberg
                  enddo
                enddo
                close(4)
c$$$
c$$$                open(4,file=trim(seed_name)//'_gfull_'
c$$$     $            //trim(kstar_ch)
c$$$     $            //'_'//trim(i_omega_ch)//'.dat', form='unformatted')
c$$$                write(4) g_out5/rydberg
c$$$                close(4)
              enddo
            enddo
          enddo

c$$$          if(nproc/=1)
c$$$     $      call mpi_allreduce_dcmplx
c$$$     $      (g_out5,num_wann**2*nqdiv*(n_omega+1),mpi_sum,
c$$$     $      mpi_comm_world);


          deallocate(g_out1)
          deallocate(g_out2)
          deallocate(g_out3)
          deallocate(g_out4)
          deallocate(g_out5)


! projection of hf hamiltonian

          allocate(g_out1(nbndf, nbndf))
          allocate(g_out2(nbndf, nbndf))
          allocate(g_out3(num_bands, num_bands))
          allocate(g_out4(num_bands, num_wann))
          allocate(g_out5(num_wann, num_wann))
          allocate(hfmat(num_wann, num_wann,nqdiv))
          hfmat=0.0d0
          
          do ind_k=1,ndim3_k(me3_k+1)
            k=n3_mpi_k(me3_k+1)+ind_k
            ndim=n_bnd(k,ispin)              
            
            g_out1=0.0d0
            g_out2=0.0d0
            do ibnd=1, nbndf
              g_out1(ibnd,ibnd)=e_bnd(ibnd,k,ispin)
            enddo
            
            call zgemm('n','n',ndim,ndim,
     $        ndim,(1.d0,0.d0),
     &        g_out1,nbndf,
     $        q_qp_2(1,1,ind_k,ispin),nbndf,(0.d0,0.d0),
     $        g_out2,nbndf)              
            g_out1=0.0d0
            call zgemm('c','n',ndim,ndim,
     $        ndim,(1.d0,0.d0),
     &        q_qp_2(1,1,ind_k,ispin), nbndf, g_out2,nbndf,
     $        (0.d0,0.d0),
     $        g_out1,nbndf) 
c$$$  endif
            
!     green's function in selected qp basis
            g_out3=0.0d0
            
            do ibnd=1, num_bands
              ii=include_bands(ibnd)
              do jbnd=1, num_bands
                jj=include_bands(jbnd)
                g_out3(ibnd,jbnd)=g_out1(ii,jj)
              enddo
            enddo
            
!     green's function in tau basis for all k in the first BZ
            
            do ikstar=1, k_star(k)
              kstar=k_list(ikstar,k)
              
c$$$  if (kstar .lt. 10) then
c$$$  write(kstar_ch, '(i1)') kstar
c$$$  elseif (kstar .lt. 100) then
c$$$  write(kstar_ch, '(i2)') kstar
c$$$  elseif (kstar .lt. 1000) then
c$$$  write(kstar_ch, '(i3)') kstar
c$$$  elseif (kstar .lt. 10000) then
c$$$  write(kstar_ch, '(i4)') kstar
c$$$  elseif (kstar .lt. 100000) then
c$$$  write(kstar_ch, '(i5)') kstar
c$$$  endif
              
              g_out4=0.0d0                
              g_out5=0.0d0
              call zgemm('n','n',num_bands,num_wann,
     $          num_bands,(1.d0,0.d0),
     $          g_out3,num_bands,
     $          v_matrix(1,1,kstar),num_bands,(0.d0,0.d0),
     $          g_out4,num_bands)
              call zgemm('c','n',num_wann,num_wann,
     $          num_bands,(1.d0,0.d0),
     $          v_matrix(1,1,kstar),num_bands,
     $          g_out4,num_bands,(0.d0,0.d0),
     $          g_out5,num_wann)
              
              hfmat(:,:,kstar)=g_out5
              
              
            enddo
          enddo
          call mpi_allreduce_dcmplx
     $      (hfmat,num_wann**2*nqdiv,mpi_sum,comm_tau)          
          if (maswrk) then
            open(unit=112, file=trim(seed_name)//'_hf.mat')
            do irk=1, nqdiv
              do jwan=1, num_wann
                do iwan=1, num_wann
                  write(112, '(3i5, 2f20.12)') iwan,jwan,irk,
     $              hfmat(iwan,jwan,irk)*rydberg
                enddo
              enddo
            enddo
            close(112)
          endif

          
c$$$          if(nproc/=1)
c$$$     $      call mpi_allreduce_dcmplx
c$$$     $      (g_out5,num_wann**2*nqdiv*(n_omega+1),mpi_sum,
c$$$     $      mpi_comm_world);


          deallocate(g_out1)
          deallocate(g_out2)
          deallocate(g_out3)
          deallocate(g_out4)
          deallocate(g_out5)
          deallocate(hfmat)
        endif


        call wan90_local_z_factor
        
c$$$        deallocate(wan_correlated)
        deallocate(ldim_wan)
        deallocate(rdim_wan)
        deallocate(lval_wan)
        deallocate(radial_coeff_wan)
        
        deallocate(tempdble)
        deallocate(tempdble2)
c$$$        deallocate(tempdcmplx)
c$$$        deallocate(tempmatrix1)
c$$$        deallocate(tempmatrix2)
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
        deallocate(atom_wan_initial)
        deallocate(atom_wan)
        deallocate(atom_wan_distance)
        deallocate(eigenvalues)
        deallocate(z_wan_bnd)
        deallocate(ev_wan_bnd)
c$$$        deallocate(z_wan_all)
c$$$        deallocate(ev_wan_all)        
        deallocate(U_matrix)
        deallocate(wan_g_overlap)
        deallocate(v_matrix)
        deallocate(U_matrix_opt)
        deallocate(lwindow)
        deallocate(wann_centers)
        deallocate(wann_spreads)
c$$$  deallocate(unitarymat)
c$$$  deallocate(kpath_point)
c$$$  deallocate(kpath)
c$$$  deallocate(tempval)

      enddo


      end


!================================!
      subroutine wigner_seitz(nrpts,ndegen, irvec)
!================================!
!     Calculates a grid of lattice vectors r that fall inside (and eventually  !
!     on the surface of) the Wigner-Seitz supercell centered on the            ! 
!     origin of the Bravais superlattice with primitive translations           !
!     ndiv(1)*a_1, ndiv(2)*a_2, and ndiv(3)*a_3                       !
!==========================================================================!
      use solid_mod
      use manager_mod
      use units_mod
      use parallel_mod
      implicit none


!     irvec(i,irpt)     The irpt-th Wigner-Seitz grid point has components
!     irvec(1:3,irpt) in the basis of the lattice vectors
!     ndegen(irpt)      Weight of the irpt-th point is 1/ndegen(irpt)
!     nrpts             number of Wigner-Seitz grid points

      integer, intent(out) :: nrpts, ndegen(8*nqdiv), irvec(3, 8*nqdiv)

      integer       :: ndiff(3),rpt_origin
      double precision :: dist(125),tot,dist_min,real_metric(3,3)
      integer       :: n1,n2,n3,i1,i2,i3,icnt,i,j,ir

      ndegen=0
      irvec=0

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

      nrpts = 0  
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
              nrpts = nrpts + 1  
c$$$  if(.not. count_pts) then
              ndegen(nrpts)=0
              do i=1,125
                if(abs(dist(i)-dist_min).lt.1.0d-7)
     $            ndegen(nrpts)=ndegen(nrpts)+1
              end do
              irvec(1, nrpts) = n1  
              irvec(2, nrpts) = n2   
              irvec(3, nrpts) = n3   
!     
!     Remember which grid point r is at the origin
!     
              if (n1==0 .and. n2==0 .and. n3==0) rpt_origin=nrpts
c$$$  endif
            end if
!     n3
          enddo
!     n2
        enddo
!     n1
      enddo
!     
c$$$  if(count_pts) then
c$$$  c$$$        if (timing_level>1.and.on_root)
c$$$  c$$$     $    call io_stopwatch('postw90_common: wigner_seitz',2)
c$$$  return
c$$$  end if

c$$$  if(iprint>=3.and.on_root) then
      if (maswrk) then
        write(iun,'(1x,i4,a,/)') nrpts,
     $    ' lattice points in Wigner-Seitz supercell:'
        do ir=1,nrpts
          write(iun,'(4x,a,3(i3,1x),a,i2)') '  vector ',irvec(1,ir),
     $      irvec(2,ir),irvec(3,ir),'  degeneracy: ',ndegen(ir)
        enddo
      endif
!     Check the "sum rule"
      tot = 0.0d0
      do ir = 1, nrpts  
!     
!     Corrects weights in Fourier sums for R-vectors on the boundary of the 
!     W-S supercell 
!     
        tot=tot+1.0d0/dble(ndegen(ir))
      enddo
      if (abs(tot-dble(ndiv(1)*ndiv(2)*ndiv(3))) > 1.0d-8)
     $  call ending
     $  ('ERROR in wigner_seitz: error in finding Wigner-Seitz points')

c$$$  if (timing_level>1.and.on_root)
c$$$  $  call io_stopwatch('postw90_common: wigner_seitz',2)
c$$$  
c$$$  return  
      end subroutine wigner_seitz


      subroutine wannier_interpolation_path(value_coarse, infilename,
     $  kpath_size,kpath)

      use wannier90_mod
      use solid_mod
      use parallel_mod
      use units_mod

      implicit none
      include 'mpif.h'


      integer,intent(in) :: kpath_size
      double precision, intent(in) :: value_coarse(num_bands, nqdiv),
     $  kpath(3,kpath_size)
      character(len=*),intent(in) :: infilename


      integer :: ndim_wann2(nproc_tau),n_mpi_wann2(nproc_tau),me1,
     $  iwan,jwan,ind_wann,wanntot,irk,iband,
     $  ndim_kpath(nproc_k),n_mpi_kpath(nproc_k),
     $  i1,j1,ik,ndim_kpath_max,
     $  ndim_nrpts(nproc_k),n_mpi_nrpts(nproc_k),ind_k,i2
      double precision :: 
     $  h_kpath_wan_val(num_wann,kpath_size)
      complex*16 :: h_k_wan_co(nqdiv),ai,
     $  h_r_wan(nrpts)

      complex*16, allocatable :: h_kpath_wan(:,:,:)



      ai=dcmplx(0.0d0, 1.0d0)

      me1=me+1

      call size_shift_par(num_wann**2,nproc_tau,ndim_wann2,n_mpi_wann2)
      call size_shift_par(kpath_size,nproc_k,ndim_kpath,n_mpi_kpath)        
      ndim_kpath_max=maxval(ndim_kpath)
      allocate(h_kpath_wan(num_wann,num_wann,ndim_kpath_max))
      h_kpath_wan=0.0d0
      call size_shift_par(nrpts,nproc_k,ndim_nrpts,n_mpi_nrpts)        

      h_r_wan=0.0d0
      h_kpath_wan=0.0d0
      h_kpath_wan_val=0.0d0
      
      do ind_wann=1, ndim_wann2(me3_tau+1)
        wanntot=n_mpi_wann2(me3_tau+1)+ind_wann
        
        iwan=(wanntot-1)/num_wann+1
        jwan=wanntot-(iwan-1)*num_wann

        h_r_wan=0.0d0
        h_k_wan_co=0.0d0
        do ind_k=1, ndim3_kk(me3_k+1)
          irk=n3_mpi_kk(me3_k+1)+ind_k
          do iband=1, num_bands                
            if (lwindow(iband,irk) .eq. .true.) then          
              h_k_wan_co(irk)=h_k_wan_co(irk)+
     $          value_coarse(iband,irk)*v_matrix(iband,jwan,irk)
     $          *dconjg(v_matrix(iband,iwan,irk))
            endif
          enddo
        enddo

        call mpi_allreduce_dcmplx(h_k_wan_co,nqdiv,mpi_sum,comm_tau)

!     fourier transform to r space

        do i2=1, ndim_nrpts(me3_k+1)
          i1=n_mpi_nrpts(me3_k+1)+i2
          do j1=1, nqdiv
            h_r_wan(i1)=h_r_wan(i1)
     $        +h_k_wan_co(j1)*
     $        cdexp(-2.0*ai*pi*
     $        sum(irvec(:,i1)*kpt_latt(:,j1)))
     $        /dble(nqdiv)
          enddo
        enddo

        call mpi_allreduce_dcmplx(h_r_wan,nrpts,mpi_sum,comm_tau)

!     band interpolation along a path
        do i2=1, ndim_kpath(me3_k+1)
          i1=n_mpi_kpath(me3_k+1)+i2
          do j1=1, nrpts
            h_kpath_wan(iwan,jwan,i2)
     $        =h_kpath_wan(iwan,jwan,i2)
     $        +h_r_wan(j1)*
     $        cdexp(2.0*ai*pi*
     $        sum(kpath(:,i1)*irvec(:,j1))
     $        )/dble(ndegen(j1))
          enddo
        enddo
      enddo

      call mpi_allreduce_dcmplx(h_kpath_wan, ndim_kpath_max*num_wann**2,
     $  mpi_sum, comm_k)

      do ik=1, ndim_kpath(me3_k+1)
        irk=n_mpi_kpath(me3_k+1)+ik
        call eig_val_solver(num_wann,num_wann,
     $    h_kpath_wan(:,:,ik),h_kpath_wan_val(:,irk))
      enddo

      call mpi_allreduce_dble(h_kpath_wan_val,kpath_size*num_wann,
     $  mpi_sum,comm_tau)

!     projected hamiltonian

      if (maswrk) then
        open(unit = 990, file = trim(infilename)//'.interpolatepath',
     $    form = 'formatted')
        do iwan=1, num_wann
          do irk=1, kpath_size
            write(990,*) irk, h_kpath_wan_val(iwan,irk)
          enddo
        enddo
        close(990)
      endif
      deallocate(h_kpath_wan)
      end

      subroutine wannier_interpolation_grid(value_coarse, 
     $  ndiv_fi, nqdiv_fi,h_k_wan_val,h_k_wan_vec,
     $  ndim_kk_fi,n_mpi_kk_fi,ndim_kk_fi_max)

      use wannier90_mod
      use parallel_mod
      use solid_mod
      use units_mod

      implicit none
      include 'mpif.h'
      integer, intent(in) :: ndiv_fi(3), nqdiv_fi,
     $  ndim_kk_fi(nproc),n_mpi_kk_fi(nproc),ndim_kk_fi_max
      double precision, intent(in) :: value_coarse(num_bands, nqdiv)
      double precision,intent(out) :: 
     $  h_k_wan_val(num_wann,ndim_kk_fi_max)
      complex*16,intent(out) :: 
     $  h_k_wan_vec(num_wann,num_wann,ndim_kk_fi_max)



      integer :: me1,iwan,jwan,ind_wann,wanntot,irk,iband,
     $  ndim_wann2(nproc),n_mpi_wann2(nproc),i1,i2,i3,j1,ik,k1,k2,k3,
     $  jj,ii
      complex*16 :: h_k_wan_co(nqdiv),ai,h_r_wan(nrpts)

      ai=dcmplx(0.0d0, 1.0d0)
      me1=me+1

      h_r_wan=0.0d0
      h_k_wan_vec=0.0d0
      h_k_wan_val=0.0d0      
      
      do iwan=1, num_wann
        do jwan=1, num_wann
          
          h_r_wan=0.0d0
          h_k_wan_co=0.0d0
          do irk=1, nqdiv
            do iband=1, num_bands                
              if (lwindow(iband,irk) .eq. .true.) then          
                h_k_wan_co(irk)=h_k_wan_co(irk)+
     $            value_coarse(iband,irk)*v_matrix(iband,jwan,irk)
     $            *dconjg(v_matrix(iband,iwan,irk))
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
            k1=(irk-1)/(ndiv_fi(3)*ndiv_fi(2))+1-1
            k2=(irk-(k1)*ndiv_fi(3)*ndiv_fi(2)-1)/ndiv_fi(3)+1-1
            k3=irk-(k1)*ndiv_fi(2)*ndiv_fi(3)-(k2)*ndiv_fi(3)-1
            
            do j1=1, nrpts
              h_k_wan_vec(iwan,jwan,ik)
     $          =h_k_wan_vec(iwan,jwan,ik)
     $          +h_r_wan(j1)*
     $          cdexp(2.0*ai*pi*
     $          sum(
     $          (/dble(k1)/dble(ndiv_fi(1)),
     $          dble(k2)/dble(ndiv_fi(2)),
     $          dble(k3)/dble(ndiv_fi(3))/)
     $          *irvec(:,j1))
     $          )/dble(ndegen(j1))
            enddo
          enddo
        enddo
      enddo

      do ik=1, ndim_kk_fi(me1)
        call eig_val_solver(num_wann,num_wann,
     $    h_k_wan_vec(:,:,ik),h_k_wan_val(:,ik))
      enddo



      end      



      subroutine cal_e_wan_all(value_coarse, e_wan_all)

      use wannier90_mod
      use parallel_mod
      use solid_mod
      use units_mod

      implicit none
      include 'mpif.h'
      double precision, intent(in) :: value_coarse(num_bands, nqdiv)
      complex*16 :: e_wan_all(num_wann)
      
      integer :: me1,iwan,jwan,ind_wann,wanntot,irk,iband,
     $  ndim_wann2(nproc),n_mpi_wann2(nproc),i1,i2,i3,j1,ik,k1,k2,k3,
     $  jj,ii
      complex*16 :: h_k_wan_co(nqdiv),ai

      ai=dcmplx(0.0d0, 1.0d0)
      me1=me+1

      e_wan_all=0.0d0

      do iwan=1, num_wann
        
        h_k_wan_co=0.0d0
        do irk=1, nqdiv
          do iband=1, num_bands                
            if (lwindow(iband,irk) .eq. .true.) then          
              h_k_wan_co(irk)=h_k_wan_co(irk)+
     $          value_coarse(iband,irk)*v_matrix(iband,iwan,irk)
     $          *dconjg(v_matrix(iband,iwan,irk))
            endif
          enddo
        enddo
        
!     fourier transform to r space
        
        do j1=1, nqdiv
          e_wan_all(iwan)=e_wan_all(iwan)
     $      +h_k_wan_co(j1)
     $      /dble(nqdiv)
        enddo

      enddo

      end      



c$$$  !=======================================!
c$$$  subroutine overwrite_wannier_chkpt
c$$$  
c$$$  use wannier90_mod
c$$$  use solid_mod 
c$$$  use units_mod
c$$$  use parallel_mod
c$$$  use manager_mod
c$$$  
c$$$  implicit none
c$$$  
c$$$  integer :: chk_unit,i,j,k,l,nkp,ierr,irk,iband1,iband2
c$$$  character*33 :: header, ch
c$$$  
c$$$  if (.not. maswrk) return
c$$$  chk_unit=980
c$$$  
c$$$  if (.not. allocated(ndimwin)) allocate(ndimwin(nqdiv))
c$$$  ndimwin=0
c$$$  open(unit=chk_unit,file=trim(seed_name)//'.chk',
c$$$  $  status='old',form='unformatted')
c$$$  
c$$$  
c$$$  
c$$$  
c$$$  read(chk_unit) header
c$$$  read(chk_unit) 
c$$$  read(chk_unit) 
c$$$  read(chk_unit) 
c$$$  read(chk_unit) 
c$$$  read(chk_unit) 
c$$$  read(chk_unit) 
c$$$  read(chk_unit) 
c$$$  read(chk_unit) 
c$$$  read(chk_unit) 
c$$$  read(chk_unit) 
c$$$  read(chk_unit) checkpoint ! checkpoint
c$$$  read(chk_unit) have_disentangled ! whether a disentanglement has been performed
c$$$  if (have_disentangled) then
c$$$  read(chk_unit) omega_invariant ! omega invariant
c$$$  read(chk_unit) 
c$$$  !     ndimwixn
c$$$  read(chk_unit) (ndimwin(nkp),nkp=1,nqdiv)
c$$$  !     U_matrix_opt
c$$$  read(chk_unit) 
c$$$  endif
c$$$  read(chk_unit) 
c$$$  !     M_matrix
c$$$  read(chk_unit)
c$$$  !     wannier_centres
c$$$  read(chk_unit) 
c$$$  read(chk_unit)
c$$$  close(chk_unit)
c$$$  
c$$$  open(unit=chk_unit,file=trim(seed_name)//'.chk',
c$$$  $  form='unformatted')
c$$$  
c$$$  write(chk_unit) header    ! Date and time from the read file
c$$$  write(chk_unit) num_bands ! Number of bands
c$$$  write(chk_unit) nexband   ! Number of excluded bands
c$$$  write(chk_unit) (exclude_bands(i),i=1,nexband) ! Excluded bands 
c$$$  write(chk_unit) ((rbas(j,i)*par*bohr,i=1,3),j=1,3) ! Real lattice
c$$$  write(chk_unit) ((gbas(j,i)*2.0d0*pi/par/bohr,i=1,3),j=1,3) ! Reciprocal lattice
c$$$  write(chk_unit) nqdiv     ! Number of k-points
c$$$  write(chk_unit) (ndiv(i),i=1,3) ! M-P grid
c$$$  write(chk_unit) ((kpt_latt(i,nkp),i=1,3),nkp=1,nqdiv) ! K-points
c$$$  write(chk_unit) nntot     ! Number of nearest k-point neighbours
c$$$  write(chk_unit) num_wann  ! Number of wannier functions
c$$$  !     Next is correct: it always print out 20 characters
c$$$  write(chk_unit) checkpoint ! Position of checkpoint
c$$$  write(chk_unit) have_disentangled ! Whether a disentanglement has been performed
c$$$  if (have_disentangled) then
c$$$  write(chk_unit) omega_invariant ! Omega invariant
c$$$  !     lwindow, ndimwin and U_matrix_opt 
c$$$  write(chk_unit) ((lwindow(i,nkp),i=1,num_bands),nkp=1,nqdiv)
c$$$  write(chk_unit) (ndimwin(nkp),nkp=1,nqdiv)
c$$$  write(chk_unit) (((u_matrix_opt(i,j,nkp),i=1,num_bands),
c$$$  $    j=1,num_wann),nkp=1,nqdiv)
c$$$  endif
c$$$  write(chk_unit) (((u_matrix(i,j,k),i=1,num_wann),j=1,num_wann),
c$$$  $  k=1,nqdiv)              ! U_matrix
c$$$  write(chk_unit) ((((m_matrix(i,j,k,l),i=1,num_wann),j=1,num_wann),
c$$$  $  k=1,nntot),l=1,nqdiv)   ! M_matrix
c$$$  write(chk_unit) ((wann_centers(i,j),i=1,3),j=1,num_wann)
c$$$  write(chk_unit) (wann_spreads(i),i=1,num_wann)
c$$$  close(chk_unit)
c$$$  
c$$$  write(iun,*) 'from overwrite routine'
c$$$  write(iun,*) checkpoint
c$$$  write(iun,*) have_disentangled
c$$$  write(iun,*) ndimwin
c$$$  write(ch, '(a, a)') trim(seed_name), '.umatrix_temp'
c$$$  open(unit = 990, file = trim(ch), form = 'formatted')
c$$$  write(990,*) 
c$$$  write(990,'(3i8)') num_wann, num_wann, nqdiv
c$$$  do irk = 1, nqdiv
c$$$  do iband1 = 1, num_wann
c$$$  do iband2 = 1, num_wann
c$$$  write(990,'(2F17.12)') u_matrix(iband1,iband2,irk)
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  close(990)
c$$$  
c$$$  end
c$$$  
c$$$  

c$$$  subroutine wannier_md(num_wann_md, ind_md, h_k_wan_md,
c$$$  $  eigenval,eigenvec)
c$$$  
c$$$  ! md: materials dependent
c$$$  use manager_mod
c$$$  use wannier90_mod
c$$$  use parallel_mod
c$$$  use solid_mod
c$$$  use units_mod
c$$$  
c$$$  implicit none
c$$$  include 'mpif.h'
c$$$  integer, intent(in) :: num_wann_md, ind_md(num_wann_md)
c$$$  double precision, intent(in) :: eigenval(num_wann)
c$$$  complex*16, intent(in) :: eigenvec(num_wann, num_wann)
c$$$  complex*16, intent(out) :: h_k_wan_md(num_wann_md,num_wann_md)
c$$$  
c$$$  
c$$$  integer :: iwan,jwan,irk,iband,ii,jj, ind1, ind2
c$$$  double precision :: 
c$$$  $  h_k_wan_val(num_wann_md), 
c$$$  
c$$$  
c$$$  
c$$$  h_k_wan_md=0.0d0
c$$$  h_k_wan_val=0.0d0
c$$$  
c$$$  do iwan=1, num_wann_md 
c$$$  do jwan=1, num_wann_md
c$$$  do iband=1, num_bands                
c$$$  h_k_wan_md(iwan, jwan)=h_k_wan_md(iwan,jwan)+
c$$$  $        eigenval(iband)
c$$$  $        *eigenvec(iband,ind_md(jwan))
c$$$  $        *dconjg(eigenvec(iband,ind_md(iwan)))
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$      
c$$$      call eig_val_solver(num_wann_md,num_wann_md,
c$$$     $  h_k_wan_md,h_k_wan_val)
c$$$      
c$$$      if (maswrk) then
c$$$        write(179,'(a, i5)')
c$$$     $    "hamiltonian and eigenvalues for "//
c$$$     $    "d_{x^2-y^2} and p orbitals", 
c$$$        do ii=1, num_wann_md
c$$$          write(179,'(f12.6, 5x, 100(f12.6, 2x))')
c$$$     $      h_k_wan_val(ii),
c$$$     $      (cdabs(h_k_wan_md(jj, ii))**2, jj=1, num_wann_md)
c$$$        enddo
c$$$     $    "getting rid of i k tau phase "//
c$$$     $    "d_{x^2-y^2} and p orbitals", 
c$$$        do ii=1, num_wann_md
c$$$          write(179,'(f12.6, 5x, 100(f12.6, 2x))')
c$$$     $      h_k_wan_val(ii),
c$$$     $      (cdabs(h_k_wan_md(jj, ii))**2, jj=1, num_wann_md)
c$$$        enddo
c$$$      endif
c$$$      
c$$$c$$$  enddo
c$$$      end      
      


c$$$      subroutine dmft_projector(wan_correlated, ispin)
c$$$      use atom_mod
c$$$      use parallel_mod
c$$$      use solid_mod
c$$$      use manager_mod
c$$$      use wannier90_mod
c$$$      implicit none
c$$$      include "mpif.h"
c$$$
c$$$      integer, intent(in) :: wan_correlated(maxel_loc, natom_loc), ispin
c$$$      integer :: ndiv_fi(3), nqdiv_fi, ndim_kk_fi(nproc),
c$$$     $  n_mpi_kk_fi(nproc),ndim_kk_fi_max, iatom_c, iatom,isort,
c$$$     $  atom_count,iwan,jwan,ik,irk,iband, ii,jj, me1,k1,k2,k3,kk
c$$$      double precision, allocatable :: 
c$$$     $  h_k_wan_val(:,:)
c$$$      complex*16, allocatable :: 
c$$$     $  h_k_wan_vec(:,:,:), e_wan_all(:), hmat(:,:,:), rhomat(:,:,:),
c$$$     $  proj(:,:,:,:)
c$$$      character*30 :: ch_kvec
c$$$
c$$$
c$$$!     calculate projector in a fine grid
c$$$      me1=me+1
c$$$      ndiv_fi=(/10,10,10/)
c$$$      nqdiv_fi=ndiv_fi(1)*ndiv_fi(2)*ndiv_fi(3)
c$$$      
c$$$      allocate(ndegen(8*nqdiv))
c$$$      ndegen=0
c$$$      allocate(irvec(3,8*nqdiv))
c$$$      irvec=0
c$$$      
c$$$      call wigner_seitz(nrpts, ndegen, irvec)
c$$$      
c$$$      if (maswrk) then
c$$$        write(iun,*) 'wigner seitz'
c$$$        do ii=1, nrpts
c$$$          write(iun,'(5i5)') ii, irvec(:,ii), ndegen(ii)
c$$$        enddo
c$$$      endif
c$$$
c$$$      call size_shift_par(nqdiv_fi,nproc,ndim_kk_fi,n_mpi_kk_fi)
c$$$      ndim_kk_fi_max=maxval(ndim_kk_fi)
c$$$      
c$$$      allocate(h_k_wan_val(num_wann,ndim_kk_fi_max))
c$$$      h_k_wan_val=0.0d0
c$$$      allocate(h_k_wan_vec(num_wann,num_wann,ndim_kk_fi_max))
c$$$      h_k_wan_vec=0.0d0
c$$$        
c$$$      call wannier_interpolation_grid(eigenvalues,ndiv_fi, nqdiv_fi,
c$$$     $  h_k_wan_val,h_k_wan_vec,ndim_kk_fi,n_mpi_kk_fi,ndim_kk_fi_max)
c$$$        
c$$$      allocate(e_wan_all(num_wann))
c$$$      e_wan_all=0.0d0
c$$$      allocate(hmat(maxel_loc,maxel_loc,natom_loc))
c$$$      hmat=0.0d0
c$$$      allocate(rhomat(maxel_loc,maxel_loc,natom_loc))
c$$$      rhomat=0.0d0
c$$$      allocate(proj(maxel_loc,natom_loc,num_wann,ndim_kk_fi_max))
c$$$      proj=0.0d0
c$$$      
c$$$      call cal_e_wan_all(eigenvalues, e_wan_all)
c$$$      
c$$$      if (maswrk) then
c$$$c$$$        if (nspin .eq. 1) then
c$$$c$$$          open(unit=889, file='e_wan_all.dat')
c$$$c$$$        else
c$$$c$$$          if (ispin .eq. 1) then
c$$$c$$$            open(unit=889, file='e_wan_all_up.dat')
c$$$c$$$          else
c$$$c$$$            open(unit=889, file='e_wan_all_dn.dat')
c$$$c$$$          endif
c$$$c$$$        endif
c$$$        do iwan=1, num_wann
c$$$          write(993,'(a, i5, 2f12.6)')
c$$$     $      'e_wan_all', iwan, e_wan_all(iwan)
c$$$        enddo
c$$$c$$$        close(889)
c$$$      endif
c$$$        
c$$$      proj=0.0d0
c$$$      
c$$$      do iatom_c=1, natom_loc
c$$$        iatom=iat_loc(iatom_c)        
c$$$        isort=is(iatom)
c$$$        atom_count=0
c$$$        do jj=1, nfun_loc(isort)
c$$$          do ik=1, ndim_kk_fi(me1)
c$$$            irk=n_mpi_kk_fi(me1)+ik
c$$$            do iband=1, num_wann
c$$$              proj(jj,iatom_c,iband,ik)=
c$$$     $          proj(jj,iatom_c,iband,ik)
c$$$     $          +h_k_wan_vec(wan_correlated(jj,iatom_c),iband,ik)
c$$$            enddo
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$      
c$$$!     projected energy at r=0
c$$$
c$$$      hmat=0.0d0
c$$$      rhomat=0.0d0
c$$$      
c$$$      do ik=1, ndim_kk_fi(me1)
c$$$        irk=n_mpi_kk_fi(me1)+ik
c$$$        do iatom_c=1, natom_loc
c$$$          iatom=iat_loc(iatom_c)        
c$$$          isort=is(iatom)
c$$$          do iwan=1, nfun_loc(isort)
c$$$            do jwan=1, nfun_loc(isort)
c$$$              do iband=1, num_wann
c$$$c$$$  
c$$$                hmat(iwan,jwan,iatom_c)=hmat(iwan,jwan,iatom_c)+
c$$$     $            h_k_wan_val(iband,ik)*proj(iwan,iatom_c,iband,ik)
c$$$     $            *dconjg(proj(jwan,iatom_c,iband,ik))/dble(nqdiv_fi)
c$$$                
c$$$                if (h_k_wan_val(iband,ik) .lt. 0.0d0) then
c$$$                  rhomat(iwan,jwan,iatom_c)=rhomat(iwan,jwan,iatom_c)+
c$$$     $              proj(iwan,iatom_c,iband,ik)
c$$$     $              *dconjg(proj(jwan,iatom_c,iband,ik))
c$$$     $              /dble(nqdiv_fi)
c$$$                endif
c$$$              enddo
c$$$            enddo
c$$$          enddo
c$$$        enddo
c$$$      enddo        
c$$$      
c$$$      call mpi_allreduce_dcmplx(rhomat,maxel_loc**2*natom_loc,
c$$$     $  mpi_sum,mpi_comm_world)
c$$$      call mpi_allreduce_dcmplx(hmat,maxel_loc**2*natom_loc,
c$$$     $  mpi_sum,mpi_comm_world)
c$$$      
c$$$      do iatom_c=1,natom_loc
c$$$        iatom=iat_loc(iatom_c)        
c$$$        isort=is(iatom)
c$$$        if (maswrk) then
c$$$          write(iun,'(a)') 'rhomat'
c$$$          do ii=1, nfun_loc(isort)
c$$$            write(iun,'(100(2f12.6, 5x))')
c$$$     $        (rhomat(ii,jj,iatom_c), jj=1,nfun_loc(isort))
c$$$          enddo
c$$$          write(iun,'(a)') 'hmat'
c$$$          do ii=1, nfun_loc(isort)
c$$$            write(iun,'(5(2f12.6, 5x))')
c$$$     $        (hmat(ii,jj,iatom_c), jj=1,5)
c$$$          enddo
c$$$        endif
c$$$      enddo
c$$$      
c$$$
c$$$c$$$      do ik=1, ndim_kk_fi(me1)
c$$$c$$$        irk=n_mpi_kk_fi(me1)+ik
c$$$c$$$        if (irk .lt. 10) then
c$$$c$$$          write(ch_kvec,'(i1)') irk
c$$$c$$$        elseif (irk .lt. 100) then
c$$$c$$$          write(ch_kvec,'(i2)') irk
c$$$c$$$        elseif (irk .lt. 1000) then
c$$$c$$$          write(ch_kvec,'(i3)') irk
c$$$c$$$        elseif (irk .lt. 10000) then
c$$$c$$$          write(ch_kvec,'(i4)') irk
c$$$c$$$        elseif (irk .lt. 100000) then
c$$$c$$$          write(ch_kvec,'(i5)') irk
c$$$c$$$        elseif (irk .lt. 1000000) then
c$$$c$$$          write(ch_kvec,'(i6)') irk
c$$$c$$$        elseif (irk .lt. 10000000) then
c$$$c$$$          write(ch_kvec,'(i7)') irk
c$$$c$$$        elseif (irk .lt. 100000000) then
c$$$c$$$          write(ch_kvec,'(i8)') irk
c$$$c$$$        endif
c$$$c$$$        if (nspin .eq. 1) then
c$$$c$$$          open(unit=998, file='projector_'//trim(ch_kvec)
c$$$c$$$     $      //'.dat')
c$$$c$$$        else
c$$$c$$$          if (ispin .eq. 1) then
c$$$c$$$            open(unit=998, file='projector_'//trim(ch_kvec)
c$$$c$$$     $        //'_up.dat')
c$$$c$$$          else
c$$$c$$$            open(unit=998, file='projector_'//trim(ch_kvec)
c$$$c$$$     $        //'_dn.dat')
c$$$c$$$          endif
c$$$c$$$        endif
c$$$c$$$        do ii=1, num_wann
c$$$c$$$          k1=(irk-1)/(ndiv_fi(3)*ndiv_fi(2))+1-1
c$$$c$$$          k2=(irk-(k1)*ndiv_fi(3)*ndiv_fi(2)-1)/ndiv_fi(3)+1-1
c$$$c$$$          k3=irk-(k1)*ndiv_fi(2)*ndiv_fi(3)-(k2)*ndiv_fi(3)-1
c$$$c$$$          write(998, '(3i5, f20.16, 3x, 500(2f20.16,1x))')
c$$$c$$$     $      k1,k2,k3,
c$$$c$$$     $      h_k_wan_val(ii,ik),
c$$$c$$$     $      ((proj(jj, kk, ii,ik),
c$$$c$$$     $      jj=1, maxel_loc), kk=1, natom_loc)
c$$$c$$$        enddo
c$$$c$$$        close(998)
c$$$c$$$      enddo
c$$$
c$$$      deallocate(h_k_wan_val)
c$$$      deallocate(h_k_wan_vec)
c$$$      deallocate(e_wan_all)
c$$$      deallocate(hmat)      
c$$$      deallocate(rhomat)      
c$$$      deallocate(proj)      
c$$$      end
c$$$
c$$$
      subroutine za_wan_for_u(ispin)

      use parallel_mod
      use manager_mod
      use wannier90_mod
      use solid_mod
      use atom_mod
      implicit none
      include "mpif.h"
      
      integer, intent(in) :: ispin
      integer :: ind_k, indtmp,iatom_c,k,
     $  iatom,isort,ii,ibnd,jbnd,k0
      complex*16 :: 
     &  znew(nfun,num_bands,nqdiv),anew(nbasmpw,num_bands,nqdiv)
      
      znew=0.0d0
      anew=0.0d0      

      do ind_k=1,ndim_kk(me+1)  ! k vector
        k=n_mpi_kk(me+1)+ind_k
        k0=i_kref(k)
        call sym_z_0(znew(1,1,k),k,z_wan_bnd(1,1,k0),
     &    num_bands,k_group(k),pnt(1,k))
        call sym_a_2(anew(1,1,k),k,ev_wan_bnd(1,1,k0),
     &    num_bands,k_group(k),k0)
      enddo

      if(nproc/=1) then
        call mpi_allreduce_dcmplx
     $    (znew,nfun*num_bands*nqdiv,
     $    mpi_sum,mpi_comm_world)
        call mpi_allreduce_dcmplx
     $    (anew,nbasmpw*num_bands*nqdiv,
     $    mpi_sum,mpi_comm_world)
      endif
      if (maswrk) then
        open(unit=334, file=trim(seed_name)//'_z_wan_band.dat',
     $    form='unformatted')
        write(334) znew
        close(334)
        open(unit=334, file=trim(seed_name)//'_a_wan_band.dat',
     $    form='unformatted')
        write(334) anew
        close(334)
      endif
      end
      
