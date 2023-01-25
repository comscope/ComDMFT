      subroutine boson_tensor(uloc,tempchar,fac)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use wannier90_mod
      implicit none
      include "mpif.h"
      complex*16, intent(in) ::
     $  uloc(num_orb_proj,num_orb_proj,num_orb_proj,num_orb_proj,
     $  nspin,nspin)
      double precision , intent(in) :: fac
      character(len=*),intent(in) :: tempchar      

      integer :: isort,l,ij,kl,i_len,iatom_c,iatom,nf,ispin,jspin,lm1,
     &  lm2,lm3,lm4,km1,km2,km3,km4,ll,ig,lm_min,ll0,llm,ii,
     &  iii,indx,icg1,icg2,jcg1,jcg2,icg,jjcg,lmi,lmj,ind_nu,
     &  i_nu,li,mudim,lf,ie,m4,m3,m2,m1,iat,
     $  i1,i2,i3,i4,kk,kk0,i1mi2,i3mi4,iwan,jwan,kwan,lwan,
     $  ind1,ind2,ind3,ind4
      
      complex*16,allocatable :: 
     $  umatrix_temp(:,:,:,:)
      character(len=5) :: catom
      
      i_len=len_trim(allfile)
      
      if (maswrk) then
        
        do iat=1,natom_cw
          iatom=iatom_cw(iat)
          if (iatom .lt. 10)then
            write(catom,'(i1)') iatom
          elseif (iatom .lt. 10)then
            write(catom,'(i2)') iatom
          elseif (iatom .lt. 100)then
            write(catom,'(i3)') iatom
          endif
          
          allocate(umatrix_temp(n_cw_atom(iatom), n_cw_atom(iatom),
     $      n_cw_atom(iatom),n_cw_atom(iatom)))
          
          open(3,file=allfile(1:i_len)
     $      //'_'//trim(adjustl(tempchar))//'_tensor_'
     $      //trim(adjustl(catom))
     $      //'.dat')
          write(3,'(a, 10000i4)')
     $      '# tensor, wannier index: ',
     $      (wan_correlated(corwan_ind_atom(iwan,iat)),
     $      iwan=1, n_cw_atom(iatom))
          do jspin=1,nspin
            do ispin=1,nspin
              if(maswrk) then
                write(3,'(a15,i1,a1,i1,a11,i2,2(a6,i2))')
     &            '# ispin/jspin = ',ispin,'/',jspin,' iatom = ',
     &            iatom
              endif
              umatrix_temp=0.0d0
              do iwan=1, n_cw_atom(iatom)
                ind1=corwan_ind_atom(iwan,iat)
                do jwan=1, n_cw_atom(iatom)
                  ind2=corwan_ind_atom(jwan,iat)
                  do kwan=1, n_cw_atom(iatom)
                    ind3=corwan_ind_atom(iwan,iat)
                    do lwan=1, n_cw_atom(iatom)
                      ind4=corwan_ind_atom(jwan,iat)
                      
                      umatrix_temp(iwan,jwan,kwan,lwan)
     $                  =uloc(ind1,ind2,ind3,ind4,ispin,jspin)
                    enddo
                  enddo
                enddo
              enddo
              write(3,'(10000000(5x,(2(e22.15,3x))))')
     &          ((((umatrix_temp(iwan,jwan,kwan,lwan)*fac,
     $          iwan=1,n_cw_atom(iatom)),
     &          jwan=1,n_cw_atom(iatom)),
     $          kwan=1,n_cw_atom(iatom)),
     &          lwan=1,n_cw_atom(iatom))
            enddo
          enddo
          deallocate(umatrix_temp)
          if(maswrk) close(3)
        enddo
      endif
      end

      

      subroutine boson_density(uloc,tempchar,fac)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use wannier90_mod
      implicit none
      include "mpif.h"
      complex*16, intent(in) ::
     $  uloc(num_orb_proj,num_orb_proj,num_orb_proj,num_orb_proj,
     $  nspin,nspin,ndim3_nu)
      double precision , intent(in) :: fac
      character(len=*),intent(in) :: tempchar      

      integer :: isort,l,ij,kl,i_len,iatom_c,iatom,nf,ispin,jspin,lm1,
     &  lm2,lm3,lm4,km1,km2,km3,km4,ll,ig,lm_min,ll0,llm,ii,
     &  iii,indx,icg1,icg2,jcg1,jcg2,icg,jjcg,lmi,lmj,ind_nu,
     &  i_nu,li,mudim,lf,ie,m4,m3,m2,m1,iat,
     $  i1,i2,i3,i4,kk,kk0,i1mi2,i3mi4,iwan,jwan,
     $  ind1,ind2
      
      complex*16,allocatable :: 
     $  umatrix_temp(:,:,:)
      character(len=5) :: catom
      

      i_len=len_trim(allfile)
      
      do iat=1,natom_cw
        iatom=iatom_cw(iat)
        if (iatom .lt. 10)then
          write(catom,'(i1)') iatom
        elseif (iatom .lt. 10)then
          write(catom,'(i2)') iatom
        elseif (iatom .lt. 100)then
          write(catom,'(i3)') iatom
        endif
        
        allocate(umatrix_temp(n_cw_atom(iatom), n_cw_atom(iatom),
     $    0:n_nu))
        
        if(maswrk) then
          open(3,file=allfile(1:i_len)
     $      //'_'//trim(adjustl(tempchar))//'_density_'
     $      //trim(adjustl(catom))
     $      //'.dat')
          write(3,'(a, 10000i4)')
     $      '# density density term, wannier index: ',
     $      (wan_correlated(corwan_ind_atom(iwan,iat)),
     $      iwan=1, n_cw_atom(iatom))
        endif
        
        do jspin=1,nspin
          do ispin=1,nspin
            if(maswrk) then
              write(3,'(a15,i1,a1,i1,a11,i2,2(a6,i2))')
     &          '# ispin/jspin = ',ispin,'/',jspin,' iatom = ',
     &          iatom
            endif
            umatrix_temp=0.0d0
            do ind_nu=1,ndim3_nu
              i_nu=me3_tau*ndim3_nu+ind_nu-1                              
              
              do iwan=1, n_cw_atom(iatom)
                ind1=corwan_ind_atom(iwan,iat)
                do jwan=1, n_cw_atom(iatom)
                  ind2=corwan_ind_atom(jwan,iat)
                  umatrix_temp(iwan,jwan,i_nu)
     $              =uloc(ind1,ind1,ind2,ind2,ispin,jspin,ind_nu)
                enddo
              enddo
            enddo
            if(nproc_tau/=1) then
              call mpi_allreduce_dcmplx
     $          (umatrix_temp,n_cw_atom(iatom)**2*(n_nu+1),
     $          mpi_sum,comm_k)
            endif
            if(maswrk) then
              do i_nu=0,n_nu
                write(3,'(100000(5x,(2(e22.15,2x))))')
     $            w_nu(i_nu)*evolt/2,
     &            ((umatrix_temp(iwan,jwan,i_nu)*fac,
     $            iwan=1,n_cw_atom(iatom)),
     &            jwan=1,n_cw_atom(iatom))
              enddo
            endif
          enddo
        enddo
        deallocate(umatrix_temp)
        if(maswrk) close(3)
      enddo       
      end

      

      subroutine boson_slater_l(uloc,tempchar,fac)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use wannier90_mod
      implicit none
      include "mpif.h"
      complex*16, intent(in) ::
     $  uloc(num_orb_proj,num_orb_proj,num_orb_proj,num_orb_proj,
     $  nspin,nspin,ndim3_nu)
      double precision , intent(in) :: fac
      character(len=*),intent(in) :: tempchar      

      integer :: isort,l,ij,kl,i_len,iatom_c,iatom,nf,ispin,jspin,lm1,
     &  lm2,lm3,lm4,km1,km2,km3,km4,ll,ig,lm_min,ll0,llm,ii,
     &  iii,indx,icg1,icg2,jcg1,jcg2,icg,jjcg,lmi,lmj,ind_nu,
     &  i_nu,li,mudim,lf,ie,m4,m3,m2,m1,iat,
     $  i1,i2,i3,i4,kk,kk0,i1mi2,i3mi4
      
      double precision :: c,cc
      complex*16,allocatable :: ff(:,:),rotmat_cmplx2real(:,:),
     $  umatrix_temp(:,:,:,:)
      character(len=5) :: lchar,catom
      

      allocate(ff(0:3,0:n_nu))
      i_len=len_trim(allfile)
      

      do iat=1,natom_cw
        iatom=iatom_cw(iat)
        isort=is(iatom)
        if (iatom .lt. 10)then
          write(catom,'(i1)') iatom
        elseif (iatom .lt. 10)then
          write(catom,'(i2)') iatom
        elseif (iatom .lt. 100)then
          write(catom,'(i3)') iatom
        endif
        
        do li=1, nrel*lmb(isort)+1
          l=li-1
          mudim=l+l+1
          write(lchar,'(i2)') l
          if ((corr(l,iatom) .eq. 2)
     $      .and.
     $      (at_bonding(l,iatom) .eq. 0)
     $      ) then
            if(maswrk) then
              open(3,file=allfile(1:i_len)
     $          //'_'//trim(adjustl(tempchar))//'_Slater_'
     $          //trim(adjustl(catom))
     $          //'_'//trim(adjustl(lchar))//'.dat')
              write(3,'(a)')
     $          '# Slater integral representation for U_matrix'
            endif

            allocate(rotmat_cmplx2real(mudim,mudim))
            rotmat_cmplx2real=0.0d0
            call cal_rotmat_cmplx2real(mudim, rotmat_cmplx2real)
            allocate(umatrix_temp(mudim,mudim,mudim,mudim))
            
            do jspin=1,nspin
              do ispin=1,nspin

                if(maswrk) then
                  write(3,'(a15,i1,a1,i1,a11,i2,2(a6,i2))')
     &              '# ispin/jspin = ',ispin,'/',jspin,' iatom = ',
     &              iatom,' L = ',l
                endif
                ff=0.d0                
                do ind_nu=1,ndim3_nu
                  i_nu=me3_tau*ndim3_nu+ind_nu-1                
                  umatrix_temp=0.0d0
                  
                  do m4=1,mudim
                    km4=mt2ciw(m4,li,iatom)
                    
                    do m3=1,mudim
                      km3=mt2ciw(m3,li,iatom)
                      
                      do m2=1,mudim
                        km2=mt2ciw(m2,li,iatom)
                        
                        do m1=1,mudim
                          km1=mt2ciw(m1,li,iatom)
                          umatrix_temp(m1,m2,m3,m4)
     $                      =umatrix_temp(m1,m2,m3,m4)
     $                      +uloc(km1,km2,km3,km4,
     $                      ispin,jspin,ind_nu)
                        enddo
                      enddo
                    enddo
                  enddo
                  call rotate_umatrix(mudim,umatrix_temp,
     $              rotmat_cmplx2real,1)
                  call cal_slater_integral(umatrix_temp,l,ff(0,i_nu))
                enddo
                
                if(nproc_tau/=1) then
                  call mpi_allreduce_dcmplx
     $              (ff,4*(n_nu+1),mpi_sum,comm_k)
                endif
                if(maswrk) then
                  do i_nu=0,n_nu
                    write(3,'(10(1x,e22.15))')w_nu(i_nu)*evolt/2,
     &                (dble(ff(ll,i_nu))*fac,ll=0,l)
                  enddo
                endif
              enddo
            enddo
            deallocate(rotmat_cmplx2real)
            deallocate(umatrix_temp)
            if(maswrk) close(3)
          endif
        enddo     
      enddo       
      deallocate(ff)
      end


      subroutine v_slater_l
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use wannier90_mod
      implicit none
      include "mpif.h"
      
      integer :: isort,l,ij,kl,i_len,iatom_c,iatom,nf,ispin,jspin,lm1,
     &  lm2,lm3,lm4,km1,km2,km3,km4,ll,ig,lm_min,ll0,llm,ii,
     &  iii,indx,icg1,icg2,jcg1,jcg2,icg,jjcg,lmi,lmj,
     &  li,mudim,lf,ie,m4,m3,m2,m1,iat,
     $  i1,i2,i3,i4,kk,kk0,i1mi2,i3mi4
      
      double precision :: c,cc
      complex*16,allocatable :: ff(:),rotmat_cmplx2real(:,:),
     $  umatrix_temp(:,:,:,:)
      character(len=5) :: lchar,catom
      
      
      allocate(ff(0:3))
      i_len=len_trim(allfile)
      
      
      do iat=1,natom_cw
        iatom=iatom_cw(iat)
        isort=is(iatom)
        if (iatom .lt. 10)then
          write(catom,'(i1)') iatom
        elseif (iatom .lt. 10)then
          write(catom,'(i2)') iatom
        elseif (iatom .lt. 100)then
          write(catom,'(i3)') iatom
        endif
        
        do li=1, nrel*lmb(isort)+1
          l=li-1
          mudim=l+l+1
          write(lchar,'(i2)') l
          if ((corr(l,iatom) .eq. 2)
     $      .and.
     $      (at_bonding(l,iatom) .eq. 0)
     $      ) then
            if(maswrk) then
              open(3,file=allfile(1:i_len)
     $          //'_v_Slater_'
     $          //trim(adjustl(catom))
     $          //'_'//trim(adjustl(lchar))//'.dat')
              write(3,'(a)')
     $          '# Slater integral representation for v_matrix'
            endif
            
            allocate(rotmat_cmplx2real(mudim,mudim))
            rotmat_cmplx2real=0.0d0
            call cal_rotmat_cmplx2real(mudim, rotmat_cmplx2real)
            allocate(umatrix_temp(mudim,mudim,mudim,mudim))
            
            do jspin=1,nspin
              do ispin=1,nspin
                
                if(maswrk) then
                  write(3,'(a15,i1,a1,i1,a11,i2,2(a6,i2))')
     &              '# ispin/jspin = ',ispin,'/',jspin,' iatom = ',
     &              iatom,' L = ',l
                endif
                ff=0.d0                
                umatrix_temp=0.0d0
                
                do m4=1,mudim
                  km4=mt2ciw(m4,li,iatom)
                  
                  do m3=1,mudim
                    km3=mt2ciw(m3,li,iatom)
                    
                    do m2=1,mudim
                      km2=mt2ciw(m2,li,iatom)
                      
                      do m1=1,mudim
                        km1=mt2ciw(m1,li,iatom)
                        umatrix_temp(m1,m2,m3,m4)
     $                    =umatrix_temp(m1,m2,m3,m4)
     $                    +v_loc(km1,km2,km3,km4,
     $                    ispin,jspin)
                      enddo
                    enddo
                  enddo
                enddo
                call rotate_umatrix(mudim,umatrix_temp,
     $            rotmat_cmplx2real,1)
                call cal_slater_integral(umatrix_temp,l,ff(0))
                
                
                if(maswrk) then
                  write(3,'(10(1x,e22.15))')
     &              (dble(ff(ll))*rydberg,ll=0,l)
                endif
              enddo
            enddo            
            deallocate(rotmat_cmplx2real)
            deallocate(umatrix_temp)
            if(maswrk) close(3)
          endif
        enddo     
      enddo       
      deallocate(ff)
      end


c$$$
c$$$      subroutine v_tensor_l
c$$$      use atom_mod
c$$$      use manager_mod
c$$$      use parallel_mod
c$$$      use solid_mod
c$$$      use units_mod
c$$$      use wannier90_mod
c$$$      implicit none
c$$$      include "mpif.h"
c$$$      
c$$$      integer :: isort,l,ij,kl,i_len,iatom_c,iatom,nf,ispin,jspin,lm1,
c$$$     &  lm2,lm3,lm4,km1,km2,km3,km4,ll,ig,lm_min,ll0,llm,ii,
c$$$     &  iii,indx,icg1,icg2,jcg1,jcg2,icg,jjcg,lmi,lmj,
c$$$     &  li,mudim,lf,ie,m4,m3,m2,m1,iat,
c$$$     $  i1,i2,i3,i4,kk,kk0,i1mi2,i3mi4
c$$$      
c$$$      double precision :: c,cc
c$$$      character(len=5) :: lchar,catom
c$$$      
c$$$      
c$$$      i_len=len_trim(allfile)
c$$$      
c$$$      
c$$$      do iat=1,natom_ciw
c$$$        iatom=iatom_ciw(iat)
c$$$        isort=is(iatom)
c$$$        if (iatom .lt. 10)then
c$$$          write(catom,'(i1)') iatom
c$$$        elseif (iatom .lt. 10)then
c$$$          write(catom,'(i2)') iatom
c$$$        elseif (iatom .lt. 100)then
c$$$          write(catom,'(i3)') iatom
c$$$        endif
c$$$        
c$$$        do li=1, nrel*lmb(isort)+1
c$$$          l=li-1
c$$$          mudim=l+l+1
c$$$          write(lchar,'(i2)') l
c$$$          if ((corr_iter(l,iatom) .eq. 2)
c$$$     $      .and.
c$$$     $      (at_bonding(l,iatom) .eq. 0)
c$$$     $      ) then
c$$$            if(maswrk) then
c$$$              open(3,file=allfile(1:i_len)
c$$$     $          //'_v_tensor_'
c$$$     $          //trim(adjustl(catom))
c$$$     $          //'_'//trim(adjustl(lchar))//'.dat')
c$$$              write(3,'(a)')
c$$$     $          '# v_tensor'
c$$$            endif
c$$$            
c$$$            do jspin=1,nspin
c$$$              do ispin=1,nspin
c$$$                
c$$$                do m4=1,mudim
c$$$                  km4=mt2ciw(m4,li,iatom)
c$$$                  
c$$$                  do m3=1,mudim
c$$$                    km3=mt2ciw(m3,li,iatom)
c$$$                    
c$$$                    do m2=1,mudim
c$$$                      km2=mt2ciw(m2,li,iatom)
c$$$                      
c$$$                      do m1=1,mudim
c$$$                        km1=mt2ciw(m1,li,iatom)
c$$$                        if(maswrk) then                        
c$$$                          write(3,'(6i5,2x,2(e22.15,2x))')
c$$$     $                      m1,m2,m3,m4,ispin,jspin,
c$$$     $                      v_loc(km1,km2,km3,km4,
c$$$     $                      ispin,jspin)*rydberg
c$$$                        endif
c$$$                      enddo
c$$$                    enddo
c$$$                  enddo
c$$$                enddo
c$$$              enddo     
c$$$            enddo
c$$$            if(maswrk) close(3)            
c$$$          endif
c$$$        enddo
c$$$      enddo
c$$$      end
c$$$


      integer function find_wan(jj)

      use wannier90_mod
      use manager_mod
      use parallel_mod
      implicit none
      integer, intent(in) :: jj
      integer :: ii
      
      find_wan=0
      do ii=1, num_orb_cor
        if (jj .eq. wan_correlated(ii)) then
          find_wan=ii
        endif
      enddo

      do ii=1, num_orb_iter
        if (jj .eq. wan_itinerant(ii)) then
          find_wan=ii+num_orb_cor
        endif
      enddo      
      
      if (find_wan .eq. 0) then
        call all_log('problem in find_wan')
        call ending_qft
      endif
      end      


