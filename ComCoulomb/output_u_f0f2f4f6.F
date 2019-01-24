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
      

      do iat=1,natom_ciw
        iatom=iatom_ciw(iat)
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
          if ((corr_iter(l,iatom) .eq. 2)
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
      
      
      do iat=1,natom_ciw
        iatom=iatom_ciw(iat)
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
          if ((corr_iter(l,iatom) .eq. 2)
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


