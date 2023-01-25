      subroutine v_slater_rel_l      
!     ff
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use wannier90_mod      
      implicit none
      include "mpif.h"
      
      integer :: isort,l,ij,kl,i_len,iatom_c,iatom,nf,ispin,jspin,lm1,
     &  lm2,lm3,lm4,km1,km2,km3,km4,ll,ig,lm_min,ll0,llm,ii,jj,
     &  iii,indx,icg1,icg2,jcg1,jcg2,icg,jjcg,lmi,lmj,
     &  li,ndim,lf,ie,m4,m3,m2,m1,n_half,n,m12,m34,tmpind,
     $  k1,k2,k3,k4,i1,i2,i3,i4,iistart,ind1,ind2,ind3,ind4,
     $  li1,li2,li3,li4,mu1,mu2,mu3,mu4,mudim1,mudim2,mudim3,mudim4,iat
      double precision :: c,cc

      complex*16, allocatable :: ff(:),rotmat_mlms2jmj(:,:),
     $  uloc_temp(:,:,:,:),
     $  uloc_temp2(:,:,:,:)

      character(len=5) :: lchar,catom


      if (nspin .ne. 1) then
        call all_log('nspin should be 1 to use boson_slater_rel_1')
        call ending
      endif

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
        
        do l=0, lmb(isort)
          
          if ((corr(l,iatom) .eq. 2)
     $      .and.
     $      (at_bonding(l,iatom) .eq. 0)
     $      ) then
            
            write(lchar,'(i1)') l            
            
            if(maswrk) then
              open(3,file=allfile(1:i_len)
     $          //'_v_Slater_'
     $          //trim(adjustl(catom))
     $          //'_'//trim(adjustl(lchar))//'.dat')
              write(3,'(a)')
     $          '# Slater integral representation for v_matrix'
              write(3,'(a15,i1,a1,i1,a11,i2,2(a6,i2))')
     &          '# ispin/jspin = ',1,'/',1,' iatom = ',
     &          iatom,' L = ',l
            endif            
            
            n=4*l+2
            n_half=n/2
            
            allocate(rotmat_mlms2jmj(n,n))
            rotmat_mlms2jmj=0.0d0
            
            call cal_rotmat_mlms2jmj(l, rotmat_mlms2jmj)
            
            allocate(uloc_temp(n,n,n,n))
            uloc_temp=0.0d0
            allocate(uloc_temp2(n_half,n_half,n_half,n_half))
            uloc_temp2=0.0d0
            
            ff=0.0d0
            
c$$$  uloc in real spherical harmonics
            
            uloc_temp=0.0d0
            uloc_temp2=0.0d0
            
            ind4=0
            do m4=1, n
              km4=mt2ciw_rel_l(m4,l+1,iatom)
              do m3=1, n
                km3=mt2ciw_rel_l(m3,l+1,iatom)                  
                do m2=1, n
                  km2=mt2ciw_rel_l(m2,l+1,iatom)                    
                  do m1=1, n
                    km1=mt2ciw_rel_l(m1,l+1,iatom)
                    uloc_temp(m1,m2,m3,m4)
     $                =v_loc(km1,km2,km3,km4,
     $                1,1)
                  enddo
                enddo
              enddo
            enddo
            call rotate_umatrix(n,uloc_temp,
     $        rotmat_mlms2jmj,1)
            
            uloc_temp2
     $        =(
     $        uloc_temp(1:n_half,1:n_half,1:n_half,1:n_half)
     $        +uloc_temp((n_half+1):n,(n_half+1):n,1:n_half,
     $        1:n_half)
     $        +uloc_temp(1:n_half,1:n_half,(n_half+1):n,
     $        (n_half+1):n)              
     $        +uloc_temp((n_half+1):n,(n_half+1):n,
     $        (n_half+1):n,(n_half+1):n)
     $        )/4.0d0
            
            call cal_slater_integral(uloc_temp2,l,ff(0))
            
            
            deallocate(rotmat_mlms2jmj)
            deallocate(uloc_temp)
            deallocate(uloc_temp2)
            
            if(maswrk) then
              write(3,'(10(1x,e22.15))')
     &          (dble(ff(ll))*rydberg,ll=0,l)
              close(3)
            endif
          endif
        enddo
      enddo  
      deallocate(ff)      
      end


c$$$      subroutine v_tensor_rel_l      
c$$$!     ff
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
c$$$     &  lm2,lm3,lm4,km1,km2,km3,km4,ll,ig,lm_min,ll0,llm,ii,jj,
c$$$     &  iii,indx,icg1,icg2,jcg1,jcg2,icg,jjcg,lmi,lmj,
c$$$     &  li,ndim,lf,ie,m4,m3,m2,m1,n_half,n,m12,m34,tmpind,
c$$$     $  k1,k2,k3,k4,i1,i2,i3,i4,iistart,ind1,ind2,ind3,ind4,
c$$$     $  li1,li2,li3,li4,mu1,mu2,mu3,mu4,mudim1,mudim2,mudim3,mudim4,iat
c$$$      double precision :: c,cc
c$$$
c$$$
c$$$      character(len=5) :: lchar,catom
c$$$
c$$$
c$$$      if (nspin .ne. 1) then
c$$$        call all_log('nspin should be 1 to use boson_slater_rel_1')
c$$$        call ending
c$$$      endif
c$$$
c$$$      i_len=len_trim(allfile)
c$$$
c$$$
c$$$      do iat=1,natom_ciw
c$$$        iatom=iatom_ciw(iat)
c$$$        isort=is(iatom)
c$$$
c$$$        if (iatom .lt. 10)then
c$$$          write(catom,'(i1)') iatom
c$$$        elseif (iatom .lt. 10)then
c$$$          write(catom,'(i2)') iatom
c$$$        elseif (iatom .lt. 100)then
c$$$          write(catom,'(i3)') iatom
c$$$        endif
c$$$        
c$$$        do l=0, lmb(isort)
c$$$          
c$$$          if ((corr_iter(l,iatom) .eq. 2)
c$$$     $      .and.
c$$$     $      (at_bonding(l,iatom) .eq. 0)
c$$$     $      ) then
c$$$            
c$$$            write(lchar,'(i1)') l            
c$$$            
c$$$            if(maswrk) then
c$$$              open(3,file=allfile(1:i_len)
c$$$     $          //'_v_tensor_'
c$$$     $          //trim(adjustl(catom))
c$$$     $          //'_'//trim(adjustl(lchar))//'.dat')
c$$$              write(3,'(a)')
c$$$     $          '# v_tensor'
c$$$            endif            
c$$$            
c$$$            ind4=0
c$$$
c$$$            n=4*l+2
c$$$            
c$$$            do m4=1, n
c$$$              km4=mt2ciw_rel_l(m4,l+1,iatom)
c$$$              do m3=1, n
c$$$                km3=mt2ciw_rel_l(m3,l+1,iatom)                  
c$$$                do m2=1, n
c$$$                  km2=mt2ciw_rel_l(m2,l+1,iatom)                    
c$$$                  do m1=1, n
c$$$                    if(maswrk) then                        
c$$$                      write(3,'(4i5,2x,2(e22.15, 2x))')                    
c$$$     $                  m1,m2,m3,m4,
c$$$     $                  v_loc(km1,km2,km3,km4,
c$$$     $                  1,1)
c$$$                    endif
c$$$                  enddo
c$$$                enddo
c$$$              enddo
c$$$            enddo
c$$$            if(maswrk) close(3)                        
c$$$          endif
c$$$        enddo
c$$$      enddo  
c$$$      end      

      subroutine boson_slater_rel_l(uloc,tempchar,fac)
      
!     ff
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
      double precision, intent(in) :: fac
      character(len=*),intent(in) :: tempchar
      
      integer :: isort,l,ij,kl,i_len,iatom_c,iatom,nf,ispin,jspin,lm1,
     &  lm2,lm3,lm4,km1,km2,km3,km4,ll,ig,lm_min,ll0,llm,ii,jj,
     &  iii,indx,icg1,icg2,jcg1,jcg2,icg,jjcg,lmi,lmj,ind_nu,
     &  i_nu,li,ndim,lf,ie,m4,m3,m2,m1,n_half,n,m12,m34,tmpind,
     $  k1,k2,k3,k4,i1,i2,i3,i4,iistart,ind1,ind2,ind3,ind4,
     $  li1,li2,li3,li4,mu1,mu2,mu3,mu4,mudim1,mudim2,mudim3,mudim4,iat
      double precision :: c,cc

      complex*16, allocatable :: ff(:,:),rotmat_mlms2jmj(:,:),
     $  uloc_temp(:,:,:,:),
     $  uloc_temp2(:,:,:,:)

      character(len=5) :: lchar,catom


      if (nspin .ne. 1) then
        call all_log('nspin should be 1 to use boson_slater_rel_1')
        call ending
      endif

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
        
        do l=0, lmb(isort)
          
          if ((corr(l,iatom) .eq. 2)
     $      .and.
     $      (at_bonding(l,iatom) .eq. 0)
     $      ) then

            write(lchar,'(i1)') l            

            if(maswrk) then
              open(3,file=allfile(1:i_len)
     $          //'_'//trim(adjustl(tempchar))//'_Slater_'
     $          //trim(adjustl(catom))
     $          //'_'//trim(adjustl(lchar))//'.dat')
              write(3,'(a)')
     $          '# Slater integral representation for U_matrix'
              write(3,'(a15,i1,a1,i1,a11,i2,2(a6,i2))')
     &          '# ispin/jspin = ',1,'/',1,' iatom = ',
     &          iatom,' L = ',l
            endif            
            
            n=4*l+2
            n_half=n/2
            
            allocate(rotmat_mlms2jmj(n,n))
            rotmat_mlms2jmj=0.0d0
            
            call cal_rotmat_mlms2jmj(l, rotmat_mlms2jmj)
            
            allocate(uloc_temp(n,n,n,n))
            uloc_temp=0.0d0
            allocate(uloc_temp2(n_half,n_half,n_half,n_half))
            uloc_temp2=0.0d0
            
            ff=0.0d0
            
c$$$  uloc in real spherical harmonics
            
            do ind_nu=1,ndim3_nu
              i_nu=me3_tau*ndim3_nu+ind_nu-1              
              uloc_temp=0.0d0
              uloc_temp2=0.0d0
              
              ind4=0
              do m4=1, n
                km4=mt2ciw_rel_l(m4,l+1,iatom)
                do m3=1, n
                  km3=mt2ciw_rel_l(m3,l+1,iatom)                  
                  do m2=1, n
                    km2=mt2ciw_rel_l(m2,l+1,iatom)                    
                    do m1=1, n
                      km1=mt2ciw_rel_l(m1,l+1,iatom)
                      uloc_temp(m1,m2,m3,m4)
     $                  =uloc(km1,km2,km3,km4,
     $                  1,1,ind_nu)                      
                    enddo
                  enddo
                enddo
              enddo
              call rotate_umatrix(n,uloc_temp,
     $          rotmat_mlms2jmj,1)
              
              uloc_temp2
     $          =(
     $          uloc_temp(1:n_half,1:n_half,1:n_half,1:n_half)
     $          +uloc_temp((n_half+1):n,(n_half+1):n,1:n_half,
     $          1:n_half)
     $          +uloc_temp(1:n_half,1:n_half,(n_half+1):n,
     $          (n_half+1):n)              
     $          +uloc_temp((n_half+1):n,(n_half+1):n,
     $          (n_half+1):n,(n_half+1):n)
     $          )/4.0d0

              call cal_slater_integral(uloc_temp2,l,ff(0,i_nu))
            enddo

            
            deallocate(rotmat_mlms2jmj)
            deallocate(uloc_temp)
            deallocate(uloc_temp2)
            
            if(nproc_tau/=1) then
              call mpi_allreduce_dcmplx
     $          (ff,4*(n_nu+1),mpi_sum,comm_k)
            endif
            if(maswrk) then
              do i_nu=0,n_nu
                write(3,'(10(1x,e22.15))')w_nu(i_nu)*rydberg,
     &            (dble(ff(ll,i_nu))*fac,ll=0,l)
              enddo
            endif
            if(maswrk) close(3)                      
          endif
        enddo
      enddo  
      deallocate(ff)      
      end
      
C$$$  subroutine output_u_f0f2f4f6_rel_j(uloc,tempchar,is_inv)
C$$$  use atom_mod
C$$$  use manager_mod
C$$$  use models_mod
C$$$  use parallel_mod
C$$$  use solid_mod
C$$$  use units_mod
C$$$  use wannier90_mod      
C$$$  implicit none
C$$$  include "mpif.h"
C$$$  integer, intent(in) :: is_inv
C$$$  complex*16, intent(in) ::
C$$$  $  uloc(nwp_max,nwp_max,nwp_max,nwp_max,
C$$$  $  natom_ciw,natom_ciw,nspin_0,nspin_0,ndim3_nu)
C$$$  character*300,intent(in) :: tempchar      
C$$$  integer :: isort,l,ij,kl,i_len,iatom_c,iatom,nf,ispin,jspin,lm1,
C$$$  &  lm2,lm3,lm4,km1,km2,km3,km4,ll,lm_min,ll0,llm,ii,icg1,
C$$$  &  icg2,jcg1,jcg2,icg,jjcg,lmi,lmj,ind_nu,i_nu,li,jj,mj1,
C$$$  &  mj,lm,l22,lm22,lget,ndim,lf,ie,m4,m3,m2,m1,iat
C$$$  double precision :: cc,n_l(0:3),fac
C$$$  double precision, allocatable :: fct(:)
C$$$  complex*16, allocatable :: ff(:,:)
C$$$  integer, external :: find_ind

C$$$  fac=evolt/2.0d0
C$$$  if (is_inv .eq. 1) then
C$$$  fac=(1.0d0/(evolt/2.0d0))
C$$$  endif
      
C$$$  allocate(fct(0:20))
C$$$  allocate(ff(0:n_nu,0:3))
C$$$  call facc(fct,20)
C$$$  i_len=len_trim(allfile)
C$$$  if(maswrk) then
C$$$  open(3,file=allfile(1:i_len)//'_F0_F2_F4_F6'
C$$$  $    //trim(adjustl(tempchar))
C$$$  $    //'_'//trim(ubi)//'_j.dat')
C$$$  write(3,*)'# F0 F2 F4 F6 representation for U_matrix'
C$$$  endif
C$$$  do jspin=1,nspin
C$$$  do ispin=1,nspin
C$$$  do iatom=1, natom
C$$$  isort=is(iatom)
C$$$  do li=1,nrel*lmb(isort)+1
C$$$  l=li/2
C$$$  ndim=((li+1)/2)*2
C$$$  lm_min=nrel*l*l              
      
C$$$  if ((corr_iter_rel(li,iatom) .eq. 2)
C$$$  $          .and.
C$$$  $          (at_bonding_rel(li,iatom) .eq. 0)
C$$$  $          ) then
      
C$$$  if(maswrk) then
C$$$  write(3,'(a15,i1,a1,i1,a11,i2,2(a6,i2))')
C$$$  &	            '# ispin/jspin = ',ispin,'/',jspin,' iatom = ',
C$$$  &              iatom,' Li = ',li
C$$$  write(3,*)'# Frequency F0 F2 F4 F6: (in eV)'
C$$$  endif
C$$$  c     -------------- Get N_l ----------------------------------------
C$$$  ii=mod(li,2)
C$$$  if(ii.eq.0) ii=-1
C$$$  jj=l+l+ii                
C$$$  if (jj .eq. 1) then
C$$$  n_l(0)=2.0d0**2
C$$$  elseif (jj .eq. 3) then
C$$$  n_l(0)=2.0d0**4
C$$$  n_l(1)=2.0d0**4/5.0d0**3
C$$$  elseif (jj .eq. 5) then
C$$$  n_l(0)=2.0d0**2*3.0d0**2
C$$$  n_l(1)=2.0d0**8*3.0d0**2/5.0d0**3/7.0d0**2
C$$$  n_l(2)=2.0d0**4/3.0d0**2/7.0d0**2
C$$$  elseif (jj .eq. 7) then
C$$$  n_l(0)=2.0d0**6
C$$$  n_l(1)=2.0d0**6/3.0d0**2*5.0d0**1/7.0d0**2
C$$$  n_l(2)=2.0d0**6*3.0d0**2/7.0d0**2/11.0d0**2
C$$$  n_l(3)=2.0d0**6/3.0d0**2
C$$$  $              *5.0d0**4/11.0d0**2/13.0d0**3
C$$$  endif
C$$$  c     ---------------------------------------------------------------	      
C$$$  ff=0.d0
C$$$  do m4=1,ndim
C$$$  km4=find_ind(mt2wan(m4,li,iatom),
C$$$  $              n_ciw_atom(iatom),ciw_ind(1,iat))
C$$$  lm4=lm_min+m4                  
C$$$  do m3=1,ndim
C$$$  km3=find_ind(mt2wan(m3,li,iatom),
C$$$  $                n_ciw_atom(iatom),ciw_ind(1,iat))                    
C$$$  lm3=lm_min+m3
C$$$  do m2=1,ndim
C$$$  km2=find_ind(mt2wan(m2,li,iatom),
C$$$  $                  n_ciw_atom(iatom),ciw_ind(1,iat))
C$$$  lm2=lm_min+m2                      
C$$$  do m1=1,ndim
C$$$  km1=find_ind(mt2wan(m1,li,iatom),
C$$$  $                    n_ciw_atom(iatom),ciw_ind(1,iat))
C$$$  lm1=lm_min+m1
C$$$  do ll=0,2*l,2
C$$$  ll0=ll/2
C$$$  cc=4.d0*pi/n_l(ll0) ! /(ll+ll+1.d0)
C$$$  llm=(ll+1)**2 !! take one of possible var-ts
C$$$  icg1=indxcgr(lm2,lm1,1,0)
C$$$  icg2=indxcgr(lm2,lm1,2,0)
C$$$  do icg = icg1, icg2
C$$$  lmi = jcgr(icg)
C$$$  if(lmi==llm) then
C$$$  jcg1=indxcgr(lm3,lm4,1,0)
C$$$  jcg2=indxcgr(lm3,lm4,2,0)
C$$$  do jjcg = jcg1, jcg2
C$$$  lmj = jcgr(jjcg)
C$$$  if(lmj==llm) then
C$$$  do ind_nu=1,ndim3_nu
C$$$  i_nu=me3_tau*ndim3_nu+ind_nu-1
C$$$  ff(i_nu,ll0)=ff(i_nu,ll0)
C$$$  &                                +cc*conjg(cgr(icg))*cgr(jjcg)
C$$$  &                                *uloc(km1,km4,km3,km2,iat,iat,
C$$$  $                                ispin,jspin,ind_nu)
C$$$  enddo
C$$$  endif
C$$$  enddo
C$$$  endif
C$$$  enddo
C$$$  enddo   !! over ll
C$$$  enddo     !! over m1
C$$$  enddo       !! over m2
C$$$  enddo         !! over m3
C$$$  enddo           !! over m4
C$$$  if(nproc_tau/=1) then
C$$$  call mpi_allreduce_dcmplx(ff,4*(n_nu+1),
C$$$  $              mpi_sum,comm_k)
C$$$  endif
C$$$  if(maswrk) then
C$$$  do i_nu=0,n_nu
C$$$  write(3,'(10(1x,e12.5))')w_nu(i_nu)*evolt/2,
C$$$  &                (ff(i_nu,ll)*fac,ll=0,l)
C$$$  enddo
C$$$  endif
C$$$  endif
C$$$  enddo
C$$$  enddo  
C$$$  enddo    
C$$$  enddo      
C$$$  if(maswrk) close(3)
C$$$  deallocate(fct,ff)
C$$$  end
      



c$$$  subroutine cal_rotmat_mlms2jmj(ll, rotmat)
c$$$  
c$$$  implicit none
c$$$  integer, intent(in) :: ll
c$$$  complex*16, intent(out) :: rotmat(4*ll+2,4*ll+2)
c$$$  
c$$$  integer:: jj2,indms,indml,indmlms,indi,indmj,indjmj,
c$$$  $  iorb,jorb
c$$$  double precision :: jj,ulmu,fac1,fac2
c$$$  
c$$$  do indms=-1,1,2      
c$$$  do indml=-ll, ll
c$$$  indmlms=indml+ll+1+(indms+1)/2*(2*ll+1)
c$$$  do indi=-1,1,2
c$$$  jj2=2*ll+indi
c$$$  do indmj=-jj2, jj2, 2
c$$$  indjmj=(indmj+jj2)/2+1+(indi+1)/2*(jj2-1)
c$$$  ulmu=dble(indmj)/2.0d0/(dble(ll)+0.5d0)
c$$$  if (indi .eq. indms) then
c$$$  fac1=1.0d0
c$$$  else
c$$$  fac1=-1.0d0
c$$$  endif
c$$$  fac2=1.0d0
c$$$  if (indi .eq. -1 .and. indms .eq. 1) then
c$$$  fac2=-1.0d0
c$$$  endif
c$$$  if (indmj .eq. indms+indml*2) then
c$$$  rotmat(indjmj,indmlms)=1.0d0/dsqrt(2.0d0)
c$$$  $            *dsqrt(1+ulmu*fac1)*fac2
c$$$  endif
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  
c$$$  
c$$$  end      
c$$$  
c$$$  
c$$$  subroutine cal_rotmat_cmplx2real(norb, rotmat)
c$$$  implicit none
c$$$  integer, intent(in) :: norb
c$$$  complex*16, intent(out) :: rotmat(norb,norb)
c$$$  
c$$$  integer :: iorb,jorb
c$$$  
c$$$  rotmat=0.0d0
c$$$  
c$$$  rotmat(norb/2+1,norb/2+1)=1.0d0       
c$$$  do iorb=1, norb/2
c$$$  rotmat(iorb+norb/2+1,iorb+norb/2+1)
c$$$  $    =1.0d0/dsqrt(2.0d0)*(-1)**iorb
c$$$  rotmat(iorb+norb/2+1,-iorb+norb/2+1)
c$$$  $    =1.0d0/dsqrt(2.0d0)
c$$$  rotmat(-iorb+norb/2+1,iorb+norb/2+1)
c$$$  $    =1.0d0/dsqrt(2.0d0)*dcmplx(0.0d0, 1.0d0)*(-1)**(iorb+1)
c$$$  rotmat(-iorb+norb/2+1,-iorb+norb/2+1)
c$$$  $    =1.0d0/dsqrt(2.0d0)*dcmplx(0.0d0, 1.0d0)
c$$$  enddo
c$$$  end      







c$$$  
c$$$  subroutine output_u_f0f2f4f6_rel(uloc)!,key)
c$$$  use atom_mod
c$$$  use manager_mod
c$$$  use models_mod
c$$$  use parallel_mod
c$$$  use solid_mod
c$$$  use units_mod
c$$$  implicit none
c$$$  include "mpif.h"
c$$$  c$$$      integer, intent(in) :: key
c$$$  complex*16, intent(in) ::
c$$$  $  uloc(maxel_loc**2,maxel_loc**2,natom_loc,
c$$$  &  nspin_0,nspin_0,ndim3_nu)
c$$$  integer :: isort,l,ij,kl,i_len,iatom_c,iatom,nf,ispin,jspin,lm1,
c$$$  &  lm2,lm3,lm4,km1,km2,km3,km4,ll,lm_min,ll0,llm,ii,icg1,
c$$$  &  icg2,jcg1,jcg2,icg,jjcg,lmi,lmj,ind_nu,i_nu,li,jj,mj1,
c$$$  &  mj,lm,l22,lm22,lget,ndim,lf,ie,m4,m3,m2,m1
c$$$  double precision :: cc,n_l(0:3)
c$$$  double precision, allocatable :: fct(:)
c$$$  complex*16, allocatable :: ff(:,:)
c$$$  character*300 :: tempchar
c$$$  
c$$$  c$$$      if (key .eq. 0) then
c$$$  c$$$        tempchar='_r0'
c$$$  c$$$      else
c$$$  tempchar=''
c$$$  c$$$      endif
c$$$  allocate(fct(0:20))
c$$$  allocate(ff(0:n_nu,0:3))
c$$$  call facc(fct,20)
c$$$  i_len=len_trim(allfile)
c$$$  if(maswrk) then
c$$$  open(3,file=allfile(1:i_len)//'_F0_F2_F4_F6'//trim(tempchar)
c$$$  $    //'_'//trim(ubi)//'.dat')
c$$$  write(3,*)' F0 F2 F4 F6 representation for U_matrix'
c$$$  endif
c$$$  do jspin=1,nspin
c$$$  do ispin=1,nspin
c$$$  do isort=1,nsort
c$$$  iatom=iat_1(isort)
c$$$  
c$$$  if(atom_loc(iatom)) then
c$$$  iatom_c=iat_loc_back(iatom)
c$$$  do lf=1,nfun0_loc(isort)
c$$$  li=li_loc(lf,isort)
c$$$  ie=ie_loc(lf,isort)
c$$$  l=li/2
c$$$  ndim=((li+1)/2)*2
c$$$  nf=nfun_loc(isort)
c$$$  if(maswrk) then
c$$$  write(3,*)' '
c$$$  write(3,'(a15,i1,a1,i1,a11,i2,2(a6,i2))')
c$$$  &	            ' ispin/jspin = ',ispin,'/',jspin,' iatom_c = ',
c$$$  &              iatom_c,' Li = ',li,' Ie = ',ie
c$$$  write(3,*)' '
c$$$  write(3,*)' Frequency F0 F2 F4 F6: (in eV)'
c$$$  write(3,*)' '
c$$$  endif
c$$$  c     -------------- Get N_l ----------------------------------------
c$$$  do ll=0,l+l,2
c$$$  ii=mod(li,2)
c$$$  if(ii.eq.0) ii=-1
c$$$  jj=l+l+ii
c$$$  n_l(ll/2)=0.d0
c$$$  do mj1=-jj,jj,2
c$$$  call getlimj(lm1,l,ii,mj1,li,1)
c$$$  do mj=-jj,jj,2
c$$$  call getlimj(lm,l,ii,mj,li,1)
c$$$  icg1=indxcgr(lm,lm1,1,0)
c$$$  icg2=indxcgr(lm,lm1,2,0)
c$$$  do icg = icg1, icg2
c$$$  lm22 = jcgr(icg)
c$$$  l22=lget(lm22)
c$$$  if(l22==ll) n_l(ll/2)=n_l(ll/2)
c$$$  &                    +abs(cgr(icg))**2
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  c     ---------------------------------------------------------------	      
c$$$  lm_min=nrel*l*l
c$$$  if(mod(li,2)==1) lm_min=lm_min+2*l
c$$$  ff=0.d0
c$$$  do m4=1,ndim
c$$$  lm4=lm_min+m4
c$$$  km4=indbasa_loc(m4,lf,isort)
c$$$  do m3=1,ndim
c$$$  lm3=lm_min+m3
c$$$  km3=indbasa_loc(m3,lf,isort)
c$$$  kl=nf*(km4-1)+km3
c$$$  do m2=1,ndim
c$$$  lm2=lm_min+m2
c$$$  km2=indbasa_loc(m2,lf,isort)
c$$$  do m1=1,ndim
c$$$  lm1=lm_min+m1
c$$$  km1=indbasa_loc(m1,lf,isort)
c$$$  ij=nf*(km2-1)+km1
c$$$  do ll=0,l+l,2
c$$$  ll0=ll/2
c$$$  cc=(ll+ll+1.d0)/4.d0/pi/n_l(ll0)**2
c$$$  llm=(ll+1)**2 !! take one of possible var-ts
c$$$  icg1=indxcgr(lm2,lm1,1,0)
c$$$  icg2=indxcgr(lm2,lm1,2,0)
c$$$  do icg = icg1, icg2
c$$$  lmi = jcgr(icg)
c$$$  if(lmi==llm) then
c$$$  jcg1=indxcgr(lm3,lm4,1,0)
c$$$  jcg2=indxcgr(lm3,lm4,2,0)
c$$$  do jjcg = jcg1, jcg2
c$$$  lmj = jcgr(jjcg)
c$$$  if(lmj==llm) then
c$$$  do ind_nu=1,ndim3_nu
c$$$  i_nu=me3_tau*ndim3_nu+ind_nu-1
c$$$  ff(i_nu,ll0)=ff(i_nu,ll0)
c$$$  &                                +cc*conjg(cgr(icg))*cgr(jjcg)
c$$$  &                                *uloc(ij,kl,iatom_c,1,1,
c$$$  &                                ind_nu)
c$$$  enddo
c$$$  endif
c$$$  enddo
c$$$  endif
c$$$  enddo
c$$$  enddo   !! over ll
c$$$  enddo     !! over m1
c$$$  enddo       !! over m2
c$$$  enddo         !! over m3
c$$$  enddo           !! over m4
c$$$  if(nproc_tau/=1) then
c$$$  call mpi_allreduce_dcmplx
c$$$  $              (ff,4*(n_nu+1),mpi_sum,comm_k)
c$$$  endif
c$$$  if(maswrk) then
c$$$  do i_nu=0,n_nu
c$$$  write(3,'(10(1x,e12.5))')w_nu(i_nu)*evolt/2,
c$$$  &                (ff(i_nu,ll)*evolt/2,ll=0,l)
c$$$  enddo
c$$$  endif
c$$$  enddo             !! over lf
c$$$  endif
c$$$  enddo                 !! over isort
c$$$  enddo                   !! over ispin
c$$$  enddo                     !! over jspin
c$$$  if(maswrk) close(3)
c$$$  deallocate(fct,ff)
c$$$  end
      
