c$$$      subroutine boson_wpb_k(w_s,matdim,ind_k,
c$$$     $  s_multiply1,s_multiply2,qvopt_multiply,
c$$$     $  boson_wpb)
c$$$      
c$$$!     need s_multiply
c$$$!     need qvopt_multiply
c$$$      use atom_mod
c$$$      use manager_mod
c$$$      use parallel_mod
c$$$      use solid_mod
c$$$      use units_mod
c$$$      use vertex_mod
c$$$      use wannier90_mod
c$$$
c$$$      implicit none
c$$$      include "mpif.h"
c$$$      integer, intent(in) :: ind_k,
c$$$     $  s_multiply1,s_multiply2,qvopt_multiply,matdim
c$$$      complex*16, intent(in) :: w_s(matdim,matdim)
c$$$      complex*16, intent(out) ::
c$$$     $  boson_wpb(nwpb_max,natom_ciw,nwpb_max,natom_ciw,k_star_max)
c$$$
c$$$      integer :: ii,jj,k,n_pbt,num_qpoint,matsize,kst,iq,
c$$$     $  ispin,iorb,jspin,jorb,iorb1,iorb2,jorb1,jorb2,ind,ind_i,ind_j,
c$$$     $  iat,jat
c$$$      complex*16 :: v_tmp(n_pbtot,n_pbtot),tmp4(n_pbtot,nvopt_max),
c$$$     $  v_tmp2(n_pbtot,n_pbtot), s(nplw_gw,nplw_gw),
c$$$     $  smat(n_pbtot,n_pbtot),tmp2(n_pbtot,n_pbtot)
c$$$      complex*16,allocatable ::v_tmp3(:,:),v_tmp4(:,:)
c$$$
c$$$
c$$$      k=n3_mpi_k(me3_k+1)+ind_k
c$$$
c$$$      v_tmp=0.0d0
c$$$
c$$$      if (qvopt_multiply .eq. 1) then
c$$$        call zgemm('n','n',n_pbtot,matdim,matdim,(1.d0,0.d0),
c$$$     $    qvopt(1,1,ind_k),n_pbtot,w_s,matdim,(0.d0,0.d0),tmp4,
c$$$     $    n_pbtot)
c$$$        call zgemm('n','c',n_pbtot,n_pbtot,matdim,(1.d0,0.d0),tmp4,
c$$$     $    n_pbtot,qvopt(1,1,ind_k),n_pbtot,(0.d0,0.d0),v_tmp,n_pbtot)
c$$$      else
c$$$        do ii=1, matdim
c$$$          do jj=1, matdim
c$$$            v_tmp(ii,jj)=w_s(ii,jj)
c$$$          enddo
c$$$        enddo
c$$$      endif
c$$$
c$$$      boson_wpb=0.0d0
c$$$
c$$$c     ----------------------------------------------------------------------
c$$$      n_pbt=n_pbmt+nplwgw(k)
c$$$
c$$$      if(.not.include_gamma_w) then
c$$$        num_qpoint=nqdiv-1
c$$$      else
c$$$        num_qpoint=nqdiv
c$$$      endif
c$$$      matsize=nwpb_max*natom_ciw
c$$$
c$$$      allocate(v_tmp3(n_pbtot,matsize))
c$$$      allocate(v_tmp4(matsize,matsize))
c$$$
c$$$      do kst=1,k_star(k)
c$$$        iq=k_list(kst,k)
c$$$        v_tmp3=0.0d0
c$$$        v_tmp4=0.0d0
c$$$        call zone1_number(pnt(:,iq),rb0,ndiv,ind)
c$$$        if((.not.include_gamma_w).and.(iq.eq.1)) cycle
c$$$c     -------- Obtaining W for the point Q from W in Q0 -----------------
c$$$        v_tmp2=v_tmp
c$$$        call sym_w(iq,k,v_tmp2)
c$$$
c$$$        if ((s_multiply1 .eq. 1) .or.(s_multiply2 .eq. 1) ) then
c$$$          call s_overlap(0,s,nplw_gw,nplwgw(k),
c$$$     &      indpw_gw(1,iq))
c$$$          smat=0.0d0
c$$$
c$$$          do ii=1, n_pbmt
c$$$            smat(ii,ii)=1.0d0
c$$$          enddo
c$$$
c$$$          do ii=1, nplwgw(k)
c$$$            do jj=1, nplwgw(k)
c$$$              smat(ii+n_pbmt,jj+n_pbmt)=s(ii,jj)
c$$$            enddo
c$$$          enddo          
c$$$
c$$$          if (s_multiply1 .eq. 1) then
c$$$            tmp2=0.0d0            
c$$$            call zgemm('n','n',n_pbt,n_pbt,n_pbt,
c$$$     $        (1.0d0,0.0d0),smat,n_pbtot,
c$$$     $        v_tmp2,n_pbtot,(0.d0,0.d0),tmp2,n_pbtot)
c$$$            v_tmp2=tmp2
c$$$          endif
c$$$          
c$$$          
c$$$          if (s_multiply2 .eq. 1) then
c$$$            tmp2=0.0d0            
c$$$            call zgemm('n','n',n_pbt,n_pbt,n_pbt,
c$$$     $        (1.0d0,0.0d0),v_tmp2,n_pbtot,
c$$$     $        smat,n_pbtot,(0.d0,0.d0),tmp2,n_pbtot)
c$$$            v_tmp2=tmp2
c$$$          endif
c$$$        endif
c$$$        
c$$$        call zgemm('n','n',n_pbt,matsize,n_pbt,
c$$$     $    (1.d0,0.d0),v_tmp2,n_pbtot,
c$$$     $    m_pi(1,1,1,ind_k,kst),n_pbtot,
c$$$     $    (0.d0,0.d0),v_tmp3,n_pbtot)
c$$$        call zgemm('c','n',matsize,matsize,n_pbt,
c$$$     $    (1.d0,0.d0),m_pi(1,1,1,ind_k,kst),n_pbtot,
c$$$     $    v_tmp3,n_pbtot,
c$$$     $    (0.d0,0.d0),v_tmp4,matsize)
c$$$        
c$$$
c$$$        do iat=1, natom_ciw
c$$$          do iorb=1,nwpb(iat)
c$$$            call indexing(matsize,2,(/nwpb_max,natom_ciw/),1,ind_i,
c$$$     $        (/iorb,iat/))
c$$$            do jat=1, natom_ciw
c$$$              do jorb=1,nwpb(jat)
c$$$                call indexing(matsize,2,(/nwpb_max,natom_ciw/),1,ind_j,
c$$$     $            (/jorb,jat/))                
c$$$                boson_wpb(iorb,iat,jorb,jat,kst)
c$$$     $            =v_tmp4(ind_i,ind_j)/dble(num_qpoint)
c$$$              enddo
c$$$            enddo
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$      
c$$$      deallocate(v_tmp3,v_tmp4)
c$$$      end
c$$$

      subroutine boson_wan_k(w_s,matdim,ind_k,
     $  s_multiply1,s_multiply2,qvopt_multiply,
     $  boson_k_den, boson_local)
      
!     need s_multiply
!     need qvopt_multiply
      use atom_mod
!     use heg_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      use wannier90_mod

      implicit none
      include "mpif.h"
      integer, intent(in) :: ind_k,
     $  s_multiply1,s_multiply2,qvopt_multiply,matdim
      complex*16, intent(in) :: w_s(matdim,matdim)
      complex*16, intent(out) ::
     $  boson_k_den(num_orb_proj,nspin,num_orb_proj,nspin,k_star_max),
     $  boson_local(num_orb_proj,num_orb_proj,
     $  num_orb_proj,num_orb_proj,nspin,nspin)



      integer :: ii,jj,k,n_pbt,num_qpoint,matsize,kst,iq,
     $  ispin,iorb,jspin,jorb,iorb1,iorb2,jorb1,jorb2,ind,ind_i,ind_j
      complex*16 :: v_tmp(n_pbtot,n_pbtot),tmp5(n_pbtot,nvopt_max),
     $  v_tmp2(n_pbtot,n_pbtot), s(nplw_gw,nplw_gw),
     $  smat(n_pbtot,n_pbtot),tmp2(n_pbtot,n_pbtot)      
      complex*16,allocatable ::v_tmp3(:,:),v_tmp4(:,:)

      
      k=n3_mpi_k(me3_k+1)+ind_k

      v_tmp=0.0d0

      if (qvopt_multiply .eq. 1) then
        call zgemm('n','n',n_pbtot,matdim,matdim,(1.d0,0.d0),
     $    qvopt(1,1,ind_k),n_pbtot,w_s,matdim,(0.d0,0.d0),tmp5,
     $    n_pbtot)
        call zgemm('n','c',n_pbtot,n_pbtot,matdim,(1.d0,0.d0),tmp5,
     $    n_pbtot,qvopt(1,1,ind_k),n_pbtot,(0.d0,0.d0),v_tmp,n_pbtot)
      else
        do ii=1, matdim
          do jj=1, matdim
            v_tmp(ii,jj)=w_s(ii,jj)
          enddo
        enddo
      endif

      boson_k_den=0.0d0
      boson_local=0.0d0

c     ----------------------------------------------------------------------
      n_pbt=n_pbmt+nplwgw(k)

      if(.not.include_gamma_w) then
        num_qpoint=nqdiv-1
      else
        num_qpoint=nqdiv
      endif
      matsize=num_orb_proj**2*nspin

      allocate(v_tmp3(n_pbtot,matsize))
      allocate(v_tmp4(matsize,matsize))
      
      do kst=1,k_star(k)
        iq=k_list(kst,k)
        v_tmp3=0.0d0
        v_tmp4=0.0d0
        call zone1_number(pnt(:,iq),rb0,ndiv,ind)
        if((.not.include_gamma_w).and.(iq.eq.1)) cycle
c     -------- Obtaining W for the point Q from W in Q0 -----------------
        v_tmp2=v_tmp
        call sym_w(iq,k,v_tmp2)


        if ((s_multiply1 .eq. 1) .or.(s_multiply2 .eq. 1) ) then
          call s_overlap(0,s,nplw_gw,nplwgw(k),
     &      indpw_gw(1,iq))
          smat=0.0d0

          do ii=1, n_pbmt
            smat(ii,ii)=1.0d0
          enddo

          do ii=1, nplwgw(k)
            do jj=1, nplwgw(k)
              smat(ii+n_pbmt,jj+n_pbmt)=s(ii,jj)
            enddo
          enddo          
          
          if (s_multiply1 .eq. 1) then
            tmp2=0.0d0            
            call zgemm('n','n',n_pbt,n_pbt,n_pbt,
     $        (1.0d0,0.0d0),smat,n_pbtot,
     $        v_tmp2,n_pbtot,(0.d0,0.d0),tmp2,n_pbtot)
            v_tmp2=tmp2
          endif
          
          if (s_multiply2 .eq. 1) then
            tmp2=0.0d0            
            call zgemm('n','n',n_pbt,n_pbt,n_pbt,
     $        (1.0d0,0.0d0),v_tmp2,n_pbtot,
     $        smat,n_pbtot,(0.d0,0.d0),tmp2,n_pbtot)
            v_tmp2=tmp2
          endif
        endif        

        call zgemm('n','n',n_pbt,matsize,n_pbt,
     $    (1.d0,0.d0),v_tmp2,n_pbtot,
     $    prod_wan2(1,1,1,1,ind_k,kst),n_pbtot,
     $    (0.d0,0.d0),v_tmp3,n_pbtot)
        call zgemm('c','n',matsize,matsize,n_pbt,
     $    (1.d0,0.d0),prod_wan2(1,1,1,1,ind_k,kst),n_pbtot,
     $    v_tmp3,n_pbtot,
     $    (0.d0,0.d0),v_tmp4,matsize)

        do ispin=1, nspin
          do iorb=1, num_orb_proj
            call indexing(matsize,3,(/num_orb_proj,num_orb_proj,nspin/),
     $        1,ind_i,(/iorb,iorb,ispin/))
            do jspin=1, nspin
              do jorb=1, num_orb_proj
                call indexing(matsize,3,
     $            (/num_orb_proj,num_orb_proj,nspin/),1,ind_j,
     $            (/jorb,jorb,jspin/))
                boson_k_den(iorb,ispin,jorb,jspin,kst)
     $            =v_tmp4(ind_i,ind_j)/dble(num_qpoint)
              enddo
            enddo
          enddo
        enddo
        
        do ispin=1, nspin
          do iorb2=1, num_orb_proj ! 2
            do iorb1=1, num_orb_proj !1
              call indexing(matsize,3,
     $          (/num_orb_proj,num_orb_proj,nspin/),
     $          1,ind_i,(/iorb2,iorb1,ispin/))
              do jspin=1, nspin
                do jorb1=1, num_orb_proj ! 3
                  do jorb2=1, num_orb_proj ! 4
                    call indexing(matsize,3,
     $                (/num_orb_proj,num_orb_proj,nspin/),1,ind_j,
     $                (/jorb1,jorb2,jspin/))
                    boson_local(iorb1,iorb2,jorb1,jorb2,ispin,jspin)
     $                =boson_local(iorb1,iorb2,jorb1,jorb2,ispin,jspin)
     $                +v_tmp4(ind_i,ind_j)/dble(num_qpoint) ! 1(r)* 2(r) 3(r')* 4(r')
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      
      deallocate(v_tmp3,v_tmp4)
      end


c$$$  subroutine chi0_projection
c$$$  
c$$$  use atom_mod
c$$$  use manager_mod
c$$$  use parallel_mod
c$$$  use solid_mod
c$$$  use units_mod
c$$$  use vertex_mod
c$$$  use wannier90_mod
c$$$  
c$$$  implicit none
c$$$  include "mpif.h"
c$$$  
c$$$  integer :: ind_nu,i_nu,iat,iorb,jorb,iorb1,iorb2,jorb1,jorb2,
c$$$  $  ind1,ind2,jnd1,jnd2,ispin,jspin,ndim,iatom,isort,
c$$$  $  ind_i,ind_j,li,l,mudim,mudim_half,iwpb,
c$$$  $  m1,m2,m3,m4,km1,km2,km3,km4,ii
c$$$  
c$$$  complex*16 :: ff(0:3,0:n_nu)
c$$$  complex*16,allocatable :: 
c$$$  $  tempmat3(:,:),tempmat4(:,:,:,:),
c$$$  
c$$$  $  chi0_wpb(:,:),
c$$$  
c$$$  $  w_wpb(:,:),w_wpb_p(:,:),chi0_wpb_p(:,:),      
c$$$  $  eps_wpb(:,:),epsinv_wpb(:,:),u_wpb(:,:),uloc_temp(:,:,:,:),
c$$$  $  uloc_temp2(:,:,:,:),
c$$$  $  rotmat(:,:),projection(:,:),chi0_wpb_nu(:,:,:)
c$$$  character(len=5) :: lchar,catom
c$$$  
c$$$  integer, external :: find_ind
c$$$  
c$$$  do iat=1, natom_ciw
c$$$  iatom=iatom_ciw(iat)
c$$$  isort=is(iatom)
c$$$  
c$$$  if (iatom .lt. 10)then
c$$$  write(catom,'(i1)') iatom
c$$$  elseif (iatom .lt. 10)then
c$$$  write(catom,'(i2)') iatom
c$$$  elseif (iatom .lt. 100)then
c$$$  write(catom,'(i3)') iatom
c$$$  endif        
c$$$  
c$$$  do l=0, lmb(isort)
c$$$  
c$$$  write(lchar,'(i2)') l          
c$$$  
c$$$  if ((corr_iter(l,iatom) .eq. 2)
c$$$  $      .and.
c$$$  $      (at_bonding(l,iatom) .eq. 0)
c$$$  $      ) then
c$$$  
c$$$  if(maswrk) then
c$$$  open(3,file=adjustl(trim(allfile))
c$$$  $          //'_chi0_low_Slater_'
c$$$  $          //trim(adjustl(txtel(isort)))//trim(adjustl(catom))
c$$$  $          //'_'//trim(adjustl(lchar))//'.dat')
c$$$  write(3,'(a)')
c$$$  $          '# Slater integral representation for U_matrix'
c$$$  endif
c$$$  
c$$$  
c$$$  
c$$$  allocate(tempmat3(n_ciw_max**2,nwpb(iat)))
c$$$  allocate(tempmat4(n_ciw_max,n_ciw_max,n_ciw_max,n_ciw_max))
c$$$  allocate(w_wpb(nwpb(iat),nwpb(iat)))
c$$$  allocate(w_wpb_p(nwpb(iat),nwpb(iat)))            
c$$$  allocate(chi0_wpb(nwpb(iat),nwpb(iat)))
c$$$  allocate(chi0_wpb_nu(nwpb(iat),nwpb(iat),ndim3_nu))            
c$$$  allocate(chi0_wpb_p(nwpb(iat),nwpb(iat)))            
c$$$  allocate(eps_wpb(nwpb(iat),nwpb(iat)))
c$$$  allocate(epsinv_wpb(nwpb(iat),nwpb(iat)))
c$$$  allocate(u_wpb(nwpb(iat),nwpb(iat)))
c$$$  
c$$$  call chi0_loc_correlated (l,iat,chi0_wpb_nu)
c$$$  
c$$$  
c$$$  if (nrel .eq. 1) then          
c$$$  li=l+1
c$$$  mudim=l+l+1
c$$$  
c$$$  
c$$$  allocate(rotmat(mudim,mudim))
c$$$  rotmat=0.0d0
c$$$  call cal_rotmat_cmplx2real(mudim, rotmat)
c$$$  else
c$$$  
c$$$  mudim=4*l+2
c$$$  mudim_half=mudim/2
c$$$  
c$$$  allocate(rotmat(mudim,mudim))
c$$$  rotmat=0.0d0
c$$$  
c$$$  call cal_rotmat_mlms2jmj(l, rotmat)
c$$$  
c$$$  allocate(uloc_temp2(mudim_half,mudim_half,
c$$$  $          mudim_half,mudim_half))
c$$$  uloc_temp2=0.0d0              
c$$$  endif
c$$$  
c$$$  allocate(uloc_temp(mudim,mudim,mudim,mudim))            
c$$$  
c$$$  allocate(projection(nwpb(iat),nwpb(iat)))
c$$$  
c$$$  
c$$$  do jspin=1,nspin
c$$$  do ispin=1,nspin
c$$$  
c$$$  if(maswrk) then
c$$$  write(3,'(a15,i1,a1,i1,a11,i2,2(a6,i2))')
c$$$  &              '# ispin/jspin = ',ispin,'/',jspin,' iatom = ',
c$$$  &              iatom,' L = ',l
c$$$  endif
c$$$  
c$$$  ff=0.d0
c$$$  
c$$$  do ind_nu=1,ndim3_nu
c$$$  i_nu=me3_tau*ndim3_nu+ind_nu-1
c$$$  
c$$$  w_wpb=0.0d0
c$$$  w_wpb_p=0.0d0                  
c$$$  
c$$$  chi0_wpb_p=0.0d0
c$$$  chi0_wpb=0.0d0
c$$$  
c$$$  
c$$$  
c$$$  do iorb=1,nwpb(iat)
c$$$  do jorb=1,nwpb(iat)
c$$$  chi0_wpb(iorb,jorb)
c$$$  $                  =chi0_wpb_nu(iorb,jorb,ind_nu)
c$$$  enddo
c$$$  enddo
c$$$  
c$$$  ndim=n_ciw_max*n_ciw_max                  
c$$$  
c$$$  call zgemm('n','n',ndim,nwpb(iat),nwpb(iat),
c$$$  $              (1.d0,0.d0),pi_ww(1,1,1,jspin,iat),ndim,chi0_wpb,
c$$$  $              nwpb(iat),(0.0d0,0.0d0),tempmat3,ndim)                  
c$$$  
c$$$  call zgemm('n','c',ndim,ndim,nwpb(iat),(1.d0,0.d0),
c$$$  $              tempmat3,ndim,pi_ww(1,1,1,ispin,iat),ndim,
c$$$  $              (0.0d0,0.0d0),tempmat4,ndim)                  
c$$$  
c$$$  uloc_temp=0.0d0
c$$$  do m4=1,mudim
c$$$  if (nrel .eq. 1) then
c$$$  km4=find_ind(mt2ciw(m4,l+1,iatom),
c$$$  $                  n_ciw_atom(iatom),ciw_ind(1:n_ciw_atom(iatom),
c$$$  $                  iat))
c$$$  else
c$$$  km4=find_ind(mt2ciw_rel_l(m4,l+1,iatom),
c$$$  $                  n_ciw_atom(iatom),ciw_ind(1:n_ciw_atom(iatom),
c$$$  $                  iat))                      
c$$$  endif
c$$$  do m3=1,mudim
c$$$  if (nrel .eq. 1) then
c$$$  km3=find_ind(mt2ciw(m3,l+1,iatom),
c$$$  $                    n_ciw_atom(iatom),ciw_ind(1:n_ciw_atom(iatom),
c$$$  $                    iat))                        
c$$$  else
c$$$  km3=find_ind(mt2ciw_rel_l(m3,l+1,iatom),
c$$$  $                    n_ciw_atom(iatom),ciw_ind(1:n_ciw_atom(iatom),
c$$$  $                    iat))                                                
c$$$  endif                      
c$$$  do m2=1,mudim
c$$$  if (nrel .eq. 1) then
c$$$  km2=find_ind(mt2ciw(m2,l+1,iatom),
c$$$  $                      n_ciw_atom(iatom),
c$$$  $                      ciw_ind(1:n_ciw_atom(iatom),iat))                                                  
c$$$  else
c$$$  km2=find_ind(mt2ciw_rel_l(m2,l+1,iatom),
c$$$  $                      n_ciw_atom(iatom),
c$$$  $                      ciw_ind(1:n_ciw_atom(iatom),iat))                                                                          
c$$$  endif                                              
c$$$  do m1=1,mudim
c$$$  if (nrel .eq. 1) then
c$$$  km1=find_ind(mt2ciw(m1,l+1,iatom),
c$$$  $                        n_ciw_atom(iatom),
c$$$  $                        ciw_ind(1:n_ciw_atom(iatom),iat))                                                    
c$$$  else
c$$$  km1=find_ind(mt2ciw_rel_l(m1,l+1,iatom),
c$$$  $                        n_ciw_atom(iatom),
c$$$  $                        ciw_ind(1:n_ciw_atom(iatom),iat))                                                                            
c$$$  endif                                              
c$$$  uloc_temp(m1,m2,m3,m4)
c$$$  $                      =uloc_temp(m1,m2,m3,m4)
c$$$  $                      +tempmat4(km2,km1,km3,km4)
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  if (nrel .eq. 1) then
c$$$  call rotate_umatrix(mudim,uloc_temp,
c$$$  $                rotmat,1)
c$$$  call cal_slater_integral(uloc_temp,l,ff(0,i_nu))
c$$$  else
c$$$  call rotate_umatrix(mudim,uloc_temp,
c$$$  $                rotmat,1)
c$$$  
c$$$  uloc_temp2=0.0d0
c$$$  uloc_temp2
c$$$  $                =(
c$$$  $                uloc_temp(1:mudim_half,1:mudim_half,
c$$$  $                1:mudim_half,1:mudim_half)
c$$$  $                +uloc_temp((mudim_half+1):mudim,
c$$$  $                (mudim_half+1):mudim,1:mudim_half,
c$$$  $                1:mudim_half)
c$$$  $                +uloc_temp(1:mudim_half,1:mudim_half,
c$$$  $                (mudim_half+1):mudim,
c$$$  $                (mudim_half+1):mudim)              
c$$$  $                +uloc_temp((mudim_half+1):mudim,
c$$$  $                (mudim_half+1):mudim,
c$$$  $                (mudim_half+1):mudim,(mudim_half+1):mudim)
c$$$  $                )/4.0d0
c$$$  
c$$$  call cal_slater_integral(uloc_temp2,l,ff(0,i_nu))
c$$$  
c$$$  endif
c$$$  enddo
c$$$  
c$$$  if(nproc_tau/=1) then
c$$$  call mpi_allreduce_dcmplx
c$$$  $              (ff,4*(n_nu+1),mpi_sum,comm_k)
c$$$  endif
c$$$  if(maswrk) then
c$$$  do i_nu=0,n_nu
c$$$  write(3,'(10(1x,e12.5))')w_nu(i_nu)*evolt/2,
c$$$  &                (dble(ff(ii,i_nu))/rydberg,ii=0,l)
c$$$  enddo
c$$$  endif                  
c$$$  enddo
c$$$  enddo                
c$$$  deallocate(tempmat3)
c$$$  deallocate(tempmat4)
c$$$  deallocate(w_wpb)
c$$$  deallocate(w_wpb_p)            
c$$$  deallocate(chi0_wpb)
c$$$  deallocate(chi0_wpb_nu)            
c$$$  deallocate(chi0_wpb_p)            
c$$$  deallocate(eps_wpb)
c$$$  deallocate(epsinv_wpb)
c$$$  deallocate(u_wpb)
c$$$  deallocate(uloc_temp)
c$$$  deallocate(rotmat)
c$$$  deallocate(projection)
c$$$  if (nrel .eq. 2) then
c$$$  deallocate(uloc_temp2)              
c$$$  endif
c$$$  if (maswrk) close(3)
c$$$  endif
c$$$  enddo
c$$$  enddo
c$$$  end      
c$$$      
c$$$
c$$$
c$$$      subroutine u_weiss
c$$$      
c$$$      use atom_mod
c$$$      use manager_mod
c$$$      use parallel_mod
c$$$      use solid_mod
c$$$      use units_mod
c$$$      use vertex_mod
c$$$      use wannier90_mod
c$$$
c$$$      implicit none
c$$$      include "mpif.h"
c$$$
c$$$      integer :: ind_nu,i_nu,iat,iorb,jorb,iorb1,iorb2,jorb1,jorb2,
c$$$     $  ind1,ind2,jnd1,jnd2,ispin,jspin,ndim,iatom,isort,
c$$$     $  ind_i,ind_j,li,l,mudim,mudim_half,iwpb,
c$$$     $  m1,m2,m3,m4,km1,km2,km3,km4,ii,
c$$$     $  iat_c,
c$$$     $  dim_pi(natom_cw)
c$$$
c$$$      complex*16 :: ff(0:3,0:n_nu),cc(0:3,0:n_nu),
c$$$     $  pi_projection(nwpb_max,nwpb_max,natom_cw),
c$$$     $  pi_ww_corr(nwpb_max,n_cw_max,n_cw_max,nspin,natom_cw),
c$$$     $  pi_ww_corr_dual(nwpb_max,n_cw_max,n_cw_max,nspin,natom_cw)      
c$$$      complex*16,allocatable ::
c$$$     $  chi0_wpb(:,:,:), 
c$$$      
c$$$     $  chi0_wpb_p(:,:) ,
c$$$     $  chi0_wpb_pp(:,:), 
c$$$
c$$$     $  w_wpb(:,:), 
c$$$     $  w_wpb_p(:,:),
c$$$     $  w_wpb_pp(:,:),
c$$$
c$$$     $  tempmat1(:,:),
c$$$     $  tempmat2(:,:), 
c$$$      
c$$$     $  tempmat3(:,:),
c$$$     $  tempmat4(:,:),
c$$$      
c$$$     $  eps_wpb_pp(:,:), 
c$$$     $  epsinv_wpb_pp(:,:),
c$$$     $  u_wpb_pp(:,:), 
c$$$     $  uloc_temp(:,:,:,:), 
c$$$     $  uloc_temp2(:,:,:,:), 
c$$$     $  rotmat(:,:)
c$$$      character(len=5) :: lchar,catom
c$$$
c$$$      integer, external :: find_ind
c$$$
c$$$      call cal_projection_and_pi_ww_corr(dim_pi, pi_projection,
c$$$     $  pi_ww_corr, pi_ww_corr_dual)
c$$$
c$$$      if (maswrk) then
c$$$        write(iun,*) 'here2'
c$$$      endif      
c$$$
c$$$      do iat_c=1, natom_cw
c$$$        iatom=iatom_cw(iat_c)
c$$$        iat=find_ind(iatom, natom_ciw, iatom_ciw)        
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
c$$$          write(lchar,'(i2)') l          
c$$$
c$$$          if ((corr(l,iatom) .eq. 2)
c$$$     $      .and.
c$$$     $      (at_bonding(l,iatom) .eq. 0)
c$$$     $      ) then
c$$$            
c$$$            if(maswrk) then
c$$$              open(3,file=adjustl(trim(allfile))
c$$$     $          //'_u_Slater_'
c$$$     $          //trim(adjustl(txtel(isort)))//trim(adjustl(catom))
c$$$     $          //'_'//trim(adjustl(lchar))//'.dat')
c$$$              write(3,'(a)')
c$$$     $          '# Slater integral representation for U_matrix'
c$$$              open(4,file=adjustl(trim(allfile))
c$$$     $          //'_chi0_Slater_'
c$$$     $          //trim(adjustl(txtel(isort)))//trim(adjustl(catom))
c$$$     $          //'_'//trim(adjustl(lchar))//'.dat')
c$$$              write(4,'(a)')
c$$$     $          '# Slater integral representation for chi0_matrix'              
c$$$            endif
c$$$
c$$$            ndim=n_cw_max*n_cw_max
c$$$
c$$$            allocate(chi0_wpb(nwpb(iat),nwpb(iat),ndim3_nu))
c$$$            allocate(chi0_wpb_p(nwpb(iat),dim_pi(iat_c)))
c$$$            allocate(chi0_wpb_pp(dim_pi(iat_c),dim_pi(iat_c)))
c$$$            
c$$$            allocate(w_wpb(nwpb(iat),nwpb(iat)))
c$$$            allocate(w_wpb_p(nwpb(iat),dim_pi(iat_c)))
c$$$            allocate(w_wpb_pp(dim_pi(iat_c),dim_pi(iat_c)))            
c$$$
c$$$            allocate(eps_wpb_pp(dim_pi(iat_c),dim_pi(iat_c)))
c$$$            allocate(epsinv_wpb_pp(dim_pi(iat_c),dim_pi(iat_c)))
c$$$            allocate(u_wpb_pp(dim_pi(iat_c),dim_pi(iat_c)))                                    
c$$$
c$$$
c$$$            allocate(tempmat3(dim_pi(iat_c),ndim))
c$$$            allocate(tempmat4(ndim,ndim))
c$$$            allocate(tempmat1(dim_pi(iat_c),ndim))
c$$$            allocate(tempmat2(ndim,ndim))            
c$$$
c$$$
c$$$            call chi0_loc_correlated (l,iat,chi0_wpb)
c$$$
c$$$            if (maswrk) then
c$$$              write(iun,*) 'here3'
c$$$            endif            
c$$$            
c$$$            if (nrel .eq. 1) then          
c$$$              li=l+1
c$$$              mudim=l+l+1
c$$$
c$$$              allocate(rotmat(mudim,mudim))
c$$$              rotmat=0.0d0
c$$$              call cal_rotmat_cmplx2real(mudim, rotmat)
c$$$            else
c$$$              mudim=4*l+2
c$$$              mudim_half=mudim/2
c$$$              
c$$$              allocate(rotmat(mudim,mudim))
c$$$              rotmat=0.0d0
c$$$              
c$$$              call cal_rotmat_mlms2jmj(l, rotmat)
c$$$              
c$$$              allocate(uloc_temp2(mudim_half,mudim_half,
c$$$     $          mudim_half,mudim_half))
c$$$              uloc_temp2=0.0d0              
c$$$            endif
c$$$            allocate(uloc_temp(mudim,mudim,mudim,mudim))                                    
c$$$
c$$$            if (maswrk) then
c$$$              write(iun,*) 'here4'
c$$$            endif
c$$$
c$$$c$$$  allocate(projection(nwpb(iat),nwpb(iat)))
c$$$c$$$  
c$$$c$$$  call projection_coulomb(iat,l,projection)
c$$$
c$$$            
c$$$            do jspin=1,nspin
c$$$              do ispin=1,nspin
c$$$                
c$$$                if(maswrk) then
c$$$                  write(3,'(a15,i1,a1,i1,a11,i2,2(a6,i2))')
c$$$     &              '# ispin/jspin = ',ispin,'/',jspin,' iatom = ',
c$$$     &              iatom,' L = ',l
c$$$                  write(4,'(a15,i1,a1,i1,a11,i2,2(a6,i2))')
c$$$     &              '# ispin/jspin = ',ispin,'/',jspin,' iatom = ',
c$$$     &              iatom,' L = ',l
c$$$                endif
c$$$                
c$$$                ff=0.0d0
c$$$                cc=0.0d0
c$$$                
c$$$                do ind_nu=1,ndim3_nu
c$$$                  i_nu=me3_tau*ndim3_nu+ind_nu-1
c$$$                  
c$$$                  w_wpb=0.0d0
c$$$                  w_wpb_p=0.0d0                  
c$$$                  
c$$$                  do iorb=1,nwpb(iat)
c$$$                    do jorb=1,nwpb(iat)
c$$$                      w_wpb(iorb,jorb)
c$$$     $                  =w_wpb_local(iorb,iat,jorb,iat,ind_nu)
c$$$                    enddo
c$$$                  enddo
c$$$
c$$$                  if (maswrk) then
c$$$                    write(iun,*) 'here5'
c$$$                  endif
c$$$                  
c$$$                  call zgemm('n','n',nwpb(iat),dim_pi(iat_c),nwpb(iat),
c$$$     $              (1.0d0,0.0d0),w_wpb,nwpb(iat),
c$$$     $              pi_projection(1,1,iat_c),
c$$$     $              nwpb_max,(0.0d0,0.0d0),w_wpb_p,nwpb(iat))
c$$$
c$$$                  call zgemm('c','n',dim_pi(iat_c),dim_pi(iat_c),
c$$$     $              nwpb(iat),
c$$$     $              (1.0d0,0.0d0),pi_projection(1,1,iat_c),nwpb_max,
c$$$     $              w_wpb_p,nwpb(iat),
c$$$     $              (0.0d0,0.0d0),w_wpb_pp,dim_pi(iat_c))
c$$$                  
c$$$                  call zgemm('n','n',nwpb(iat),dim_pi(iat_c),nwpb(iat),
c$$$     $              (1.0d0,0.0d0),chi0_wpb(1,1,ind_nu),nwpb(iat),
c$$$     $              pi_projection(1,1,iat_c),nwpb_max,(0.0d0,0.0d0),
c$$$     $              chi0_wpb_p,nwpb(iat))
c$$$
c$$$                  call zgemm('c','n',dim_pi(iat_c),dim_pi(iat_c),
c$$$     $              nwpb(iat),
c$$$     $              (1.0d0,0.0d0),pi_projection(1,1,iat_c),nwpb_max,
c$$$     $              chi0_wpb_p,nwpb(iat),
c$$$     $              (0.0d0,0.0d0),chi0_wpb_pp,dim_pi(iat_c))
c$$$c$$$
c$$$c$$$
c$$$                  call zgemm('n','n',dim_pi(iat_c),ndim,dim_pi(iat_c),
c$$$     $              (1.d0,0.d0),chi0_wpb_pp,dim_pi(iat_c),
c$$$     $              pi_ww_corr_dual(1,1,1,ispin,iat_c),nwpb_max,
c$$$     $              (0.0d0,0.0d0),tempmat1,dim_pi(iat_c))                  
c$$$                  
c$$$                  call zgemm('c','n',ndim,ndim,dim_pi(iat_c),
c$$$     $              (1.d0,0.d0),pi_ww_corr_dual(1,1,1,jspin,iat_c),
c$$$     $              nwpb_max,
c$$$     $              tempmat1,dim_pi(iat_c),
c$$$     $              (0.0d0,0.0d0),tempmat2,ndim)                                    
c$$$
c$$$
c$$$                  if (maswrk) then
c$$$                    write(iun,*) 'here6'
c$$$                  endif
c$$$                  eps_wpb_pp=0.0d0
c$$$                  do iwpb=1, dim_pi(iat_c)
c$$$                    eps_wpb_pp(iwpb,iwpb)=1.0d0
c$$$                  enddo
c$$$                  
c$$$                  call zgemm('n','n',
c$$$     $              dim_pi(iat_c),dim_pi(iat_c),dim_pi(iat_c),
c$$$     $              (1.0d0,0.0d0),w_wpb_pp,dim_pi(iat_c),chi0_wpb_pp,
c$$$     $              dim_pi(iat_c),(1.0d0,0.0d0),
c$$$     $              eps_wpb_pp,dim_pi(iat_c))
c$$$                  
c$$$                  epsinv_wpb_pp=0.0d0
c$$$                  
c$$$                  call dcmplx_matinv(eps_wpb_pp, epsinv_wpb_pp,
c$$$     $              dim_pi(iat_c),dim_pi(iat_c))
c$$$                  
c$$$                  u_wpb_pp=0.0d0
c$$$
c$$$                  if (maswrk) then
c$$$                    write(iun,*) 'here7'
c$$$                  endif                  
c$$$                  
c$$$                  call zgemm('n','n',
c$$$     $              dim_pi(iat_c),dim_pi(iat_c),dim_pi(iat_c),                  
c$$$     $              (1.0d0,0.0d0),epsinv_wpb_pp,dim_pi(iat_c),w_wpb_pp,
c$$$     $              dim_pi(iat_c),(0.0d0,0.0d0),u_wpb_pp,dim_pi(iat_c))
c$$$
c$$$                  call zgemm('n','n',dim_pi(iat_c),ndim,dim_pi(iat_c),
c$$$     $              (1.d0,0.d0),u_wpb_pp,dim_pi(iat_c),
c$$$     $              pi_ww_corr(1,1,1,ispin,iat_c),nwpb_max,
c$$$     $              (0.0d0,0.0d0),tempmat3,dim_pi(iat_c))                  
c$$$                  
c$$$                  call zgemm('c','n',ndim,ndim,dim_pi(iat_c),
c$$$     $              (1.d0,0.d0),pi_ww_corr(1,1,1,jspin,iat_c),nwpb_max,
c$$$     $              tempmat3,dim_pi(iat_c),
c$$$     $              (0.0d0,0.0d0),tempmat4,ndim)                  
c$$$                  
c$$$                  uloc_temp=0.0d0
c$$$                  do m4=1,mudim
c$$$                    do m3=1,mudim
c$$$                      call indexing(ndim,3,
c$$$     $                  (/n_cw_max,n_cw_max,nspin/),
c$$$     $                  1,ind1,(/m3,m4,ispin/))
c$$$                      do m2=1,mudim
c$$$                        do m1=1,mudim
c$$$                          call indexing(ndim,3,
c$$$     $                      (/n_cw_max,n_cw_max,nspin/),
c$$$     $                      1,ind2,(/m1,m2,jspin/))                            
c$$$                          uloc_temp(m2,m1,m3,m4)
c$$$     $                      =uloc_temp(m2,m1,m3,m4)
c$$$     $                      +tempmat4(ind2,ind1)
c$$$                        enddo
c$$$                      enddo
c$$$                    enddo
c$$$                  enddo
c$$$
c$$$                  if (maswrk) then
c$$$                    write(iun,*) 'here8'
c$$$                  endif                 
c$$$!!!!!!from here
c$$$                  if (nrel .eq. 1) then
c$$$                    call rotate_umatrix(mudim,uloc_temp,
c$$$     $                rotmat,1)
c$$$                    call cal_slater_integral(uloc_temp,l,ff(0,i_nu))
c$$$                  else
c$$$                    call rotate_umatrix(mudim,uloc_temp,
c$$$     $                rotmat,1)
c$$$
c$$$                    uloc_temp2=0.0d0
c$$$                    uloc_temp2
c$$$     $                =(
c$$$     $                uloc_temp(1:mudim_half,1:mudim_half,
c$$$     $                1:mudim_half,1:mudim_half)
c$$$     $                +uloc_temp((mudim_half+1):mudim,
c$$$     $                (mudim_half+1):mudim,1:mudim_half,
c$$$     $                1:mudim_half)
c$$$     $                +uloc_temp(1:mudim_half,1:mudim_half,
c$$$     $                (mudim_half+1):mudim,
c$$$     $                (mudim_half+1):mudim)              
c$$$     $                +uloc_temp((mudim_half+1):mudim,
c$$$     $                (mudim_half+1):mudim,
c$$$     $                (mudim_half+1):mudim,(mudim_half+1):mudim)
c$$$     $                )/4.0d0
c$$$                    
c$$$                    call cal_slater_integral(uloc_temp2,l,ff(0,i_nu))
c$$$
c$$$                    if (maswrk) then
c$$$                      write(iun,*) 'here9'
c$$$                    endif                    
c$$$                    
c$$$                  endif
c$$$
c$$$
c$$$                  uloc_temp=0.0d0
c$$$                  do m4=1,mudim
c$$$                    do m3=1,mudim
c$$$                        call indexing(ndim,3,
c$$$     $                    (/n_cw_max,n_cw_max,nspin/),
c$$$     $                    1,ind1,(/m3,m4,ispin/))
c$$$                      do m2=1,mudim
c$$$                        do m1=1,mudim
c$$$                            call indexing(ndim,3,
c$$$     $                        (/n_cw_max,n_cw_max,nspin/),
c$$$     $                        1,ind2,(/m1,m2,jspin/))                            
c$$$                          uloc_temp(m2,m1,m3,m4)
c$$$     $                      =uloc_temp(m2,m1,m3,m4)
c$$$     $                      +tempmat2(ind2,ind1)
c$$$                        enddo
c$$$                      enddo
c$$$                    enddo
c$$$                  enddo
c$$$  
c$$$!!!!!!from here
c$$$                  if (nrel .eq. 1) then
c$$$                    call rotate_umatrix(mudim,uloc_temp,
c$$$     $                rotmat,1)
c$$$                    call cal_slater_integral(uloc_temp,l,cc(0,i_nu))
c$$$                  else
c$$$                    call rotate_umatrix(mudim,uloc_temp,
c$$$     $                rotmat,1)
c$$$                    
c$$$                    uloc_temp2=0.0d0
c$$$                    uloc_temp2
c$$$     $                =(
c$$$     $                uloc_temp(1:mudim_half,1:mudim_half,
c$$$     $                1:mudim_half,1:mudim_half)
c$$$     $                +uloc_temp((mudim_half+1):mudim,
c$$$     $                (mudim_half+1):mudim,1:mudim_half,
c$$$     $                1:mudim_half)
c$$$     $                +uloc_temp(1:mudim_half,1:mudim_half,
c$$$     $                (mudim_half+1):mudim,
c$$$     $                (mudim_half+1):mudim)              
c$$$     $                +uloc_temp((mudim_half+1):mudim,
c$$$     $                (mudim_half+1):mudim,
c$$$     $                (mudim_half+1):mudim,(mudim_half+1):mudim)
c$$$     $                )/4.0d0
c$$$                    
c$$$                    call cal_slater_integral(uloc_temp2,l,cc(0,i_nu))
c$$$                    
c$$$                    if (maswrk) then
c$$$                      write(iun,*) 'here9'
c$$$                    endif                    
c$$$                    
c$$$                  endif
c$$$                  
c$$$                enddo
c$$$                
c$$$                if(nproc_tau/=1) then
c$$$                  call mpi_allreduce_dcmplx
c$$$     $              (ff(0,0),4*(n_nu+1),mpi_sum,comm_k)
c$$$                  call mpi_allreduce_dcmplx
c$$$     $              (cc,4*(n_nu+1),mpi_sum,comm_k)                  
c$$$                endif
c$$$                if(maswrk) then
c$$$                  do i_nu=0,n_nu
c$$$                    write(3,'(10(1x,e12.5))')w_nu(i_nu)*evolt/2,
c$$$     &                (dble(ff(ii,i_nu))*rydberg,ii=0,l)
c$$$                    write(4,'(10(1x,e12.5))')w_nu(i_nu)*evolt/2,
c$$$     &                (dble(cc(ii,i_nu))/rydberg,ii=0,l)
c$$$                  enddo
c$$$                endif                  
c$$$              enddo
c$$$            enddo
c$$$
c$$$            deallocate(chi0_wpb)
c$$$            
c$$$            deallocate(tempmat3)
c$$$            deallocate(tempmat4)
c$$$            deallocate(tempmat1)
c$$$            deallocate(tempmat2)            
c$$$            deallocate(w_wpb)
c$$$            deallocate(w_wpb_p)
c$$$            deallocate(w_wpb_pp)                        
c$$$
c$$$            deallocate(chi0_wpb_p)
c$$$            deallocate(chi0_wpb_pp)                        
c$$$            deallocate(eps_wpb_pp)
c$$$            deallocate(epsinv_wpb_pp)
c$$$            deallocate(u_wpb_pp)
c$$$            deallocate(uloc_temp)
c$$$            deallocate(rotmat)
c$$$            if (nrel .eq. 2) then
c$$$              deallocate(uloc_temp2)              
c$$$            endif
c$$$            if (maswrk) close(3)
c$$$            if (maswrk) close(4)            
c$$$          endif
c$$$        enddo
c$$$      enddo
c$$$      end
c$$$
c$$$
      

c$$$  subroutine projection_test
c$$$  
c$$$  use atom_mod
c$$$  use manager_mod
c$$$  use parallel_mod
c$$$  use solid_mod
c$$$  use units_mod
c$$$  use vertex_mod
c$$$  use wannier90_mod
c$$$  
c$$$  implicit none
c$$$  include "mpif.h"
c$$$  
c$$$  integer :: ind_nu,i_nu,iat,iorb,jorb,iorb1,iorb2,jorb1,jorb2,
c$$$  $  ind1,ind2,jnd1,jnd2,ispin,jspin,ndim,iatom,isort,
c$$$  $  ind_i,ind_j,li,l,mudim,mudim_half,iwpb,
c$$$  $  m1,m2,m3,m4,km1,km2,km3,km4,ii
c$$$  
c$$$  complex*16 :: ff(0:3,0:n_nu)
c$$$  complex*16,allocatable :: tempmat1(:,:),
c$$$  $  tempmat2(:,:),tempmat3(:,:),tempmat4(:,:,:,:),
c$$$  $  w_wpb(:,:), chi0_wpb(:,:),eps_wpb(:,:),
c$$$  $  epsinv_wpb(:,:),u_wpb(:,:),uloc_temp(:,:,:,:),
c$$$  $  uloc_temp2(:,:,:,:),
c$$$  $  rotmat(:,:),projection(:,:)
c$$$  character(len=5) :: lchar,catom
c$$$  
c$$$  integer, external :: find_ind
c$$$  
c$$$  do iat=1, natom_ciw
c$$$  iatom=iatom_ciw(iat)
c$$$  isort=is(iatom)
c$$$  write(iun,*) 'projection_test', iat
c$$$  
c$$$  call projection_wpb(iat)
c$$$  
c$$$  enddo
c$$$  end



c$$$  subroutine projection_coulomb(iat,l,proj_coulomb)
c$$$  use wannier90_mod
c$$$  use manager_mod
c$$$  use parallel_mod
c$$$  implicit none
c$$$  
c$$$  integer, intent(in) :: iat
c$$$  integer, intent(in) :: l
c$$$  complex*16, intent(out) :: proj_coulomb(nwpb(iat),nwpb(iat))
c$$$  
c$$$  integer :: iorb,jorb,korb,lorb,iwpb,jwpb,m1,m2,ispin,jspin,iatom,
c$$$  $  ndim
c$$$  c$$$      complex*16 :: proj_coulomb(nwpb(iat),nwpb(iat))
c$$$  complex*16, allocatable ::
c$$$  $  pi_ww_at_dual(:,:,:,:),
c$$$  $  pi_ww_at(:,:,:,:)      
c$$$  c$$$     $  identity_test(:,:,:,:,:,:),
c$$$  c$$$     $  resolution(:,:)
c$$$  
c$$$  
c$$$  integer, external :: find_ind
c$$$  
c$$$  iatom=iatom_ciw(iat)
c$$$  
c$$$  
c$$$  allocate(pi_ww_at(nwpb(iatom),
c$$$  $  n_ciw_atom(iatom),n_ciw_atom(iatom),nspin))
c$$$  pi_ww_at=0.0d0
c$$$  allocate(pi_ww_at_dual(nwpb(iatom),
c$$$  $  n_ciw_atom(iatom),n_ciw_atom(iatom),nspin))
c$$$  pi_ww_at_dual=0.0d0      
c$$$  
c$$$  
c$$$  
c$$$  c$$$      allocate(identity_test(n_ciw_atom(iatom),n_ciw_atom(iatom),nspin,
c$$$  c$$$     $  n_ciw_atom(iatom),n_ciw_atom(iatom),nspin))
c$$$  c$$$      allocate(resolution(nwpb(iat),nwpb(iat)))
c$$$  
c$$$  
c$$$  c$$$  projection=0.0d0
c$$$  pi_ww_at=0.0d0      
c$$$  pi_ww_at_dual=0.0d0
c$$$  c$$$      identity_test=0.0d0
c$$$  c$$$  resolution=0.0d0
c$$$  proj_coulomb=0.0d0
c$$$  
c$$$  do ispin=1, nspin
c$$$  do iorb=1, n_ciw_atom(iatom)
c$$$  do jorb=1, n_ciw_atom(iatom)
c$$$  do iwpb=1, nwpb(iat)
c$$$  pi_ww_at(iwpb,iorb,jorb,ispin)
c$$$  $          =dconjg(pi_ww(iorb,jorb,iwpb,ispin,iat))
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  enddo
c$$$  ndim=n_ciw_atom(iatom)**2*nspin
c$$$  
c$$$  if (maswrk) then
c$$$  write(iun,'(a, 2i5)') 'dual vector for iatom, l: ', iat, l
c$$$  endif
c$$$  call dual_vector(nwpb(iat),ndim,pi_ww_at,pi_ww_at_dual)
c$$$  
c$$$  c$$$      call zgemm('c','n',n_ciw_atom(iatom)**2*nspin,
c$$$  c$$$     $  n_ciw_atom(iatom)**2*nspin,
c$$$  c$$$     $  nwpb(iat),(1.0d0,0.0d0),pi_ww_at,nwpb(iat),
c$$$  c$$$     $  pi_ww_at_dual,nwpb(iat),(0.0d0,0.0d0),
c$$$  c$$$     $  identity_test,n_ciw_atom(iatom)**2*nspin)
c$$$  
c$$$  c$$$      call zgemm('n','c',nwpb(iat),nwpb(iat),n_ciw_atom(iatom)**2*nspin,
c$$$  c$$$     $  (1.0d0,0.0d0),pi_ww_at,nwpb(iat),
c$$$  c$$$     $  pi_ww_at_dual,nwpb(iat),(0.0d0,0.0d0),
c$$$  c$$$     $  resolution,nwpb(iat))
c$$$  c$$$
c$$$  c$$$!     resolution test
c$$$  c$$$      do iwpb=1, nwpb(iat)
c$$$  c$$$        do jwpb=1, nwpb(iat)
c$$$  c$$$          
c$$$  c$$$          if (iwpb .eq. jwpb) then
c$$$  c$$$            if (cdabs(resolution(iwpb,jwpb)-1.0d0) .gt. 1.0d-6) then
c$$$  c$$$              write(iun,'(a, 2i5, 2f12.6)')
c$$$  c$$$     $          'resolution matrix fail', iwpb,jwpb,
c$$$  c$$$     $          resolution(iwpb,jwpb)
c$$$  c$$$            endif
c$$$  c$$$          else
c$$$  c$$$            if (cdabs(resolution(iwpb,jwpb)) .gt. 1.0d-6) then            
c$$$  c$$$              write(iun,'(a, 2i5, 2f12.6)')
c$$$  c$$$     $          'resolution matrix fail', iwpb,jwpb,
c$$$  c$$$     $          resolution(iwpb,jwpb)
c$$$  c$$$            endif
c$$$  c$$$          endif
c$$$  c$$$        enddo
c$$$  c$$$      enddo
c$$$  c$$$
c$$$  c$$$
c$$$  c$$$      
c$$$  !     identity operator
c$$$  c$$$
c$$$  c$$$      do ispin=1, nspin
c$$$  c$$$        do iorb=1, n_ciw_atom(iatom)
c$$$  c$$$          do jorb=1, n_ciw_atom(iatom)
c$$$  c$$$            do jspin=1, nspin
c$$$  c$$$              do korb=1, n_ciw_atom(iatom)
c$$$  c$$$                do lorb=1, n_ciw_atom(iatom)
c$$$  c$$$                  
c$$$  c$$$                  if ((ispin .eq. jspin) .and.
c$$$  c$$$     $              (iorb .eq. korb) .and.
c$$$  c$$$     $              (jorb .eq. lorb)) then
c$$$  c$$$                    if (cdabs(identity_test
c$$$  c$$$     $                (iorb,jorb,ispin,korb,lorb,jspin)-1.0d0)
c$$$  c$$$     $                .gt. 1.0d-6) then
c$$$  c$$$                      write(iun,'(a, 6i5, 2f12.6)')
c$$$  c$$$     $                  'identity_test fail',
c$$$  c$$$     $                  iorb,jorb,ispin,korb,lorb,jspin,
c$$$  c$$$     $                  identity_test(iorb,jorb,ispin,korb,lorb,jspin)
c$$$  c$$$                    endif
c$$$  c$$$                  else
c$$$  c$$$                    if (cdabs(identity_test
c$$$  c$$$     $                (iorb,jorb,ispin,korb,lorb,jspin))
c$$$  c$$$     $                .gt. 1.0d-6) then
c$$$  c$$$                      write(iun,'(a, 6i5, 2f12.6)')
c$$$  c$$$     $                  'identity_test fail',
c$$$  c$$$     $                  iorb,jorb,ispin,korb,lorb,jspin,
c$$$  c$$$     $                  identity_test(iorb,jorb,ispin,korb,lorb,jspin)
c$$$  c$$$                    endif                      
c$$$  c$$$                  endif
c$$$  c$$$                enddo
c$$$  c$$$              enddo
c$$$  c$$$            enddo
c$$$  c$$$          enddo
c$$$  c$$$        enddo
c$$$  c$$$      enddo      
c$$$  
c$$$  
c$$$  do iwpb=1, nwpb(iat)
c$$$  do jwpb=1, nwpb(iat)
c$$$  
c$$$  do ispin=1, nspin
c$$$  if (nrel .eq. 1) then            
c$$$  do m1=1, 2*l+1
c$$$  
c$$$  iorb=find_ind
c$$$  $            (mt2ciw(m1,l+1,iatom),n_ciw_atom(iatom),
c$$$  $            ciw_ind(1:n_ciw_atom(iatom),iat))
c$$$  
c$$$  do m2=1, 2*l+1
c$$$  jorb=find_ind
c$$$  $              (mt2ciw(m2,l+1,iatom),n_ciw_atom(iatom),
c$$$  $              ciw_ind(1:n_ciw_atom(iatom),iat))
c$$$  
c$$$  proj_coulomb(iwpb,jwpb)
c$$$  $              = proj_coulomb(iwpb,jwpb)
c$$$  $              +pi_ww_at(iwpb,iorb,jorb,ispin)
c$$$  $              *dconjg(pi_ww_at_dual(jwpb,iorb,jorb,ispin))
c$$$  enddo
c$$$  enddo
c$$$  else
c$$$  do m1=1,4*l+2
c$$$  
c$$$  iorb=find_ind
c$$$  $            (mt2ciw_rel_l(m1,l+1,iatom),n_ciw_atom(iatom),
c$$$  $            ciw_ind(1:n_ciw_atom(iatom),iat))
c$$$  
c$$$  do m2=1,4*l+2
c$$$  jorb=find_ind
c$$$  $              (mt2ciw_rel_l(m2,l+1,iatom),n_ciw_atom(iatom),
c$$$  $              ciw_ind(1:n_ciw_atom(iatom),iat))
c$$$  
c$$$  proj_coulomb(iwpb,jwpb)
c$$$  $              = proj_coulomb(iwpb,jwpb)
c$$$  $              +pi_ww_at(iwpb,iorb,jorb,ispin)
c$$$  $              *dconjg(pi_ww_at_dual(jwpb,iorb,jorb,ispin))                  
c$$$  enddo
c$$$  enddo
c$$$  endif
c$$$  enddo
c$$$  
c$$$  c$$$          if (iwpb .eq. jwpb) then
c$$$  c$$$            if (cdabs(proj_coulomb(iwpb,jwpb)-1.0d0) .gt. 1.0d-6) then
c$$$  c$$$              write(iun,'(a, 2i5, 2f12.6)')
c$$$  c$$$     $          'proj_coulomb matrix fail', iwpb,jwpb,
c$$$  c$$$     $          proj_coulomb(iwpb,jwpb)
c$$$  c$$$            endif
c$$$  c$$$          else
c$$$  c$$$            if (cdabs(proj_coulomb(iwpb,jwpb)) .gt. 1.0d-6) then            
c$$$  c$$$              write(iun,'(a, 2i5, 2f12.6)')
c$$$  c$$$     $          'proj_coulomb matrix fail', iwpb,jwpb,
c$$$  c$$$     $          proj_coulomb(iwpb,jwpb)
c$$$  c$$$            endif
c$$$  c$$$          endif          
c$$$  
c$$$  enddo
c$$$  enddo
c$$$  
c$$$  deallocate(pi_ww_at)      
c$$$  deallocate(pi_ww_at_dual)
c$$$  c$$$      deallocate(identity_test)
c$$$  c$$$      deallocate(resolution)
c$$$  
c$$$  end




      subroutine cal_projection_and_pi_ww_corr(dim_pi, pi_projection,
     $  pi_ww_corr, pi_ww_corr_dual)
      use wannier90_mod
      use manager_mod
      use parallel_mod
      use solid_mod      
      implicit none

      integer, intent(out) :: dim_pi(natom_cw)
      complex*16, intent(out) ::
     $  pi_projection(nwpb_max,nwpb_max,natom_cw),
     $  pi_ww_corr(nwpb_max,n_cw_max,n_cw_max,nspin,natom_cw),
     $  pi_ww_corr_dual(nwpb_max,n_cw_max,n_cw_max,nspin,natom_cw)
      
      integer :: iat_c,iat,iatom,iorb,jorb,korb,lorb,iwpb,jwpb,m1,m2,
     $  ispin,jspin,ndim, iorb_i, jorb_j, ndim2,
     $  dim_pi_max,ii,jj,kk,
     $  maxdim,mindim,info,lwork

      complex*16 ::
     $  pi_projection_temp(nwpb_max,nwpb_max,natom_cw),
     $  pi_ww_corr_temp(nwpb_max,n_cw_max,n_cw_max,nspin,natom_cw),
     $  pi_ww_corr_dual_temp(nwpb_max,n_cw_max,n_cw_max,nspin,natom_cw)      
      
      double precision, allocatable :: ss(:),rwork(:)
      complex*16,allocatable :: work(:),worktemp(:),matin_temp(:,:),
     $  uu(:,:),vt(:,:)
      complex*16, allocatable ::
     $  pi_ww_at(:,:,:,:)

      complex*16, allocatable :: testmat(:,:)
      complex*16 :: testval
      integer, external :: find_ind

      pi_projection_temp=0.0d0
      pi_ww_corr_temp=0.0d0
      pi_ww_corr_dual_temp=0.0d0
      dim_pi=0
      
      do iat_c=1, natom_cw
        
        iatom=iatom_cw(iat_c)
        iat=find_ind(iatom, natom_ciw, iatom_ciw)
        
        ndim=n_cw_atom(iatom)

        allocate(pi_ww_at
     $    (nwpb(iat),ndim,ndim,nspin))
        pi_ww_at=0.0d0
        
        do ispin=1, nspin
          do iorb=1, ndim
            iorb_i=find_ind(cw_ind(iorb, iat_c),
     $        n_ciw_atom(iatom), ciw_ind(1,iat))
            do jorb=1, ndim
              jorb_j=find_ind(cw_ind(jorb, iat_c),
     $          n_ciw_atom(iatom), ciw_ind(1,iat))            
              do iwpb=1, nwpb(iat)
                pi_ww_at(iwpb,iorb,jorb,ispin) ! iwpb and iorb are with *
     $            =dconjg(pi_ww(iorb_i,jorb_j,iwpb,ispin,iat))
              enddo
            enddo
          enddo
        enddo

        ndim2=ndim**2*nspin
        
        maxdim=max(nwpb(iat),ndim2)
        mindim=min(nwpb(iat),ndim2)
        
        allocate(ss(mindim))
        ss=0.0d0
        allocate(rwork(5*mindim))
        rwork=0.0d0

        allocate(matin_temp(nwpb(iat),ndim2))
        matin_temp=0.0d0

        do ispin=1, nspin
          do iorb=1, ndim
            do jorb=1, ndim
              call indexing(ndim2,3,(/ndim,ndim,nspin/),1,ii,
     $          (/jorb,iorb,ispin/))                          
              do iwpb=1, nwpb(iat)        
                matin_temp(iwpb,ii)=pi_ww_at(iwpb,jorb,iorb,ispin)
              enddo
            enddo
          enddo
        enddo
        
        allocate(worktemp(3*maxdim))
        worktemp=0.0d0
        allocate(uu(nwpb(iat),nwpb(iat)))
        uu=0.0d0
        allocate(vt(ndim2,ndim2))
        vt=0.0d0        
        
        call zgesvd('A','A',nwpb(iat),ndim2,matin_temp,nwpb(iat),
     $    ss,uu,nwpb(iat),vt,ndim2,worktemp,-1,rwork,info)
        
        lwork=max(1, int(worktemp(1)))
        allocate(work(lwork))
        work=0.0d0


        do ispin=1, nspin
          do iorb=1, ndim
            do jorb=1, ndim
              call indexing(ndim2,3,(/ndim,ndim,nspin/),1,ii,
     $          (/jorb,iorb,ispin/))                          
              do iwpb=1, nwpb(iat)        
                matin_temp(iwpb,ii)=pi_ww_at(iwpb,jorb,iorb,ispin)
              enddo
            enddo
          enddo
        enddo
        
        call zgesvd('A','A',nwpb(iat),ndim2,matin_temp,nwpb(iat),
     $    ss,uu,nwpb(iat),vt,ndim2,work,lwork,rwork,info)
        deallocate(work)
        deallocate(worktemp)      
        deallocate(rwork)
        deallocate(matin_temp)
        
c$$$  if (maswrk) then
c$$$  write(iun,'(a, 10000f20.12)')
c$$$  $      'singular value product basis', ss
c$$$  endif
        
        do ii=1, mindim
          if (ss(ii)>0.00001d0) then
            dim_pi(iat_c)=dim_pi(iat_c)+1
          endif
        enddo
        
        do ii=1, dim_pi(iat_c)
          pi_projection_temp(:,ii,iat_c)=uu(:,ii)
        enddo
        
        do ii=1, dim_pi(iat_c)
          do iorb=1, n_cw_atom(iatom)
            do jorb=1, n_cw_atom(iatom)
              do ispin=1, nspin
                call indexing(ndim2,3,
     $            (/n_cw_atom(iatom),n_cw_atom(iatom),nspin/),1,jj,
     $            (/iorb,jorb,ispin/))            
                pi_ww_corr_temp(ii,iorb,jorb,ispin,iat_c)
     $            =pi_ww_corr_temp(ii,iorb,jorb,ispin,iat_c)
     $            +ss(ii)*vt(ii,jj)
                pi_ww_corr_dual_temp(ii,iorb,jorb,ispin,iat_c)
     $            =pi_ww_corr_dual_temp(ii,iorb,jorb,ispin,iat_c)
     $            +1.0d0/ss(ii)*vt(ii,jj)            
              enddo
            enddo
          enddo
        enddo
        deallocate(uu)
        deallocate(vt)
        deallocate(ss)
        deallocate(pi_ww_at)
      enddo

      dim_pi_max=maxval(dim_pi)

      pi_projection=0.0d0
      pi_ww_corr=0.0d0

      pi_ww_corr_dual=0.0d0

      do iat_c=1, natom_cw
        iatom=iatom_cw(iat_c)
        iat=find_ind(iatom, natom_ciw, iatom_ciw)
        do ii=1, dim_pi(iat_c)        
          do jj=1,nwpb(iat)
            pi_projection(jj,ii,iat_c)
     $        =pi_projection(jj,ii,iat_c)
     $        +pi_projection_temp(jj,ii,iat_c)
          enddo
        enddo

      enddo

      do iat_c=1, natom_cw
        iatom=iatom_cw(iat_c)
        iat=find_ind(iatom, natom_ciw, iatom_ciw)
        do ispin=1, nspin
          do ii=1, n_cw_atom(iatom)
            do jj=1, n_cw_atom(iatom)
              do kk=1, dim_pi(iat_c)
                pi_ww_corr(kk,jj,ii,ispin,iat_c)
     $            =pi_ww_corr(kk,jj,ii,ispin,iat_c)
     $            +pi_ww_corr_temp(kk,jj,ii,ispin,iat_c)
                pi_ww_corr_dual(kk,jj,ii,ispin,iat_c)
     $            =pi_ww_corr_dual(kk,jj,ii,ispin,iat_c)
     $            +pi_ww_corr_dual_temp(kk,jj,ii,ispin,iat_c)
              enddo
            enddo
          enddo
        enddo

      enddo
      
      end



      


c$$$  subroutine dual_vector(dim1,dim2,matin,matout)
c$$$  use manager_mod
c$$$  use parallel_mod
c$$$  implicit none
c$$$  integer, intent(in) :: dim1,dim2
c$$$  complex*16, intent(in) :: matin(dim1,dim2)
c$$$  complex*16, intent(out) :: matout(dim1,dim2)
c$$$  
c$$$  integer :: info,lwork,maxdim,mindim,ii
c$$$  complex*16 :: matin_temp(dim1,dim2),uu(dim1,dim1),vt(dim2,dim2),
c$$$  $  tempmat(dim1,dim2),ss_inv_mat(dim1,dim2)
c$$$  double precision, allocatable :: ss(:),rwork(:)
c$$$  complex*16,allocatable :: work(:),worktemp(:)
c$$$  
c$$$  
c$$$  maxdim=max(dim1,dim2)
c$$$  mindim=min(dim1,dim2)      
c$$$  allocate(ss(mindim))
c$$$  allocate(rwork(5*mindim))
c$$$  
c$$$  matin_temp=matin
c$$$  
c$$$  allocate(worktemp(3*maxdim))
c$$$  
c$$$  call zgesvd('A','A',dim1,dim2,matin_temp,dim1,
c$$$  $  ss,uu,dim1,vt,dim2,worktemp,-1,rwork,info)
c$$$  
c$$$  lwork=max(1, int(worktemp(1)))
c$$$  allocate(work(lwork))
c$$$  
c$$$  call zgesvd('A','A',dim1,dim2,matin_temp,dim1,
c$$$  $  ss,uu,dim1,vt,dim2,work,lwork,rwork,info)
c$$$  deallocate(work)
c$$$  deallocate(worktemp)      
c$$$  
c$$$  deallocate(rwork)
c$$$  
c$$$  if (maswrk) then
c$$$  write(iun,'(a, 10000f20.12)')
c$$$  $    'singular value product basis', ss
c$$$  endif
c$$$  
c$$$  ss_inv_mat=0.0d0
c$$$  
c$$$  do ii=1, mindim
c$$$  ss_inv_mat(ii,ii)=1.0d0/ss(ii)
c$$$  enddo
c$$$  
c$$$  deallocate(ss)
c$$$  
c$$$  call zgemm('n','n',dim1,dim2,dim1,
c$$$  $  (1.0d0,0.0d0),uu,dim1,ss_inv_mat,dim1,
c$$$  $  (0.0d0,0.0d0),tempmat,dim1)
c$$$  call zgemm('n','n',dim1,dim2,dim2,
c$$$  $  (1.0d0,0.0d0),tempmat,dim1,vt,dim2,
c$$$  $  (0.0d0,0.0d0),matout,dim1)      
c$$$  
c$$$  end      


      subroutine chi0_loc_correlated(llval,iat, pola_wpb_nu)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      use wannier90_mod      
      implicit none
      include 'mpif.h'

      integer, intent(in) ::  llval,iat
      complex*16, intent(out) ::
     $  pola_wpb_nu(nwpb(iat),nwpb(iat),ndim3_nu)
      

      integer :: mudim,iatom,num_kpoint,iwan1,iwan2,jwan1,jwan2,itau,
     $  ispin,ind_k,ik,ibnd,m1,m2,m3,m4,iorb,jorb,iwpb,ind_tau,
     $  mtau,tau1,tau2,tau3,tau4,ind_nu

      double precision :: de

      double precision ::
     $  pola_wpb_tau_packed(nwpb(iat),nwpb(iat),ndim3_tau),
     $  pola_wpb_nu_packed(nwpb(iat),nwpb(iat),ndim3_nu)
      complex*16 :: pola_wpb_tau(nwpb(iat),nwpb(iat),ndim3_tau)


      complex*16, allocatable :: g_tau(:,:,:,:),pi_ww_at(:,:,:,:),
     $  pola_tau(:,:,:,:),tempmat(:,:)
      
      double precision, external :: green0
      
      integer, external :: find_ind

      if (nrel .eq. 1) then
        mudim=2*llval+1
      else
        mudim=4*llval+2
      endif

      pola_wpb_tau_packed=0.0d0
      pola_wpb_nu_packed=0.0d0
      pola_wpb_tau=0.0d0
      pola_wpb_nu=0.0d0      

      allocate(g_tau(mudim,mudim,0:n_tau,nspin))
      g_tau=0.0d0
      allocate(pi_ww_at(nwpb(iat),mudim,mudim,nspin))
      pi_ww_at=0.0d0
      allocate(pola_tau(mudim,mudim,mudim,mudim))
      pola_tau=0.0d0
      allocate(tempmat(nwpb(iat),mudim**2))
      tempmat=0.0d0

      pola_wpb_tau=0.0d0
      pola_wpb_tau_packed=0.0d0

      iatom=iatom_ciw(iat)

!     g_loc

      num_kpoint=nqdiv

      do ispin=1, nspin
        do ind_k=1,ndim_kk(me+1) ! k vector
          ik=n_mpi_kk(me+1)+ind_k
          do ibnd=1, num_bands
            de=eigenvalues(ibnd,ik)/rydberg

            do iwan1=1,mudim
              if (nrel .eq. 1) then              
                iwan2=mt2ciw(iwan1,llval+1,iatom)
              else
                iwan2=mt2ciw_rel_l(iwan1,llval+1,iatom)
              endif
              do jwan1=1, mudim
                if (nrel .eq. 1) then              
                  jwan2=mt2ciw(jwan1,llval+1,iatom)
                else
                  jwan2=mt2ciw_rel_l(jwan1,llval+1,iatom)
                endif
c$$$  
                do itau=0,n_tau
c$$$  write(179,'(9i5)')
c$$$  $              ispin,ind_k,ik,ibnd,iwan1,iwan2,jwan1,jwan2,itau
c$$$  c$$$                  write(179,*) g_tau(iwan1,jwan1,itau,ispin) ! good
c$$$  write(179,*) green0(de,tau_mesh(itau))                              
c$$$  c$$$                  write(179,*) proj_renorm(ibnd,iwan2,ind_k,ispin) ! good
c$$$  c$$$                  write(179,*) proj_renorm(ibnd,jwan2,ind_k,ispin)!
                  g_tau(iwan1,jwan1,itau,ispin)                  
     $              =g_tau(iwan1,jwan1,itau,ispin)
     $              +dconjg(proj_renorm(ibnd,iwan2,ind_k,ispin))
     $              *green0(de,tau_mesh(itau))
     $              *proj_renorm(ibnd,jwan2,ind_k,ispin)
     $              /dble(num_kpoint)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo

      call mpi_allreduce_dcmplx(g_tau,mudim**2*(n_tau+1)*nspin,mpi_sum,
     $  mpi_comm_world)

      do ispin=1, nspin
        do m1=1, mudim
          if (nrel .eq. 1) then
            if (maswrk) then
              write(iun,'(1000i5)')
     $          mt2ciw(m1,llval+1,iatom), n_ciw_atom(iatom),
     $          ciw_ind(1:n_ciw_atom(iatom),iat)
            endif
            iorb=find_ind
     $        (mt2ciw(m1,llval+1,iatom),n_ciw_atom(iatom),
     $        ciw_ind(1:n_ciw_atom(iatom),iat))
          else
            iorb=find_ind
     $        (mt2ciw_rel_l(m1,llval+1,iatom),n_ciw_atom(iatom),
     $        ciw_ind(1:n_ciw_atom(iatom),iat))
          endif
          do m2=1, mudim
            if (nrel .eq. 1) then
              if (maswrk) then
                write(iun,'(1000i5)')
     $            mt2ciw(m2,llval+1,iatom), n_ciw_atom(iatom),
     $            ciw_ind(1:n_ciw_atom(iatom),iat)
              endif              
              jorb=find_ind
     $          (mt2ciw(m2,llval+1,iatom),n_ciw_atom(iatom),
     $          ciw_ind(1:n_ciw_atom(iatom),iat))
            else
              jorb=find_ind
     $          (mt2ciw_rel_l(m2,llval+1,iatom),n_ciw_atom(iatom),
     $          ciw_ind(1:n_ciw_atom(iatom),iat))
            endif
            do iwpb=1, nwpb(iat)
              pi_ww_at(iwpb,m1,m2,ispin)
     $          =dconjg(pi_ww(iorb,jorb,iwpb,ispin,iat))
            enddo
          enddo
        enddo
      enddo      

      do ind_tau=1,ndim3_tau
        itau=me3_tau*ndim3_tau+ind_tau-1
        mtau=n_tau-itau
        pola_tau=0.0d0
        do ispin=1, nspin
          do tau1=1, mudim
            do tau2=1, mudim
              do tau3=1, mudim
                do tau4=1, mudim
                  pola_tau(tau4,tau1,tau3,tau2)
     $              =pola_tau(tau4,tau1,tau3,tau2)                  
     $              -g_tau(tau1,tau2,itau,ispin)
     $              *g_tau(tau3,tau4,mtau,ispin)
     $              *2.0d0/dble(nspin)/dble(nrel)
                enddo
              enddo
            enddo
          enddo
          call zgemm('n','n',nwpb(iat),mudim**2,mudim**2,(1.d0,0.d0),
     $      pi_ww_at(1,1,1,ispin),nwpb(iat),pola_tau,mudim**2,
     $      (0.d0,0.d0),tempmat,nwpb(iat))
          call zgemm('n','c',nwpb(iat),nwpb(iat),mudim**2,(1.d0,0.d0),
     $      tempmat,nwpb(iat),pi_ww_at(1,1,1,ispin),nwpb(iat),
     $      (0.d0,0.d0),pola_wpb_tau(1,1,ind_tau),nwpb(iat))
        enddo
      enddo
      do ind_tau=1,ndim3_tau      
        call boson_pack_tau(1.0d0,pola_wpb_tau(1,1,ind_tau),nwpb(iat),
     $    0.0d0,pola_wpb_tau_packed(1,1,ind_tau),nwpb(iat),nwpb(iat))
      enddo
      
      call from_tau_to_nu_ba(pola_wpb_nu_packed,nwpb(iat),
     $  pola_wpb_tau_packed,nwpb(iat),nwpb(iat))
      
      do ind_nu=1,ndim3_nu
        call unpack_hermit(pola_wpb_nu(1,1,ind_nu),
     $    pola_wpb_nu_packed(1,1,ind_nu),
     $    nwpb(iat),nwpb(iat),nwpb(iat),
     &    (0.d0,0.d0),(1.d0,0.d0))
      enddo
      deallocate(g_tau)
      deallocate(pi_ww_at)
      deallocate(pola_tau)
      deallocate(tempmat)
      end


      subroutine chi0_from_g_correlated(pola_tau_vop_packed)
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      use wannier90_mod      
      implicit none
      include 'mpif.h'


      double precision, intent(out) ::
     $  pola_tau_vop_packed(nvopt_max,nvopt_max,ndim3_tau,
     $  ndim3_k(me3_k+1))

      integer :: ispin,ind_k,ik,ibnd,iwan1,jwan1,itau,ind_tau,mtau,
     $  ir1,ir2,ir3,mr1,mr2,mr3,ir_ind,mr_ind,tau1,tau2,tau3,tau4,k,
     $  ind,iq,ii,iorb,jorb

      double precision :: de,fac

      complex*16 ::
     $  g_tau(nqdiv,num_orb_cor,num_orb_cor,0:n_tau,nspin),
     $  pola_tau(nqdiv,num_orb_cor,num_orb_cor,
     $  num_orb_cor,num_orb_cor,nspin),
     $  tempmat1(num_orb_cor,num_orb_cor,num_orb_cor,num_orb_cor),
     $  s(nplw_gw,nplw_gw),
     $  vop_wan2(nvopt_max,num_orb_proj,num_orb_proj),
     $  vop_wan2_cor(nvopt_max,num_orb_cor,num_orb_cor),
     $  tmp1(nvopt_max,nplw_gw),
     $  tempmat2(nvopt_max,num_orb_cor,num_orb_cor),
     $  pola_tau_vop(nvopt_max,nvopt_max,ndim3_k(me3_k+1))

      double precision, external :: green0
      
      g_tau=0.0d0

      do ispin=1, nspin
        do ind_k=1,ndim_kk(me+1) ! k vector
          ik=n_mpi_kk(me+1)+ind_k
          call zone1_number(pnt(:,ik),rb0,ndiv,ind)          
          do ibnd=1, num_bands
            de=eigenvalues(ibnd,ik)/rydberg
            
            do iwan1=1,num_orb_cor
              do jwan1=1,num_orb_cor
                do itau=0,n_tau                                  
                  g_tau(ind,iwan1,jwan1,itau,ispin)
     $              =g_tau(ind,iwan1,jwan1,itau,ispin)
     $              +dconjg(proj_renorm(ibnd,iwan1,ind_k,ispin))
     $              *green0(de,tau_mesh(itau))
     $              *proj_renorm(ibnd,jwan1,ind_k,ispin)
     $              /dble(nqdiv)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      
      if(nproc/=1) call mpi_allreduce_dcmplx
     $  (g_tau,nqdiv*num_orb_cor**2*nspin*(n_tau+1),
     $  mpi_sum,mpi_comm_world)

      call fft3(ndiv(1),ndiv(2),ndiv(3),num_orb_cor**2*nspin*(n_tau+1),
     $  g_tau,-1)


      pola_tau_vop_packed=0.0d0
      
      do ind_tau=1,ndim3_tau
        itau=me3_tau*ndim3_tau+ind_tau-1
        mtau=n_tau-itau

        pola_tau=0.0d0

        pola_tau_vop=0.0d0
        
        do ispin=1, nspin                
          
          do ir1=0,ndiv(1)-1
            do ir2=0,ndiv(2)-1
              do ir3=0,ndiv(3)-1
                mr1=ndiv(1)-ir1
                if (mr1 .eq. ndiv(1)) mr1=0
                mr2=ndiv(2)-ir2
                if (mr2 .eq. ndiv(2)) mr2=0                        
                mr3=ndiv(3)-ir3
                if (mr3 .eq. ndiv(3)) mr3=0                                                
                
                call indexing
     $            (nqdiv,3,ndiv,1,ir_ind,(/ir1+1,ir2+1,ir3+1/))
                call indexing
     $            (nqdiv,3,ndiv,1,mr_ind,(/mr1+1,mr2+1,mr3+1/))
                do tau1=1, num_orb_cor
                  do tau2=1, num_orb_cor
                    do tau3=1, num_orb_cor 
                      do tau4=1, num_orb_cor
                        pola_tau(ir_ind,tau4,tau1,tau3,tau2,ispin)
     $                    =pola_tau(ir_ind,tau4,tau1,tau3,tau2,ispin)
     $                    -g_tau(ir_ind,tau1,tau2,itau,ispin)
     $                    *g_tau(mr_ind,tau3,tau4,mtau,ispin)
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo

        call fft3(ndiv(1),ndiv(2),ndiv(3),num_orb_cor**4*nspin,
     $    pola_tau,1)        
        
        do ind_k=1,ndim3_k(me3_k+1)
          k=n3_mpi_k(me3_k+1)+ind_k
          
          call zone1_number(pnt(:,k),rb0,ndiv,ind)

          do ispin=1,nspin
            tempmat1=0.0d0            
            do tau1=1, num_orb_cor
              do tau2=1, num_orb_cor
                do tau3=1, num_orb_cor 
                  do tau4=1, num_orb_cor          
                    tempmat1(tau1,tau2,tau3,tau4)
     $                =pola_tau(ind,tau1,tau2,tau3,tau4,ispin)
                  enddo
                enddo          
              enddo
            enddo

            iq=k_list(1,k)
!     mt part
            vop_wan2=0.0d0
            call zgemm('c','n',nvopt(k),num_orb_proj**2,n_pbmt,
     $        (1.d0,0.d0),qvopt(1,1,ind_k),n_pbtot,
     $        prod_wan2(1,1,1,ispin,ind_k,1),n_pbtot,
     &        (1.d0,0.d0),vop_wan2,nvopt_max)                        

!     it part
            s=0.0d0
            call s_overlap(0,s,nplw_gw,nplwgw(k),indpw_gw(1,iq))

            tmp1=0.0d0
            call zgemm('c','n',nvopt(k),nplwgw(k),nplwgw(k),(1.d0,0.d0),
     &        qvopt(n_pbmt+1,1,ind_k),n_pbtot,s,nplw_gw,
     &        (0.d0,0.d0),tmp1,nvopt_max)
            call zgemm('n','n',nvopt(k),num_orb_proj**2,nplwgw(k),
     &        (1.d0,0.d0),tmp1,nvopt_max,
     $        prod_wan2(n_pbmt+1,1,1,ispin,ind_k,1),
     $        n_pbtot,(1.d0,0.d0),vop_wan2,nvopt_max)

            vop_wan2_cor=0.0d0
            
            do ii=1, nvopt(k)
              do iorb=1,num_orb_cor
                do jorb=1,num_orb_cor
                  vop_wan2_cor(ii,iorb,jorb)
     $              =vop_wan2(ii,iorb,jorb)
                enddo
              enddo
            enddo

            tempmat2=0.0d0
            
            call zgemm('n','n',nvopt(k),num_orb_cor**2,num_orb_cor**2,
     $        (1.d0,0.d0),vop_wan2_cor,nvopt_max,tempmat1,
     $        num_orb_cor**2,(0.d0,0.d0),tempmat2,nvopt_max)

            fac=dcmplx(2.0d0,0.0d0)/dble(nspin*nrel)
            call zgemm('n','c',nvopt(k),nvopt(k),num_orb_cor**2,
     $        fac,tempmat2,nvopt_max,vop_wan2_cor,nvopt_max,
     $        (1.d0,0.d0),pola_tau_vop(1,1,ind_k),nvopt_max)
          enddo
        enddo


        do ind_k=1,ndim3_k(me3_k+1)
          k=n3_mpi_k(me3_k+1)+ind_k          
          call boson_pack_tau(1.0d0,pola_tau_vop(1,1,ind_k),nvopt_max,
     $      0.0d0,pola_tau_vop_packed(1,1,ind_tau,ind_k),nvopt_max,
     $      nvopt(k))
        enddo
      enddo
      end      
      

      

      
      

      


      
