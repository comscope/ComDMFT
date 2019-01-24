      subroutine boson_wpb_k(w_s,matdim,ind_k,
     $  s_multiply1,s_multiply2,qvopt_multiply,
     $  boson_wpb)
      
!     need s_multiply
!     need qvopt_multiply
      use atom_mod
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
     $  boson_wpb(nwpb_max,natom_ciw,nwpb_max,natom_ciw,k_star_max)

      integer :: ii,jj,k,n_pbt,num_qpoint,matsize,kst,iq,
     $  ispin,iorb,jspin,jorb,iorb1,iorb2,jorb1,jorb2,ind,ind_i,ind_j,
     $  iat,jat
      complex*16 :: v_tmp(n_pbtot,n_pbtot),tmp4(n_pbtot,nvopt_max),
     $  v_tmp2(n_pbtot,n_pbtot), s(nplw_gw,nplw_gw),
     $  tmp1(nplw_gw,nplw_gw),tmp2(nplw_gw,nplw_gw),
     $  tmp3(nplw_gw,nplw_gw)
      complex*16,allocatable ::v_tmp3(:,:),v_tmp4(:,:)


      k=n3_mpi_k(me3_k+1)+ind_k

      v_tmp=0.0d0

      if (qvopt_multiply .eq. 1) then
        call zgemm('n','n',n_pbtot,matdim,matdim,(1.d0,0.d0),
     $    qvopt(1,1,ind_k),n_pbtot,w_s,matdim,(0.d0,0.d0),tmp4,
     $    n_pbtot)
        call zgemm('n','c',n_pbtot,n_pbtot,matdim,(1.d0,0.d0),tmp4,
     $    n_pbtot,qvopt(1,1,ind_k),n_pbtot,(0.d0,0.d0),v_tmp,n_pbtot)
      else
        do ii=1, matdim
          do jj=1, matdim
            v_tmp(ii,jj)=w_s(ii,jj)
          enddo
        enddo
      endif

      boson_wpb=0.0d0

c     ----------------------------------------------------------------------
      n_pbt=n_pbmt+nplwgw(k)

      if(.not.include_gamma_w) then
        num_qpoint=nqdiv-1
      else
        num_qpoint=nqdiv
      endif
      matsize=nwpb_max*natom_ciw

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
          tmp1=0.0d0
          tmp2=0.0d0
          tmp3=0.0d0
          do ii=1, nplwgw(k)
            do jj=1, nplwgw(k)
              tmp1(ii,jj)=v_tmp2(ii+n_pbmt,jj+n_pbmt)
            enddo
          enddo
          
          if (s_multiply2 .eq. 1) then
            call zgemm('n','n',nplwgw(k),nplwgw(k),nplwgw(k),
     $        (1.0d0,0.0d0),tmp1,nplw_gw,
     $        s,nplw_gw,(0.d0,0.d0),tmp2,nplw_gw)
          else
            tmp2=tmp1
          endif
          if (s_multiply1 .eq. 1) then
            call zgemm('n','n',nplwgw(k),nplwgw(k),nplwgw(k),
     $        (1.d0,0.d0),s,nplw_gw,tmp2,nplw_gw,
     $        (0.d0,0.d0),tmp3,nplw_gw)
          else
            tmp3=tmp2
          endif

          do ii=1, nplwgw(k)
            do jj=1, nplwgw(k)
              v_tmp2(ii+n_pbmt,jj+n_pbmt)=tmp3(ii,jj)
            enddo
          enddo
        endif
        
        call zgemm('n','n',n_pbt,matsize,n_pbt,
     $    (1.d0,0.d0),v_tmp2,n_pbtot,
     $    m_pi(1,1,1,ind_k,kst),n_pbtot,
     $    (0.d0,0.d0),v_tmp3,n_pbtot)
        call zgemm('c','n',matsize,matsize,n_pbt,
     $    (1.d0,0.d0),m_pi(1,1,1,ind_k,kst),n_pbtot,
     $    v_tmp3,n_pbtot,
     $    (0.d0,0.d0),v_tmp4,matsize)
        

        do iat=1, natom_ciw
          do iorb=1,nwpb(iat)
            call indexing(matsize,2,(/nwpb_max,natom_ciw/),1,ind_i,
     $        (/iorb,iat/))
            do jat=1, natom_ciw
              do jorb=1,nwpb(jat)
                call indexing(matsize,2,(/nwpb_max,natom_ciw/),1,ind_j,
     $            (/jorb,jat/))                
                boson_wpb(iorb,iat,jorb,jat,kst)
     $            =v_tmp4(ind_i,ind_j)/dble(num_qpoint)
              enddo
            enddo
          enddo
        enddo
      enddo
      
      deallocate(v_tmp3,v_tmp4)
      end


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
      complex*16 :: v_tmp(n_pbtot,n_pbtot),tmp4(n_pbtot,nvopt_max),
     $  v_tmp2(n_pbtot,n_pbtot), s(nplw_gw,nplw_gw),
     $  tmp1(nplw_gw,nplw_gw),tmp2(nplw_gw,nplw_gw),
     $  tmp3(nplw_gw,nplw_gw)
      complex*16,allocatable ::v_tmp3(:,:),v_tmp4(:,:)



      
      k=n3_mpi_k(me3_k+1)+ind_k

      v_tmp=0.0d0

      if (qvopt_multiply .eq. 1) then
        call zgemm('n','n',n_pbtot,matdim,matdim,(1.d0,0.d0),
     $    qvopt(1,1,ind_k),n_pbtot,w_s,matdim,(0.d0,0.d0),tmp4,
     $    n_pbtot)
        call zgemm('n','c',n_pbtot,n_pbtot,matdim,(1.d0,0.d0),tmp4,
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
          tmp1=0.0d0
          tmp2=0.0d0
          tmp3=0.0d0
          do ii=1, nplwgw(k)
            do jj=1, nplwgw(k)
              tmp1(ii,jj)=v_tmp2(ii+n_pbmt,jj+n_pbmt)
            enddo
          enddo
          
          if (s_multiply2 .eq. 1) then
            call zgemm('n','n',nplwgw(k),nplwgw(k),nplwgw(k),
     $        (1.0d0,0.0d0),tmp1,nplw_gw,
     $        s,nplw_gw,(0.d0,0.d0),tmp2,nplw_gw)
          else
            tmp2=tmp1
          endif
          if (s_multiply1 .eq. 1) then
            call zgemm('n','n',nplwgw(k),nplwgw(k),nplwgw(k),
     $        (1.d0,0.d0),s,nplw_gw,tmp2,nplw_gw,
     $        (0.d0,0.d0),tmp3,nplw_gw)
          else
            tmp3=tmp2
          endif

          do ii=1, nplwgw(k)
            do jj=1, nplwgw(k)
              v_tmp2(ii+n_pbmt,jj+n_pbmt)=tmp3(ii,jj)
            enddo
          enddo
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


      subroutine chi0_projection
      
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      use wannier90_mod

      implicit none
      include "mpif.h"

      integer :: ind_nu,i_nu,iat,iorb,jorb,iorb1,iorb2,jorb1,jorb2,
     $  ind1,ind2,jnd1,jnd2,ispin,jspin,ndim,iatom,isort,
     $  ind_i,ind_j,li,l,mudim,mudim_half,iwpb,
     $  m1,m2,m3,m4,km1,km2,km3,km4,ii

      complex*16 :: ff(0:3,0:n_nu)
      complex*16,allocatable :: 
     $  tempmat3(:,:),tempmat4(:,:,:,:),
     $  w_wpb(:,:),w_wpb_p(:,:),chi0_wpb_p(:,:),chi0_wpb(:,:),
     $  eps_wpb(:,:),epsinv_wpb(:,:),u_wpb(:,:),uloc_temp(:,:,:,:),
     $  uloc_temp2(:,:,:,:),
     $  rotmat(:,:),projection(:,:),chi0_wpb_nu(:,:,:)
      character(len=5) :: lchar,catom

      integer, external :: find_ind

      do iat=1, natom_ciw
        iatom=iatom_ciw(iat)
        isort=is(iatom)

        if (iatom .lt. 10)then
          write(catom,'(i1)') iatom
        elseif (iatom .lt. 10)then
          write(catom,'(i2)') iatom
        elseif (iatom .lt. 100)then
          write(catom,'(i3)') iatom
        endif        

        do l=0, lmb(isort)

          write(lchar,'(i2)') l          

          if ((corr_iter(l,iatom) .eq. 2)
     $      .and.
     $      (at_bonding(l,iatom) .eq. 0)
     $      ) then
            
            if(maswrk) then
              open(3,file=adjustl(trim(allfile))
     $          //'_chi0_low_Slater_'
     $          //trim(adjustl(txtel(isort)))//trim(adjustl(catom))
     $          //'_'//trim(adjustl(lchar))//'.dat')
              write(3,'(a)')
     $          '# Slater integral representation for U_matrix'
            endif



            allocate(tempmat3(n_ciw_max**2,nwpb(iat)))
            allocate(tempmat4(n_ciw_max,n_ciw_max,n_ciw_max,n_ciw_max))
            allocate(w_wpb(nwpb(iat),nwpb(iat)))
            allocate(w_wpb_p(nwpb(iat),nwpb(iat)))            
            allocate(chi0_wpb(nwpb(iat),nwpb(iat)))
            allocate(chi0_wpb_nu(nwpb(iat),nwpb(iat),ndim3_nu))            
            allocate(chi0_wpb_p(nwpb(iat),nwpb(iat)))            
            allocate(eps_wpb(nwpb(iat),nwpb(iat)))
            allocate(epsinv_wpb(nwpb(iat),nwpb(iat)))
            allocate(u_wpb(nwpb(iat),nwpb(iat)))

            call chi0_loc_correlated (l,iat,chi0_wpb_nu)

            
            if (nrel .eq. 1) then          
              li=l+1
              mudim=l+l+1

              
              allocate(rotmat(mudim,mudim))
              rotmat=0.0d0
              call cal_rotmat_cmplx2real(mudim, rotmat)
            else

              mudim=4*l+2
              mudim_half=mudim/2
              
              allocate(rotmat(mudim,mudim))
              rotmat=0.0d0
              
              call cal_rotmat_mlms2jmj(l, rotmat)
              
              allocate(uloc_temp2(mudim_half,mudim_half,
     $          mudim_half,mudim_half))
              uloc_temp2=0.0d0              
            endif

            allocate(uloc_temp(mudim,mudim,mudim,mudim))            
            
            allocate(projection(nwpb(iat),nwpb(iat)))
            
            
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
                  
                  w_wpb=0.0d0
                  w_wpb_p=0.0d0                  
                  
                  chi0_wpb_p=0.0d0
                  chi0_wpb=0.0d0

                  
                  
                  do iorb=1,nwpb(iat)
                    do jorb=1,nwpb(iat)
                      chi0_wpb(iorb,jorb)
     $                  =chi0_wpb_nu(iorb,jorb,ind_nu)
                    enddo
                  enddo

                  ndim=n_ciw_max*n_ciw_max                  
                  
                  call zgemm('n','n',ndim,nwpb(iat),nwpb(iat),
     $              (1.d0,0.d0),pi_ww(1,1,1,jspin,iat),ndim,chi0_wpb,
     $              nwpb(iat),(0.0d0,0.0d0),tempmat3,ndim)                  
                  
                  call zgemm('n','c',ndim,ndim,nwpb(iat),(1.d0,0.d0),
     $              tempmat3,ndim,pi_ww(1,1,1,ispin,iat),ndim,
     $              (0.0d0,0.0d0),tempmat4,ndim)                  
                  
                  uloc_temp=0.0d0
                  do m4=1,mudim
                    if (nrel .eq. 1) then
                      km4=find_ind(mt2ciw(m4,l+1,iatom),
     $                  n_ciw_atom(iatom),ciw_ind(1:n_ciw_atom(iatom),
     $                  iat))
                    else
                      km4=find_ind(mt2ciw_rel_l(m4,l+1,iatom),
     $                  n_ciw_atom(iatom),ciw_ind(1:n_ciw_atom(iatom),
     $                  iat))                      
                    endif
                    do m3=1,mudim
                      if (nrel .eq. 1) then
                        km3=find_ind(mt2ciw(m3,l+1,iatom),
     $                    n_ciw_atom(iatom),ciw_ind(1:n_ciw_atom(iatom),
     $                    iat))                        
                      else
                        km3=find_ind(mt2ciw_rel_l(m3,l+1,iatom),
     $                    n_ciw_atom(iatom),ciw_ind(1:n_ciw_atom(iatom),
     $                    iat))                                                
                      endif                      
                      do m2=1,mudim
                        if (nrel .eq. 1) then
                          km2=find_ind(mt2ciw(m2,l+1,iatom),
     $                      n_ciw_atom(iatom),
     $                      ciw_ind(1:n_ciw_atom(iatom),iat))                                                  
                        else
                          km2=find_ind(mt2ciw_rel_l(m2,l+1,iatom),
     $                      n_ciw_atom(iatom),
     $                      ciw_ind(1:n_ciw_atom(iatom),iat))                                                                          
                        endif                                              
                        do m1=1,mudim
                          if (nrel .eq. 1) then
                            km1=find_ind(mt2ciw(m1,l+1,iatom),
     $                        n_ciw_atom(iatom),
     $                        ciw_ind(1:n_ciw_atom(iatom),iat))                                                    
                          else
                            km1=find_ind(mt2ciw_rel_l(m1,l+1,iatom),
     $                        n_ciw_atom(iatom),
     $                        ciw_ind(1:n_ciw_atom(iatom),iat))                                                                            
                          endif                                              
                          uloc_temp(m1,m2,m3,m4)
     $                      =uloc_temp(m1,m2,m3,m4)
     $                      +tempmat4(km2,km1,km3,km4)
                        enddo
                      enddo
                    enddo
                  enddo
                  if (nrel .eq. 1) then
                    call rotate_umatrix(mudim,uloc_temp,
     $                rotmat,1)
                    call cal_slater_integral(uloc_temp,l,ff(0,i_nu))
                  else
                    call rotate_umatrix(mudim,uloc_temp,
     $                rotmat,1)

                    uloc_temp2=0.0d0
                    uloc_temp2
     $                =(
     $                uloc_temp(1:mudim_half,1:mudim_half,
     $                1:mudim_half,1:mudim_half)
     $                +uloc_temp((mudim_half+1):mudim,
     $                (mudim_half+1):mudim,1:mudim_half,
     $                1:mudim_half)
     $                +uloc_temp(1:mudim_half,1:mudim_half,
     $                (mudim_half+1):mudim,
     $                (mudim_half+1):mudim)              
     $                +uloc_temp((mudim_half+1):mudim,
     $                (mudim_half+1):mudim,
     $                (mudim_half+1):mudim,(mudim_half+1):mudim)
     $                )/4.0d0
                    
                    call cal_slater_integral(uloc_temp2,l,ff(0,i_nu))
                    
                  endif
                enddo
                
                if(nproc_tau/=1) then
                  call mpi_allreduce_dcmplx
     $              (ff,4*(n_nu+1),mpi_sum,comm_k)
                endif
                if(maswrk) then
                  do i_nu=0,n_nu
                    write(3,'(10(1x,e12.5))')w_nu(i_nu)*evolt/2,
     &                (dble(ff(ii,i_nu))/rydberg,ii=0,l)
                  enddo
                endif                  
              enddo
            enddo                
            deallocate(tempmat3)
            deallocate(tempmat4)
            deallocate(w_wpb)
            deallocate(w_wpb_p)            
            deallocate(chi0_wpb)
            deallocate(chi0_wpb_nu)            
            deallocate(chi0_wpb_p)            
            deallocate(eps_wpb)
            deallocate(epsinv_wpb)
            deallocate(u_wpb)
            deallocate(uloc_temp)
            deallocate(rotmat)
            deallocate(projection)
            if (nrel .eq. 2) then
              deallocate(uloc_temp2)              
            endif
            if (maswrk) close(3)
          endif
        enddo
      enddo
      end      
      


      subroutine u_weiss
      
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      use wannier90_mod

      implicit none
      include "mpif.h"

      integer :: ind_nu,i_nu,iat,iorb,jorb,iorb1,iorb2,jorb1,jorb2,
     $  ind1,ind2,jnd1,jnd2,ispin,jspin,ndim,iatom,isort,
     $  ind_i,ind_j,li,l,mudim,mudim_half,iwpb,
     $  m1,m2,m3,m4,km1,km2,km3,km4,ii

      complex*16 :: ff(0:3,0:n_nu)
      complex*16,allocatable :: 
     $  tempmat3(:,:),tempmat4(:,:,:,:),
     $  w_wpb(:,:),w_wpb_p(:,:),chi0_wpb_p(:,:),chi0_wpb(:,:,:),
     $  eps_wpb(:,:),epsinv_wpb(:,:),u_wpb(:,:),uloc_temp(:,:,:,:),
     $  uloc_temp2(:,:,:,:),
     $  rotmat(:,:),projection(:,:)
      character(len=5) :: lchar,catom

      integer, external :: find_ind

      do iat=1, natom_ciw
        iatom=iatom_ciw(iat)
        isort=is(iatom)

        if (iatom .lt. 10)then
          write(catom,'(i1)') iatom
        elseif (iatom .lt. 10)then
          write(catom,'(i2)') iatom
        elseif (iatom .lt. 100)then
          write(catom,'(i3)') iatom
        endif        

        do l=0, lmb(isort)

          write(lchar,'(i2)') l          

          if ((corr_iter(l,iatom) .eq. 2)
     $      .and.
     $      (at_bonding(l,iatom) .eq. 0)
     $      ) then
            
            if(maswrk) then
              open(3,file=adjustl(trim(allfile))
     $          //'_u_Slater_'
     $          //trim(adjustl(txtel(isort)))//trim(adjustl(catom))
     $          //'_'//trim(adjustl(lchar))//'.dat')
              write(3,'(a)')
     $          '# Slater integral representation for U_matrix'
            endif

            allocate(tempmat3(n_ciw_max**2,nwpb(iat)))
            allocate(tempmat4(n_ciw_max,n_ciw_max,n_ciw_max,n_ciw_max))
            allocate(w_wpb(nwpb(iat),nwpb(iat)))
            allocate(w_wpb_p(nwpb(iat),nwpb(iat)))            
            allocate(chi0_wpb(nwpb(iat),nwpb(iat),ndim3_nu))
            allocate(chi0_wpb_p(nwpb(iat),nwpb(iat)))            
            allocate(eps_wpb(nwpb(iat),nwpb(iat)))
            allocate(epsinv_wpb(nwpb(iat),nwpb(iat)))
            allocate(u_wpb(nwpb(iat),nwpb(iat)))

            call chi0_loc_correlated (l,iat,chi0_wpb)            

            
            if (nrel .eq. 1) then          
              li=l+1
              mudim=l+l+1

              
              allocate(rotmat(mudim,mudim))
              rotmat=0.0d0
              call cal_rotmat_cmplx2real(mudim, rotmat)
            else

              mudim=4*l+2
              mudim_half=mudim/2
              
              allocate(rotmat(mudim,mudim))
              rotmat=0.0d0
              
              call cal_rotmat_mlms2jmj(l, rotmat)
              
              allocate(uloc_temp2(mudim_half,mudim_half,
     $          mudim_half,mudim_half))
              uloc_temp2=0.0d0              
            endif

            allocate(uloc_temp(mudim,mudim,mudim,mudim))            
            
            allocate(projection(nwpb(iat),nwpb(iat)))
            
            call projection_coulomb(iat,l,projection)

            
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
                  
                  w_wpb=0.0d0
                  w_wpb_p=0.0d0                  
                  
                  do iorb=1,nwpb(iat)
                    do jorb=1,nwpb(iat)
                      w_wpb(iorb,jorb)
     $                  =w_wpb_local(iorb,iat,jorb,iat,ind_nu)
                    enddo
                  enddo
                  
                  call zgemm('n','n',nwpb(iat),nwpb(iat),nwpb(iat),
     $              (1.0d0,0.0d0),w_wpb,nwpb(iat),projection,
     $              nwpb(iat),(0.0d0,0.0d0),w_wpb_p,nwpb(iat))
                  
                  
                  chi0_wpb_p=0.0d0
                  
                  call zgemm('n','c',nwpb(iat),nwpb(iat),nwpb(iat),
     $              (1.0d0,0.0d0),chi0_wpb(1,1,ind_nu),nwpb(iat),
     $              projection,nwpb(iat),(0.0d0,0.0d0),
     $              chi0_wpb_p,nwpb(iat))
                  
                  eps_wpb=0.0d0
                  do iwpb=1, nwpb(iat)
                    eps_wpb(iwpb,iwpb)=1.0d0
                  enddo
                  
                  call zgemm('n','n',nwpb(iat),nwpb(iat),nwpb(iat),
     $              (1.0d0,0.0d0),w_wpb_p,nwpb(iat),chi0_wpb_p,
     $              nwpb(iat),(1.0d0,0.0d0),eps_wpb,nwpb(iat))
                  
                  epsinv_wpb=0.0d0
                  
                  call dcmplx_matinv(eps_wpb, epsinv_wpb,nwpb(iat),
     $              nwpb(iat))
                  
                  u_wpb=0.0d0
                  call zgemm('n','n',nwpb(iat),nwpb(iat),nwpb(iat),
     $              (1.0d0,0.0d0),epsinv_wpb,nwpb(iat),w_wpb,
     $              nwpb(iat),(0.0d0,0.0d0),u_wpb,nwpb(iat))
                  
                  ndim=n_ciw_max*n_ciw_max

                  call zgemm('n','n',ndim,nwpb(iat),nwpb(iat),
     $              (1.d0,0.d0),pi_ww(1,1,1,jspin,iat),ndim,u_wpb,
     $              nwpb(iat),(0.0d0,0.0d0),tempmat3,ndim)                  
                  
                  call zgemm('n','c',ndim,ndim,nwpb(iat),(1.d0,0.d0),
     $              tempmat3,ndim,pi_ww(1,1,1,ispin,iat),ndim,
     $              (0.0d0,0.0d0),tempmat4,ndim)                  
                  
                  uloc_temp=0.0d0
                  do m4=1,mudim
                    if (nrel .eq. 1) then
                      km4=find_ind(mt2ciw(m4,l+1,iatom),
     $                  n_ciw_atom(iatom),ciw_ind(1:n_ciw_atom(iatom),
     $                  iat))
                    else
                      km4=find_ind(mt2ciw_rel_l(m4,l+1,iatom),
     $                  n_ciw_atom(iatom),ciw_ind(1:n_ciw_atom(iatom),
     $                  iat))                      
                    endif
                    do m3=1,mudim
                      if (nrel .eq. 1) then
                        km3=find_ind(mt2ciw(m3,l+1,iatom),
     $                    n_ciw_atom(iatom),ciw_ind(1:n_ciw_atom(iatom),
     $                    iat))                        
                      else
                        km3=find_ind(mt2ciw_rel_l(m3,l+1,iatom),
     $                    n_ciw_atom(iatom),ciw_ind(1:n_ciw_atom(iatom),
     $                    iat))                                                
                      endif                      
                      do m2=1,mudim
                        if (nrel .eq. 1) then
                          km2=find_ind(mt2ciw(m2,l+1,iatom),
     $                      n_ciw_atom(iatom),
     $                      ciw_ind(1:n_ciw_atom(iatom),iat))                                                  
                        else
                          km2=find_ind(mt2ciw_rel_l(m2,l+1,iatom),
     $                      n_ciw_atom(iatom),
     $                      ciw_ind(1:n_ciw_atom(iatom),iat))                                                                          
                        endif                                              
                        do m1=1,mudim
                          if (nrel .eq. 1) then
                            km1=find_ind(mt2ciw(m1,l+1,iatom),
     $                        n_ciw_atom(iatom),
     $                        ciw_ind(1:n_ciw_atom(iatom),iat))                                                    
                          else
                            km1=find_ind(mt2ciw_rel_l(m1,l+1,iatom),
     $                        n_ciw_atom(iatom),
     $                        ciw_ind(1:n_ciw_atom(iatom),iat))                                                                            
                          endif                                              
                          uloc_temp(m1,m2,m3,m4)
     $                      =uloc_temp(m1,m2,m3,m4)
     $                      +tempmat4(km2,km1,km3,km4)
                        enddo
                      enddo
                    enddo
                  enddo
                  if (nrel .eq. 1) then
                    call rotate_umatrix(mudim,uloc_temp,
     $                rotmat,1)
                    call cal_slater_integral(uloc_temp,l,ff(0,i_nu))
                  else
                    call rotate_umatrix(mudim,uloc_temp,
     $                rotmat,1)

                    uloc_temp2=0.0d0
                    uloc_temp2
     $                =(
     $                uloc_temp(1:mudim_half,1:mudim_half,
     $                1:mudim_half,1:mudim_half)
     $                +uloc_temp((mudim_half+1):mudim,
     $                (mudim_half+1):mudim,1:mudim_half,
     $                1:mudim_half)
     $                +uloc_temp(1:mudim_half,1:mudim_half,
     $                (mudim_half+1):mudim,
     $                (mudim_half+1):mudim)              
     $                +uloc_temp((mudim_half+1):mudim,
     $                (mudim_half+1):mudim,
     $                (mudim_half+1):mudim,(mudim_half+1):mudim)
     $                )/4.0d0
                    
                    call cal_slater_integral(uloc_temp2,l,ff(0,i_nu))
                    
                  endif
                enddo
                
                if(nproc_tau/=1) then
                  call mpi_allreduce_dcmplx
     $              (ff,4*(n_nu+1),mpi_sum,comm_k)
                endif
                if(maswrk) then
                  do i_nu=0,n_nu
                    write(3,'(10(1x,e12.5))')w_nu(i_nu)*evolt/2,
     &                (dble(ff(ii,i_nu))*rydberg,ii=0,l)
                  enddo
                endif                  
              enddo
            enddo                
            deallocate(tempmat3)
            deallocate(tempmat4)
            deallocate(w_wpb)
            deallocate(w_wpb_p)            
            deallocate(chi0_wpb)
            deallocate(chi0_wpb_p)            
            deallocate(eps_wpb)
            deallocate(epsinv_wpb)
            deallocate(u_wpb)
            deallocate(uloc_temp)
            deallocate(rotmat)
            deallocate(projection)
            if (nrel .eq. 2) then
              deallocate(uloc_temp2)              
            endif
            if (maswrk) close(3)
          endif
        enddo
      enddo
      end


      
c$$$      subroutine projection_test
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
c$$$     $  m1,m2,m3,m4,km1,km2,km3,km4,ii
c$$$
c$$$      complex*16 :: ff(0:3,0:n_nu)
c$$$      complex*16,allocatable :: tempmat1(:,:),
c$$$     $  tempmat2(:,:),tempmat3(:,:),tempmat4(:,:,:,:),
c$$$     $  w_wpb(:,:), chi0_wpb(:,:),eps_wpb(:,:),
c$$$     $  epsinv_wpb(:,:),u_wpb(:,:),uloc_temp(:,:,:,:),
c$$$     $  uloc_temp2(:,:,:,:),
c$$$     $  rotmat(:,:),projection(:,:)
c$$$      character(len=5) :: lchar,catom
c$$$
c$$$      integer, external :: find_ind
c$$$
c$$$      do iat=1, natom_ciw
c$$$        iatom=iatom_ciw(iat)
c$$$        isort=is(iatom)
c$$$        write(iun,*) 'projection_test', iat
c$$$
c$$$        call projection_wpb(iat)
c$$$
c$$$      enddo
c$$$      end



      subroutine projection_coulomb(iat,l,proj_coulomb)
      use wannier90_mod
      use manager_mod
      use parallel_mod
      implicit none
      
      integer, intent(in) :: iat
      integer, intent(in) :: l
      complex*16, intent(out) :: proj_coulomb(nwpb(iat),nwpb(iat))
      
      integer :: iorb,jorb,korb,lorb,iwpb,jwpb,m1,m2,ispin,jspin,iatom,
     $  ndim
c$$$      complex*16 :: proj_coulomb(nwpb(iat),nwpb(iat))
      complex*16, allocatable ::
     $  pi_ww_at_dual(:,:,:,:),
     $  pi_ww_at(:,:,:,:)      
c$$$     $  identity_test(:,:,:,:,:,:),
c$$$     $  resolution(:,:)
      

      integer, external :: find_ind

      iatom=iatom_ciw(iat)
      

      allocate(pi_ww_at(nwpb(iatom),
     $  n_ciw_atom(iatom),n_ciw_atom(iatom),nspin))
      pi_ww_at=0.0d0
      allocate(pi_ww_at_dual(nwpb(iatom),
     $  n_ciw_atom(iatom),n_ciw_atom(iatom),nspin))
      pi_ww_at_dual=0.0d0      
      

      
c$$$      allocate(identity_test(n_ciw_atom(iatom),n_ciw_atom(iatom),nspin,
c$$$     $  n_ciw_atom(iatom),n_ciw_atom(iatom),nspin))
c$$$      allocate(resolution(nwpb(iat),nwpb(iat)))
      
       
c$$$  projection=0.0d0
      pi_ww_at=0.0d0      
      pi_ww_at_dual=0.0d0
c$$$      identity_test=0.0d0
c$$$  resolution=0.0d0
      proj_coulomb=0.0d0

      do ispin=1, nspin
        do iorb=1, n_ciw_atom(iatom)
          do jorb=1, n_ciw_atom(iatom)
            do iwpb=1, nwpb(iat)
              pi_ww_at(iwpb,iorb,jorb,ispin)
     $          =dconjg(pi_ww(iorb,jorb,iwpb,ispin,iat))
            enddo
          enddo
        enddo
      enddo
      ndim=n_ciw_atom(iatom)**2*nspin

      if (maswrk) then
        write(iun,'(a, 2i5)') 'dual vector for iatom, l: ', iat, l
      endif
      call dual_vector(nwpb(iat),ndim,pi_ww_at,pi_ww_at_dual)
              
c$$$      call zgemm('c','n',n_ciw_atom(iatom)**2*nspin,
c$$$     $  n_ciw_atom(iatom)**2*nspin,
c$$$     $  nwpb(iat),(1.0d0,0.0d0),pi_ww_at,nwpb(iat),
c$$$     $  pi_ww_at_dual,nwpb(iat),(0.0d0,0.0d0),
c$$$     $  identity_test,n_ciw_atom(iatom)**2*nspin)

c$$$      call zgemm('n','c',nwpb(iat),nwpb(iat),n_ciw_atom(iatom)**2*nspin,
c$$$     $  (1.0d0,0.0d0),pi_ww_at,nwpb(iat),
c$$$     $  pi_ww_at_dual,nwpb(iat),(0.0d0,0.0d0),
c$$$     $  resolution,nwpb(iat))
c$$$
c$$$!     resolution test
c$$$      do iwpb=1, nwpb(iat)
c$$$        do jwpb=1, nwpb(iat)
c$$$          
c$$$          if (iwpb .eq. jwpb) then
c$$$            if (cdabs(resolution(iwpb,jwpb)-1.0d0) .gt. 1.0d-6) then
c$$$              write(iun,'(a, 2i5, 2f12.6)')
c$$$     $          'resolution matrix fail', iwpb,jwpb,
c$$$     $          resolution(iwpb,jwpb)
c$$$            endif
c$$$          else
c$$$            if (cdabs(resolution(iwpb,jwpb)) .gt. 1.0d-6) then            
c$$$              write(iun,'(a, 2i5, 2f12.6)')
c$$$     $          'resolution matrix fail', iwpb,jwpb,
c$$$     $          resolution(iwpb,jwpb)
c$$$            endif
c$$$          endif
c$$$        enddo
c$$$      enddo
c$$$
c$$$
c$$$      
!     identity operator
c$$$
c$$$      do ispin=1, nspin
c$$$        do iorb=1, n_ciw_atom(iatom)
c$$$          do jorb=1, n_ciw_atom(iatom)
c$$$            do jspin=1, nspin
c$$$              do korb=1, n_ciw_atom(iatom)
c$$$                do lorb=1, n_ciw_atom(iatom)
c$$$                  
c$$$                  if ((ispin .eq. jspin) .and.
c$$$     $              (iorb .eq. korb) .and.
c$$$     $              (jorb .eq. lorb)) then
c$$$                    if (cdabs(identity_test
c$$$     $                (iorb,jorb,ispin,korb,lorb,jspin)-1.0d0)
c$$$     $                .gt. 1.0d-6) then
c$$$                      write(iun,'(a, 6i5, 2f12.6)')
c$$$     $                  'identity_test fail',
c$$$     $                  iorb,jorb,ispin,korb,lorb,jspin,
c$$$     $                  identity_test(iorb,jorb,ispin,korb,lorb,jspin)
c$$$                    endif
c$$$                  else
c$$$                    if (cdabs(identity_test
c$$$     $                (iorb,jorb,ispin,korb,lorb,jspin))
c$$$     $                .gt. 1.0d-6) then
c$$$                      write(iun,'(a, 6i5, 2f12.6)')
c$$$     $                  'identity_test fail',
c$$$     $                  iorb,jorb,ispin,korb,lorb,jspin,
c$$$     $                  identity_test(iorb,jorb,ispin,korb,lorb,jspin)
c$$$                    endif                      
c$$$                  endif
c$$$                enddo
c$$$              enddo
c$$$            enddo
c$$$          enddo
c$$$        enddo
c$$$      enddo      


      do iwpb=1, nwpb(iat)
        do jwpb=1, nwpb(iat)

          do ispin=1, nspin
            if (nrel .eq. 1) then            
              do m1=1, 2*l+1
                
                iorb=find_ind
     $            (mt2ciw(m1,l+1,iatom),n_ciw_atom(iatom),
     $            ciw_ind(1:n_ciw_atom(iatom),iat))
                
                do m2=1, 2*l+1
                  jorb=find_ind
     $              (mt2ciw(m2,l+1,iatom),n_ciw_atom(iatom),
     $              ciw_ind(1:n_ciw_atom(iatom),iat))
                  
                  proj_coulomb(iwpb,jwpb)
     $              = proj_coulomb(iwpb,jwpb)
     $              +pi_ww_at(iwpb,iorb,jorb,ispin)
     $              *dconjg(pi_ww_at_dual(jwpb,iorb,jorb,ispin))
                enddo
              enddo
            else
              do m1=1,4*l+2
                
                iorb=find_ind
     $            (mt2ciw_rel_l(m1,l+1,iatom),n_ciw_atom(iatom),
     $            ciw_ind(1:n_ciw_atom(iatom),iat))
                
                do m2=1,4*l+2
                  jorb=find_ind
     $              (mt2ciw_rel_l(m2,l+1,iatom),n_ciw_atom(iatom),
     $              ciw_ind(1:n_ciw_atom(iatom),iat))

                  proj_coulomb(iwpb,jwpb)
     $              = proj_coulomb(iwpb,jwpb)
     $              +pi_ww_at(iwpb,iorb,jorb,ispin)
     $              *dconjg(pi_ww_at_dual(jwpb,iorb,jorb,ispin))                  
                enddo
              enddo
            endif
          enddo

c$$$          if (iwpb .eq. jwpb) then
c$$$            if (cdabs(proj_coulomb(iwpb,jwpb)-1.0d0) .gt. 1.0d-6) then
c$$$              write(iun,'(a, 2i5, 2f12.6)')
c$$$     $          'proj_coulomb matrix fail', iwpb,jwpb,
c$$$     $          proj_coulomb(iwpb,jwpb)
c$$$            endif
c$$$          else
c$$$            if (cdabs(proj_coulomb(iwpb,jwpb)) .gt. 1.0d-6) then            
c$$$              write(iun,'(a, 2i5, 2f12.6)')
c$$$     $          'proj_coulomb matrix fail', iwpb,jwpb,
c$$$     $          proj_coulomb(iwpb,jwpb)
c$$$            endif
c$$$          endif          
          
        enddo
      enddo

      deallocate(pi_ww_at)      
      deallocate(pi_ww_at_dual)
c$$$      deallocate(identity_test)
c$$$      deallocate(resolution)
        
      end



      subroutine dual_vector(dim1,dim2,matin,matout)
      use manager_mod
      use parallel_mod
      implicit none
      integer, intent(in) :: dim1,dim2
      complex*16, intent(in) :: matin(dim1,dim2)
      complex*16, intent(out) :: matout(dim1,dim2)
      
      integer :: info,lwork,maxdim,mindim,ii
      complex*16 :: matin_temp(dim1,dim2),uu(dim1,dim1),vt(dim2,dim2),
     $  tempmat(dim1,dim2),ss_inv_mat(dim1,dim2)
      double precision, allocatable :: ss(:),rwork(:)
      complex*16,allocatable :: work(:),worktemp(:)
      

      maxdim=max(dim1,dim2)
      mindim=min(dim1,dim2)      
      allocate(ss(mindim))
      allocate(rwork(5*mindim))

      matin_temp=matin

      allocate(worktemp(3*maxdim))

      call zgesvd('A','A',dim1,dim2,matin_temp,dim1,
     $  ss,uu,dim1,vt,dim2,worktemp,-1,rwork,info)
      
      lwork=max(1, int(worktemp(1)))
      allocate(work(lwork))
      
      call zgesvd('A','A',dim1,dim2,matin_temp,dim1,
     $  ss,uu,dim1,vt,dim2,work,lwork,rwork,info)
      deallocate(work)
      deallocate(worktemp)      

      deallocate(rwork)

      if (maswrk) then
        write(iun,'(a, 10000f20.12)')
     $    'singular value product basis', ss
      endif

      ss_inv_mat=0.0d0
      
      do ii=1, mindim
        ss_inv_mat(ii,ii)=1.0d0/ss(ii)
      enddo
      
      deallocate(ss)
      
      call zgemm('n','n',dim1,dim2,dim1,
     $  (1.0d0,0.0d0),uu,dim1,ss_inv_mat,dim1,
     $  (0.0d0,0.0d0),tempmat,dim1)
      call zgemm('n','n',dim1,dim2,dim2,
     $  (1.0d0,0.0d0),tempmat,dim1,vt,dim2,
     $  (0.0d0,0.0d0),matout,dim1)      

      end      


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

      ! g_loc

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
c$$$                  write(179,'(9i5)')
c$$$     $              ispin,ind_k,ik,ibnd,iwan1,iwan2,jwan1,jwan2,itau
c$$$c$$$                  write(179,*) g_tau(iwan1,jwan1,itau,ispin) ! good
c$$$                  write(179,*) green0(de,tau_mesh(itau))                              
c$$$c$$$                  write(179,*) proj_renorm(ibnd,iwan2,ind_k,ispin) ! good
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
     $             =pola_tau(tau4,tau1,tau3,tau2)                  
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
      

      

      
      

        


      
