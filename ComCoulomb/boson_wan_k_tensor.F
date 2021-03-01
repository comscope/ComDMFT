      subroutine boson_wan_k_tensor(w_s,matdim,ind_k,
     $  s_multiply1,s_multiply2,qvopt_multiply)
      
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

      integer :: ii,jj,k,n_pbt,num_qpoint,matsize,kst,iq,
     $  ispin,iorb,jspin,jorb,iorb1,iorb2,jorb1,jorb2,ind,ind_i,ind_j
      complex*16 :: v_tmp(n_pbtot,n_pbtot),tmp4(n_pbtot,nvopt_max),
     $  v_tmp2(n_pbtot,n_pbtot), s(nplw_gw,nplw_gw),
     $  tmp1(nplw_gw,nplw_gw),tmp2(nplw_gw,nplw_gw),
     $  tmp3(nplw_gw,nplw_gw)
      complex*16,allocatable ::v_tmp3(:,:),v_tmp4(:,:)

      complex*16::
     $  boson_k_tensor(num_orb_proj,num_orb_proj,
     $  num_orb_proj,num_orb_proj,nspin,nspin),


      
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

      boson_k_tensor=0.0d0

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

        boson_k_tensor=0.0d0

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
                    boson_k_tensor(iorb1,iorb2,jorb1,jorb2,ispin,jspin)
     $                =boson_k_tensor(iorb1,iorb2,jorb1,jorb2,
     $                ispin,jspin)
     $                +v_tmp4(ind_i,ind_j)/dble(num_qpoint) ! 1(r)* 2(r) 3(r')* 4(r')
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
!!! print here at each k point
!!!   for each atom and l value (li=l+1 for the system without spinorbit coupling)
c$$$  do m4=1,mudim
c$$$  km4=mt2ciw(m4,li,iatom)
c$$$  
c$$$  do m3=1,mudim
c$$$  km3=mt2ciw(m3,li,iatom)
c$$$  
c$$$  do m2=1,mudim
c$$$  km2=mt2ciw(m2,li,iatom)
c$$$  
c$$$  do m1=1,mudim
c$$$  km1=mt2ciw(m1,li,iatom)
c$$$  boson_k_tensor(km1,km2,km3,km4, ispin, jspin)=> tensor you want

      enddo
      
      deallocate(v_tmp3,v_tmp4)
      end      
