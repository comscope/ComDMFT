      subroutine v_matr
      use atom_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: igg,nd0(3),ia,ib,ic,nqdiv0,k,ind_k,numpw
      double precision :: alfa,w0,qg(3),qg2,pi2a,p(3)
      complex*16, allocatable :: str_con(:,:)


      call timel('****    start v_matr    ************')
      pi2a=(pi+pi)/par
      alfa=3.d0
c     ---------------------------------------------------------------------------
      nd0=ndiv                  !*2
      nqdiv0=nd0(1)*nd0(2)*nd0(3)
      w0=1.d0/dfloat(nqdiv0)
      v_gg_00=amega/4.d0/pi/sqrt(alfa*pi)/w0 !+alfa
      do ic=0,nd0(3)-1
        do ib=0,nd0(2)-1
          do ia=0,nd0(1)-1
            p=dfloat(ia)*gbas(:,1)/nd0(1)+dfloat(ib)*gbas(:,2)/nd0(2)
     &        +dfloat(ic)*gbas(:,3)/nd0(3)
            do igg=1,numplw
              qg=(gbs(:,igg)+p)*pi2a
              qg2=dot_product(qg,qg)
              if(qg2>1.d-9) v_gg_00=v_gg_00-exp(-alfa*qg2)/qg2
            enddo
          enddo
        enddo
      enddo
c     ------- Calculate V_MM V_II for k=1,npnt ----------------------------
      v_mm=(0.d0,0.d0) ! nonlocal part of MT-MT coulmb matrix
      call timel('****    start set_lattices    ******')
      call set_lattices
      call timel('****   finish set_lattices    ******')
      allocate(str_con(nlb2,nlb2))
      do ind_k=1,ndim3_k(me3_k+1)
        k=n3_mpi_k(me3_k+1)+ind_k
c     ----- generate STR_CON ----------------------------------------------
        call timel('**** start get_str_con, ind_k   ****')
        call get_str_con(str_con,pnt(1,k))
c$$$        write(179,*) sum(str_con)
        call timel('**** finish get_str_con, ind_k  ****')
c     -------- Mt-Mt -----------------------------------------------------
        v_mm(:,:,ind_k)=4.d0*pi*str_con
          call timel('** start v_coul_mi_old, ind_k   ****')
          call v_coul_mi_old(ind_k)
          call timel('** finish v_coul_mi_old, ind_k  ****')
c     --------- Int-Int ---------------------------------------------------
        call timel('** start v_coul_ii_old, ind_k   ****')
        call v_coul_ii_gg(ind_k,numpw)
        call timel('** finish v_coul_ii_old, ind_k  ****')
      enddo                     !! over k
      call timel('***    start v_spectrum         ****')
      call v_spectrum
      call timel('***   finish v_spectrum         ****')
c     --------------------------------------------------------------------
      deallocate(str_con)
c     ------- Calculate S_TT_R0 for R=1,npnt ------------------------------
      call timel('******    start v_mm_r_from_k   ****')
      call v_mm_r_from_k
      call timel('*******  V_MATR finished ***********')
      end
