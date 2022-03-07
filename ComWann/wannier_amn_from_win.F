      subroutine wannier_orbital(ispin)
      use comwann_mod
      implicit none
      include "mpif.h"            
      integer, intent(in) :: ispin
      integer ::jj,ii, 
     $  isort,l,iival,jjval,iatom,
     $  li,ist,mdim,lval,i_index,
     $  is1,is2,mm,
     $  ind_k,k,k0,ibnd,ind0,llval,mjval,lival,lldim,ie1,inn1,in1,
     $  radind1,m1,lm1,km1,mjind,jjdim,irk,iomega,i_len,mt0,ir
      double precision :: weight, occ,
     $  llval_dble,jjval_dble,omega
c$$$  $  ,pdos(-500:500,0:maxb,nrel,nsort)
      complex*16 :: znew(nfun,num_bands),
     $  proj(0:maxb,nrel,2*maxb+2,natom,num_bands,ndim_kk(me+1))

      complex*16, allocatable :: wfn_mt0(:,:), wfn_mt1(:,:),
     $  rot_matrix(:,:)

!     radial function determination
c$$$  dimension determination`

      if (ubi_num==4) then      
        call t_t1_x_all
      endif
      
c$$$      weight_in_froz_win=0.0d0


      
      if (irel.le.1) then
        do iatom=1, natom
          isort=is(iatom)
          do l=0,lmb(isort)     ! l
            li=l+1              !li
            call cal_overlapmax_bandprojection_to_mtorb_low_rmode
     $        (iatom,l,iival,jjval,li,ispin,
     $        dim_radial(l,1,iatom),coeff_radial(1,l,1,iatom),
     $        weight,occ)
            weight_in_froz_win(l,1,isort)
     $        =weight_in_froz_win(l,1,isort)
     $        +weight/sort_degen(isort)
            occ_iwan(l,1,isort)
     $        =occ_iwan(l,1,isort)
     $        +occ/sort_degen(isort)            
          enddo
        enddo
      else
        do iatom=1, natom 
          isort=is(iatom)
          li=0
          do l=0,lmb(isort)     ! l
            ist=-1
            if(l.eq.0) ist=1
            do iival=ist,1,2    ! i
              i_index=(iival+1)/2+1
              li=li+1           !li
              jjval=l+l+iival   ! j
              call cal_overlapmax_bandprojection_to_mtorb_low_rmode
     $          (iatom,l,iival,jjval,li,ispin,
     $          dim_radial(l,i_index,iatom),
     $          coeff_radial(1,l,i_index,iatom),weight,occ)
              weight_in_froz_win(l,i_index,isort)
     $          =weight_in_froz_win(l,i_index,isort)
     $          +weight/sort_degen(isort)
              occ_iwan(l,i_index,isort)
     $          =occ_iwan(l,i_index,isort)
     $          +occ/sort_degen(isort)              
            enddo
          enddo
        enddo
      endif


c$$$      if (me .eq. 0) then
c$$$        i_len=len_trim(allfile)        
c$$$        open(3,file=allfile(1:i_len)//'_v_mt_h_'
c$$$     $    //trim(ubi_0)//'_2.txt')
c$$$        do isort=1,nsort
c$$$          mt0=indmt(1,isort,1)          
c$$$          write(3,'(i5, 10000(f20.12,2x))'), isort,
c$$$     $      (v_mt_h(mt0+ir), ir=0, nrad(isort))
c$$$        enddo
c$$$        close(3)
c$$$c$$$      elseif (me .eq. 5) then
c$$$c$$$        open(3,file=allfile(1:i_len)//'_v_mt_h_'
c$$$c$$$     $    //trim(ubi_0)//'_2.txt')
c$$$c$$$        do isort=1,nsort
c$$$c$$$          mt0=indmt(1,isort,1)          
c$$$c$$$          write(3,'(i5, 10000(f20.12,2x))'), isort,
c$$$c$$$     $      (v_mt_h(mt0+ir), ir=0, nrad(isort))
c$$$c$$$        enddo
c$$$c$$$        close(3)
c$$$      endif      

      call print_radial_function(ispin)


c$$$      if (me .eq. 0) then
c$$$        i_len=len_trim(allfile)        
c$$$        open(3,file=allfile(1:i_len)//'_v_mt_h_'
c$$$     $    //trim(ubi_0)//'_3.txt')
c$$$        do isort=1,nsort
c$$$          mt0=indmt(1,isort,1)          
c$$$          write(3,'(i5, 10000(f20.12,2x))'), isort,
c$$$     $      (v_mt_h(mt0+ir), ir=0, nrad(isort))
c$$$        enddo
c$$$        close(3)
c$$$c$$$      elseif (me .eq. 5) then
c$$$c$$$        open(3,file=allfile(1:i_len)//'_v_mt_h_'
c$$$c$$$     $    //trim(ubi_0)//'_2.txt')
c$$$c$$$        do isort=1,nsort
c$$$c$$$          mt0=indmt(1,isort,1)          
c$$$c$$$          write(3,'(i5, 10000(f20.12,2x))'), isort,
c$$$c$$$     $      (v_mt_h(mt0+ir), ir=0, nrad(isort))
c$$$c$$$        enddo
c$$$c$$$        close(3)
c$$$      endif      
      

c$$$  if (maswrk) then
c$$$  write(iun,*) 'dim_radial', dim_radial
c$$$  endif
      
!     calculate levels
      
      proj=0.0d0
      do ind_k=1,ndim_kk(me+1)  ! k vector
        k=n_mpi_kk(me+1)+ind_k
        k0=i_kref(k)
        
        call sym_z_0(znew,k,z_wan_bnd(1,1,k0),
     &    num_bands,k_group(k),pnt(1,k))
        
        do ibnd=1,num_bands

          do iatom=1, natom
            isort=is(iatom)            
            do llval=0,lmb(isort) ! l          
              
              ind0=io_lem(iatom)-1
              
              if (irel .le. 1) then
                do mjval=-llval, llval

                  llval_dble=dble(llval)                  
                  lival=llval+1
                  lldim=2*llval+1

                  allocate(rot_matrix(lldim,lldim))
                  rot_matrix=0.0d0
                  
                  call rotmat_wan_proj(atom_xaxis(1,iatom),
     $              atom_zaxis(1,iatom),llval_dble,lldim,1,
     $              rot_matrix)

                  allocate(wfn_mt0(lldim,dim_radial(llval,1,iatom)))
                  wfn_mt0=0.0d0
                  allocate(wfn_mt1(lldim,dim_radial(llval,1,iatom)))
                  wfn_mt1=0.0d0                  
                  radind1=0
                  do ie1=1,ntle(llval,isort)
                    in1=1
                    if(augm(ie1,llval,isort)/='LOC') in1=2
                    do inn1=1,in1
                      radind1=radind1+1
                      do m1=-llval, llval
                        lm1=llval*(llval+1)+m1+1
                        km1=ind0+indbasa(inn1,ie1,lm1,isort)
                        mjind=m1+llval+1
                        wfn_mt0(mjind,radind1)
     $                    =wfn_mt0(mjind,radind1)
     $                    +znew(km1,ibnd)
                      enddo                    
                    enddo
                  enddo
                  wfn_mt1=
     $              matmul(rot_matrix,wfn_mt0)
                  mjind=mjval+llval+1
                  do radind1=1, dim_radial(llval,1,iatom)
                    proj(llval,1,mjind,iatom,ibnd,ind_k)
     $                =proj(llval,1,mjind,iatom,ibnd,ind_k)
     $                +dconjg(coeff_radial(radind1,llval,1,iatom))
     $                *wfn_mt1(mjind,radind1)
                  enddo
                  deallocate(rot_matrix)
                  deallocate(wfn_mt0)
                  deallocate(wfn_mt1)
                enddo
              else
                do iival=-1,1, 2
                  if ((llval .eq. 0) .and. (iival .eq. -1)) cycle
                  i_index=(iival+1)/2+1
                  jjval=llval*2+iival
                  jjval_dble=llval+iival/2.0d0
                  jjdim=nint(2*jjval_dble+1)
                  lival=llval+llval+(iival+1)/2
                  do mjval=-jjval, jjval, 2
                    
                    allocate(rot_matrix(jjdim,jjdim))
                    rot_matrix=0.0d0
                    allocate(wfn_mt0(jjdim,dim_radial
     $                (llval,i_index,iatom)))
                    wfn_mt0=0.0d0
                    allocate(wfn_mt1(jjdim,dim_radial
     $                (llval,i_index,iatom)))
                    wfn_mt1=0.0d0                                  
                    
                    call rotmat_wan_proj(atom_xaxis(1,iatom),
     $                atom_zaxis(1,iatom),jjval_dble,jjdim,0,
     $                rot_matrix)
                    
                    radind1=0
                    do ie1=1,ntle(llval,isort)
                      in1=1
                      if(augm(ie1,llval,isort)/='LOC') in1=2
                      do inn1=1,in1
                        radind1=radind1+1
                        do m1=-jjval, jjval, 2
                          call getlimj(lm1,llval,iival,m1,lival,1)                      
                          km1=ind0+indbasa(inn1,ie1,lm1,isort)
                          mjind=(m1+jjval+2)/2
c$$$  write(179,'(6i5)') llval, iival,
c$$$  $                      mjind, radind1,jjdim,dim_radial
c$$$  $                      (llval,i_index,iatom)
                          wfn_mt0(mjind,radind1)
     $                      =wfn_mt0(mjind,radind1)
     $                      +znew(km1,ibnd)
                        enddo                    
                      enddo
                    enddo
                    wfn_mt1=
     $                matmul(rot_matrix,wfn_mt0)
                    
                    mjind=(mjval+jjval+2)/2
                    do radind1=1, dim_radial(llval,i_index,iatom)
                      proj(llval,i_index,mjind,iatom,ibnd,ind_k)
     $                  =proj(llval,i_index,mjind,iatom,ibnd,ind_k)
     $                  +dconjg(coeff_radial(radind1,llval,
     $                  i_index,iatom))
     $                  *wfn_mt1(mjind,radind1)
                    enddo
                    
                    deallocate(rot_matrix)
                    deallocate(wfn_mt0)
                    deallocate(wfn_mt1)

                  enddo
                enddo
              endif
            enddo
          enddo
        enddo
      enddo      


c$$$  write(179,*) 'sum(proj)', sum(proj)
      
c$$$      weight_low=0.0d0
c$$$      weight_total=0.0d0      
c$$$  pdos=0.0d0
      
      do iatom=1, natom
        isort=is(iatom)            
        do llval=0,lmb(isort)
          if (nrel .eq. 1) then
            lival=llval+1
            do mjval=-llval, llval
              mjind=mjval+llval+1
              do ind_k=1,ndim_kk(me+1) ! k vector
                irk=n_mpi_kk(me+1)+ind_k
                k0=i_kref(irk)                
                do ibnd=1, num_bands
c$$$  do iomega=-500,500
c$$$  omega=iomega*0.1d0
c$$$  pdos(iomega,llval,1,isort)
c$$$  $                =pdos(iomega,llval,1,isort)
c$$$  $                +1.0d0/dsqrt(2.0d0*pi)/0.1d0
c$$$  $                *dexp(-(omega-eigenvalues(ibnd,irk))**2/0.02d0)
c$$$  $                *cdabs(proj(llval,1,mjind,iatom,ibnd,ind_k))**2
c$$$  $                /dble(nqdiv)/sort_degen(isort)
c$$$  enddo
                  weight_total(llval,1,isort)
     $              =weight_total(llval,1,isort)
     $              +cdabs(proj(llval,1,mjind,iatom,ibnd,ind_k))**2
     $              /dble(nqdiv)/sort_degen(isort)                  
                  if (eigenvalues(ibnd,irk) .lt. dis_froz_min) then
                    weight_low(llval,1,isort)
     $                =weight_low(llval,1,isort)
     $                +cdabs(proj(llval,1,mjind,iatom,ibnd,ind_k))**2
     $                /dble(nqdiv)/sort_degen(isort)
                  endif
                enddo
              enddo
            enddo                
          else
            do iival=-1,1,2
              if ((llval .eq. 0) .and. (iival .eq. -1)) cycle
              i_index=(iival+1)/2+1              
              jjval=llval*2+iival
              jjval_dble=llval+iival/2.0d0
              jjdim=nint(2*jjval_dble+1)
              lival=llval+llval+(iival+1)/2

              do mjval=-jjval, jjval, 2
                mjind=(mjval+jjval+2)/2
c$$$  write(179,'(4i5)') iatom,llval,iival,mjval
                
                do ind_k=1,ndim_kk(me+1) ! k vector
                  irk=n_mpi_kk(me+1)+ind_k
                  k0=i_kref(irk)                                
                  do ibnd=1, num_bands
c$$$  do iomega=-500,500
c$$$  omega=iomega*0.1d0
c$$$  pdos(iomega,llval,i_index,isort)
c$$$  $                  =pdos(iomega,llval,i_index,isort)
c$$$  $                  +1.0d0/dsqrt(2.0d0*pi)/0.1d0
c$$$  $                  *dexp(-(omega-eigenvalues(ibnd,irk))**2/0.02d0)
c$$$  $                  *cdabs(proj(llval,i_index,mjind,
c$$$  $                  iatom,ibnd,ind_k))**2
c$$$  $                  /dble(nqdiv)/sort_degen(isort)
c$$$  enddo
                    weight_total(llval,i_index,isort)
     $                =weight_total(llval,i_index,isort)
     $                +cdabs(proj(llval,i_index,mjind,
     $                iatom,ibnd,ind_k))**2
     $                /dble(nqdiv)/sort_degen(isort)                    
                    if (eigenvalues(ibnd,irk) .lt. dis_froz_min) then
c$$$  write(179,'(6i5)')
c$$$  $                  iatom,llval,iival,mjval,irk,ibnd                      
                      weight_low(llval,i_index,isort)
     $                  =weight_low(llval,i_index,isort)
     $                  +cdabs(proj(llval,i_index,mjind,
     $                  iatom,ibnd,ind_k))**2
     $                  /dble(nqdiv)/sort_degen(isort)
                    endif
                  enddo
                enddo
              enddo                
            enddo
          endif            
        enddo
      enddo

c$$$  write(179,*) 'weight_low', sum(weight_low)      


      if(nproc/=1) then
        call mpi_allreduce_dble
     $    (weight_low(0,1,1),(maxb+1)*nrel*nsort,
     $    mpi_sum,mpi_comm_world)
        call mpi_allreduce_dble
     $    (weight_total(0,1,1),(maxb+1)*nrel*nsort,
     $    mpi_sum,mpi_comm_world)        
c$$$  call mpi_allreduce_dble
c$$$  $    (pdos(-500,0,1,1),1001*(maxb+1)*nrel*nsort,
c$$$  $    mpi_sum,mpi_comm_world)        
      endif

c$$$  write(179,*) 'weight_low', sum(weight_low)            
      
      if (maswrk) then
c$$$  open(unit=100, file='pdos.out')
c$$$  do iomega=-500,500
c$$$  omega=iomega*0.1d0
c$$$  write(100,'(10000f12.6)') omega,
c$$$  $      (((pdos(iomega,llval,i_index,isort),i_index=1, nrel),
c$$$  $      llval=0,maxb), isort=1, nsort)
c$$$  enddo
c$$$  close(100)

        open(unit=41, file='weight_trial.dat')
        if (irel .le. 1) then        
          do isort=1, nsort
            do l=0,lmb(isort)   ! l
              li=l+1            !li              
              write(41, '(2(a, i4, 5x), 100f12.6)')
     $          'isort', isort, 'li', li,
     $          weight_in_froz_win(l,1,isort),weight_low(l,1,isort),
     $          weight_total(l,1,isort)
            enddo
          enddo
        else
          do isort=1, nsort
            do l=0,lmb(isort)   ! l
              ist=-1              
              if(l.eq.0) ist=1
              do iival=ist,1,2  ! i
                i_index=(iival+1)/2+1
                li=l*2+(iival+1)/2
                jjval=l+l+iival ! j              
                write(41, '(2(a, i4, 5x), 100f12.6)')
     $            'isort', isort, 'li', li,
     $            weight_in_froz_win(l,i_index,isort),
     $            weight_low(l,i_index,isort),
     $            weight_total(l,i_index,isort)                
              enddo
            enddo
          enddo
        endif
        close(41)

        open(unit=41, file='occ_trial.dat')
        if (irel .le. 1) then        
          do isort=1, nsort
            do l=0,lmb(isort)   ! l
              li=l+1            !li              
              write(41, '(2(a, i4, 5x), 100f12.6)')
     $          'isort', isort, 'li', li,
     $          occ_iwan(l,1,isort)
            enddo
          enddo
        else
          do isort=1, nsort
            do l=0,lmb(isort)   ! l
              ist=-1              
              if(l.eq.0) ist=1
              do iival=ist,1,2  ! i
                i_index=(iival+1)/2+1
                li=l*2+(iival+1)/2
                jjval=l+l+iival ! j              
                write(41, '(2(a, i4, 5x), 100f12.6)')
     $            'isort', isort, 'li', li,
     $            occ_iwan(l,i_index,isort)
              enddo
            enddo
          enddo
        endif
        close(41)
        
        

        write(iun,*) 'weight summary'
        if (irel .le. 1) then        
          do isort=1, nsort
            do l=0,lmb(isort)   ! l
              li=l+1            !li              
              write(iun, '(4(a, i4, 5x), 100f12.6)')
     $          'isort', isort, 'l', l, 'mu_dim', 2*l+1, 'li', li,
     $          weight_in_froz_win(l,1,isort),weight_low(l,1,isort),
     $          weight_total(l,1,isort)
            enddo
          enddo
        else
          do isort=1, nsort
            do l=0,lmb(isort)   ! l
              ist=-1              
              if(l.eq.0) ist=1
              do iival=ist,1,2  ! i
                i_index=(iival+1)/2+1
                li=l*2+(iival+1)/2
                jjval=l+l+iival ! j              
                write(iun, '(5(a, i4, 5x), 100f12.6)')
     $            'isort', isort, 'l', l, 'i',
     $            iival,'mu_dim', jjval+1, 'li', li,
     $            weight_in_froz_win(l,i_index,isort),
     $            weight_low(l,i_index,isort),
     $            weight_total(l,i_index,isort)                
              enddo
            enddo
          enddo
        endif


        write(iun,*) 'occ summary'
        if (irel .le. 1) then        
          do isort=1, nsort
            do l=0,lmb(isort)   ! l
              li=l+1            !li              
              write(iun, '(4(a, i4, 5x), 100f12.6)')
     $          'isort', isort, 'l', l, 'mu_dim', 2*l+1, 'li', li,
     $          occ_iwan(l,1,isort)
            enddo
          enddo
        else
          do isort=1, nsort
            do l=0,lmb(isort)   ! l
              ist=-1              
              if(l.eq.0) ist=1
              do iival=ist,1,2  ! i
                i_index=(iival+1)/2+1
                li=l*2+(iival+1)/2
                jjval=l+l+iival ! j              
                write(iun, '(5(a, i4, 5x), 100f12.6)')
     $            'isort', isort, 'l', l, 'i',
     $            iival,'mu_dim', jjval+1, 'li', li,
     $            occ_iwan(l,i_index,isort)
              enddo
            enddo
          enddo
        endif
        
        
        num_wann=0
        open(unit=42, file='orb_for_froz_win.dat')
c$$$        write(42,*) 'nbnd_in_froz', nbnd_in_froz
        open(unit=43, file=adjustl(trim(seed_name))//'.inip')

        do iatom=1, natom
          isort=is(iatom)
          do lval=0, lmb(isort)
            mdim=lval*2+1
            if (irel .le. 1) then              
              if (
     $          (weight_in_froz_win(lval,1,isort).ge.cut_froz)
     $          .and.
     $          (weight_low(lval,1,isort).le.cut_low)
     $          )then                      
                num_wann=num_wann+mdim
                do mm=-lval, lval
                  write(43,'(3i5,6f12.6)') iatom, lval, mm,
     $              atom_xaxis(:,iatom),atom_zaxis(:,iatom)
                enddo
                write(42,'(3i5, f12.6)') iatom, isort, lval,
     $            weight_in_froz_win(lval,1,isort)
              endif
            else
              ist=-1
              if(lval.eq.0) ist=1
              do iival=ist,1,2  ! 2i
                i_index=(iival+1)/2+1
                jjval=2*lval+iival ! 2j              
                mdim=jjval+1
                if (
     $            (weight_in_froz_win(lval,i_index,isort).ge.cut_froz)
     $            .and.
     $            (weight_low(lval,i_index,isort).le.cut_low)
     $            )then                                      
                  num_wann=num_wann+mdim
                  do mm=-jjval, jjval, 2
                    write(43,'(2i5,2x,2(f4.1, 2x),6f12.6)')
     $                iatom, lval, iival/2.0d0,mm/2.0d0,
     $                atom_xaxis(:,iatom),atom_zaxis(:,iatom)
                  enddo
                  write(42,'(3i5, 2x,f4.1, 2x, 1f12.6)')iatom,isort,
     $              lval,iival/2.0d0,
     $              weight_in_froz_win(lval,i_index,isort)
                endif
              enddo
            endif
          enddo
        enddo
        write(42, *) 'num_wann', num_wann
        close(42)        
        close(43)          
      endif

      call mpi_bcast(num_wann,1,mpi_integer,0,mpi_comm_world,mpierr)

      end

      subroutine count_nbnd_in_froz_win
      use comwann_mod
      implicit none

      integer :: ik,tempint,ibnd, ind(1)
      


      if (maswrk) then
        open(unit=112, file='band_in_froz_win.dat')
      endif      

      do ik=1, nqdiv


        ind=minloc(eigenvalues(:,ik),
     $    mask=eigenvalues(:,ik) .gt. dis_froz_min)
        band_ind_froz(1,ik)=ind(1)                
        
        band_eig_froz(1,ik)        
     $    =eigenvalues(band_ind_froz(1,ik),ik)


        ind=maxloc(eigenvalues(:,ik),
     $    mask=eigenvalues(:,ik) .lt. dis_froz_max)
        band_ind_froz(2,ik)=ind(1)                

        band_eig_froz(2,ik)        
     $    =eigenvalues(band_ind_froz(2,ik),ik)        

        band_ind_froz(3,ik)
     $    = band_ind_froz(2,ik)-band_ind_froz(1,ik)+1

        if (maswrk) then
          write(112,
     $      '(a,2x,i5,2x,a,2x,2(i5,2x),a,2x,2(f20.12,2x),a,2x,i5)')
     $      'ik=', ik, 'index:',
     $      band_ind_froz(1:2,ik),  'energy:',
     $      band_eig_froz(1:2,ik),
     $      'total number:' ,
     $      band_ind_froz(3,ik)
        endif
      enddo
      if(maswrk) then
        close(112)
      endif
      end


      subroutine reset_dis_froz_max
      use comwann_mod
      implicit none
c$$$      integer :: ik, ind
c$$$      double precision :: tempval

c$$$      tempval=0.0d0
c$$$      do ik=1, nqdiv
c$$$        ind=band_ind_froz(1,ik)+band_ind_froz(3,ik)
c$$$        if 
      dis_froz_max=dis_froz_max*dis_froz_max_fac
      if (maswrk) then
        write(iun,*) 
        write(iun,*) 'dis_froz_max has been reset to ', dis_froz_max
        write(iun,*)         
      endif
      end
      
      

      subroutine poorman_wan_proj(ispin)
      use comwann_mod
      implicit none
      include "mpif.h"      
      integer, intent(in) :: ispin
      integer :: ie1,inn1,in1,
     $  ii,jj,k,k0,ibnd,radind1,m1,lm1,
     $  km1,isort,iwan,
     $  iival,jjval,llval,mjval,lival,iatom,ind0,
     $  ist,ind_k,km,irk,iband,iorb, lldim,jjdim,mjind,i_index
      integer, allocatable :: ipiv(:)
      double precision :: llval_dble,jjval_dble,weight
      complex*16 :: znew(nfun,num_bands),
     $  proj_poor_old(num_wann,num_bands),
     $  proj_poor_new(num_wann,num_bands,nqdiv)
      complex*16, allocatable :: wfn_mt0(:,:), wfn_mt1(:,:),
     $  rot_matrix(:,:)

!     radial function determination
c$$$  dimension determination`
      proj_poor_new=0.0d0
      proj_poor_old=0.0d0


      do ind_k=1,ndim_kk(me+1)  ! k vector
        k=n_mpi_kk(me+1)+ind_k
        k0=i_kref(k)
        
        call sym_z_0(znew,k,z_wan_bnd(1,1,k0),
     &    num_bands,k_group(k),pnt(1,k))
        
        proj_poor_old=0.0d0        
        do ibnd=1,num_bands            
          
          do iwan=1, num_wann
            do iorb=1, nalimj(iwan)
              iatom=alimj(1,iorb,iwan)
              isort=is(iatom)
              ind0=io_lem(iatom)-1
              llval=alimj(2,iorb,iwan)
              
              if (irel .le. 1) then
                mjval=alimj(3,iorb,iwan)        
                lival=llval+1
                i_index=1                
                
                allocate(rot_matrix(2*llval+1,2*llval+1))
                rot_matrix=0.0d0

                llval_dble=dble(llval)
                lldim=2*llval+1
                call rotmat_wan_proj(alimj_xaxis(1,iorb,iwan),
     $            alimj_zaxis(1,iorb,iwan),llval_dble,lldim,1,
     $            rot_matrix)
                allocate(wfn_mt0(2*llval+1,
     $            dim_radial(llval,i_index,iatom)))
                wfn_mt0=0.0d0
                allocate(wfn_mt1(2*llval+1,
     $            dim_radial(llval,i_index,iatom)))
                wfn_mt1=0.0d0                  
                radind1=0
                do ie1=1,ntle(llval,isort)
                  in1=1
                  if(augm(ie1,llval,isort)/='LOC') in1=2
                  do inn1=1,in1
                    radind1=radind1+1
                    do m1=-llval, llval
                      lm1=llval*(llval+1)+m1+1
                      km1=ind0+indbasa(inn1,ie1,lm1,isort)
                      wfn_mt0(m1+llval+1,radind1)
     $                  =wfn_mt0(m1+llval+1,radind1)
     $                  +znew(km1,ibnd)
                    enddo                    
                  enddo
                enddo
                wfn_mt1=
     $            matmul(rot_matrix,wfn_mt0)
                do radind1=1, dim_radial(llval,i_index,iatom)
                  proj_poor_old(iwan,ibnd)
     $              =proj_poor_old(iwan,ibnd)
     $              +alimj_coeff(iorb,iwan)
     $              *dconjg(coeff_radial(radind1,llval,i_index,iatom))
     $              *wfn_mt1(mjval,radind1)
                enddo
                deallocate(rot_matrix)
                deallocate(wfn_mt0)
                deallocate(wfn_mt1)
              else
                iival=alimj(3,iorb,iwan)
                mjval=alimj(4,iorb,iwan)                                        
                jjval=llval*2+iival
                jjval_dble=llval+iival/2.0d0
                jjdim=nint(2*jjval_dble+1)
                lival=llval+llval+(iival+1)/2
                i_index=(iival+1)/2+1                                        


                allocate(rot_matrix(jjdim,jjdim))
                rot_matrix=0.0d0
                allocate(wfn_mt0(jjdim,dim_radial(llval,i_index,iatom)))
                wfn_mt0=0.0d0
                allocate(wfn_mt1(jjdim,dim_radial(llval,i_index,iatom)))
                wfn_mt1=0.0d0                                  
                
                call rotmat_wan_proj(alimj_xaxis(1,iorb,iwan),
     $            alimj_zaxis(1,iorb,iwan),jjval_dble,jjdim,0,
     $            rot_matrix)

                radind1=0
                do ie1=1,ntle(llval,isort)
                  in1=1
                  if(augm(ie1,llval,isort)/='LOC') in1=2
                  do inn1=1,in1
                    radind1=radind1+1
                    do m1=-jjval, jjval, 2
                      call getlimj(lm1,llval,iival,m1,lival,1)                      
                      km1=ind0+indbasa(inn1,ie1,lm1,isort)
                      mjind=(m1+jjval+2)/2
                      wfn_mt0(mjind,radind1)
     $                  =wfn_mt0(mjind,radind1)
     $                  +znew(km1,ibnd)
                    enddo                    
                  enddo
                enddo
                wfn_mt1=
     $            matmul(rot_matrix,wfn_mt0)
                
                do radind1=1, dim_radial(llval,i_index,iatom)
                  mjind=(mjval+jjval+2)/2                  
                  proj_poor_old(iwan,ibnd)
     $              =proj_poor_old(iwan,ibnd)
     $              +alimj_coeff(iorb,iwan)
     $              *dconjg(coeff_radial(radind1,llval,i_index,iatom))                  
     $              *wfn_mt1(mjind,radind1)
                enddo                

                deallocate(rot_matrix)
                deallocate(wfn_mt0)
                deallocate(wfn_mt1)
                
              endif
            enddo
          enddo
        enddo
c     write(179,'(a, 2i5)')
c     $    'symmetric_orthonormalization in poorman_wan_proj',
c     $    ind_k, k
        call symmetric_orthonormalization(proj_poor_old,
     $    proj_poor_new(1,1,k),num_wann,num_bands)
      enddo


      if(nproc/=1) then
        call mpi_allreduce_dcmplx
     $    (proj_poor_new,num_bands*num_wann*nqdiv,
     $    mpi_sum,mpi_comm_world)
      endif

      a_matrix=0.0d0
      do irk = 1, nqdiv
        do iwan = 1, num_wann
          do iband = 1, num_bands
            a_matrix(iband,iwan,irk)
     $        =dconjg(proj_poor_new(iwan,iband,irk))
          enddo
        enddo
      enddo

      end      
      

      subroutine symmetric_orthonormalization(proj_old,proj_new,n1,n2)

      implicit none
      integer, intent(in) :: n1,n2
      complex*16, intent(in) :: proj_old(n1,n2)
      complex*16, intent(out) :: proj_new(n1,n2)

      integer :: j,ii,jj
      double precision :: diag(n1)
      complex*16 :: ov(n1,n1),t0(n1,n1),t(n1,n1)


      proj_new=0.0d0
      call zgemm('n','c',n1,n1,n2,(1.d0,0.d0),proj_old,
     &  n1,proj_old,n1,(0.d0,0.d0),ov,
     &  n1)

c     -------- Spectrum of Overlap matrix -------------------------------
      call hermitianeigen_cmplxdouble(n1,diag,ov)
c     write(179,'(a, f12.6, a, f12.6)')
c     $  'maxval', diag(n1), '     minval', diag(1)
c     write(179,'(10000f12.6)') diag
      
      
c     ------ Here we get T = O^(-1/2) -----------------------------------
      t0=0.0d0
      do j=1,n1
        if (diag(j) .gt. 0.0001d0) then
          t0(:,j)=ov(:,j)/sqrt(diag(j))
        endif
      enddo
      call zgemm('n','c',n1,n1,n1,
     &  (1.d0,0.d0),t0,n1,ov,n1,(0.d0,0.d0),
     &  t,n1)
      
c$$$  --------- Orthonormalized projector ----------------------
      call zgemm('c','n',n1,n2,n1,(1.d0,0.d0),t,
     &  n1,proj_old,n1,(0.d0,0.d0),proj_new,n1)
      end



      subroutine cal_overlapmax_bandprojection_to_mtorb_low_rmode
     $  (iatom,l,iival,jjval,li,ispin,dim_radial_orb,
     $  coeff_radial_orb,weight,occ)


      use comwann_mod
      implicit none
      include "mpif.h"
      integer, intent(in) :: iatom,l,iival,jjval,li,ispin! if irel.le.1 then iival and jjval is dummy variable

      integer, intent(out) :: dim_radial_orb
      double precision, intent(out) :: weight,occ
      complex*16, intent(out) :: coeff_radial_orb(100)

      if (rmode .eq. 0) then
        
        call cal_overlapmax_bandprojection_to_mtorb_low
     $    (iatom,l,iival,jjval,li,ispin,dim_radial_orb,
     $    coeff_radial_orb,weight,occ)

      elseif (rmode .eq. 1) then
        
        call cal_overlapmax_bandprojection_to_mtorb_low_0
     $    (iatom,l,iival,jjval,li,ispin,dim_radial_orb,
     $    coeff_radial_orb,weight,occ)
      elseif (rmode .eq. 2) then
        call cal_overlapmax_bandprojection_to_mtorb_low_new
     $    (iatom,l,iival,jjval,li,ispin,dim_radial_orb,
     $    coeff_radial_orb,weight,occ)
      endif
      end
      

      subroutine cal_overlapmax_bandprojection_to_mtorb_low
     $  (iatom,l,iival,jjval,li,ispin,dim_radial_orb,
     $  coeff_radial_orb,weight,occ)
      
      use comwann_mod
      implicit none
      include "mpif.h"
      integer, intent(in) :: iatom,l,iival,jjval,li,ispin ! if irel.le.1 then iival and jjval is dummy variable

      integer, intent(out) :: dim_radial_orb
      double precision, intent(out) :: weight,occ
      complex*16, intent(out) :: coeff_radial_orb(100)
      
      integer :: ie1,inn1,inn2,ie2,in2,in1,ind1,ind2,
     $  ii,jj,k,k0,ibnd,radind1,radind2,radind3,radind4,m1,lm1,
     $  km1,km2,isort,ind0,ind_k,
     $  mj,ival, cnt, cnt1,ir,nres,mtw,maxval_loc(1)
      double precision :: tempemin, tempemax           
      complex*16 :: znew(nfun,num_bands)
      double precision, allocatable :: overlap(:,:),t0(:,:),shalf(:,:),
     $  partial_proj_eigval(:),diag(:),smhalf(:,:),smat(:,:),sinv(:,:)
      complex*16, allocatable :: partial_proj(:,:)

!     radial function determination
c$$$  dimension determination`


      coeff_radial_orb=0.0d0

      isort=is(iatom)
      ind0=io_lem(iatom)-1
      
      dim_radial_orb=0
      do ie1=1,ntle(l,isort)
        in1=1
        if(augm(ie1,l,isort)/='LOC') in1=2
        do inn1=1,in1
          dim_radial_orb=dim_radial_orb+1
        enddo
      enddo

      
      allocate(overlap(dim_radial_orb,dim_radial_orb))
      overlap=0.0d0
      allocate(t0(dim_radial_orb,dim_radial_orb))
      t0=0.0d0
      allocate(shalf(dim_radial_orb,dim_radial_orb))
      shalf=0.0d0
      allocate(smhalf(dim_radial_orb,dim_radial_orb))
      smhalf=0.0d0
      allocate(smat(dim_radial_orb,dim_radial_orb))
      smat=0.0d0
      allocate(sinv(dim_radial_orb,dim_radial_orb))
      smat=0.0d0                              
      allocate(partial_proj(dim_radial_orb,dim_radial_orb))
      partial_proj=0.0d0
      allocate(partial_proj_eigval(dim_radial_orb))
      partial_proj_eigval=0.0d0
      allocate(diag(dim_radial_orb))
      diag=0.0d0
      
      ind1=0

      do ie1=1,ntle(l,isort)
        in1=1
        if(augm(ie1,l,isort)/='LOC') in1=2
        do inn1=1,in1
          ind1=ind1+1
          ind2=0
          do ie2=1,ntle(l,isort)
            in2=1
            if(augm(ie2,l,isort)/='LOC') in2=2
            do inn2=1,in2
              ind2=ind2+1
              overlap(ind1, ind2)
     $          =ffsmt(inn1,inn2,ie1,ie2,li,isort,ispin)
            enddo
          enddo
        enddo
      enddo
      smat=overlap

      

!     symmetric orthonormalization with neglection of linear contribution
      call symmetriceigen_double(dim_radial_orb,diag,overlap)
      t0=0.0d0
      do jj=1,dim_radial_orb
        t0(:,jj)=overlap(:,jj)*dsqrt(diag(jj)) 
      enddo
      call dgemm('n','c',dim_radial_orb,dim_radial_orb,dim_radial_orb,
     &  1.0d0,t0,dim_radial_orb,overlap,dim_radial_orb,0.0d0,
     &  shalf,dim_radial_orb)

      t0=0.0d0
      do jj=1,dim_radial_orb
        t0(:,jj)=overlap(:,jj)/dsqrt(diag(jj)) 
      enddo
      call dgemm('n','c',dim_radial_orb,dim_radial_orb,dim_radial_orb,
     &  1.0d0,t0,dim_radial_orb,overlap,dim_radial_orb,0.0d0,
     &  smhalf,dim_radial_orb)

      t0=0.0d0
      do jj=1,dim_radial_orb
        t0(:,jj)=overlap(:,jj)/diag(jj)
      enddo
      call dgemm('n','c',dim_radial_orb,dim_radial_orb,dim_radial_orb,
     &  1.0d0,t0,dim_radial_orb,overlap,dim_radial_orb,0.0d0,
     &  sinv,dim_radial_orb)
      
      
      partial_proj=0.0d0

c$$$  if (l .eq. 3)then
c$$$  tempemin=dis_froz_min
c$$$  tempemax=0.0d0
c$$$  else
      tempemin=dis_froz_min
      tempemax=dis_froz_max
c$$$  endif
c$$$

      do ind_k=1,ndim_kk(me+1)  ! k vector
        k=n_mpi_kk(me+1)+ind_k
        k0=i_kref(k)
        
        call sym_z_0(znew,k,z_wan_bnd(1,1,k0),
     &    num_bands,k_group(k),pnt(1,k))
        
        do ibnd=1,num_bands            
          if ((eigenvalues(ibnd,k) .gt. tempemin) .and.
     $      (eigenvalues(ibnd,k) .lt. tempemax)) then
            radind1=0
            do ie1=1,ntle(l,isort)
              in1=1
              if(augm(ie1,l,isort)/='LOC') in1=2
              do inn1=1,in1
                radind1=radind1+1
                
                radind2=0
                do ie2=1,ntle(l,isort)
                  in2=1
                  if(augm(ie2,l,isort)/='LOC') in2=2
                  do inn2=1,in2
                    radind2=radind2+1
                    if (irel .le. 1) then
                      do m1=-l, l
                        lm1=l*(l+1)+m1+1
                        km1=ind0+indbasa(inn1,ie1,lm1,isort)
                        km2=ind0+indbasa(inn2,ie2,lm1,isort)
                        do radind3=1,dim_radial_orb
                          do radind4=1,dim_radial_orb
                            partial_proj(radind3,radind4)=
     $                        partial_proj(radind3,radind4)+
     $                        znew(km1,ibnd)*dconjg(znew(km2,ibnd))
     $                        *shalf(radind2,radind4)
     $                        *shalf(radind3,radind1)
     $                        /dble(nqdiv)
     $                        /(exp(eigenvalues(ibnd,k)
     $                        /(8.617332262145d0*0.9d0*10.0d0**-2))
     $                        +1.0d0)
                          enddo
                        enddo
                      enddo
                    else
                      do mj=-jjval, jjval,2
                        call getlimj(lm1,l,iival,mj,li,1)
                        km1=ind0+indbasa(inn1,ie1,lm1,isort)
                        km2=ind0+indbasa(inn2,ie2,lm1,isort)
                        do radind3=1,dim_radial_orb
                          do radind4=1,dim_radial_orb
                            partial_proj(radind3,radind4)=
     $                        partial_proj(radind3,radind4)+
     $                        znew(km1,ibnd)*dconjg(znew(km2,ibnd))
     $                        *shalf(radind2,radind4)
     $                        *shalf(radind3,radind1)
     $                        /dble(nqdiv)
     $                        /(exp(eigenvalues(ibnd,k)
     $                        /(8.617332262145d0*0.9d0*10.0d0**-2))
     $                        +1.0d0)
                          enddo
                        enddo
                      enddo
                    endif
                  enddo
                enddo
              enddo
            enddo
          endif
        enddo
      enddo
      
      
      
      if(nproc/=1) then
        call mpi_allreduce_dcmplx
     $    (partial_proj,dim_radial_orb**2,mpi_sum,mpi_comm_world)        
      endif


      call hermitianeigen_cmplxdouble
     $  (dim_radial_orb,partial_proj_eigval,partial_proj)

      coeff_radial_orb=0.0d0
      do ii=1, dim_radial_orb
        do jj=1, dim_radial_orb
          coeff_radial_orb(ii)=coeff_radial_orb(ii)+
     $      shalf(ii,jj)
     $      *partial_proj(jj,dim_radial_orb)
        enddo
      enddo

      ind1=0
      do ie1=1,ntle(l,isort)
        in1=1
        if(augm(ie1,l,isort)/='LOC') in1=2
        do inn1=1,in1
          ind1=ind1+1
          bound_radfun_overlap(inn1,ie1,li,iatom,ispin)
     $      =coeff_radial_orb(ind1)
        enddo
      enddo

      ind1=0
      do ie1=1,ntle(l,isort)
        in1=1
        if(augm(ie1,l,isort)/='LOC') in1=2
        do inn1=1,in1
          ind1=ind1+1
          
          ind2=0
          do ie2=1,ntle(l,isort)
            in2=1
            if(augm(ie2,l,isort)/='LOC') in2=2
            do inn2=1,in2
              ind2=ind2+1                
              
              bound_radfun_coeff(inn1,ie1,li,iatom,ispin)
     $          =bound_radfun_coeff(inn1,ie1,li,iatom,ispin)
     $          +sinv(ind1,ind2)
     $          *bound_radfun_overlap(inn2,ie2,li,iatom,ispin)
            enddo
          enddo
        enddo
      enddo


      if (maswrk) then
        write(iun, '(a, i5, a, i3, a, i5)')
     $    'bound radial function overlap and coeff, iatom',
     $    iatom, 'li', li, 'ispin', ispin
        do ie1=1,ntle(l,isort)
          in1=1
          if(augm(ie1,l,isort)/='LOC') in1=2
          do inn1=1,in1          
            write(iun,'(a,i5,a,i5,a,2f12.6)')
     $        'ie', ie1, 'jn', inn1, 'overlap', 
     $        bound_radfun_overlap(inn1,ie1,li,iatom,ispin),
     $        bound_radfun_coeff(inn1,ie1,li,iatom,ispin)
          enddo
        enddo
      endif      

      do ir=0, nrad(isort)
        do ie1=1,ntle(l,isort)
          in1=1
          if(augm(ie1,l,isort)/='LOC') in1=2
          do inn1=1,in1
            nres=indfun0(inn1,ie1,li,isort)
            mtw=ind_wf(nres,isort)
            bound_radfun(ir,1,li,iatom,ispin)
     $        =bound_radfun(ir,1,li,iatom,ispin)
     $        +bound_radfun_coeff(inn1,ie1,li,iatom,ispin)
     $        *gfun(mtw+ir,ispin)
            if (irel .le. 1) then
              bound_radfun(ir,2,li,iatom,ispin)
     $          =bound_radfun(ir,2,li,iatom,ispin)
     $          +bound_radfun_coeff(inn1,ie1,li,iatom,ispin)
     $          *gfund(mtw+ir,ispin)
            endif
          enddo
        enddo
      enddo


            

        

c$$$  if (l .eq. 3)then
c$$$  coeff_radial_orb(2)=coeff_radial_orb(2)*2.0
c$$$  endif
      weight=0.0d0

      occ=0.0d0

      do ind_k=1,ndim_kk(me+1)  ! k vector
        k=n_mpi_kk(me+1)+ind_k
        k0=i_kref(k)
        
        call sym_z_0(znew,k,z_wan_bnd(1,1,k0),
     &    num_bands,k_group(k),pnt(1,k))
        
        do ibnd=1,num_bands            
          radind1=0
          do ie1=1,ntle(l,isort)
            in1=1
            if(augm(ie1,l,isort)/='LOC') in1=2
            do inn1=1,in1
              radind1=radind1+1
              
              radind2=0
              do ie2=1,ntle(l,isort)
                in2=1
                if(augm(ie2,l,isort)/='LOC') in2=2
                do inn2=1,in2
                  radind2=radind2+1
                  
                  if (irel .le. 1) then
                    do m1=-l, l
                      lm1=l*(l+1)+m1+1
                      km1=ind0+indbasa(inn1,ie1,lm1,isort)
                      km2=ind0+indbasa(inn2,ie2,lm1,isort)
                      occ=occ
     $                  +dble(znew(km1,ibnd)*dconjg(znew(km2,ibnd))
     $                  *dconjg(coeff_radial_orb(radind1))
     $                  *coeff_radial_orb(radind2))
     $                  /dble(nqdiv)
     $                  /(dexp(eigenvalues(ibnd,k)
     $                  /(8.617332262145d-5*0.9d3))
     $                  +1.0d0)
                      
                      if ((eigenvalues(ibnd,k) .gt. tempemin) .and.
     $                  (eigenvalues(ibnd,k) .lt. tempemax)) then
                        weight=weight
     $                    +dble(znew(km1,ibnd)*dconjg(znew(km2,ibnd))
     $                    *dconjg(coeff_radial_orb(radind1))
     $                    *coeff_radial_orb(radind2))                        
     $                    /dble(nqdiv)
                      endif
                    enddo
                  else
                    do mj=-jjval, jjval,2
                      call getlimj(lm1,l,iival,mj,li,1)
                      km1=ind0+indbasa(inn1,ie1,lm1,isort)
                      km2=ind0+indbasa(inn2,ie2,lm1,isort)
                      occ=occ
     $                  +dble(znew(km1,ibnd)*dconjg(znew(km2,ibnd))
     $                  *dconjg(coeff_radial_orb(radind1))
     $                  *coeff_radial_orb(radind2))                      
     $                  /dble(nqdiv)
     $                  /(exp(eigenvalues(ibnd,k)
     $                  /(8.617332262145d-5*0.9d3))
     $                  +1.0d0)
                      
                      if ((eigenvalues(ibnd,k) .gt. tempemin) .and.
     $                  (eigenvalues(ibnd,k) .lt. tempemax)) then
                        weight=weight
     $                    +dble(znew(km1,ibnd)*dconjg(znew(km2,ibnd))
     $                    *dconjg(coeff_radial_orb(radind1))
     $                    *coeff_radial_orb(radind2))                        
     $                    /dble(nqdiv)
                      endif
                    enddo
                  endif
                enddo
              enddo
            enddo
          enddo
c$$$  endif
        enddo
      enddo
      
      
      
      if(nproc/=1) then
        call mpi_allreduce_dble
     $    (weight,1,mpi_sum,mpi_comm_world)
        call mpi_allreduce_dcmplx
     $    (occ,1,mpi_sum,mpi_comm_world)                
      endif
      
      if (maswrk) then
        write(iun, '(3(a, i4, 5x), 100e20.12)')
     $    'occ, weight, iatom', iatom, 'l', l, 'li', li,
     $    occ,weight
      endif
      
      deallocate(overlap)      
      deallocate(shalf)
      deallocate(smhalf)
      deallocate(smat)
      deallocate(sinv)                  
      deallocate(t0)
      deallocate(diag)            
      deallocate(partial_proj)
      deallocate(partial_proj_eigval)

      end



      subroutine cal_overlapmax_bandprojection_to_mtorb_low_0
     $  (iatom,l,iival,jjval,li,ispin,dim_radial_orb,
     $  coeff_radial_orb,weight,occ)
      
      use comwann_mod
      implicit none
      include "mpif.h"
      integer, intent(in) :: iatom,l,iival,jjval,li,ispin ! if irel.le.1 then iival and jjval is dummy variable

      integer, intent(out) :: dim_radial_orb
      double precision, intent(out) :: weight,occ
      complex*16, intent(out) :: coeff_radial_orb(100)
      
      integer :: ie1,inn1,inn2,ie2,in2,in1,ind1,ind2,
     $  ii,jj,k,k0,ibnd,radind1,radind2,radind3,radind4,m1,lm1,
     $  km1,km2,isort,ind0,ind_k,
     $  mj,ival, cnt, cnt1,ir,nres,mtw,maxval_loc(1)
      double precision :: tempemin, tempemax           
      complex*16 :: znew(nfun,num_bands)
      double precision, allocatable :: overlap(:,:),t0(:,:),shalf(:,:),
     $  partial_proj_eigval(:),diag(:),smhalf(:,:),smat(:,:),sinv(:,:)
      complex*16, allocatable :: partial_proj(:,:)

!     radial function determination
c$$$  dimension determination`


      coeff_radial_orb=0.0d0

      isort=is(iatom)
      ind0=io_lem(iatom)-1
      
      dim_radial_orb=0
      do ie1=1,ntle(l,isort)
        in1=1
        if(augm(ie1,l,isort)/='LOC') in1=2
        do inn1=1,in1
          dim_radial_orb=dim_radial_orb+1
        enddo
      enddo

      
      allocate(overlap(dim_radial_orb,dim_radial_orb))
      overlap=0.0d0
      allocate(t0(dim_radial_orb,dim_radial_orb))
      t0=0.0d0
      allocate(shalf(dim_radial_orb,dim_radial_orb))
      shalf=0.0d0
      allocate(smhalf(dim_radial_orb,dim_radial_orb))
      smhalf=0.0d0
      allocate(smat(dim_radial_orb,dim_radial_orb))
      smat=0.0d0
      allocate(sinv(dim_radial_orb,dim_radial_orb))
      sinv=0.0d0                              
      allocate(partial_proj(dim_radial_orb,dim_radial_orb))
      partial_proj=0.0d0
      allocate(partial_proj_eigval(dim_radial_orb))
      partial_proj_eigval=0.0d0
      allocate(diag(dim_radial_orb))
      diag=0.0d0
      
      ind1=0

      do ie1=1,ntle(l,isort)
        in1=1
        if(augm(ie1,l,isort)/='LOC') in1=2
        do inn1=1,in1
          ind1=ind1+1
          ind2=0
          do ie2=1,ntle(l,isort)
            in2=1
            if(augm(ie2,l,isort)/='LOC') in2=2
            do inn2=1,in2
              ind2=ind2+1
              overlap(ind1, ind2)
     $          =ffsmt(inn1,inn2,ie1,ie2,li,isort,ispin)
            enddo
          enddo
        enddo
      enddo
      smat=overlap

      

!     symmetric orthonormalization with neglection of linear contribution
      call symmetriceigen_double(dim_radial_orb,diag,overlap)
      t0=0.0d0
      do jj=1,dim_radial_orb
        t0(:,jj)=overlap(:,jj)*dsqrt(diag(jj)) 
      enddo
      call dgemm('n','c',dim_radial_orb,dim_radial_orb,dim_radial_orb,
     &  1.0d0,t0,dim_radial_orb,overlap,dim_radial_orb,0.0d0,
     &  shalf,dim_radial_orb)

      t0=0.0d0
      do jj=1,dim_radial_orb
        t0(:,jj)=overlap(:,jj)/dsqrt(diag(jj)) 
      enddo
      call dgemm('n','c',dim_radial_orb,dim_radial_orb,dim_radial_orb,
     &  1.0d0,t0,dim_radial_orb,overlap,dim_radial_orb,0.0d0,
     &  smhalf,dim_radial_orb)

      t0=0.0d0
      do jj=1,dim_radial_orb
        t0(:,jj)=overlap(:,jj)/diag(jj)
      enddo
      call dgemm('n','c',dim_radial_orb,dim_radial_orb,dim_radial_orb,
     &  1.0d0,t0,dim_radial_orb,overlap,dim_radial_orb,0.0d0,
     &  sinv,dim_radial_orb)      


      tempemin=dis_froz_min
      tempemax=dis_froz_max
      
      partial_proj=0.0d0

      coeff_radial_orb=0.0d0

      ind1=0
      do ie1=1,ntle(l,isort)
        in1=1
        if(augm(ie1,l,isort)/='LOC') in1=2
        do inn1=1,in1
          ind1=ind1+1
          if (inn1 .eq. 2) then
            do ii=1, dim_radial_orb
              coeff_radial_orb(ii)=smat(ii,inn1)
     $          /dsqrt(smat(inn1,inn1))
            enddo
          endif
        enddo
      enddo


      ind1=0
      do ie1=1,ntle(l,isort)
        in1=1
        if(augm(ie1,l,isort)/='LOC') in1=2
        do inn1=1,in1
          ind1=ind1+1
          bound_radfun_overlap(inn1,ie1,li,iatom,ispin)
     $      =coeff_radial_orb(ind1)
        enddo
      enddo

      ind1=0
      do ie1=1,ntle(l,isort)
        in1=1
        if(augm(ie1,l,isort)/='LOC') in1=2
        do inn1=1,in1
          ind1=ind1+1
          
          ind2=0
          do ie2=1,ntle(l,isort)
            in2=1
            if(augm(ie2,l,isort)/='LOC') in2=2
            do inn2=1,in2
              ind2=ind2+1                
              
              bound_radfun_coeff(inn1,ie1,li,iatom,ispin)
     $          =bound_radfun_coeff(inn1,ie1,li,iatom,ispin)
     $          +sinv(ind1,ind2)
     $          *bound_radfun_overlap(inn2,ie2,li,iatom,ispin)
            enddo
          enddo
        enddo
      enddo


      if (maswrk) then
        write(iun, '(a, i5, a, i3, a, i5)')
     $    'bound radial function overlap and coeff, iatom',
     $    iatom, 'li', li, 'ispin', ispin
        do ie1=1,ntle(l,isort)
          in1=1
          if(augm(ie1,l,isort)/='LOC') in1=2
          do inn1=1,in1          
            write(iun,'(a,i5,a,i5,a,2f12.6)')
     $        'ie', ie1, 'jn', inn1, 'overlap', 
     $        bound_radfun_overlap(inn1,ie1,li,iatom,ispin),
     $        bound_radfun_coeff(inn1,ie1,li,iatom,ispin)
          enddo
        enddo
      endif

      do ir=0, nrad(isort)
        do ie1=1,ntle(l,isort)
          in1=1
          if(augm(ie1,l,isort)/='LOC') in1=2
          do inn1=1,in1
            nres=indfun0(inn1,ie1,li,isort)
            mtw=ind_wf(nres,isort)
            bound_radfun(ir,1,li,iatom,ispin)
     $        =bound_radfun(ir,1,li,iatom,ispin)
     $        +bound_radfun_coeff(inn1,ie1,li,iatom,ispin)
     $        *gfun(mtw+ir,ispin)
            if (irel .le. 1) then
              bound_radfun(ir,2,li,iatom,ispin)
     $          =bound_radfun(ir,2,li,iatom,ispin)
     $          +bound_radfun_coeff(inn1,ie1,li,iatom,ispin)
     $          *gfund(mtw+ir,ispin)
            endif
          enddo
        enddo
      enddo
            
      weight=0.0d0

      occ=0.0d0

      do ind_k=1,ndim_kk(me+1)  ! k vector
        k=n_mpi_kk(me+1)+ind_k
        k0=i_kref(k)
        
        call sym_z_0(znew,k,z_wan_bnd(1,1,k0),
     &    num_bands,k_group(k),pnt(1,k))
        
        do ibnd=1,num_bands            
          radind1=0
          do ie1=1,ntle(l,isort)
            in1=1
            if(augm(ie1,l,isort)/='LOC') in1=2
            do inn1=1,in1
              radind1=radind1+1
              
              radind2=0
              do ie2=1,ntle(l,isort)
                in2=1
                if(augm(ie2,l,isort)/='LOC') in2=2
                do inn2=1,in2
                  radind2=radind2+1
                  
                  if (irel .le. 1) then
                    do m1=-l, l
                      lm1=l*(l+1)+m1+1
                      km1=ind0+indbasa(inn1,ie1,lm1,isort)
                      km2=ind0+indbasa(inn2,ie2,lm1,isort)
                      occ=occ
     $                  +dble(znew(km1,ibnd)*dconjg(znew(km2,ibnd))
     $                  *dconjg(coeff_radial_orb(radind1))
     $                  *coeff_radial_orb(radind2))
     $                  /dble(nqdiv)
     $                  /(dexp(eigenvalues(ibnd,k)
     $                  /(8.617332262145d-5*0.9d3))
     $                  +1.0d0)
                      
                      if ((eigenvalues(ibnd,k) .gt. tempemin) .and.
     $                  (eigenvalues(ibnd,k) .lt. tempemax)) then
                        weight=weight
     $                    +dble(znew(km1,ibnd)*dconjg(znew(km2,ibnd))
     $                    *dconjg(coeff_radial_orb(radind1))
     $                    *coeff_radial_orb(radind2))                        
     $                    /dble(nqdiv)
                      endif
                    enddo
                  else
                    do mj=-jjval, jjval,2
                      call getlimj(lm1,l,iival,mj,li,1)
                      km1=ind0+indbasa(inn1,ie1,lm1,isort)
                      km2=ind0+indbasa(inn2,ie2,lm1,isort)
                      occ=occ
     $                  +dble(znew(km1,ibnd)*dconjg(znew(km2,ibnd))
     $                  *dconjg(coeff_radial_orb(radind1))
     $                  *coeff_radial_orb(radind2))                      
     $                  /dble(nqdiv)
     $                  /(exp(eigenvalues(ibnd,k)
     $                  /(8.617332262145d-5*0.9d3))
     $                  +1.0d0)
                      
                      if ((eigenvalues(ibnd,k) .gt. tempemin) .and.
     $                  (eigenvalues(ibnd,k) .lt. tempemax)) then
                        weight=weight
     $                    +dble(znew(km1,ibnd)*dconjg(znew(km2,ibnd))
     $                    *dconjg(coeff_radial_orb(radind1))
     $                    *coeff_radial_orb(radind2))                        
     $                    /dble(nqdiv)
                      endif
                    enddo
                  endif
                enddo
              enddo
            enddo
          enddo
c$$$  endif
        enddo
      enddo
      
      
      
      if(nproc/=1) then
        call mpi_allreduce_dble
     $    (weight,1,mpi_sum,mpi_comm_world)
        call mpi_allreduce_dcmplx
     $    (occ,1,mpi_sum,mpi_comm_world)                
      endif
      
      if (maswrk) then
        write(iun, '(3(a, i4, 5x), 100e20.12)')
     $    'occ, weight, iatom', iatom, 'l', l, 'li', li,
     $    occ,weight
      endif
      
c$$$
c$$$c$$$  if (printflag .eq. 1) then
c$$$c$$$  if (irel .le. 1) then
c$$$c$$$  ival=1
c$$$c$$$  else
c$$$c$$$  ival=(iival+1)/2+1
c$$$c$$$  endif
c$$$c$$$  weight_in_froz_win(l,ival,iatom)
c$$$c$$$  $  =weight_in_froz_win(l,ival,iatom)
c$$$c$$$  $  +partial_proj_eigval(dim_radial_orb)
c$$$      
c$$$      if (maswrk) then
c$$$        if (irel .le. 1) then
c$$$
c$$$          cnt1=0
c$$$          cnt=0
c$$$          do ie1=1,ntle(l,isort)
c$$$            in1=1
c$$$            if(augm(ie1,l,isort)/='LOC') in1=2
c$$$            do inn1=1,in1
c$$$              cnt=cnt+1              
c$$$              if (inn1 .eq. 2) then
c$$$                cnt1=cnt
c$$$              endif
c$$$            enddo
c$$$          enddo
c$$$
c$$$          write(iun, '(5(a, i4, 5x), 100f12.6)')
c$$$     $      'iatom', iatom, 'l', l, 'mu_dim', 2*l+1, 'li', li,
c$$$     $      'r2index', cnt1, 
c$$$     $      partial_proj_eigval(dim_radial_orb)
c$$$          
c$$$          cnt=0
c$$$          do ie1=1,ntle(l,isort)
c$$$            in1=1
c$$$            if(augm(ie1,l,isort)/='LOC') in1=2
c$$$            do inn1=1,in1
c$$$              cnt=cnt+1
c$$$              write(iun,'(3i5, 4f12.6)')
c$$$     $          ie1, inn1, cnt, shalf(cnt, cnt1),
c$$$     $          smat(cnt,cnt1)/dsqrt(smat(cnt1,cnt1)),
c$$$     $          coeff_radial_orb(cnt)
c$$$            enddo
c$$$          enddo
c$$$          
c$$$        else
c$$$          cnt1=0
c$$$          cnt=0
c$$$          do ie1=1,ntle(l,isort)
c$$$            in1=1
c$$$            if(augm(ie1,l,isort)/='LOC') in1=2
c$$$            do inn1=1,in1
c$$$              cnt=cnt+1              
c$$$              if (inn1 .eq. 2) then
c$$$                cnt1=cnt
c$$$              endif
c$$$            enddo
c$$$          enddo
c$$$          
c$$$          write(iun, '(6(a, i4, 5x), 100f12.6)')
c$$$     $      'iatom', iatom, 'l', l, 'i',
c$$$     $      iival,'mu_dim', jjval+1, 'li', li,
c$$$     $      'r2index', cnt1, 
c$$$     $      partial_proj_eigval(dim_radial_orb)
c$$$
c$$$          cnt=0
c$$$          do ie1=1,ntle(l,isort)
c$$$            in1=1
c$$$            if(augm(ie1,l,isort)/='LOC') in1=2
c$$$            do inn1=1,in1
c$$$              cnt=cnt+1
c$$$              write(iun,'(3i5, 4f12.6)')
c$$$     $          ie1, inn1, cnt, shalf(cnt, cnt1),
c$$$     $          smat(cnt,cnt1)/dsqrt(smat(cnt1,cnt1)),              
c$$$     $          coeff_radial_orb(cnt)
c$$$            enddo
c$$$          enddo          
c$$$          
c$$$        endif
c$$$      endif
c$$$  endif
      
      deallocate(overlap)      
      deallocate(shalf)
      deallocate(smhalf)
      deallocate(smat)
      deallocate(sinv)                  
      deallocate(t0)
      deallocate(diag)            
      deallocate(partial_proj)
      deallocate(partial_proj_eigval)

      end      



      subroutine cal_overlapmax_bandprojection_to_mtorb_low_new
     $  (iatom,l,iival,jjval,li,ispin,dim_radial_orb,
     $  coeff_radial_orb,weight,occ)
      
      use comwann_mod
      implicit none
      include "mpif.h"
      integer, intent(in) :: iatom,l,iival,jjval,li,ispin ! if irel.le.1 then iival and jjval is dummy variable

      integer, intent(out) :: dim_radial_orb
      double precision, intent(out) :: weight, occ
      complex*16, intent(out) :: coeff_radial_orb(100)
      
      integer :: ie1,inn1,inn2,ie2,in2,in1,ind1,ind2,
     $  ii,jj,k,k0,ibnd,radind1,radind2,radind3,radind4,m1,lm1,
     $  km1,km2,isort,ind0,ind_k,ie,in,jn,ir,
     $  mj,ival, cnt, cnt1,key11,nrad_bound(1),node_in,kc,nitec,
     $  node_out,nres,mtw,mt0,nre,mt,maxval_loc(1)
      double precision :: tempemin, tempemax,e,fifi,fi,dfi,convc,dl,acc,
     $  norm,c2
      complex*16 :: znew(nfun,num_bands)
      double precision, allocatable :: overlap(:,:),t0(:,:),shalf(:,:),
     $  diag(:),smhalf(:,:),smat(:,:),sinv(:,:),
     $  pbound(:),qbound(:),f_ex(:),f_exq(:),work(:),fr(:),v0(:)
      double precision, external :: dqdall

      
!     radial function determination
      c2=clight**2

      coeff_radial_orb=0.0d0
      
      isort=is(iatom)
      ind0=io_lem(iatom)-1
      
      dim_radial_orb=0
      do ie1=1,ntle(l,isort)
        in1=1
        if(augm(ie1,l,isort)/='LOC') in1=2
        do inn1=1,in1
          dim_radial_orb=dim_radial_orb+1
        enddo
      enddo

      
      allocate(overlap(dim_radial_orb,dim_radial_orb))
      overlap=0.0d0
      allocate(t0(dim_radial_orb,dim_radial_orb))
      t0=0.0d0
      allocate(shalf(dim_radial_orb,dim_radial_orb))
      shalf=0.0d0
      allocate(smhalf(dim_radial_orb,dim_radial_orb))
      smhalf=0.0d0
      allocate(smat(dim_radial_orb,dim_radial_orb))
      smat=0.0d0
      allocate(sinv(dim_radial_orb,dim_radial_orb))
      sinv=0.0d0                              
      allocate(diag(dim_radial_orb))
      diag=0.0d0
      
      ind1=0

      do ie1=1,ntle(l,isort)
        in1=1
        if(augm(ie1,l,isort)/='LOC') in1=2
        do inn1=1,in1
          ind1=ind1+1
          ind2=0
          do ie2=1,ntle(l,isort)
            in2=1
            if(augm(ie2,l,isort)/='LOC') in2=2
            do inn2=1,in2
              ind2=ind2+1
              overlap(ind1, ind2)
     $          =ffsmt(inn1,inn2,ie1,ie2,li,isort,ispin)
            enddo
          enddo
        enddo
      enddo
      smat=overlap


      call symmetriceigen_double(dim_radial_orb,diag,overlap)
      t0=0.0d0
      do jj=1,dim_radial_orb
        t0(:,jj)=overlap(:,jj)*dsqrt(diag(jj)) 
      enddo
      call dgemm('n','c',dim_radial_orb,dim_radial_orb,dim_radial_orb,
     &  1.0d0,t0,dim_radial_orb,overlap,dim_radial_orb,0.0d0,
     &  shalf,dim_radial_orb)

      t0=0.0d0
      do jj=1,dim_radial_orb
        t0(:,jj)=overlap(:,jj)/dsqrt(diag(jj)) 
      enddo
      call dgemm('n','c',dim_radial_orb,dim_radial_orb,dim_radial_orb,
     &  1.0d0,t0,dim_radial_orb,overlap,dim_radial_orb,0.0d0,
     &  smhalf,dim_radial_orb)

      t0=0.0d0
      do jj=1,dim_radial_orb
        t0(:,jj)=overlap(:,jj)/diag(jj)
      enddo
      call dgemm('n','c',dim_radial_orb,dim_radial_orb,dim_radial_orb,
     &  1.0d0,t0,dim_radial_orb,overlap,dim_radial_orb,0.0d0,
     &  sinv,dim_radial_orb)      
      

      allocate(pbound(0:nrad(isort)))
      allocate(qbound(0:nrad(isort)))
      allocate(f_ex(0:nrad(isort)))
      allocate(f_exq(0:nrad(isort)))
      allocate(work(0:nrad(isort)))
      allocate(fr(0:nrad(isort)))
      allocate(v0(0:nrad(isort)))      
      pbound=0.0d0
      qbound=0.0d0
      f_ex=0.0d0
      f_exq=0.0d0
      work=0.0d0
      fr=0.0d0
      v0=0.0d0

      occ=0.0d0
      weight=0.0d0

      do ie=1,ntle(l,isort)
        
        if(augm(ie,l,isort).ne.'LOC') then          
          node_in=int(ptnl(ie,l,isort,ispin))-l-1           
c$$$  do it=1,nit
          
c     -------------- First we solve for Fi -----------------------------         
          e=eny(ie,li,isort,ispin)
c     --------- Solve for the Bound solution (P=0), and (P'=0) ------
          if (ubi_num==1) then
            key11=0
          elseif (ubi_num==4) then
            key11=1
          endif
          nrad_bound
     $      =minloc(dabs(r(0:nrad(isort),isort)-(smt(isort)*mt_fac)))-1
          if (maswrk) then
            write(iun,'(a, 2i5)') 'nrad_bound', iatom, li
            write(iun,'(3(f12.6,2x), i5)') smt(isort)*mt_fac,
     $        smt(isort), r(nrad(isort),isort), nrad_bound(1)
          endif
          pbound=p_f(0:nrad(isort),ie,li,isort,ispin)
          qbound=q_f(0:nrad(isort),ie,li,isort,ispin)
          
          mt0=indmt(1,isort,1)

          mt=indmt(1,isort,ispin)          

          nre=nrad_bound(1)          
          if (key11 .eq. 1) then
            call f_ex_new(f_ex,f_exq,isort,pbound,qbound,ispin,
     $        li,nrad(isort))
            v0(0:nrad(isort))=v_mt_h(mt0:mt0+nrad(isort))
     $        /dsqrt(4.0d0*pi)
          else
            v0(0:nrad(isort))=
     $        (v_mt_h(mt0:mt0+nrad(isort))+v_mt_xc(mt:mt+nrad(isort)))
     $        /dsqrt(4.0d0*pi)                                  
          endif

          pbound=p_f(0:nrad(isort),ie,li,isort,ispin)
          qbound=q_f(0:nrad(isort),ie,li,isort,ispin)          

c$$$          if (maswrk) then
c$$$            write(iun,*)
c$$$            write(iun,*)
c$$$            write(iun,'(2(a, i5, 2x))') 'isort', isort, 'li', li                           
c$$$            write(iun, '(a,6(i5,2x),6(f12.6,2x))')
c$$$     $        'rad_eq in',irel,nrad_bound(1),node_in,key11,l,iival,
c$$$     $        h(isort),z(isort),e,e-3,e+3,clight
c$$$          endif
          
          if (me .eq. 0) then
            write(iun,*) 'input test'
            write(iun,'(6e20.12)') v0(0:5)
            write(iun,'(6e20.12)') r(0:5,isort)
            write(iun,'(6e20.12)') pbound(0:5)
            write(iun,'(6e20.12)') qbound(0:5)
            write(iun,'(6e20.12)') f_ex(0:5)                            
            write(iun,*) nre
            write(iun,'(6e20.12)') f_exq(0:5)
            write(iun,'(6e20.12)') dr(0:5,isort)
            write(iun,'(6e20.12)') work(0:5)                                                        
          endif          

          call rad_eq(e-10,e+10,e,1.d-14,z(isort),l,iival,node_in,0.d0,
     &      -1.d-20,v0(0),fifi,r(0,isort),nrad_bound(1),
     $      pbound(0),qbound(0),
     &      f_ex(0),fi,dfi,kc,nitec,convc,nre,
     &      f_exq(0),dr(0,isort),key11,h(isort),
     &      dl,work(0),acc,node_out,irel,clight)

c$$$          if (maswrk) then          
c$$$            write(iun, '(a,i5,2x,f12.6)')
c$$$     $        'rad_eq inout',nre,e
c$$$            write(iun, '(a,3(i5,2x),6(f12.6,2x))')
c$$$     $        'rad_eq out',node_out,kc,nitec,acc,dl,convc,dfi,fi,fifi
c$$$          endif          
          
c$$$          if (maswrk) then
c$$$            write(iun, '(a,f12.6,2x,7(i5,2x),f12.6)')
c$$$     $        'bound radial function',e,
c$$$     $        isort, li, nrad_bound(1), node_in, node_out, nitec,
c$$$     $        2000,acc
c$$$  endif
          
          if (maswrk) then
            write(iun, '(a,f12.6,2x,7(i5,2x),2(f12.6,2x))')
     $        'bound radial function',e,
     $        isort, li, nrad_bound(1), node_in, node_out, nitec,
     $        5000,acc,convc
          endif          
          
          do ir=nrad_bound(1)+1, nrad(isort)
            pbound(ir)=0.0d0
            qbound(ir)=0.0d0
          enddo
          
          fr=0.0d0
          do ir=0,nrad(isort)
            fr(ir)=pbound(ir)**2
            if(irel.ge.1) fr(ir)=fr(ir)+qbound(ir)**2/c2
            fr(ir)=fr(ir)*dr(ir,isort)
          enddo
          norm=dsqrt(dqdall(h(isort),fr,nrad(isort)))
          
          do ir=1, nrad(isort)
            bound_radfun(ir,1,li,iatom,ispin)
     $        =pbound(ir)/r(ir,isort)/norm
            if (irel .le. 1) then            
              bound_radfun(ir,2,li,iatom,ispin)
     $          =qbound(ir)/r(ir,isort)/norm
            endif
          enddo
          call fit_zero(bound_radfun(0,1,li,iatom,ispin),r(0,isort))
          call fit_zero(bound_radfun(0,2,li,iatom,ispin),r(0,isort))
        endif
      enddo
      
      do ie=1,ntle(l,isort)        
        in=1
        if(augm(ie,l,isort)/='LOC') in=2
        do jn=1,in
          
          fr=0.0d0

          nres=indfun0(jn,ie,li,isort)
          mtw=ind_wf(nres,isort)          
          
          do ir=0,nrad(isort)
            fr(ir)=gfun(mtw+ir,ispin)*bound_radfun(ir,1,li,iatom,ispin)
     $        *r(ir,isort)**2
            if(irel.ge.1) fr(ir)=fr(ir)+gfund(mtw+ir,ispin)
     &        *bound_radfun(ir,2,li,iatom,ispin)*r(ir,isort)**2/c2
            fr(ir)=fr(ir)*dr(ir,isort)
          enddo
          bound_radfun_overlap(jn,ie,li,iatom,ispin)=
     &      dqdall(h(isort),fr,nrad(isort))
        enddo
      enddo

      ind1=0
      do ie1=1,ntle(l,isort)
        in1=1
        if(augm(ie1,l,isort)/='LOC') in1=2
        do inn1=1,in1
          ind1=ind1+1
          
          ind2=0
          do ie2=1,ntle(l,isort)
            in2=1
            if(augm(ie2,l,isort)/='LOC') in2=2
            do inn2=1,in2
              ind2=ind2+1                
              
              bound_radfun_coeff(inn1,ie1,li,iatom,ispin)
     $          =bound_radfun_coeff(inn1,ie1,li,iatom,ispin)
     $          +sinv(ind1,ind2)
     $          *bound_radfun_overlap(inn2,ie2,li,iatom,ispin)
            enddo
          enddo
        enddo
      enddo        
      
      
      if (maswrk) then
        write(iun, '(a, i5, a, i3, a, i5)')
     $    'bound radial function overlap and coeff, iatom',
     $    iatom, 'li', li, 'ispin', ispin
        do ie=1,ntle(l,isort)
          in=1
          if(augm(ie,l,isort)/='LOC') in=2
          do jn=1,in          
            write(iun,'(a,i5,a,i5,a,2f12.6)')
     $        'ie', ie, 'jn', jn, 'overlap', 
     $        bound_radfun_overlap(jn,ie,li,iatom,ispin),
     $        bound_radfun_coeff(jn,ie,li,iatom,ispin)
          enddo
        enddo
      endif

      
      tempemin=dis_froz_min
      tempemax=dis_froz_max
      
      do ind_k=1,ndim_kk(me+1)  ! k vector
        k=n_mpi_kk(me+1)+ind_k
        k0=i_kref(k)
        
        call sym_z_0(znew,k,z_wan_bnd(1,1,k0),
     &    num_bands,k_group(k),pnt(1,k))
        
        do ibnd=1,num_bands            
c$$$  if ((eigenvalues(ibnd,k) .gt. tempemin) .and.
c$$$  $      (eigenvalues(ibnd,k) .lt. tempemax)) then
          radind1=0
          do ie1=1,ntle(l,isort)
            in1=1
            if(augm(ie1,l,isort)/='LOC') in1=2
            do inn1=1,in1
              radind1=radind1+1
              
              radind2=0
              do ie2=1,ntle(l,isort)
                in2=1
                if(augm(ie2,l,isort)/='LOC') in2=2
                do inn2=1,in2
                  radind2=radind2+1
                  
                  if (irel .le. 1) then
                    do m1=-l, l
                      lm1=l*(l+1)+m1+1
                      km1=ind0+indbasa(inn1,ie1,lm1,isort)
                      km2=ind0+indbasa(inn2,ie2,lm1,isort)
                      occ=occ
     $                  +dble(znew(km1,ibnd)*dconjg(znew(km2,ibnd)))
     $                  *bound_radfun_overlap(inn1,ie1,li,iatom,ispin)
     $                  *bound_radfun_overlap(inn2,ie2,li,iatom,ispin)                        
     $                  /dble(nqdiv)
     $                  /(dexp(eigenvalues(ibnd,k)
     $                  /(8.617332262145d-5*0.9d3))
     $                  +1.0d0)
                      
                      if ((eigenvalues(ibnd,k) .gt. tempemin) .and.
     $                  (eigenvalues(ibnd,k) .lt. tempemax)) then
                        weight=weight
     $                    +dble(znew(km1,ibnd)*dconjg(znew(km2,ibnd)))
     $                    *bound_radfun_overlap(inn1,ie1,li,iatom,ispin)
     $                    *bound_radfun_overlap(inn2,ie2,li,iatom,ispin)                        
     $                    /dble(nqdiv)
                      endif
                    enddo
                  else
                    do mj=-jjval, jjval,2
                      call getlimj(lm1,l,iival,mj,li,1)
                      km1=ind0+indbasa(inn1,ie1,lm1,isort)
                      km2=ind0+indbasa(inn2,ie2,lm1,isort)
                      occ=occ
     $                  + dble(znew(km1,ibnd)*dconjg(znew(km2,ibnd)))
     $                  *bound_radfun_overlap(inn1,ie1,li,iatom,ispin)
     $                  *bound_radfun_overlap(inn2,ie2,li,iatom,ispin)                        
     $                  /dble(nqdiv)
     $                  /(exp(eigenvalues(ibnd,k)
     $                  /(8.617332262145d-5*0.9d3))
     $                  +1.0d0)
                      
                      if ((eigenvalues(ibnd,k) .gt. tempemin) .and.
     $                  (eigenvalues(ibnd,k) .lt. tempemax)) then
                        weight=weight
     $                    + dble(znew(km1,ibnd)*dconjg(znew(km2,ibnd)))
     $                    *bound_radfun_overlap(inn1,ie1,li,iatom,ispin)
     $                    *bound_radfun_overlap(inn2,ie2,li,iatom,ispin)                        
     $                    /dble(nqdiv)
                      endif
                    enddo
                  endif
                enddo
              enddo
            enddo
          enddo
c$$$  endif
        enddo
      enddo
      
      
      
      if(nproc/=1) then
        call mpi_allreduce_dble
     $    (weight,1,mpi_sum,mpi_comm_world)
        call mpi_allreduce_dcmplx
     $    (occ,1,mpi_sum,mpi_comm_world)                
      endif
      
      if (maswrk) then
        write(iun, '(3(a, i4, 5x), 100e20.12)')
     $    'occ, weight, iatom', iatom, 'l', l, 'li', li,
     $    occ,weight
        
      endif

      ind1=0
      do ie1=1,ntle(l,isort)
        in1=1
        if(augm(ie1,l,isort)/='LOC') in1=2
        do inn1=1,in1
          ind1=ind1+1          
          coeff_radial_orb(ind1)
     $      =bound_radfun_overlap(inn1,ie1,li,iatom,ispin)
        enddo
      enddo
      
      deallocate(pbound)
      deallocate(qbound)
      deallocate(f_ex)
      deallocate(f_exq)
      deallocate(work)
      deallocate(fr)
      deallocate(v0)        
      deallocate(overlap)      
      deallocate(shalf)
      deallocate(smhalf)
      deallocate(smat)
      deallocate(sinv)                  
      deallocate(t0)
      deallocate(diag)            
      end



      subroutine cal_occ_weight


      use comwann_mod
      implicit none
      include "mpif.h"
      integer :: ind_k,k,k0,ibnd,iwan,iorb,iatom,isort,llval,mjval,
     $  lival,iival,i_index,l,li,jjval,ist
      double precision :: tempemin, tempemax

      tempemin=dis_froz_min
      tempemax=dis_froz_max      

      do ind_k=1,ndim_kk(me+1)  ! k vector
        k=n_mpi_kk(me+1)+ind_k
        k0=i_kref(k)
        
        do ibnd=1,num_bands            
          do iwan=1, num_wann
            do iorb=1, nalimj(iwan)
              iatom=alimj(1,iorb,iwan)
              isort=is(iatom)
              llval=alimj(2,iorb,iwan)
              
              if (irel .le. 1) then
                mjval=alimj(3,iorb,iwan)        
                lival=llval+1

                occ_fwan(llval,1,isort)
     $            =occ_fwan(llval,1,isort)                
     $            +cdabs(v_matrix(ibnd,iwan,k))**2
     $            /dble(nqdiv)
     $            /(dexp(eigenvalues(ibnd,k)
     $            /(8.617332262145d-5*0.9d3))
     $            +1.0d0)
     $            /sort_degen(isort)                
                if ((eigenvalues(ibnd,k) .gt. tempemin) .and.
     $            (eigenvalues(ibnd,k) .lt. tempemax)) then
                  weight_fwan(llval,1,isort)
     $              =weight_fwan(llval,1,isort)                
     $              +cdabs(v_matrix(ibnd,iwan,k))**2
     $              /dble(nqdiv)
     $              /sort_degen(isort)                                  
                endif
              else
                iival=alimj(3,iorb,iwan)
                i_index=(iival+1)/2+1                
                mjval=alimj(4,iorb,iwan)                                        
                
                occ_fwan(llval,i_index,isort)
     $            =occ_fwan(llval,i_index,isort)
     $            +cdabs(v_matrix(ibnd,iwan,k))**2
     $            /dble(nqdiv)
     $            /(dexp(eigenvalues(ibnd,k)
     $            /(8.617332262145d-5*0.9d3))
     $            +1.0d0)
     $            /sort_degen(isort)                                
                if ((eigenvalues(ibnd,k) .gt. tempemin) .and.
     $            (eigenvalues(ibnd,k) .lt. tempemax)) then
                  weight_fwan(llval,i_index,isort)
     $              =weight_fwan(llval,i_index,isort)                
     $              +cdabs(v_matrix(ibnd,iwan,k))**2
     $              /dble(nqdiv)
     $              /sort_degen(isort)                                  
                endif                
              endif
            enddo
          enddo
        enddo
      enddo


      if(nproc/=1) then
        call mpi_allreduce_dble
     $    (occ_fwan, (maxb+1)*nrel*nsort,
     $    mpi_sum,mpi_comm_world)
        call mpi_allreduce_dble
     $    (weight_fwan, (maxb+1)*nrel*nsort,
     $    mpi_sum,mpi_comm_world)        
      endif


      if (maswrk) then
        open(unit=41, file='weight_final.dat')
        if (irel .le. 1) then        
          do isort=1, nsort
            do l=0,lmb(isort)   ! l
              li=l+1            !li              
              write(41, '(2(a, i4, 5x), 100f12.6)')
     $          'isort', isort, 'li', li,
     $          weight_fwan(l,1,isort)
            enddo
          enddo
        else
          do isort=1, nsort
            do l=0,lmb(isort)   ! l
              ist=-1              
              if(l.eq.0) ist=1
              do iival=ist,1,2  ! i
                i_index=(iival+1)/2+1
                li=l*2+(iival+1)/2
                jjval=l+l+iival ! j              
                write(41, '(2(a, i4, 5x), 100f12.6)')
     $            'isort', isort, 'li', li,
     $            weight_fwan(l,i_index,isort)
              enddo
            enddo
          enddo
        endif
        close(41)

        open(unit=41, file='occ_final.dat')
        if (irel .le. 1) then        
          do isort=1, nsort
            do l=0,lmb(isort)   ! l
              li=l+1            !li              
              write(41, '(2(a, i4, 5x), 100f12.6)')
     $          'isort', isort, 'li', li,
     $          occ_fwan(l,1,isort)
            enddo
          enddo
        else
          do isort=1, nsort
            do l=0,lmb(isort)   ! l
              ist=-1              
              if(l.eq.0) ist=1
              do iival=ist,1,2  ! i
                i_index=(iival+1)/2+1
                li=l*2+(iival+1)/2
                jjval=l+l+iival ! j              
                write(41, '(2(a, i4, 5x), 100f12.6)')
     $            'isort', isort, 'li', li,
     $            occ_fwan(l,i_index,isort)
              enddo
            enddo
          enddo
        endif
        close(41)      
      endif
      end            

