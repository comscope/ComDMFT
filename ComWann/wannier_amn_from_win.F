      subroutine wannier_orbital(ispin)
      use comwann_mod
      implicit none
      include "mpif.h"            
      integer, intent(in) :: ispin
      integer ::jj,ii, 
     $  isort,l,dim_radial(0:maxb,nrel,natom),iival,jjval,iatom,
     $  li,ist,mdim,lval,i_index,
     $  is1,is2,mm,sort_degen(nsort),
     $  ind_k,k,k0,ibnd,ind0,llval,mjval,lival,lldim,ie1,inn1,in1,
     $  radind1,m1,lm1,km1,mjind,jjdim,irk,iomega
      double precision :: weight, weight_in_froz_win(0:maxb,nrel,nsort),
     $  weight_low(0:maxb,nrel,nsort),
     $  weight_total(0:maxb,nrel,nsort),
     $  llval_dble,jjval_dble,omega
c$$$     $  ,pdos(-500:500,0:maxb,nrel,nsort)
      complex*16 :: znew(nfun,num_bands),
     $  coeff_radial(100,0:maxb,nrel,natom),
     $  proj(0:maxb,nrel,2*maxb+2,natom,num_bands,ndim_kk(me+1))

      complex*16, allocatable :: wfn_mt0(:,:), wfn_mt1(:,:),
     $  rot_matrix(:,:)

!     radial function determination
c$$$  dimension determination`


      weight_in_froz_win=0.0d0

      sort_degen=0

      do iatom=1, natom
        isort=is(iatom)
        sort_degen(isort)=sort_degen(isort)+1
      enddo

      dim_radial=0
      coeff_radial=0.0d0
      
      if (irel.le.1) then
        do iatom=1, natom
          isort=is(iatom)
          do l=0,lmb(isort)     ! l
            li=l+1              !li
            call cal_overlapmax_bandprojection_to_mtorb_low
     $        (iatom,l,iival,jjval,li,ispin,
     $        dim_radial(l,1,iatom),coeff_radial(1,l,1,iatom),weight)
            weight_in_froz_win(l,1,isort)
     $        =weight_in_froz_win(l,1,isort)
     $        +weight/sort_degen(isort)
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
              call cal_overlapmax_bandprojection_to_mtorb_low
     $          (iatom,l,iival,jjval,li,ispin,
     $          dim_radial(l,i_index,iatom),
     $          coeff_radial(1,l,i_index,iatom),weight)
              weight_in_froz_win(l,i_index,isort)
     $          =weight_in_froz_win(l,i_index,isort)
     $          +weight/sort_degen(isort)
            enddo
          enddo
        enddo
      endif

c$$$      if (maswrk) then
c$$$        write(iun,*) 'dim_radial', dim_radial
c$$$      endif
      
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
c$$$                          write(179,'(6i5)') llval, iival,
c$$$     $                      mjind, radind1,jjdim,dim_radial
c$$$     $                      (llval,i_index,iatom)
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


c$$$      write(179,*) 'sum(proj)', sum(proj)
      
      weight_low=0.0d0
      weight_total=0.0d0      
c$$$      pdos=0.0d0
      
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
c$$$                  do iomega=-500,500
c$$$                    omega=iomega*0.1d0
c$$$                    pdos(iomega,llval,1,isort)
c$$$     $                =pdos(iomega,llval,1,isort)
c$$$     $                +1.0d0/dsqrt(2.0d0*pi)/0.1d0
c$$$     $                *dexp(-(omega-eigenvalues(ibnd,irk))**2/0.02d0)
c$$$     $                *cdabs(proj(llval,1,mjind,iatom,ibnd,ind_k))**2
c$$$     $                /dble(nqdiv)/sort_degen(isort)
c$$$                  enddo
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
c$$$                write(179,'(4i5)') iatom,llval,iival,mjval
                
                do ind_k=1,ndim_kk(me+1) ! k vector
                  irk=n_mpi_kk(me+1)+ind_k
                  k0=i_kref(irk)                                
                  do ibnd=1, num_bands
c$$$                    do iomega=-500,500
c$$$                      omega=iomega*0.1d0
c$$$                      pdos(iomega,llval,i_index,isort)
c$$$     $                  =pdos(iomega,llval,i_index,isort)
c$$$     $                  +1.0d0/dsqrt(2.0d0*pi)/0.1d0
c$$$     $                  *dexp(-(omega-eigenvalues(ibnd,irk))**2/0.02d0)
c$$$     $                  *cdabs(proj(llval,i_index,mjind,
c$$$     $                  iatom,ibnd,ind_k))**2
c$$$     $                  /dble(nqdiv)/sort_degen(isort)
c$$$                    enddo
                    weight_total(llval,i_index,isort)
     $                =weight_total(llval,i_index,isort)
     $                  +cdabs(proj(llval,i_index,mjind,
     $                  iatom,ibnd,ind_k))**2
     $                  /dble(nqdiv)/sort_degen(isort)                    
                    if (eigenvalues(ibnd,irk) .lt. dis_froz_min) then
c$$$                      write(179,'(6i5)')
c$$$     $                  iatom,llval,iival,mjval,irk,ibnd                      
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

c$$$      write(179,*) 'weight_low', sum(weight_low)      


      if(nproc/=1) then
        call mpi_allreduce_dble
     $    (weight_low(0,1,1),(maxb+1)*nrel*nsort,
     $    mpi_sum,mpi_comm_world)
        call mpi_allreduce_dble
     $    (weight_total(0,1,1),(maxb+1)*nrel*nsort,
     $    mpi_sum,mpi_comm_world)        
c$$$        call mpi_allreduce_dble
c$$$     $    (pdos(-500,0,1,1),1001*(maxb+1)*nrel*nsort,
c$$$     $    mpi_sum,mpi_comm_world)        
      endif

c$$$      write(179,*) 'weight_low', sum(weight_low)            
      
      if (maswrk) then
c$$$        open(unit=100, file='pdos.out')
c$$$        do iomega=-500,500
c$$$          omega=iomega*0.1d0
c$$$          write(100,'(10000f12.6)') omega,
c$$$     $      (((pdos(iomega,llval,i_index,isort),i_index=1, nrel),
c$$$     $      llval=0,maxb), isort=1, nsort)
c$$$        enddo
c$$$        close(100)
          
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
        num_wann=0
        open(unit=42, file='orb_for_froz_win.dat')
        write(42,*) 'nbnd_in_froz', nbnd_in_froz
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

      integer :: ik,tempint,ibnd

      nbnd_in_froz=0        
      do ik=1, nqdiv
        tempint=0
        do ibnd=1, num_bands
          if (
     $      (eigenvalues(ibnd,ik) .ge. dis_froz_min)
     $      .and. 
     $      (eigenvalues(ibnd,ik) .le. dis_froz_max)
     $      ) then
            tempint=tempint+1
          endif
        enddo
        if (tempint .gt. nbnd_in_froz) then
          nbnd_in_froz=tempint
        endif
      enddo
      if (maswrk) then
        write(iun, *) 'nbnd_in_froz', nbnd_in_froz
      endif
      
      end      
      
      

      subroutine poorman_wan_proj(ispin)
      use comwann_mod
      implicit none
      include "mpif.h"      
      integer, intent(in) :: ispin
      integer :: ie1,inn1,in1,
     $  ii,jj,k,k0,ibnd,radind1,m1,lm1,
     $  km1,isort,iwan,dim_radial(max_nalimj,num_wann),
     $  iival,jjval,llval,mjval,lival,iatom,ind0,
     $  ist,ind_k,km,irk,iband,iorb, lldim,jjdim,mjind
      integer, allocatable :: ipiv(:)
      double precision :: llval_dble,jjval_dble,weight
      complex*16 :: znew(nfun,num_bands),
     $  coeff_radial(100,max_nalimj,num_wann),
     $  proj_poor_old(num_wann,num_bands),
     $  proj_poor_new(num_wann,num_bands,nqdiv)
      complex*16, allocatable :: wfn_mt0(:,:), wfn_mt1(:,:),
     $  rot_matrix(:,:)

!     radial function determination
c$$$  dimension determination`
      dim_radial=0
      coeff_radial=0.0d0
      proj_poor_new=0.0d0
      proj_poor_old=0.0d0

      do iwan=1, num_wann
        
        do iorb=1, nalimj(iwan)
          
          iatom=alimj(1,iorb,iwan)
          isort=is(iatom)
          ind0=io_lem(iatom)-1
          llval=alimj(2,iorb,iwan)
          
          if (irel .le. 1) then
            mjval=alimj(3,iorb,iwan)        
            lival=llval+1
          else
            iival=alimj(3,iorb,iwan)
            mjval=alimj(4,iorb,iwan)
            jjval=llval+llval+iival ! j
            lival=llval+llval+(iival+1)/2
          endif

          call cal_overlapmax_bandprojection_to_mtorb_low(iatom,llval,
     $      iival,jjval,lival,ispin,dim_radial(iorb,iwan),
     $      coeff_radial(1,iorb,iwan),weight)
        enddo
      enddo


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
                
                allocate(rot_matrix(2*llval+1,2*llval+1))
                rot_matrix=0.0d0

                llval_dble=dble(llval)
                lldim=2*llval+1
                call rotmat_wan_proj(alimj_xaxis(1,iorb,iwan),
     $            alimj_zaxis(1,iorb,iwan),llval_dble,lldim,1,
     $            rot_matrix)
                allocate(wfn_mt0(2*llval+1,dim_radial(iorb,iwan)))
                wfn_mt0=0.0d0
                allocate(wfn_mt1(2*llval+1,dim_radial(iorb,iwan)))
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
                do radind1=1, dim_radial(iorb,iwan)
                  proj_poor_old(iwan,ibnd)
     $              =proj_poor_old(iwan,ibnd)
     $              +alimj_coeff(iorb,iwan)
     $              *dconjg(coeff_radial(radind1,iorb,iwan))
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


                allocate(rot_matrix(jjdim,jjdim))
                rot_matrix=0.0d0
                allocate(wfn_mt0(jjdim,dim_radial(iorb,iwan)))
                wfn_mt0=0.0d0
                allocate(wfn_mt1(jjdim,dim_radial(iorb,iwan)))
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
                
                do radind1=1, dim_radial(iorb,iwan)
                  mjind=(mjval+jjval+2)/2                  
                  proj_poor_old(iwan,ibnd)
     $              =proj_poor_old(iwan,ibnd)
     $              +alimj_coeff(iorb,iwan)
     $              *dconjg(coeff_radial(radind1,iorb,iwan))
     $              *wfn_mt1(mjind,radind1)
                enddo                


                deallocate(rot_matrix)
                deallocate(wfn_mt0)
                deallocate(wfn_mt1)
                
              endif
            enddo
          enddo
        enddo
c        write(179,'(a, 2i5)')
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
c      write(179,'(a, f12.6, a, f12.6)')
c     $  'maxval', diag(n1), '     minval', diag(1)
c      write(179,'(10000f12.6)') diag
      
      
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


      subroutine cal_overlapmax_bandprojection_to_mtorb_low(iatom,l,
     $  iival,jjval,li,ispin,dim_radial,coeff_radial,weight)
      
      use comwann_mod
      implicit none
      include "mpif.h"
      integer, intent(in) :: iatom,l,iival,jjval,li,ispin ! if irel.le.1 then iival and jjval is dummy variable

      integer, intent(out) :: dim_radial
      double precision, intent(out) :: weight
      complex*16, intent(out) :: coeff_radial(100)
      
      integer :: ie1,inn1,inn2,ie2,in2,in1,ind1,ind2,
     $  ii,jj,k,k0,ibnd,radind1,radind2,radind3,radind4,m1,lm1,
     $  km1,km2,isort,ind0,ind_k,
     $  mj,ival
      double precision :: tempemin, tempemax            
      complex*16 :: znew(nfun,num_bands)
      double precision, allocatable :: overlap(:,:),t0(:,:),sxmat(:,:),
     $  partial_proj_eigval(:),diag(:)
      complex*16, allocatable :: partial_proj(:,:)

!     radial function determination
c$$$  dimension determination`


      coeff_radial=0.0d0

      isort=is(iatom)
      ind0=io_lem(iatom)-1
      
      dim_radial=0
      do ie1=1,ntle(l,isort)
        in1=1
        if(augm(ie1,l,isort)/='LOC') in1=2
        do inn1=1,in1
          dim_radial=dim_radial+1
        enddo
      enddo

      
      allocate(overlap(dim_radial,dim_radial))
      overlap=0.0d0
      allocate(t0(dim_radial,dim_radial))
      t0=0.0d0
      allocate(sxmat(dim_radial,dim_radial))
      sxmat=0.0d0            
      allocate(partial_proj(dim_radial,dim_radial))
      partial_proj=0.0d0
      allocate(partial_proj_eigval(dim_radial))
      partial_proj_eigval=0.0d0
      allocate(diag(dim_radial))
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


!     symmetric orthonormalization with neglection of linear contribution
      call symmetriceigen_double(dim_radial,diag,overlap)
      do jj=1,dim_radial
        t0(:,jj)=overlap(:,jj)*dsqrt(diag(jj)) ! to calculate S*X
      enddo
      call dgemm('n','c',dim_radial,dim_radial,dim_radial,
     &  1.0d0,t0,dim_radial,overlap,dim_radial,0.0d0,
     &  sxmat,dim_radial)      
      
      partial_proj=0.0d0

      tempemin=dis_froz_min
      tempemax=dis_froz_max

      
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
                        do radind3=1,dim_radial
                          do radind4=1,dim_radial
                            partial_proj(radind3,radind4)=
     $                        partial_proj(radind3,radind4)+
     $                        znew(km1,ibnd)*dconjg(znew(km2,ibnd))
     $                        *sxmat(radind2,radind4)
     $                        *sxmat(radind1,radind3)
     $                        /dble(nqdiv)
                          enddo
                        enddo
                      enddo
                    else
                      do mj=-jjval, jjval,2
                        call getlimj(lm1,l,iival,mj,li,1)
                        km1=ind0+indbasa(inn1,ie1,lm1,isort)
                        km2=ind0+indbasa(inn2,ie2,lm1,isort)
                        do radind3=1,dim_radial
                          do radind4=1,dim_radial
                            partial_proj(radind3,radind4)=
     $                        partial_proj(radind3,radind4)+
     $                        znew(km1,ibnd)*dconjg(znew(km2,ibnd))
     $                        *sxmat(radind2,radind4)
     $                        *sxmat(radind1,radind3)
     $                        /dble(nqdiv)
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
     $    (partial_proj,dim_radial**2,mpi_sum,mpi_comm_world)        
      endif


      call hermitianeigen_cmplxdouble
     $  (dim_radial,partial_proj_eigval,partial_proj)

      
      coeff_radial=0.0d0
      do ii=1, dim_radial
        do jj=1, dim_radial
          coeff_radial(ii)=coeff_radial(ii)+
     $      sxmat(ii,jj)
     $      *partial_proj(jj,dim_radial)
        enddo
      enddo
      weight=partial_proj_eigval(dim_radial)

c$$$  if (printflag .eq. 1) then
c$$$  if (irel .le. 1) then
c$$$  ival=1
c$$$  else
c$$$  ival=(iival+1)/2+1
c$$$  endif
c$$$  weight_in_froz_win(l,ival,iatom)
c$$$  $    =weight_in_froz_win(l,ival,iatom)
c$$$  $    +partial_proj_eigval(dim_radial)
c$$$  if (maswrk) then
c$$$  if (irel .le. 1) then
c$$$  write(iun, '(4(a, i4, 5x), 100f12.6)')
c$$$  $        'iatom', iatom, 'l', l, 'mu_dim', 2*l+1, 'li', li,
c$$$  $        partial_proj_eigval(dim_radial)          
c$$$  else
c$$$  write(iun, '(5(a, i4, 5x), 100f12.6)')
c$$$  $        'iatom', iatom, 'l', l, 'i',
c$$$  $        iival,'mu_dim', jjval+1, 'li', li,
c$$$  $        partial_proj_eigval(dim_radial)
c$$$  endif
c$$$  endif
c$$$  endif
      
      deallocate(overlap)      
      deallocate(sxmat)
      deallocate(t0)
      deallocate(diag)            
      deallocate(partial_proj)
      deallocate(partial_proj_eigval)

      end


