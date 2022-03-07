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

