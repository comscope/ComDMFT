      subroutine cor_and_iter_wan
      use atom_mod
      use manager_mod
      use solid_mod
      use parallel_mod
      use vertex_mod
      use wannier90_mod
      implicit none
      include 'mpif.h'
      integer :: iatom,isort,l1,i1,ii,i,m,lm,n,lm0,ndim,li,i2,ii2,
     $  jj,in,i3,lf,iatom_c,i0,kk,l,ind,ind_c,ie1,in1,inn1,inttmp,
     $  awan,iistart,iorb,ival,lival,iwan,jval,ll,mm,mudim,mval,
     $  lval,is_bonding,cnt,iat,mind,m1,jatom
      integer, allocatable :: tempmap(:,:,:), tempind(:)
      integer, external :: find_ind

!     write(100,*)'alimj',alimj,'at_bonding',at_bonding,'at_orb_ciw',
!     $ at_orb_ciw,'ciw_ind',ciw_ind,'corr_iter',corr_iter,
!     $ 'corr_iter_rel',corr_iter_rel,'iatom_ciw',iatom_ciw,'irel',irel,
!     $ 'is',is,'lmb',lmb,'maxb',maxb,'mt2ciw',mt2ciw,'mt2ciw_rel_l',
!     $ mt2ciw_rel_l,'mt2allwan',mt2allwan,'n_ciw_atom',n_ciw_atom,'n_ciw_max',      
!     $ n_ciw_max,'nalimj',nalimj,'natom',natom,'natom_ciw',natom_ciw,
!     $ 'natom_cor',natom_cor,'natom_iter',natom_iter,'nrel',nrel,
!     $ 'num_orb_cor',num_orb_cor,'num_orb_iter',num_orb_iter,
!     $ 'num_orb_proj',num_orb_proj,'seed_name',seed_name,
!     $ 'wan_correlated',wan_correlated,'wan_itinerant',wan_itinerant
!     
!     call mpi_bcast(nproc_k,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)

      write(seed_name,'(a)') 'wannier'      
      call read_wann_inip



!     call mpi_bcast(nproc_k,1,mpi_integer,0,MPI_COMM_WORLD,mpierr)


     

      allocate(natom_cor(natom))
      natom_cor=.false.
      allocate(natom_iter(natom))
      natom_iter=.false.      
      allocate(n_ciw_atom(natom))
      n_ciw_atom=0
      allocate(n_cw_atom(natom))
      n_cw_atom=0      

      
      allocate(mt2allwan(2*maxb+nrel,nrel*maxb+1,natom)) ! m, li, atom
      mt2allwan=0

      allocate(mt2corwan(2*maxb+nrel,nrel*maxb+1,natom)) ! m, li, atom
      mt2corwan=0      

      
      allocate(mt2ciw(2*maxb+nrel,nrel*maxb+1,natom)) ! m, li, atom
      mt2ciw=0
      allocate(mt2cw(2*maxb+nrel,nrel*maxb+1,natom)) ! m, li, atom
      mt2cw=0            

      if (irel .eq. 2) then
        allocate(mt2ciw_rel_l(4*maxb+2,maxb+1,natom)) ! m, li, atom
        mt2ciw_rel_l=0
        allocate(mt2cw_rel_l(4*maxb+2,maxb+1,natom)) ! m, li, atom
        mt2cw_rel_l=0                
      endif
      
      allocate(corr_iter(0:maxb,natom)) ! l, atom
      corr_iter=0

      allocate(corr(0:maxb,natom)) ! l, atom
      corr=0            
      allocate(at_bonding(0:maxb,natom)) ! l, atom
      at_bonding=0      
      if (irel .eq. 2) then
        allocate(corr_iter_rel(nrel*maxb+1,natom)) ! m, li, atom
        corr_iter_rel=0
        allocate(corr_rel(nrel*maxb+1,natom)) ! m, li, atom
        corr_rel=0                
        allocate(at_bonding_rel(nrel*maxb+1,natom)) ! m, li, atom
        at_bonding_rel=0        
      endif




      

      do awan=1, num_orb_cor
        iwan=wan_correlated(awan)
        do iorb=1, nalimj(iwan)
          if (irel .le. 1) then
            iatom=alimj(1,iorb,iwan)
            lval=alimj(2,iorb,iwan)
            lival=lval+1
            mval=alimj(3,iorb,iwan) ! 1,2,3,4,5
          else
            iatom=alimj(1,iorb,iwan)
            lval=alimj(2,iorb,iwan)
            ival=alimj(3,iorb,iwan)           
            jval=lval+lval+ival ! j
            lival=lval+lval+(ival+1)/2
            mval=((alimj(4,iorb,iwan))+lival)/2+1 !from -5, -3, -1, 1, 3, 5 to 1, 2, 3, 4, 5, 6          
          endif
          natom_cor(iatom)=.true.
          mt2allwan(mval,lival,iatom)=iwan
          mt2corwan(mval,lival,iatom)=iwan                    
          mt2ciw(mval,lival,iatom)=awan
          mt2cw(mval,lival,iatom)=awan                    
          if (irel .eq. 2) then
            if (lval .eq. 0) then
              mind=mval
            else
              mind=mval+(ival+1)/2*(2*lval)
            endif
            mt2ciw_rel_l(mind,lval+1,iatom)=awan
            mt2cw_rel_l(mind,lval+1,iatom)=awan                        
          endif
        enddo
        n_ciw_atom(iatom)=n_ciw_atom(iatom)+1
        n_cw_atom(iatom)=n_cw_atom(iatom)+1                
      enddo




      do awan=1, num_orb_iter
        iwan=wan_itinerant(awan)
        do iorb=1, nalimj(iwan)
          if (irel .le. 1) then
            iatom=alimj(1,iorb,iwan)
            lval=alimj(2,iorb,iwan)
            lival=lval+1
            mval=alimj(3,iorb,iwan) ! 1,2,3,4,5
          else
            iatom=alimj(1,iorb,iwan)
            lval=alimj(2,iorb,iwan)
            ival=alimj(3,iorb,iwan)           
            jval=lval+lval+ival ! j
            lival=lval+lval+(ival+1)/2
            mval=((alimj(4,iorb,iwan))+lival)/2+1 !from -5, -3, -1, 1, 3, 5 to 1, 2, 3, 4, 5, 6          
          endif
          natom_iter(iatom)=.true.
          mt2allwan(mval,lival,iatom)=iwan          
          mt2ciw(mval,lival,iatom)=awan+num_orb_cor
          if (irel .eq. 2) then
            if (lval .eq. 0) then
              mind=mval
            else
              mind=mval+(ival+1)/2*(2*lval)
            endif
            mt2ciw_rel_l(mind,lval+1,iatom)=awan+num_orb_cor
          endif          
        enddo
        n_ciw_atom(iatom)=n_ciw_atom(iatom)+1        
      enddo


      do iatom=1, natom
        isort=is(iatom)
        do ll=0, lmb(isort)
          if (irel .eq. 2) then
            iistart=-1
            if (ll .eq. 0) iistart=1
            do ii=iistart, 1, 2
              li=2*ll+(ii+1)/2
              mudim=2*ll+ii+1
              inttmp=0
              is_bonding=0
              do mm=1, mudim
                if (mt2allwan(mm,li,iatom) .ne. 0) then
                  inttmp=inttmp+1
                  if (nalimj(mt2allwan(mm,li,iatom)) .ne. 1) then
                    is_bonding=is_bonding+1
                  endif
                endif
              enddo
              if (inttmp .eq. mudim) then
                corr_iter_rel(li,iatom)=2
              elseif (inttmp .ne. 0) then
                corr_iter_rel(li,iatom)=1
              endif
              if (is_bonding.eq.0) then
                at_bonding_rel(li,iatom)=0
              else
                at_bonding_rel(li,iatom)=1
              endif
            enddo
            if (ll .eq. 0) then
              corr_iter(0,iatom)=corr_iter_rel(1,iatom)
              at_bonding(0,iatom)=at_bonding_rel(1,iatom)              
            else
              if ((corr_iter_rel(2*ll,iatom).eq.2)
     $          .and.
     $          (corr_iter_rel(2*ll+1,iatom).eq.2)) then
                corr_iter(ll,iatom)=2
              elseif ((corr_iter_rel(2*ll,iatom).eq.0)
     $            .and.
     $            (corr_iter_rel(2*ll+1,iatom).eq.0)) then
                corr_iter(ll,iatom)=0
              else
                corr_iter(ll,iatom)=1
              endif
              if ((at_bonding_rel(2*ll,iatom).eq.0)
     $          .and.
     $          (at_bonding_rel(2*ll+1,iatom).eq.0)) then
                at_bonding(ll,iatom)=0
              else
                at_bonding(ll,iatom)=1
              endif                          
            endif
          else
            li=ll+1
            mudim=2*ll+1
            inttmp=0
            is_bonding=0            
            do mm=1, mudim
              if (mt2allwan(mm,li,iatom) .ne. 0) then
                inttmp=inttmp+1
                if (nalimj(mt2allwan(mm,li,iatom)) .ne. 1) then
                  is_bonding=is_bonding+1
                endif                
              endif
            enddo
            if (inttmp .eq. mudim) then
              corr_iter(ll,iatom)=2
            elseif (inttmp .ne. 0) then
              corr_iter(ll,iatom)=1
            endif
            if (is_bonding .eq.0) then
              at_bonding(ll,iatom)=0
            else
              at_bonding(ll,iatom)=1
            endif            
          endif
        enddo
      enddo


      do iatom=1, natom
        isort=is(iatom)
        do ll=0, lmb(isort)
          if (irel .eq. 2) then
            iistart=-1
            if (ll .eq. 0) iistart=1
            do ii=iistart, 1, 2
              li=2*ll+(ii+1)/2
              mudim=2*ll+ii+1
              inttmp=0
              is_bonding=0
              do mm=1, mudim
                if (mt2corwan(mm,li,iatom) .ne. 0) then
                  inttmp=inttmp+1
                  if (nalimj(mt2corwan(mm,li,iatom)) .ne. 1) then
                    is_bonding=is_bonding+1
                  endif
                endif
              enddo
              if (inttmp .eq. mudim) then
                corr_rel(li,iatom)=2
              elseif (inttmp .ne. 0) then
                corr_rel(li,iatom)=1
              endif
              if (is_bonding.eq.0) then
                at_bonding_rel(li,iatom)=0
              else
                at_bonding_rel(li,iatom)=1
              endif
            enddo
            if (ll .eq. 0) then
              corr_rel(0,iatom)=corr_rel(1,iatom)
              at_bonding(0,iatom)=at_bonding_rel(1,iatom)              
            else
              if ((corr_rel(2*ll,iatom).eq.2)
     $          .and.
     $          (corr_rel(2*ll+1,iatom).eq.2)) then
                corr(ll,iatom)=2
              elseif ((corr_rel(2*ll,iatom).eq.0)
     $            .and.
     $            (corr_rel(2*ll+1,iatom).eq.0)) then
                corr(ll,iatom)=0                  
              else
                corr(ll,iatom)=1
              endif
              if ((at_bonding_rel(2*ll,iatom).eq.0)
     $          .and.
     $          (at_bonding_rel(2*ll+1,iatom).eq.0)) then
                at_bonding(ll,iatom)=0
              else
                at_bonding(ll,iatom)=1
              endif                          
            endif
          else
            li=ll+1
            mudim=2*ll+1
            inttmp=0
            is_bonding=0            
            do mm=1, mudim
              if (mt2corwan(mm,li,iatom) .ne. 0) then
                inttmp=inttmp+1
                if (nalimj(mt2corwan(mm,li,iatom)) .ne. 1) then
                  is_bonding=is_bonding+1
                endif                
              endif
            enddo
            if (inttmp .eq. mudim) then
              corr(ll,iatom)=2                
            elseif (inttmp .ne. 0) then
              corr(ll,iatom)=1                
            endif
            if (is_bonding .eq.0) then
              at_bonding(ll,iatom)=0
            else
              at_bonding(ll,iatom)=1
            endif            
          endif
        enddo
      enddo      
      
      
      

      natom_ciw=0
      natom_cw=0
      do iatom=1, natom
        if (n_ciw_atom(iatom) .gt. 0) then
          natom_ciw=natom_ciw+1
        endif
        if (n_cw_atom(iatom) .gt. 0) then
          natom_cw=natom_cw+1
        endif        
      enddo
      allocate(iatom_ciw(natom_ciw))
      iatom_ciw=0
      cnt=1
      do iatom=1, natom
        if (n_ciw_atom(iatom) .gt. 0) then
          iatom_ciw(cnt)=iatom
          cnt=cnt+1
        endif
      enddo



      allocate(iatom_cw(natom_cw))
      iatom_cw=0
      cnt=1
      do iatom=1, natom
        if (n_cw_atom(iatom) .gt. 0) then
          iatom_cw(cnt)=iatom
          cnt=cnt+1
        endif
      enddo            


c$$$  allocate(is_ciw(natom))
c$$$  is_ciw=0
c$$$  
c$$$  do iat=1, natom_ciw
c$$$  is_ciw(iatom_ciw(iat))=iat
c$$$  enddo      

      n_ciw_max=maxval(n_ciw_atom)
      n_cw_max=maxval(n_cw_atom)      
      allocate(ciw_ind(n_ciw_max,natom_ciw))
      ciw_ind=0
      allocate(tempind(natom_ciw))
      tempind=0

      

      do iwan=1, num_orb_cor
        iatom=alimj(1,1,wan_correlated(iwan))
        iat=find_ind(iatom, natom_ciw, iatom_ciw)
        tempind(iat)=tempind(iat)+1
        ciw_ind(tempind(iat),iat)=iwan
      enddo

      if (num_orb_iter .gt. 0) then
        do iwan=1, num_orb_iter
          iatom=alimj(1,1,wan_itinerant(iwan))
          iat=find_ind(iatom, natom_ciw, iatom_ciw)
          tempind(iat)=tempind(iat)+1
          ciw_ind(tempind(iat),iat)=iwan+num_orb_cor
        enddo
      endif
      deallocate(tempind)


      allocate(cw_ind(n_cw_max,natom_cw))
      cw_ind=0

      allocate(tempind(natom_cw))
      tempind=0
      do iwan=1, num_orb_cor
        iatom=alimj(1,1,wan_correlated(iwan))
        iat=find_ind(iatom, natom_cw, iatom_cw)
        tempind(iat)=tempind(iat)+1
        cw_ind(tempind(iat),iat)=iwan
      enddo
      

      deallocate(tempind)

      

      allocate(at_orb_ciw(2,num_orb_proj)) ! inverse of ciw_ind
      at_orb_ciw=0.0d0
      do iwan=1, num_orb_cor
        iatom=alimj(1,1,wan_correlated(iwan))
        iat=find_ind(iatom, natom_ciw, iatom_ciw)        
        at_orb_ciw(1,iwan)=iat
        iat=find_ind(iatom, natom_ciw, iatom_ciw)        
        iorb=find_ind(iwan, n_ciw_max, ciw_ind(1,iat))        
        at_orb_ciw(2,iwan)=iorb        
      enddo

      if (num_orb_iter .gt. 0) then
        do iwan=1, num_orb_iter
          iatom=alimj(1,1,wan_itinerant(iwan))
          iat=find_ind(iatom, natom_ciw, iatom_ciw)                  
          at_orb_ciw(1,iwan+num_orb_cor)=iat
          iat=find_ind(iatom, natom_ciw, iatom_ciw)          
          iorb=find_ind(iwan+num_orb_cor, n_ciw_max, ciw_ind(1,iat))          
          at_orb_ciw(2,iwan+num_orb_cor)=iorb        
        enddo        
      endif

      allocate(at_orb_cw(2,num_orb_cor)) ! inverse of ciw_ind
      at_orb_cw=0.0d0
      do iwan=1, num_orb_cor
        iatom=alimj(1,1,wan_correlated(iwan))
        at_orb_cw(1,iwan)=iatom
        iat=find_ind(iatom, natom_cw, iatom_cw)                  
        iorb=find_ind(iwan, n_cw_max, cw_ind(1,iat))
        at_orb_cw(2,iwan)=iorb        
      enddo


      allocate(corwan_ind_atom(num_orb_cor,natom_cw)) ! l, atom      

      corwan_ind_atom=0

      do iat=1, natom_cw
        jatom= iatom_cw(iat)
        ii=0
        do awan=1, num_orb_cor
          iwan=wan_correlated(awan)
          do iorb=1, nalimj(iwan)
            iatom=alimj(1,iorb,iwan)

            if (iatom .eq. jatom) then
              ii=ii+1              
              corwan_ind_atom(ii, iat)=awan
            endif
          enddo
        enddo
      enddo
        
      
      


      if (maswrk) then
        write(iun,*) 'for correlated and itinerate orbitals'

        write(iun, *) 'corwan_ind_atom'
        do iat=1,natom_cw
          iatom=iatom_cw(iat)          
          write(iun, '(1000i5)') iatom, (corwan_ind_atom(ii, iat),
     $      ii=1, n_cw_atom(iatom))
        enddo

        if (irel .le.1) then
          
          write(iun,*) 'mt2allwan'
          do iatom=1, natom
            isort=is(iatom)            
            do li=1, lmb(isort)+1
              write(iun,'(2i4,4x,1000i4)') iatom,li,
     $          (mt2allwan(mm,li,iatom),mm=1,2*li-1)
            enddo
          enddo
          write(iun,*) 'mt2ciw'
          do iatom=1, natom
            isort=is(iatom)            
            do li=1, lmb(isort)+1
              write(iun,'(2i4,4x,1000i4)') iatom,li,
     $          (mt2ciw(mm,li,iatom),mm=1,2*li-1)
            enddo
          enddo
          write(iun,*) 'corr_iter'
          do iatom=1, natom
            isort=is(iatom)            
!            write(iun,'(i4,4x,1000(2i4,2x))') iatom,
!     $        ((corr_iter(ll,iatom), at_bonding(ll,iatom)),
!     $        ll=0, lmb(isort))
          enddo
        else
          write(iun,*) 'mt2allwan'
          do iatom=1, natom
            isort=is(iatom)            
            do li=1, 2*lmb(isort)+1
              write(iun,'(2i4,4x,1000i4)') iatom,li,
     $          (mt2allwan(mm,li,iatom),mm=1,((li+1)/2)*2)
            enddo
          enddo
          write(iun,*) 'mt2ciw'
          do iatom=1, natom
            isort=is(iatom)            
            do li=1, 2*lmb(isort)+1
              write(iun,'(2i4,4x,1000i4)') iatom,li,
     $          (mt2ciw(mm,li,iatom),mm=1,((li+1)/2)*2)
            enddo
          enddo

          write(iun,*) 'mt2ciw_rel_l'
          do iatom=1, natom
            isort=is(iatom)            
            do li=1, lmb(isort)+1
              write(iun,'(2i4,4x,1000i4)') iatom,li,
     $          (mt2ciw_rel_l(mm,li,iatom),mm=1,4*lmb(isort)+2)
            enddo
          enddo          
          write(iun,*) 'corr_iter_rel'
          do iatom=1, natom
            isort=is(iatom)            
!            write(iun,'(i4,4x,1000(2i4,2x))') iatom,
!     $        ((corr_iter_rel(li,iatom),at_bonding_rel(li,iatom)),
!     $        li=1, 2*lmb(isort)+1)
          enddo
          write(iun,*) 'corr_iter'
          do iatom=1, natom
            isort=is(iatom)            
!            write(iun,'(i4,4x,1000(2i4,2x))') iatom,
!     $        ((corr_iter(ll,iatom),at_bonding(ll,iatom)),
!     $        ll=0, lmb(isort))
          enddo

        endif
        write(iun,*) 'natom_ciw', natom_ciw
        write(iun,*) 'n_ciw_max', n_ciw_max
        write(iun,*) 'n_ciw_atom'
        do iatom=1, natom
          write(iun,'(2i4)') iatom,n_ciw_atom(iatom)
        enddo        
        write(iun,*) 'iatom_ciw'
        do iat=1, natom_ciw
          write(iun,'(2i4)') iat,iatom_ciw(iat)
        enddo

        write(iun,*) 'ciw_ind'
        do iat=1, natom_ciw
          iatom=iatom_ciw(iat)
          write(iun,'(i4,4x,1000i4)') iat,(ciw_ind(ii,iat),
     $      ii=1, n_ciw_atom(iatom))
        enddo

        write(iun,*) 'at_orb_ciw'
        do iwan=1,num_orb_proj
          write(iun,'(3i4)') iwan,
     $      (at_orb_ciw(ii, iwan), ii=1, 2)
        enddo


        write(iun,*) 'for correlated orbitals'        

        if (irel .le.1) then

          write(iun,*) 'mt2corwan'
          do iatom=1, natom
            isort=is(iatom)            
            do li=1, lmb(isort)+1
              write(iun,'(2i4,4x,1000i4)') iatom,li,
     $          (mt2corwan(mm,li,iatom),mm=1,2*li-1)
            enddo
          enddo          
          
          write(iun,*) 'mt2cw'
          do iatom=1, natom
            isort=is(iatom)            
            do li=1, lmb(isort)+1
              write(iun,'(2i4,4x,1000i4)') iatom,li,
     $          (mt2cw(mm,li,iatom),mm=1,2*li-1)
            enddo
          enddo          
          write(iun,*) 'corr'
          do iatom=1, natom
            isort=is(iatom)            
!            write(iun,'(i4,4x,1000(2i4,2x))') iatom,
!     $        ((corr(ll,iatom), at_bonding(ll,iatom)),
!     $        ll=0, lmb(isort))
          enddo          
        else

          write(iun,*) 'mt2corwan'
          do iatom=1, natom
            isort=is(iatom)            
            do li=1, 2*lmb(isort)+1
              write(iun,'(2i4,4x,1000i4)') iatom,li,
     $          (mt2corwan(mm,li,iatom),mm=1,((li+1)/2)*2)
            enddo
          enddo
          
          write(iun,*) 'mt2cw'
          do iatom=1, natom
            isort=is(iatom)            
            do li=1, 2*lmb(isort)+1
              write(iun,'(2i4,4x,1000i4)') iatom,li,
     $          (mt2cw(mm,li,iatom),mm=1,((li+1)/2)*2)
            enddo
          enddo          

          write(iun,*) 'mt2cw_rel_l'
          do iatom=1, natom
            isort=is(iatom)            
            do li=1, lmb(isort)+1
              write(iun,'(2i4,4x,1000i4)') iatom,li,
     $          (mt2cw_rel_l(mm,li,iatom),mm=1,4*lmb(isort)+2)
            enddo
          enddo
          
          write(iun,*) 'corr_rel'
          do iatom=1, natom
            isort=is(iatom)            
!            write(iun,'(i4,4x,1000(2i4,2x))') iatom,
!     $        ((corr_rel(li,iatom),at_bonding_rel(li,iatom)),
!     $        li=1, 2*lmb(isort)+1)
          enddo
          
          write(iun,*) 'corr'
          do iatom=1, natom
            isort=is(iatom)            
!            write(iun,'(i4,4x,1000(2i4,2x))') iatom,
!     $        ((corr(ll,iatom),at_bonding(ll,iatom)),
!     $        ll=0, lmb(isort))
          enddo          

        endif
        write(iun,*) 'natom_cw', natom_cw        
        write(iun,*) 'n_cw_max', n_cw_max        

        write(iun,*) 'n_cw_atom'
        do iatom=1, natom
          write(iun,'(2i4)') iatom,n_cw_atom(iatom)
        enddo
        
        write(iun,*) 'iatom_cw'
        do iat=1, natom_cw
          write(iun,'(2i4)') iat,iatom_cw(iat)
        enddo        

        write(iun,*) 'cw_ind'
        do iat=1, natom_cw
          iatom=iatom_cw(iat)          
          write(iun,'(i4,4x,1000i4)') iat,(cw_ind(ii,iat),
     $      ii=1, n_cw_atom(iatom))
        enddo        

        write(iun,*) 'at_orb_cw'
        do iwan=1,num_orb_cor
          write(iun,'(3i4)') iwan,
     $      (at_orb_cw(ii, iwan), ii=1, 2)
        enddo
        
        do iat=1, natom_cw
          iatom=iatom_cw(iat)
          isort=is(iatom)          
          
          do l=0, lmb(isort)
            
            if ((corr(l,iatom) .eq. 2)
     $        .and.
     $        (at_bonding(l,iatom) .eq. 0)
     $        ) then
              
              if (nrel .eq. 1) then
                mudim=l+l+1
                do m1=1, mudim
                  iorb=find_ind(mt2cw(m1,l+1,iatom),
     $              n_cw_atom(iatom),cw_ind(1:n_cw_atom(iatom),iat))
                  write(iun,'(5i5)') iat, l, m1, mt2cw(m1,l+1,iatom),
     $              iorb
                enddo
              else
                mudim=4*l+2                
                do m1=1, mudim
                  iorb=find_ind(mt2cw_rel_l(m1,l+1,iatom),
     $              n_cw_atom(iatom),cw_ind(1:n_cw_atom(iatom),iat))                  
                  write(iun,'(5i5)') iat, l, m1,
     $              mt2cw_rel_l(m1,l+1,iatom),iorb
                enddo
              endif
            endif
          enddo
        enddo                

        do iat=1, natom_ciw
          iatom=iatom_ciw(iat)
          isort=is(iatom)          
          
          do l=0, lmb(isort)
            
            if ((corr_iter(l,iatom) .eq. 2)
     $        .and.
     $        (at_bonding(l,iatom) .eq. 0)
     $        ) then
              
              if (nrel .eq. 1) then
                mudim=l+l+1
                do m1=1, mudim
                  iorb=find_ind(mt2ciw(m1,l+1,iatom),
     $              n_ciw_atom(iatom),ciw_ind(1:n_ciw_atom(iatom),iat))
                  write(iun,'(5i5)') iat, l, m1, mt2ciw(m1,l+1,iatom),
     $              iorb
                enddo
              else
                mudim=4*l+2                
                do m1=1, mudim
                  iorb=find_ind(mt2ciw_rel_l(m1,l+1,iatom),
     $              n_ciw_atom(iatom),ciw_ind(1:n_ciw_atom(iatom),iat))                  
                  write(iun,'(5i5)') iat, l, m1,
     $              mt2ciw_rel_l(m1,l+1,iatom),iorb
                enddo
              endif
            endif
          enddo
        enddo

      endif

      
      end


      subroutine read_wann_inip


      use wannier90_mod
      use parallel_mod      
      use manager_mod      
      use solid_mod      
      use atom_mod      
      use units_mod   

      implicit none
      include "mpif.h"      
      integer :: iwan,iorb,ii,jorb,iatom,lval,mval1
      double precision :: tempvec1(3),tempvec2(3),mval2,ival

      if (maswrk) then

        
        if (irel .le. 1) then
          num_alimj=3
        else
          num_alimj=4
        endif

        open(unit=111, file=adjustl(trim(wanpath))
     $    //'/'//trim(seed_name)//'.inip') ! initial projector
        allocate(nalimj(num_wann))
        do iwan=1, num_wann
          nalimj(iwan)=1
        enddo
        max_nalimj=maxval(nalimj)              
        allocate(alimj(num_alimj,max_nalimj,num_wann))
        alimj=0.0d0
        allocate(alimj_coeff(max_nalimj,num_wann))
        alimj_coeff=0.0d0
        allocate(alimj_xaxis(3,max_nalimj,num_wann))
        alimj_coeff=0.0d0
        allocate(alimj_zaxis(3,max_nalimj,num_wann))
        alimj_coeff=0.0d0                        

        do iwan=1, num_wann
          if (irel .le. 1) then
            read(111, *) 
     $        iatom, lval, mval1,
     $        (alimj_xaxis(ii, 1,iwan), ii=1, 3),
     $        (alimj_zaxis(ii, 1,iwan), ii=1, 3)
            alimj(1,1,iwan)=iatom
            alimj(2,1,iwan)=lval
            alimj(3,1,iwan)=mval1+lval+1
          else
            read(111, *) 
     $        iatom, lval, ival,mval2,
     $        (alimj_xaxis(ii, 1,iwan), ii=1, 3),
     $        (alimj_zaxis(ii, 1,iwan), ii=1, 3)
            alimj(1,1,iwan)=iatom
            alimj(2,1,iwan)=lval
            alimj(3,1,iwan)=nint(ival*2.0d0)
            alimj(4,1,iwan)=nint(mval2*2.0d0)
          endif          
          do iorb=1, nalimj(iwan)
            alimj_coeff(iorb,iwan)=dcmplx(1.0d0, 0.0d0)
          enddo
        enddo
        close(111)

      endif

      call mpi_bcast(num_alimj,1,mpi_integer,
     $  0,mpi_comm_world,mpierr)

      if (.not. maswrk) allocate(nalimj(num_wann))      
      call mpi_bcast(nalimj,num_wann,mpi_integer,
     $  0,mpi_comm_world,mpierr)
      
      if (.not. maswrk) then

        max_nalimj=maxval(nalimj)            
        allocate(alimj(num_alimj,max_nalimj,num_wann))
        alimj=0.0d0
        allocate(alimj_coeff(max_nalimj,num_wann))
        alimj_coeff=0.0d0
        allocate(alimj_xaxis(3,max_nalimj,num_wann))
        alimj_xaxis=0.0d0
        allocate(alimj_zaxis(3,max_nalimj,num_wann))
        alimj_zaxis=0.0d0                        
      endif
      call mpi_bcast(num_wann,1,mpi_integer,
     $  0,mpi_comm_world,mpierr)
      call mpi_bcast(alimj,num_alimj*max_nalimj*num_wann,mpi_integer,
     $  0,mpi_comm_world,mpierr)
      call mpi_bcast(alimj_coeff,max_nalimj*num_wann,
     $  mpi_double_complex,0,mpi_comm_world,mpierr)
      call mpi_bcast(alimj_xaxis,3*max_nalimj*num_wann,
     $  mpi_double_precision,0,mpi_comm_world,mpierr)
      call mpi_bcast(alimj_zaxis,3*max_nalimj*num_wann,
     $  mpi_double_precision,0,mpi_comm_world,mpierr)

      end
