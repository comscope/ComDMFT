      subroutine wannier_mmn(ispin)
      use comwann_mod

      implicit none
      include "mpif.h"
      integer, intent(in) :: ispin
      integer :: k,k0,kq,iq0,
     &  kq0,ind_k,iatom
!     $  ndimb
      double precision :: 
     $  bvec(3),fifi_j(maxlfun,maxlfun,0:2*maxb,nsort)
      complex*16 :: 
     &  znew(nfun,num_bands),znew1(nfun,num_bands),
     &  anew1(nbasmpw,num_bands),anew(nbasmpw,num_bands),
     &  ff(maxel,maxel,natom),
     $  m_matrix_debug_0(num_bands,num_bands)
      integer :: irk, iband1,iband2


      
c     ----------------------------------------------------	
      fifi_j=0.0d0
c     --------------------------------------------------------------------
      m_matrix=0.0d0


      do ind_k=1,ndim_kk(me+1)  ! k vector
        m_matrix_debug_0=0.0d0        
        k=n_mpi_kk(me+1)+ind_k
        k0=i_kref(k)
        do iq0=0,nntot          ! b vector
          if (iq0 .eq. 0) then
            bvec=0.0d0
          else
            bvec=matmul(gbas,
     $        (kpt_latt(:,nnlist(k,iq0))+nncell(:,k,iq0)
     $        -kpt_latt(:,k)))
          endif
          znew=0.0d0
          anew=0.0d0
          znew1=0.0d0
          anew1=0.0d0
          ff=0.0d0
          call fifi_j_prepare_onepoint(fifi_j,bvec,ispin)
          
          call integral_band_pw_band_mt_wan(fifi_j,ff,
     &      bvec)

          call sym_z_0(znew,k,z_wan_bnd(1,1,k0),
     &      num_bands,k_group(k),pnt(1,k))
c$$$  if (iq0 .eq. 0) then
c$$$  write(179,'(a, i5, 1000000f12.6)')'new', k, znew
c$$$  write(179,'(a, i5, 1000000f12.6)')'orig',k,
c$$$  $        z_wan_bnd(:,:,k0)
c$$$  endif                    

          call sym_a_2(anew,k,ev_wan_bnd(1,1,k0),
     &      num_bands,k_group(k),k0)
          if (iq0 .eq. 0) then
            kq=k
          else
            kq=nnlist(k,iq0)
          endif
          kq0=i_kref(kq)
          call sym_z_0(znew1,kq,
     $      z_wan_bnd(1,1,kq0),num_bands,k_group(kq),pnt(1,kq))
          call sym_a_2(anew1,kq,
     $      ev_wan_bnd(1,1,kq0),num_bands,k_group(kq),kq0)          
          if (iq0 .eq. 0) then
            call integral_band_pw_band_wan(ff,znew,anew,znew1,anew1,
     &        m_matrix_debug_0,k,bvec,kq,
     $        (/0,0,0/))
          else
            call integral_band_pw_band_wan(ff,znew,anew,znew1,anew1,
     &        m_matrix(1,1,iq0,k),k,bvec,kq,
     $        nncell(:,k,iq0))
          endif
        enddo
        call is_identity_dcmplx
     $    (m_matrix_debug_0,num_bands,num_bands)
      enddo
      call mpi_allreduce_dcmplx(m_matrix,num_bands**2*nqdiv*nntot,
     $  mpi_sum,mpi_comm_world)


      end
      


      subroutine wan_r0_fft(pw_k,ispin,zwan,awan,wan_fft)

      use comwann_mod
      implicit none

      integer, intent(in) :: ispin
      double precision, intent(in) :: pw_k(3)
      complex*16,intent(in) :: 
     &  zwan(nfun,num_wann,nqdiv),
     $  awan(nbasmpw,num_wann,nqdiv)
      complex*16,intent(out) ::wan_fft(num_wann,nrel)

      integer :: iq,ii,iwan,irl,ibas1,q0
      double precision :: zerovec(3),kshift(3)
      complex*16 :: outvec_mt(nfun,nrel), outvec_it(nbasmpw)
      
      wan_fft=0.0d0
      outvec_mt=0.0d0

      zerovec=0.0d0
      
      call mt_pw_projection(pw_k,ispin,zerovec,outvec_mt)      
      
      do iq=1, nqdiv
        q0=i_kref(iq)
        kshift=matmul(transpose(rbas),pnt(:,iq))-pw_k
        outvec_it=0.0d0              


        if (
     $    (dabs(kshift(1)-nint(kshift(1))).gt.1.0d0/dble(ndiv(1))/100d0)
     $    .or.
     $    (dabs(kshift(2)-nint(kshift(2))).gt.1.0d0/dble(ndiv(2))/100d0)
     $    .or.
     $    (dabs(kshift(3)-nint(kshift(3))).gt.1.0d0/dble(ndiv(3))/100d0)
     $    ) then
          cycle
        endif
        
        call it_pw_projection(iq,pw_k,ispin,zerovec,
     $    outvec_it)

        do iwan=1, num_wann        
          do ii=1, nfun
            do irl=1, nrel
              wan_fft(iwan,irl)
     $          = wan_fft(iwan,irl)
     $          + outvec_mt(ii,irl)*zwan(ii,iwan,iq)
            enddo
          enddo
        enddo
        
        do iwan=1, num_wann                
          do irl=1, nrel
            ibas1=nbask(q0)/nrel*(irl-1)            
            do ii=1, nbask(q0)/nrel
              wan_fft(iwan,irl)
     $          = wan_fft(iwan,irl)
     $          + outvec_it(ibas1+ii)*awan(ibas1+ii,iwan,iq)
            enddo
          enddo
        enddo
      enddo
      end
      
      
      subroutine mt_pw_projection(pw_k,ispin,bvec,outvec)
      
      use comwann_mod
      implicit none
      integer, intent(in) :: ispin
      double precision, intent(in) :: pw_k(3), bvec(3)
      complex*16 :: outvec(nfun,nrel)
      integer :: iatom,ind0,isort,ndimb,ir,
     $  l1,i1,mj1,li1,nm1,in1,km1,km2,lf1,mt,
     $  limu1,ie1,jn1,mdim,m1,m2,ss,lm
      double precision :: rv,
     $  bess(0:maxb,0:maxnrad),bessd(0:maxb),norm_k_b_c,
     $  unit_k_b(3),work(0:maxnrad),ylm((maxb+1)**2),
     $  radint,pw_k_c(3),bvec_c(3),theta,phi,
     $  jval,mval,j1val,m1val,j2val,m2val
      complex*16 :: ai,c1
      integer, external :: ind2jm, jm2ind,lmindex,lget
      double precision, external :: dqdall, clebschg
      complex*16,external :: spharmonics

!     bvec in the unit of reciprocal lattice vector
!     pw_k in the unit of reciprocal lattice vector
      
      if (magn .eq. 2)then
c        write(179,*) 'magn .eq. 2 has not been tested yet'
c        call flush(179)
        if (maswrk) then
          write(iun,*) 'magn .eq. 2 has not been tested yet'
          call flush(iun)        
        endif
        call ending_qft
      endif

      outvec=0.0d0
      
      ai=dcmplx(0.0d0,1.0d0)
!     pi=datan2(1.0d0,1.0d0)*4.0d0


      ylm=0.0d0
c$$$      bvec=0.0d0
      
      pw_k_c=matmul(gbas,pw_k)*2.0d0*pi/par      
      bvec_c=matmul(gbas,bvec)*2.0d0*pi/par      

      norm_k_b_c=dsqrt(sum((pw_k_c+bvec_c)**2))
      
      if(norm_k_b_c.gt.1.0d-9) then
        unit_k_b=(pw_k_c+bvec_c)/norm_k_b_c
        call sphharm(unit_k_b(1),unit_k_b(2),unit_k_b(3),
     $    maxb,ylm)
      else
        ylm=0.d0
        ylm(1)=1.0d0/dsqrt(4.0d0*pi)
      endif

      
      if(norm_k_b_c.gt.1.0d-9) then
        theta= dacos(unit_k_b(3))
        phi=datan2(unit_k_b(2),unit_k_b(1))        
      else
        theta=0.0d0
        phi=0.0d0        
      endif
      
      
      outvec=0.0d0      
      do iatom=1,natom
        ind0=io_lem(iatom)-1    ! the starting value of the {n,l,i,mu} for each atom
        isort=is(iatom)
!     nl=lfunm(isort)
        ndimb=nrel*(lmb(isort)+1)**2 ! dimemsion of the angular part
        c1=4.0d0*pi*cdexp(-ai*sum((pw_k_c+bvec_c)*tau(:,iatom)*par))
     $    /dsqrt(amega)/dsqrt(dble(nqdiv))
        do ir=0,nrad(isort)     ! bessel function
          rv=r(ir,isort)
          call BESSR(norm_k_b_c,rv,bess(0,ir),bessd(0),
     $      lmb(isort))
        enddo
        do limu1=1,ndimb        ! {l,i,mu}
          if(irel<=1) then
            l1=lget(limu1)
          else
            call getlimj(limu1,l1,i1,mj1,li1,0)
          endif
          nm1=ntle(l1,isort)    ! radial function index dimension
          do ie1=1,nm1          ! {ie1,jn1} index for radial functions
            in1=1
            if(augm(ie1,l1,isort)/='LOC') in1=2
            do jn1=1,in1
              km1=indbasa(jn1,ie1,limu1,isort)
              km2=ind0+km1
              lf1=lf_isz(km1,isort) ! radial function index
              mt=ind_wf(lf1,isort) ! starting index for the radial grid              

              work=0.0d0
              do ir=0,nrad(isort)
                rv=r(ir,isort)                  
                work(ir)=work(ir)
     $            +gfun(mt+ir,ispin)*bess(l1,ir)
     $            *rv**2*dr(ir,isort)
              enddo
              radint=dqdall(h(isort),work,nrad(isort))                              
              if (irel .le. 1) then
                outvec(km2,1)=outvec(km2,1)
     $            +c1*(-1.0d0*ai)**l1*radint*ylm(limu1)
              else
                do ss=1, 2
                  jval=l1+i1/2.0d0
                  mval=mj1/2.0d0

                  j2val=1.0d0/2.0d0
                  m2val=((ss-1)*2-1)/2.0d0

                  lm=nint(mval-m2val)
                  j1val=dble(l1)
                  m1val=dble(lm)

c$$$                  write(*,'(6f12.6)') jval,mval,j1val,m1val,j2val,m2val

                  if (abs(lm) .le. l1) then
                    outvec(km2,ss)=outvec(km2,ss)
     $                +c1*(-1.0d0*ai)**l1*radint
     $                *spharmonics(l1,lm,theta,phi)
     $                *clebschG(jval,mval,j1val,m1val,j2val,m2val)
                  endif
                enddo
              endif
            enddo              
          enddo
        enddo
      enddo
      end

      subroutine it_pw_projection(kk,pw_k,ispin,bvec,outvec)
      
      use comwann_mod
      implicit none
      integer, intent(in) :: ispin, kk
      double precision, intent(in) :: pw_k(3), bvec(3)
      complex*16 :: outvec(nbasmpw)
      
      integer :: k0,ibas,ind,ia,ib,ic,j,irl,ibas0

      double precision :: kshift(3)
!     complex*16 :: ai,c1

      if (magn .eq. 2) then
c        write(179,*) 'magn .eq. 2 has not been tested yet'
c        call flush(179)
        if (maswrk) then
          write(iun,*) 'magn .eq. 2 has not been tested yet'
          call flush(iun)        
        endif
        call ending_qft
      endif

!     bvec in the unit of reciprocal lattice vector
!     pw_k in the unit of reciprocal lattice vector

!     normalization of the plane wave
      outvec=0.0d0

      k0=i_kref(kk)
      
      kshift=matmul(transpose(rbas),pnt(:,kk))-(pw_k+bvec)
      

      do irl=1, nrel
        ibas0=nbask(k0)/nrel*(irl-1)                  
        do ibas=1,nbask(k0)/nrel
          ind=indgb(ibas,kk)
          ia=igbs(1,ind)+nint(kshift(1))
          ib=igbs(2,ind)+nint(kshift(2))
          ic=igbs(3,ind)+nint(kshift(3))
          j=indplw(ia,ib,ic)        
          if (
     $      (ia.lt.-maxia) .or. (ia.gt.maxia).or.
     $      (ib.lt.-maxib) .or. (ib.gt.maxib).or.
     $      (ic.lt.-maxic) .or. (ic.gt.maxic).or.
     $      (j .le. 0)
     $      ) then
            cycle
          endif
          
          if(complex_ro) then
            outvec(ibas0+ibas)=dcmplx(sovr(j),sovi(j))
     $        /dsqrt(dble(nqdiv))
          else
            outvec(ibas0+ibas)=dcmplx(sovr(j),0.0d0)
     $        /dsqrt(dble(nqdiv))
          endif
        enddo
      enddo
      
      end      
      
c     integer function ind2jm(j2,ind)
c     
c     
c     implicit none
c     integer, intent(in) :: j2,ind
c     
c     ind2jm = (ind-1)*2-j2
c     end


c     integer function jm2ind(j2,jm)
c     
c     !     example:
c     !     j2=3
c     !     -3,-1,1,3 to 1,2,3,4 
c     
c     implicit none
c     integer, intent(in) :: j2,jm
c     jm2ind = (jm+j2)/2+1
c     end

      integer function lm2ind(ll,lm)
      
!     example:
!     j=2
!     -2, -1, 0, 1, 2 to 1, 2, 3, 4, 5
      
      implicit none
      integer, intent(in) :: ll,lm
      lm2ind=lm+ll+1
      end      

c     integer function lmindex(ll,mm)
c     
c     !     example:
c     !     ll=2, mm=0 -> 1+3+2=6
c     !     
c     
c     implicit none
c     integer, intent(in) :: ll,mm
c     integer, external :: lm2ind
c     
c     lmindex=ll*2+lm2ind(ll,mm)
c     
c     end


      subroutine index_shift(vec,dimmat, outvec)

C$$$  from 1,2,3,4,5,6-> 0,1,2,-3,-2,-1
C$$$  from 1,2,3,4,5 -> 0,1,2,-2,-1
      implicit none

      integer, intent(in) :: vec(3), dimmat(3)
      integer, intent(out) :: outvec(3)

      integer :: n_crit,ii

      outvec=0

      do ii=1, 3
        n_crit=dimmat(ii)/2
        if ((dimmat(ii)-n_crit*2) .ne. 0) then
          n_crit=n_crit+1
        endif
        if (vec(ii) .le. n_crit) then
          outvec(ii)=vec(ii)-1
        else
          outvec(ii)=vec(ii)-dimmat(ii)-1
        endif
      enddo
      end




      subroutine fftshift3d(Nfft,fftboxx,fftboxxout)
      implicit none
      integer, intent(in) :: Nfft(3)
      complex*16, dimension(Nfft(1),Nfft(2),Nfft(3)), intent(in)
     $  :: fftboxx
      complex*16, dimension(Nfft(1),Nfft(2),Nfft(3)), intent(out)
     $  :: fftboxxout
      integer :: ii, jj, kk
      
      fftboxxout=0.0d0
      
      do ii=1, Nfft(1)
        do jj=1, Nfft(2)
          do kk=1, Nfft(3)
            if (ii .le. (Nfft(1)+1)/2) then
              if (jj .le. (Nfft(2)+1)/2) then
                if (kk .le. (Nfft(3)+1)/2) then
                  fftboxxout(ii+nfft(1)/2, jj+nfft(2)/2, kk+nfft(3)/2)
     $              =fftboxx(ii, jj, kk)
                else
                  fftboxxout(ii+nfft(1)/2, jj+nfft(2)/2,
     $              kk-(nfft(3)+1)/2)
     $              =fftboxx(ii, jj, kk)
                end if
              else
                if (kk .le. (Nfft(3)+1)/2) then
                  fftboxxout(ii+nfft(1)/2, jj-(nfft(2)+1)/2,
     $              kk+nfft(3)/2)
     $              =fftboxx(ii, jj, kk)
                else
                  fftboxxout(ii+nfft(1)/2, jj-(nfft(2)+1)/2,
     $              kk-(nfft(3)+1)/2)
     $              =fftboxx(ii, jj, kk)
                end if
              end if
            else
              if (jj .le. (Nfft(2)+1)/2) then
                if (kk .le. (Nfft(3)+1)/2) then
                  fftboxxout(ii-(nfft(1)+1)/2, jj+nfft(2)/2,
     $              kk+nfft(3)/2)
     $              =fftboxx(ii, jj, kk)
                else
                  fftboxxout(ii-(nfft(1)+1)/2, jj+nfft(2)/2,
     $              kk-(nfft(3)+1)/2)
     $              =fftboxx(ii, jj, kk)
                end if
              else
                if (kk .le. (Nfft(3)+1)/2) then
                  fftboxxout(ii-(nfft(1)+1)/2,
     $              jj-(nfft(2)+1)/2, kk+nfft(3)/2)
     $              =fftboxx(ii, jj, kk)
                else
                  fftboxxout(ii-(nfft(1)+1)/2, jj-(nfft(2)+1)/2,
     $              kk-(nfft(3)+1)/2)=fftboxx(ii, jj, kk)
                end if
              end if
            end if
          end do
        end do
      end do
      
      end subroutine fftshift3d




      double precision function f_coefficient(lval,mval,indnum)
      
      implicit none
      integer,intent(in) :: lval,mval,indnum
      
      f_coefficient=0.0d0
      if (indnum .eq. 1) then
        f_coefficient
     $    =-dsqrt(
     $    dble((lval+mval+1)*(lval+mval+2))
     $    /dble((2*lval+1)*(2*lval+3))
     $    )
      elseif(indnum .eq. 2) then
        f_coefficient
     $    =dsqrt(
     $    dble((lval-mval)*(lval-mval-1))
     $    /dble((2*lval-1)*(2*lval+1))
     $    )
        
      elseif(indnum .eq. 3) then
        f_coefficient
     $    =dsqrt(
     $    dble((lval-mval+1)*(lval-mval+2))
     $    /dble((2*lval+1)*(2*lval+3))
     $    )
        
      elseif(indnum .eq. 4) then
        f_coefficient
     $    =-dsqrt(
     $    dble((lval+mval)*(lval+mval-1))
     $    /dble((2*lval-1)*(2*lval+1))
     $    )
        
      elseif(indnum .eq. 5) then
        f_coefficient
     $    =dsqrt(
     $    dble((lval-mval+1)*(lval+mval+1))
     $    /dble((2*lval+1)*(2*lval+3))
     $    )
        
      elseif(indnum .eq. 6) then
        f_coefficient
     $    =dsqrt(
     $    dble((lval-mval)*(lval+mval))
     $    /dble((2*lval-1)*(2*lval+1))
     $    )
      endif
      end


      
      double precision function momentum_angular_spherical
     $  (l1,m1,l2,m2,idir)
      implicit none
      integer,intent(in) :: l1,m1,l2,m2, idir
      double precision, external :: f_coefficient

      momentum_angular_spherical=0.0d0
      if (idir .eq. 1) then
        if ((l1.eq.l2+1) .and. (m1.eq.m2+1)) then
          momentum_angular_spherical=f_coefficient(l2,m2,1)
        elseif ((l2.eq.l1+1) .and. (m2.eq.m1-1)) then
          momentum_angular_spherical=f_coefficient(l2,m2,2)
        endif
      elseif (idir .eq. -1) then
        if ((l1.eq.l2+1) .and. (m1.eq.m2-1)) then
          momentum_angular_spherical=f_coefficient(l2,m2,3)
        elseif ((l2.eq.l1+1) .and. (m2.eq.m1+1)) then
          momentum_angular_spherical=f_coefficient(l2,m2,4)
        endif
      elseif (idir .eq. 0) then
        if ((l1.eq.l2+1) .and. (m1.eq.m2)) then
          momentum_angular_spherical=f_coefficient(l2,m2,5)
        elseif ((l2.eq.l1+1) .and. (m2.eq.m1)) then
          momentum_angular_spherical=f_coefficient(l2,m2,6)
        endif
      endif
      end


      subroutine momentum_angular_cartesian_nospinorbit
     $  (lmbmax,lmbmax_mdim,angular_mommat)
      implicit none
      integer,intent(in) :: lmbmax,lmbmax_mdim
      complex*16,intent(out) ::
     $  angular_mommat(0:lmbmax,lmbmax_mdim,0:lmbmax,lmbmax_mdim,3)

      integer :: idir,l1,m1,l2,m2,idirp,m1p,m2p,
     $  idir_ind,m1_ind,m2_ind,idirp_ind,m1p_ind,m2p_ind
      complex*16 :: rotmat_dir(3,3), ai,
     $  angular_mommat1(0:lmbmax,lmbmax_mdim,0:lmbmax,lmbmax_mdim,3),
     $  angular_mommat2(0:lmbmax,lmbmax_mdim,0:lmbmax,lmbmax_mdim,3)

      complex*16, allocatable :: rotmat_cmplx2real(:,:)
      
      double precision, external :: momentum_angular_spherical

      ai=dcmplx(0.0d0,1.0d0)
      
      rotmat_dir(:,1)=(/1.0d0,0.0d0,1.0d0/)/2.0d0 ! x
      rotmat_dir(:,2)=(/1.0d0,0.0d0,-1.0d0/)/2.0d0/ai ! y
      rotmat_dir(:,3)=(/0.0d0,1.0d0,0.0d0/) ! z

!     rotation of basis from spherical harmonics to real harmonics
      
      
      angular_mommat1=0.0d0
      angular_mommat2=0.0d0
      angular_mommat=0.0d0            
      
!     rotation of the basis
      
      do l2=0,lmbmax
        allocate(rotmat_cmplx2real(2*l2+1,2*l2+1))
        call cal_rotmat_cmplx2real(2*l2+1, rotmat_cmplx2real)
        do l1=0,lmbmax
          if (abs(l1-l2) .eq. 1) then
            
            do idir=-1,1
              idir_ind=idir+2
              
              do m1=-l1,l1
                m1_ind=m1+l1+1        
                do m2=-l2,l2
                  m2_ind=m2+l2+1
                  do m2p=-l2,l2
                    m2p_ind=m2p+l2+1
                    
                    angular_mommat1(l1,m1_ind,l2,m2_ind,idir_ind)
     $                =angular_mommat1(l1,m1_ind,l2,m2_ind,idir_ind)
     $                +momentum_angular_spherical(l1,m1,l2,m2p,idir)
     $                *dconjg(rotmat_cmplx2real(m2_ind,m2p_ind))
                  enddo
                enddo
              enddo
            enddo
          endif
        enddo
        deallocate(rotmat_cmplx2real)
      enddo


!     rotation of the basis
      do l1=0,lmbmax
        allocate(rotmat_cmplx2real(2*l1+1,2*l1+1))        
        call cal_rotmat_cmplx2real(2*l1+1,rotmat_cmplx2real)        
        do l2=0,lmbmax
          if (abs(l1-l2) .eq. 1) then          
            do idir=-1,1
              idir_ind=idir+2
              do m1=-l1,l1
                m1_ind=m1+l1+1
                do m1p=-l1,l1
                  m1p_ind=m1p+l1+1            
                  do m2=-l2,l2
                    m2_ind=m2+l2+1
                    angular_mommat2(l1,m1_ind,l2,m2_ind,idir_ind)
     $                =angular_mommat2(l1,m1_ind,l2,m2_ind,idir_ind)
     $                +angular_mommat1(l1,m1p_ind,l2,m2_ind,idir_ind)
     $                *rotmat_cmplx2real(m1_ind,m1p_ind)
                  enddo
                enddo
              enddo
            enddo
          endif
        enddo
        deallocate(rotmat_cmplx2real)
      enddo


!     rotation of the basis
      do idir=-1,1
        idir_ind=idir+2
        do idirp=-1,1
          idirp_ind=idirp+2          
          do l1=0,lmbmax
            do m1=-l1,l1
              m1_ind=m1+l1+1
              do l2=0,lmbmax
                do m2=-l2,l2
                  m2_ind=m2+l2+1
                  angular_mommat(l1,m1_ind,l2,m2_ind,idir_ind)
     $              =angular_mommat(l1,m1_ind,l2,m2_ind,idir_ind)
     $              +angular_mommat2(l1,m1_ind,l2,m2_ind,idirp_ind)
     $              *rotmat_dir(idirp_ind,idir_ind)
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo


      
      end


      subroutine momentum_mt_nospinorbit(ispin,momentum_mt)
      
      use comwann_mod
      implicit none
      integer,intent(in) :: ispin
      complex*16,intent(out) :: momentum_mt(maxel,maxel,natom,3)

      integer :: lmbmax,lmbmax_mdim,ndimb,iatom,isort,lm1,lm2,nm1,nm2,
     $  ie1,in1,jn1,km1,lf1,ie2,in2,jn2,km2,lf2,idir,l1,m1,l2,m2
      double precision :: radial_mt(maxlfun,maxlfun,nsort)
      complex*16 :: ai
      complex*16, allocatable  :: angular_mt(:,:,:,:,:)
      
!     MT region
      ai=dcmplx(0.0d0,1.0d0)

      call momentum_radial_nospinorbit(ispin,radial_mt)
c      write(179,*) 'momentum_mt_radial done'
c      call flush(179)      

      lmbmax=maxval(lmb)
      lmbmax_mdim=2*lmbmax+1
      
      allocate(angular_mt(0:lmbmax,lmbmax_mdim,0:lmbmax,lmbmax_mdim,3))
      angular_mt=0.0d0
      
      call momentum_angular_cartesian_nospinorbit
     $  (lmbmax,lmbmax_mdim,angular_mt(0,1,0,1,1))
c      write(179,*) 'momentum_mt_angular done'
c      call flush(179)

      

      momentum_mt=0.0d0
      
      do iatom=1,natom
        isort=is(iatom)
        ndimb=nrel*(lmb(isort)+1)**2
        do lm1=1,ndimb
          do lm2=1,ndimb
            call getlm(lm1,l1,m1,0)
            call getlm(lm2,l2,m2,0)              
            nm1=ntle(l1,isort)
            nm2=ntle(l2,isort)            
            if ((l1 .eq. l2+1) .or. (l1 .eq. l2-1)) then        
              do ie1=1,nm1
                in1=1
                if(augm(ie1,l1,isort)/='LOC') in1=2
                do jn1=1,in1
                  km1=indbasa(jn1,ie1,lm1,isort)
                  lf1=lf_isz(km1,isort)
                  do ie2=1,nm2
                    in2=1
                    if(augm(ie2,l2,isort)/='LOC') in2=2
                    do jn2=1,in2
                      km2=indbasa(jn2,ie2,lm2,isort)
                      lf2=lf_isz(km2,isort)
                      do idir=1, 3
                        momentum_mt(km1,km2,iatom,idir)
     $                    =momentum_mt(km1,km2,iatom,idir)
     $                    -ai*angular_mt(l1,m1,l2,m2,idir)
     $                    *radial_mt(lf1,lf2,isort)
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            endif
          enddo
        enddo
      enddo
      deallocate(angular_mt)
c      write(179,*) 'momentum_mt_radial_angular done'
c      call flush(179)

      end

      
      subroutine moment_operator_nospinorbit(ispin,mom_mat)
      
      use comwann_mod      
      implicit none


      integer, intent(in) :: ispin
      complex*16, intent(out) :: mom_mat(num_bands,num_bands,nqdiv,3)
      integer :: idir,iq,k0,iatom,isort,ndimb,ind0,
     $  jbas,jnd,ibas,ind,ia1,ib1,ic1,j,ii
      double precision :: k_g2
      complex*16 :: mom_mat_mt(maxel,maxel,natom,3),
     $  mom_mat1(num_bands,num_bands),
     $  mom_mat2(num_bands,num_bands),
     $  tmp1(maxel,num_bands),tmp(nbndtm,num_bands),
     $  s(nbndtm,nbndtm)



!     pi=4.0d0*datan2(1.0d0,1.0d0)

      mom_mat=0.0d0

      if (irel .gt. 1) return      

c      write(179,*) 'momentum_mt started'
c      call flush(179)      
      
      call momentum_mt_nospinorbit(ispin,mom_mat_mt)
c      write(179,*) 'momentum_mt done'
c      call flush(179)
      
      do idir=1, 3
        do iq=1, nqdiv
          mom_mat1=0.0d0
          mom_mat2=0.0d0          
          k0=i_kref(iq)          
!     mt contribution
          do iatom=1, natom
            isort=is(iatom)
            ndimb=lfunm(isort)
            ind0=io_lem(iatom)
            tmp1=0.0d0          
            call zgemm('n','n',ndimb,num_bands,ndimb,
     &        (1.d0,0.d0),mom_mat_mt(1,1,iatom,idir),maxel,
     &        lte_bnd(ind0,1,iq),nfun,(0.d0,0.d0),tmp1,maxel)
            call zgemm('c','n',num_bands,num_bands,ndimb,
     &        (1.d0,0.d0),lte_bnd(ind0,1,iq),nfun,
     &        tmp1,maxel,(1.d0,0.d0),mom_mat1,
     &        num_bands)
          enddo
          
!     it contribution
          s=0.0d0
          do jbas=1,nbask(k0)/nrel !G1
            jnd=indgb(jbas,iq)
            do ibas=1,nbask(k0)/nrel !G2
              ind=indgb(ibas,iq)
              
              k_g2=0.0d0
              do ii=1, 3
                k_g2=k_g2+gbas(idir,ii)*(kpt_latt(ii,iq)+igbs(ii,ind))
     $            *2.0d0*pi/par
              enddo

              
              ia1=-igbs(1,jnd)+igbs(1,ind)
              ib1=-igbs(2,jnd)+igbs(2,ind)
              ic1=-igbs(3,jnd)+igbs(3,ind)
              j=indplw(ia1,ib1,ic1)              
              
              if(complex_ro) then
                s(jbas,ibas)=dcmplx(sovr(j),sovi(j))*k_g2
              else
                s(jbas,ibas)=dcmplx(sovr(j),0.d0)*k_g2
              endif
            enddo
          enddo
          
          tmp=0.0d0
          call zgemm('n','n',nbask(k0),num_bands,nbask(k0),
     &      (1.d0,0.d0),s,nbndtm,gks_bnd(1,1,iq),nbasmpw,
     &      (0.d0,0.d0),tmp,nbndtm)
          call zgemm('c','n',num_bands,num_bands,
     &      nbask(k0),(1.d0,0.d0),gks_bnd(1,1,iq),nbasmpw,
     &      tmp,nbndtm,(0.d0,0.d0),mom_mat2,
     &      num_bands)
          mom_mat(:,:,iq,idir)
     $      =((mom_mat1+mom_mat2)+transpose(dconjg(mom_mat1+mom_mat2)))
     $      /2.0d0

        enddo
      enddo
c      write(179,*) 'momentum_it done'
c      call flush(179)      
      
      end      

      

      subroutine momentum_radial_nospinorbit(ispin,radial)
      
      use comwann_mod      
      implicit none
      
      integer, intent(in) :: ispin
      double precision,intent(out) ::radial(maxlfun,maxlfun,nsort)
      
      integer :: isort,ir,lf1,mt1,l1,lf2,mt2,l2
!     double precision :: betta
!     complex*16 :: ai

      double precision, allocatable :: r1dr(:),r2dr(:),
     $  fun1(:),fun2(:),fun3(:),fun4(:)
      
      double precision, external :: dqdall

!     MT region
c      write(179,*) 'momentum_mt_radial started'
c      call flush(179)      

!     ai=dcmplx(0.0d0, 1.0d0)
      
      radial=0.0d0
      if (irel .gt. 1) return
      
      do isort=1,nsort
c        write(179,*) isort
c        call flush(179)        
!     ndimb=(lmb(isort)+1)**2
c        write(179,*) isort
c        call flush(179)                
        
        allocate(r2dr(0:nrad(isort)))
        allocate(r1dr(0:nrad(isort)))
        
        r2dr=0.0d0
        r1dr=0.0d0

c        write(179,*) isort
c        call flush(179)                
        do ir=0,nrad(isort)
          r2dr(ir)=r(ir,isort)**2*dr(ir,isort)
          r1dr(ir)=r(ir,isort)*dr(ir,isort)
        enddo
c        write(179,*) isort
c        call flush(179)                

        allocate(fun1(0:nrad(isort)))
        allocate(fun2(0:nrad(isort)))
        allocate(fun3(0:nrad(isort)))
        allocate(fun4(0:nrad(isort)))
c        write(179,*) isort
c        call flush(179)
        
        do lf1=1,lfun(isort)
          mt1=ind_wf(lf1,isort)
          l1=li_val(lf1,isort)-1
          do lf2=1,lfun(isort)
            mt2=ind_wf(lf2,isort)
            l2=li_val(lf2,isort)-1
            
            
c            write(179,'(3i5,f12.6)') isort,lf1,lf2,
c     $        r(nrad(isort),isort)
c            call flush(179)
            
            fun1=0.0d0
            fun2=0.0d0
            fun3=0.0d0
            fun4=0.0d0
            do ir=0,nrad(isort)
              fun1(ir)=gfun(mt1+ir,ispin)
              fun2(ir)=gfun(mt2+ir,ispin)
            enddo
            call fderiv_dble(1,nrad(isort)+1,r(0,isort),fun2(0),fun3(0))
            if (l1 .eq. l2+1) then
              fun4=fun1*fun3*r2dr-l2*fun1*fun2*r1dr
            elseif (l1 .eq. l2-1) then
              fun4=fun1*fun3*r2dr+(l1+2)*fun1*fun2*r1dr
            endif
            radial(lf1,lf2,isort)
     $        =radial(lf1,lf2,isort)
     $        +dqdall(h(isort),fun4,nrad(isort))
            if (irel .eq. 1) then

              fun1=0.0d0
              fun2=0.0d0
              fun3=0.0d0
              fun4=0.0d0
              do ir=0,nrad(isort)
                fun1(ir)=gfund(mt1+ir,ispin)/clight
                fun2(ir)=gfund(mt2+ir,ispin)/clight
              enddo
              call fderiv_dble(1,nrad(isort)+1,
     $          r(0,isort),fun2(0),fun3(0))
              if (l1 .eq. l2+1) then
                fun4=fun1*fun3*r2dr-l2*fun1*fun2*r1dr
              elseif (l1 .eq. l2-1) then
                fun4=fun1*fun3*r2dr+(l1+2)*fun1*fun2*r1dr
              endif
              radial(lf1,lf2,isort)
     $          =radial(lf1,lf2,isort)
     $          -dqdall(h(isort),fun4,nrad(isort))              
            endif
          enddo
        enddo

        deallocate(r1dr,r2dr)
        deallocate(fun1,fun2,fun3,fun4)
c        write(179,'(i5, a)') isort,'done'
c        call flush(179)        
      enddo

      end      

      


      subroutine wannier_realgrid(rrdiv,ispin)
      use comwann_mod            

      implicit none
      include "mpif.h"

      integer, intent(in) :: ispin,rrdiv(3)
      integer :: irk,ndim_fft(nproc),
     $  n_mpi_fft(nproc),ind_g,ind_gg,
     $  ind3(3),ind3_shifted(3),
     $  iwan,i1,i2,i3,ntot,dim3(3),ii,jj,ir
      double precision :: pw_k(3),fac,
     $  wan_norm_r,
     $  wan_center_r(3),
!     $  wan_spreading_r(num_wann),
     $  wan_r2,r2,rvec(3)

      complex*16 :: 
     $  wan_tmp(num_wann,nrel),
     $  wan_out(rrdiv(1)*ndiv(1),rrdiv(2)*ndiv(2),rrdiv(3)*ndiv(3),
     $  nrel),
     $  wan_kgrid(rrdiv(1)*ndiv(1),rrdiv(2)*ndiv(2),rrdiv(3)*ndiv(3),
     $  nrel,num_wann),
     $  wan_rgrid(rrdiv(1)*ndiv(1),rrdiv(2)*ndiv(2),rrdiv(3)*ndiv(3),
     $  nrel,num_wann)
      character*10 :: ch, ch_spin


      ntot=ndiv(1)*ndiv(2)*ndiv(3)*rrdiv(1)*rrdiv(2)*rrdiv(3)
      dim3=(/ndiv(1)*rrdiv(1),ndiv(2)*rrdiv(2),ndiv(3)*rrdiv(3)/)
      
      call size_shift_par(ntot,nproc,ndim_fft,n_mpi_fft)

      fac=1.0/dsqrt(amega)/dsqrt(dble(nqdiv))
      wan_rgrid=0.0d0
      wan_kgrid=0.0d0      
      do ind_g=1, ndim_fft(me+1)
        ind_gg=n_mpi_fft(me+1)+ind_g

        call indexing(ntot,3,dim3,-1,ind_gg,ind3)
        call index_shift(ind3,dim3,ind3_shifted)
        pw_k=dble(ind3_shifted)/dble(ndiv)
        wan_tmp=0.0d0
        call wan_r0_fft(pw_k,ispin,lte_wan,gks_wan,wan_tmp)
        do iwan=1, num_wann
          do ir=1, nrel
            wan_kgrid(ind3(1),ind3(2),ind3(3),ir,iwan)
     $        =wan_tmp(iwan,ir)*fac
          enddo
        enddo
      enddo


      
      call mpi_allreduce_dcmplx(wan_kgrid,ntot*num_wann*nrel,
     $  mpi_sum,mpi_comm_world)
      
      wan_rgrid=wan_kgrid
      
      call fft3(dim3(1),dim3(2),dim3(3),num_wann*nrel,wan_rgrid,1)
      
      if (maswrk) then
        
        write(iun,*) 'wan_info, norm, center, and spreading'
        


!     wan_spreading_r=0.0d0                    
        do iwan=1, num_wann
          wan_norm_r=0.0d0
          wan_center_r=0.0d0          
          wan_r2=0.0d0                  
          do ir=1, nrel
            call fftshift3d(dim3,wan_rgrid(1,1,1,ir,iwan),
     $        wan_out(1,1,1,ir))
!     normalization test
            wan_norm_r=wan_norm_r
     $        +sum(cdabs(wan_rgrid(:,:,:,ir,iwan))**2)
     $        *amega/dble(rrdiv(1)*rrdiv(2)*rrdiv(3))
!     center

            do i1=1, dim3(1)
              do i2=1, dim3(2)
                do i3=1, dim3(3)
                  ind3=(/i1,i2,i3/)
                  call index_shift(ind3,dim3,ind3_shifted)
                  rvec=matmul(rbas, ind3_shifted)/dble(rrdiv)*par*bohr
                  r2=sum((rvec-wann_centers(:,iwan))**2)
                  
                  do ii=1,3
                    wan_center_r(ii)=wan_center_r(ii)
     $                +rvec(ii)*cdabs(wan_rgrid(i1,i2,i3,ir,iwan))**2
     $                *amega/dble(rrdiv(1)*rrdiv(2)*rrdiv(3))
                  enddo
                  wan_r2=wan_r2
     $              +r2*cdabs(wan_rgrid(i1,i2,i3,ir,iwan))**2
     $              *amega/dble(rrdiv(1)*rrdiv(2)*rrdiv(3))
                enddo
              enddo
            enddo
          enddo
!     spreading
          
          write(iun,'(i5,f12.6,3x,3f12.6,3x,f12.6)')
     $      iwan, wan_norm_r,
     $      (wan_center_r(ii),ii=1, 3),
     $      wan_r2
          
          if (iwan .lt. 10) then
            write(ch,'(i1)') iwan
          elseif (iwan .lt. 100) then
            write(ch,'(i2)') iwan
          elseif (iwan .lt. 1000) then
            write(ch,'(i3)') iwan
          endif

          
          do ir=1, nrel
            if (nrel .eq. 1) then
              write(ch_spin,'(a)') 'real'              
            else
              if (ir .eq. 1) then
                write(ch_spin,'(a)') 'dn_real'
              else
                write(ch_spin,'(a)') 'up_real'              
              endif
            endif
            
            open(unit=992, file=trim(seed_name)//'_'//trim(ch)//'_'
     $        //trim(ch_spin)//'.xsf')
            write(992,'(a)') "CRYSTAL"
            write(992,'(a)') "PRIMVEC"
            write(992,'(3(f12.7,3x))') (rbas(ii,1)*par*bohr, ii=1,3)
            write(992,'(3(f12.7,3x))') (rbas(ii,2)*par*bohr, ii=1,3)
            write(992,'(3(f12.7,3x))') (rbas(ii,3)*par*bohr, ii=1,3)
            write(992,'(a)') "CONVVEC"
            write(992,'(3(f12.7,3x))') (rbas(ii,1)*par*bohr, ii=1,3)
            write(992,'(3(f12.7,3x))') (rbas(ii,2)*par*bohr, ii=1,3)
            write(992,'(3(f12.7,3x))') (rbas(ii,3)*par*bohr, ii=1,3)
            write(992,'(a)') "PRIMCOORD"
            write(992,'(2(i5, 2x))') natom, 1
            do jj=1, natom
              write(992,'(i5, 2x, 3(f12.6,3x))') nint(z(is(jj))),
     $          (tau(ii,jj)*par*bohr, ii=1,3)
            enddo
            
            write(992,'(a)') "BEGIN_BLOCK_DATAGRID_3D"                        
            write(992,'(a)') "wannier wavefunction"
            write(992,'(a)') "BEGIN_DATAGRID_3D"
            write(992,'(3i5)') dim3
            write(992,'(3(f20.12,2x))')
     $        -rbas(:,1)*par*bohr*(ndiv(1)/2.0d0+0.5d0/dim3(1))
     $        -rbas(:,2)*par*bohr*(ndiv(2)/2.0d0+0.5d0/dim3(2))
     $        -rbas(:,3)*par*bohr*(ndiv(3)/2.0d0+0.5d0/dim3(3))
            write(992,'(3(f20.12,2x))')
     $        (rbas(ii,1)*par*bohr*ndiv(1),ii=1,3) 
            write(992,'(3(f20.12,2x))') 
     $        (rbas(ii,2)*par*bohr*ndiv(2),ii=1,3)              
            write(992,'(3(f20.12,2x))')
     $        (rbas(ii,3)*par*bohr*ndiv(3),ii=1,3)
            write(992,*)
     $        (((dble(wan_out(i1,i2,i3,ir)),i1=1,dim3(1)),i2=1,dim3(2))
     $        ,i3=1,dim3(3))
            
            write(992,'(a)') "END_DATAGRID_3D"
            
            write(992,'(a)') "END_BLOCK_DATAGRID_3D"            
            close(992)

            if (nrel .eq. 1) then
              write(ch_spin,'(a)') 'imag'              
            else
              if (ir .eq. 1) then
                write(ch_spin,'(a)') 'dn_imag'
              else
                write(ch_spin,'(a)') 'up_imag'              
              endif
            endif            
            
            
            open(unit=992, file=trim(seed_name)//'_'//trim(ch)//'_'
     $        //trim(ch_spin)//'.xsf')
            write(992,'(a)') "CRYSTAL"
            write(992,'(a)') "PRIMVEC"
            write(992,'(3(f12.7,3x))') (rbas(ii,1)*par*bohr, ii=1,3)
            write(992,'(3(f12.7,3x))') (rbas(ii,2)*par*bohr, ii=1,3)
            write(992,'(3(f12.7,3x))') (rbas(ii,3)*par*bohr, ii=1,3)
            write(992,'(a)') "CONVVEC"
            write(992,'(3(f12.7,3x))') (rbas(ii,1)*par*bohr, ii=1,3)
            write(992,'(3(f12.7,3x))') (rbas(ii,2)*par*bohr, ii=1,3)
            write(992,'(3(f12.7,3x))') (rbas(ii,3)*par*bohr, ii=1,3)
            write(992,'(a)') "PRIMCOORD"
            write(992,'(2(i5, 2x))') natom, 1
            do jj=1, natom
              write(992,'(i5, 2x, 3(f12.6,3x))') nint(z(is(jj))),
     $          (tau(ii,jj)*par*bohr, ii=1,3)
            enddo
            
            write(992,'(a)') "BEGIN_BLOCK_DATAGRID_3D"                        
            write(992,'(a)') "wannier wavefunction"
            write(992,'(a)') "BEGIN_DATAGRID_3D"
            write(992,'(3i5)') dim3
            write(992,'(3(f20.12,2x))')
     $        -rbas(:,1)*par*bohr*(ndiv(1)/2.0d0+0.5d0/dim3(1))
     $        -rbas(:,2)*par*bohr*(ndiv(2)/2.0d0+0.5d0/dim3(2))
     $        -rbas(:,3)*par*bohr*(ndiv(3)/2.0d0+0.5d0/dim3(3))
            write(992,'(3(f20.12,2x))')
     $        (rbas(ii,1)*par*bohr*ndiv(1),ii=1,3) 
            write(992,'(3(f20.12,2x))') 
     $        (rbas(ii,2)*par*bohr*ndiv(2),ii=1,3)              
            write(992,'(3(f20.12,2x))')
     $        (rbas(ii,3)*par*bohr*ndiv(3),ii=1,3)
            write(992,*)
     $        (((dimag(wan_out(i1,i2,i3,ir)),i1=1,dim3(1)),i2=1,dim3(2))
     $        ,i3=1,dim3(3))
            
            write(992,'(a)') "END_DATAGRID_3D"
            
            write(992,'(a)') "END_BLOCK_DATAGRID_3D"            
            close(992)
          enddo
        enddo
      endif
      

      end


      
      subroutine cal_pw_wan(numk,ndim_kpath,n_mpi_kpath,
     $  ndim_kpath_max,pw_vec,wan_rgrid,pw_wan)

      use comwann_mod

      implicit none
      include "mpif.h"

      integer,intent(in) :: numk,ndim_kpath(nproc),n_mpi_kpath(nproc),
     $  ndim_kpath_max
      double precision, intent(in) :: pw_vec(3,numk) ! in gbas
      complex*16,intent(in) :: wan_rgrid(nrdiv(1)*ndiv(1),
     $  nrdiv(2)*ndiv(2),nrdiv(3)*ndiv(3),num_wann)
      complex*16,intent(out) :: pw_wan(3,num_wann,ndim_kpath_max)
      
      integer :: ind3(3),dim3(3),ind3_shifted(3),ind_k,k,i1,i2,i3,
     $  idir,iwan
      double precision :: pw_vec_c(3)
      complex*16 :: ai


      ai=dcmplx(0.0d0,1.0d0)
!     pi=datan2(1.0d0,1.0d0)*4.0d0
      
      pw_wan=0.0d0

      dim3=ndiv*nrdiv

      do ind_k=1,ndim_kpath(me+1) ! k vector
        k=n_mpi_kpath(me+1)+ind_k          
        pw_vec_c=matmul(gbas,pw_vec(:,k))*2.0d0*pi/par
        
        
        do i1=1,dim3(1)
          do i2=1,dim3(2)
            do i3=1,dim3(3)
              ind3=(/i1,i2,i3/)
              call index_shift(ind3,dim3,ind3_shifted)                    

c     do idir=1,3
c     rvec(idir)      ! in atomic unit
c     $            =(
c     $            rbas(idir,1)/dble(nrdiv(1))*ind3_shifted(1)
c     $            +rbas(idir,2)/dble(nrdiv(2))*ind3_shifted(2)
c     $            +rbas(idir,3)/dble(nrdiv(3))*ind3_shifted(3)
c     $            )*par
c     enddo
              
              do iwan=1,num_wann 
                do idir=1,3
                  pw_wan(idir,iwan,ind_k)
     $              =pw_wan(idir,iwan,ind_k)
     $              +wan_rgrid(i1,i2,i3,iwan)
     $              *pw_vec_c(idir)
     $              *cdexp(-ai*2.0d0*pi
     $              *sum(pw_vec(:,k)/dble(nrdiv)*dble(ind3_shifted))
     $              )
     $              /dsqrt(dble(nqdiv)) !! normalization check
                enddo
              enddo      
            enddo
          enddo
        enddo
      enddo
      
      end



      subroutine is_identity_dcmplx(matin,matsize,matsize0)

      use comwann_mod

      implicit none
      
      integer,intent(in) :: matsize, matsize0
      complex*16,intent(in) :: matin(matsize0,matsize0)
      integer :: ii, jj
      double precision :: diagerr,offdiagerr


      diagerr=1.0d-3
      offdiagerr=1.0d-3      
      do ii=1, matsize
        do jj=1, matsize        
          if (ii .eq. jj) then
            if (cdabs(matin(ii,jj)-1.0d0) .gt. diagerr) then
              if (maswrk) then
                write(iun,'(a, 2i5, 2f12.6)') 'identity_check fail',
     $            ii,jj,cdabs(matin(ii,jj)),cdabs(matin(ii,jj)-1.0d0)
              endif
c              write(179,'(a, 2i5, 2f12.6)') 'identity_check fail',
c     $          ii, jj,cdabs(matin(ii,jj)),cdabs(matin(ii,jj)-1.0d0)
              call ending
            endif
          else
            if (cdabs(matin(ii,jj)) .gt. offdiagerr) then
              if (maswrk) then
                write(iun,'(a, 2i5, f12.6)') 'identity_check fail',
     $            ii, jj, cdabs(matin(ii,jj))
              endif
c              write(179,'(a, 2i5, f12.6)') 'identity_check fail',
c     $          ii, jj, cdabs(matin(ii,jj))
              call ending              
            endif
          endif
        enddo
      enddo
      end



      complex*16 function spharmonics(l, m, theta, phi)
!     with Condon-Shortley phase
      implicit none
      integer, intent(in) :: l, m
      double precision, intent(in) :: theta, phi
      integer :: lmm, lpm
      double precision :: pi, costheta
      double precision, external :: asso_lp,ifactorial
      complex*16 :: ai

      pi=4.0d0*datan2(1.0d0, 1.0d0)
      ai=dcmplx(0.0d0, 1.0d0)

c$$$      theta= dacos(vec(3)/dsqrt(sum(vec**2)))
c$$$      phi=datan2(vec(2), vec(1))
c$$$      write(*,*) theta, phi

      lmm=l-m
      lpm=l+m
      costheta=dcos(theta)
      spharmonics=dsqrt((2*l+1.0d0)/(4.0d0*pi)
     $  *ifactorial(lmm)/ifactorial(lpm))
     $  *asso_lp(l,m,costheta)*cdexp(ai*m*phi)
      
      end


      double precision FUNCTION asso_lp(l,m,x)

      implicit none
      integer, intent(in) :: l, m
      double precision, intent(in) :: x
      integer :: absm,lmm,lpm
      double precision, external :: plgndr, ifactorial

      if (m .ge. 0) then
        asso_lp=plgndr(l,m,x)
      else
        absm=-m
        lmm=l-absm
        lpm=l+absm
        asso_lp=(-1)**absm*ifactorial(lmm)/ifactorial(lpm)
     $    *plgndr(l,absm,x)
      endif
      end
      

      double precision FUNCTION plgndr(l,m,x)
!     from numerical recipe
      
      implicit none
      INTEGER, intent(in) ::  l,m
      double precision, intent(in) :: x      
      INTEGER :: i,ll
      double precision  :: fact,pll,pmm,pmmp1,somx2
      pmm=1.0d0 !Compute P mm .
      if(m.gt.0) then
        somx2=dsqrt((1.0d0-x)*(1.0d0+x))
        fact=1.0d0
        do i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.0d0
        enddo
      endif
      if(l.eq.m) then
        plgndr=pmm
      else
        pmmp1=x*(2*m+1)*pmm     !Compute P m m+1.
        if(l.eq.m+1) then
          plgndr=pmmp1
        else                    !Compute P m l , l>m + 1.
          pll=0.0d0
          do ll=m+2,l
            pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
            pmm=pmmp1
            pmmp1=pll
          enddo
          plgndr=pll
        endif
      endif
      END      


